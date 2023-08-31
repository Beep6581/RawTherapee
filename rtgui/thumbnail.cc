/*
 *  This file is part of RawTherapee.
 *
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifdef _WIN32
#include <windows.h>
#endif

#include "cachemanager.h"
#include "multilangmgr.h"
#include "thumbnail.h"
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include "../rtengine/colortemp.h"
#include "../rtengine/imagedata.h"
#include "../rtengine/procparams.h"
#include "../rtengine/rtthumbnail.h"
#include <glib/gstdio.h>
#include <glibmm/timezone.h>

#include "../rtengine/dynamicprofile.h"
#include "../rtengine/metadata.h"
#include "../rtengine/profilestore.h"
#include "../rtengine/settings.h"
#include "guiutils.h"
#include "batchqueue.h"
#include "extprog.h"
#include "md5helper.h"
#include "pathutils.h"
#include "paramsedited.h"
#include "ppversion.h"
#include "procparamchangers.h"
#include "version.h"


namespace {

bool CPBDump(
    const Glib::ustring& commFName,
    const Glib::ustring& imageFName,
    const Glib::ustring& profileFName,
    const Glib::ustring& defaultPParams,
    const CacheImageData* cfs,
    bool flagMode
)
{
    const std::unique_ptr<Glib::KeyFile> kf(new Glib::KeyFile);

    if (!kf) {
        return false;
    }

    // open the file in write mode
    const std::unique_ptr<FILE, decltype(&std::fclose)> f(g_fopen(commFName.c_str (), "wt"), &std::fclose);

    if (!f) {
        printf ("CPBDump(\"%s\") >>> Error: unable to open file with write access!\n", commFName.c_str());
        return false;
    }

    try {
        kf->set_string ("RT General", "CachePath", options.cacheBaseDir);
        kf->set_string ("RT General", "AppVersion", RTVERSION);
        kf->set_integer ("RT General", "ProcParamsVersion", PPVERSION);
        kf->set_string ("RT General", "ImageFileName", imageFName);
        kf->set_string ("RT General", "OutputProfileFileName", profileFName);
        kf->set_string ("RT General", "DefaultProcParams", defaultPParams);
        kf->set_boolean ("RT General", "FlaggingMode", flagMode);

        kf->set_integer ("Common Data", "FrameCount", cfs->frameCount);
        kf->set_integer ("Common Data", "SampleFormat", cfs->sampleFormat);
        kf->set_boolean ("Common Data", "IsHDR", cfs->isHDR);
        kf->set_boolean ("Common Data", "IsPixelShift", cfs->isPixelShift);
        kf->set_double ("Common Data", "FNumber", cfs->fnumber);
        kf->set_double ("Common Data", "Shutter", cfs->shutter);
        kf->set_double ("Common Data", "FocalLength", cfs->focalLen);
        kf->set_integer ("Common Data", "ISO", cfs->iso);
        kf->set_string ("Common Data", "Lens", cfs->lens);
        kf->set_string ("Common Data", "Make", cfs->camMake);
        kf->set_string ("Common Data", "Model", cfs->camModel);

    } catch (const Glib::KeyFileError&) {
    }

    try {
        fprintf (f.get(), "%s", kf->to_data().c_str());
    } catch (const Glib::KeyFileError&) {
    }

    return true;
}

} // namespace

using namespace rtengine::procparams;

Thumbnail::Thumbnail(CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf) :
    fname(fname),
    cfs(*cf),
    cachemgr(cm),
    ref(1),
    enqueueNumber(0),
    tpp(nullptr),
    pparams(new ProcParams),
    pparamsValid(false),
    imageLoading(false),
    lastImg(nullptr),
    lastW(0),
    lastH(0),
    lastScale(0),
    initial_(false)
{

    loadProcParams ();

    // should be safe to use the unprotected version of loadThumbnail, since we are in the constructor
    _loadThumbnail ();
    generateExifDateTimeStrings ();

    if (cfs.rankOld >= 0) {
        // rank and inTrash were found in cache (old style), move them over to pparams

        // try to load the last saved parameters from the cache or from the paramfile file
        createProcParamsForUpdate(false, false); // this can execute customprofilebuilder to generate param file

        // TODO? should we call notifylisterners_procParamsChanged here?

        setRank(cfs.rankOld);
        setStage(cfs.inTrashOld);
    }

    delete tpp;
    tpp = nullptr;
}

Thumbnail::Thumbnail(CacheManager* cm, const Glib::ustring& fname, const std::string& md5, const std::string &xmpSidecarMd5) :
    fname(fname),
    cachemgr(cm),
    ref(1),
    enqueueNumber(0),
    tpp(nullptr),
    pparams(new ProcParams),
    pparamsValid(false),
    imageLoading(false),
    lastImg(nullptr),
    lastW(0),
    lastH(0),
    lastScale(0.0),
    initial_(true)
{


    cfs.md5 = md5;
    cfs.xmpSidecarMd5 = xmpSidecarMd5;
    loadProcParams ();
    _generateThumbnailImage ();
    cfs.recentlySaved = false;

    initial_ = false;

    delete tpp;
    tpp = nullptr;
}

Glib::ustring Thumbnail::xmpSidecarPath(const Glib::ustring &imagePath)
{
    return rtengine::Exiv2Metadata::xmpSidecarPath(imagePath);
}

void Thumbnail::_generateThumbnailImage ()
{

    //  delete everything loaded into memory
    delete tpp;
    tpp = nullptr;
    delete [] lastImg;
    lastImg = nullptr;
    tw = options.maxThumbnailWidth;
    th = options.maxThumbnailHeight;
    imgRatio = -1.;

    // generate thumbnail image
    const std::string ext = getExtension(fname).lowercase();

    if (ext.empty()) {
        return;
    }

    cfs.supported = false;
    cfs.exifValid = false;
    cfs.timeValid = false;

    if (ext == "jpg" || ext == "jpeg") {
        infoFromImage (fname);
        tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, -1, pparams->wb.equal, pparams->wb.observer);

        if (tpp) {
            cfs.format = FT_Jpeg;
        }
    } else if (ext == "png") {
        tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, -1, pparams->wb.equal, pparams->wb.observer);

        if (tpp) {
            cfs.format = FT_Png;
        }
    } else if (ext == "tif" || ext == "tiff") {
        infoFromImage (fname);
        tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, -1, pparams->wb.equal, pparams->wb.observer);

        if (tpp) {
            cfs.format = FT_Tiff;
        }
    } else {
        // RAW works like this:
        //  1. if we are here it's because we aren't in the cache so load the JPG
        //     image out of the RAW. Mark as "quick".
        //  2. if we don't find that then just grab the real image.
        bool quick = false;

        rtengine::eSensorType sensorType = rtengine::ST_NONE;
        if ( initial_ && options.internalThumbIfUntouched) {
            quick = true;
            tpp = rtengine::Thumbnail::loadQuickFromRaw (fname, sensorType, tw, th, 1, TRUE);
        }

        if ( tpp == nullptr ) {
            quick = false;
            tpp = rtengine::Thumbnail::loadFromRaw (fname, sensorType, tw, th, 1, pparams->wb.equal, pparams->wb.observer, TRUE);
        }

        cfs.sensortype = sensorType;
        if (tpp) {
            cfs.format = FT_Raw;
            cfs.thumbImgType = quick ? CacheImageData::QUICK_THUMBNAIL : CacheImageData::FULL_THUMBNAIL;
            infoFromImage (fname);
            if (!quick) {
                cfs.width = tpp->full_width;
                cfs.height = tpp->full_height;
            }
        }
    }

    if (tpp) {
        tpp->getAutoWBMultipliers(cfs.redAWBMul, cfs.greenAWBMul, cfs.blueAWBMul);
        _saveThumbnail ();
        cfs.supported = true;

        cfs.save (getCacheFileName ("data", ".txt"));

        generateExifDateTimeStrings ();
    }
}

bool Thumbnail::isSupported () const
{
    return cfs.supported;
}

const ProcParams& Thumbnail::getProcParams ()
{
    MyMutex::MyLock lock(mutex);
    return getProcParamsU();
}

// Unprotected version of getProcParams, when
const ProcParams& Thumbnail::getProcParamsU ()
{
    if (pparamsValid) {
        return *pparams;
    } else {
        *pparams = *(ProfileStore::getInstance()->getDefaultProcParams (getType() == FT_Raw));

        if (pparams->wb.method == "Camera") {
            double ct;
            getCamWB (ct, pparams->wb.green, pparams->wb.observer);
            pparams->wb.temperature = ct;
        } else if (pparams->wb.method == "autold") {
            double ct;
            getAutoWB (ct, pparams->wb.green, pparams->wb.equal, pparams->wb.observer, pparams->wb.tempBias);
            pparams->wb.temperature = ct;
        }
    }

    return *pparams; // there is no valid pp to return, but we have to return something
}

/** @brief  Create default params on demand and returns a new updatable object
 *
 *  The loaded profile may be partial, but it return a complete ProcParams (i.e. without ParamsEdited)
 *
 *  @param returnParams Ask to return a pointer to a ProcParams object if true
 *  @param force True if the profile has to be re-generated even if it already exists
 *  @param flaggingMode True if the ProcParams will be created because the file browser is being flagging an image
 *                      (rang, to trash, color labels). This parameter is passed to the CPB.
 *
 *  @return Return a pointer to a ProcPamas structure to be updated if returnParams is true and if everything went fine, NULL otherwise.
 */
rtengine::procparams::ProcParams* Thumbnail::createProcParamsForUpdate(bool returnParams, bool force, bool flaggingMode)
{
    // try to load the last saved parameters from the cache or from the paramfile file
    ProcParams* ldprof = nullptr;

    Glib::ustring defProf = getType() == FT_Raw ? options.defProfRaw : options.defProfImg;

    const CacheImageData* cfs = getCacheImageData();
    Glib::ustring defaultPparamsPath = options.findProfilePath(defProf);
    const bool create = (!hasProcParams() || force);
    const bool run_cpb = !options.CPBPath.empty() && !defaultPparamsPath.empty() && cfs && cfs->exifValid && create;

    const Glib::ustring outFName =
        (options.paramsLoadLocation == PLL_Input && options.saveParamsFile) ?
        fname + paramFileExtension :
        getCacheFileName("profiles", paramFileExtension);

    if (!run_cpb) {
        if (defProf == DEFPROFILE_DYNAMIC && create && cfs && cfs->exifValid) {
            const auto pp_deleter =
                [](PartialProfile* pp)
                {
                    pp->deleteInstance();
                    delete pp;
                };
            const std::unique_ptr<const rtengine::FramesMetaData> imageMetaData(rtengine::FramesMetaData::fromFile(fname));
            const std::unique_ptr<PartialProfile, decltype(pp_deleter)> pp(
                imageMetaData
                    ? ProfileStore::getInstance()->loadDynamicProfile(imageMetaData.get(), fname)
                    : nullptr,
                pp_deleter
            );
            if (pp && !pp->pparams->save(outFName)) {
                loadProcParams();
            }
        } else if (create && defProf != DEFPROFILE_DYNAMIC) {
            const PartialProfile* const p = ProfileStore::getInstance()->getProfile(defProf);
            if (p && !p->pparams->save(outFName)) {
                loadProcParams();
            }
        }
    } else {
        // First generate the communication file, with general values and EXIF metadata
        static int index = 0; // Will act as unique identifier during the session
        Glib::ustring tmpFileName( Glib::build_filename(options.cacheBaseDir, Glib::ustring::compose("CPB_temp_%1.txt", index++)) );

        CPBDump(tmpFileName, fname, outFName,
                defaultPparamsPath == DEFPROFILE_INTERNAL ? DEFPROFILE_INTERNAL : Glib::build_filename(defaultPparamsPath, Glib::path_get_basename(defProf) + paramFileExtension), cfs, flaggingMode);

        // For the filename etc. do NOT use streams, since they are not UTF8 safe
        Glib::ustring cmdLine = options.CPBPath + Glib::ustring(" \"") + tmpFileName + Glib::ustring("\"");

        if (rtengine::settings->verbose) {
            printf("Custom profile builder's command line: %s\n", Glib::ustring(cmdLine).c_str());
        }

        bool success = ExtProgStore::spawnCommandSync (cmdLine);

        // Now they SHOULD be there (and potentially "partial"), so try to load them and store it as a full procparam
        if (success) {
            loadProcParams();
        }

        g_remove (tmpFileName.c_str ());
    }

    if (returnParams && hasProcParams()) {
        ldprof = new ProcParams ();
        *ldprof = getProcParams ();
    }

    return ldprof;
}

void Thumbnail::notifylisterners_procParamsChanged(int whoChangedIt)
{
    for (size_t i = 0; i < listeners.size(); i++) {
        listeners[i]->procParamsChanged (this, whoChangedIt);
    }
}

/*
 * Load the procparams from the cache or from the sidecar file (priority set in
 * the Preferences).
 *
 * The result is a complete ProcParams with default values merged with the values
 * from the loaded ProcParams (sidecar or cache file).
*/
void Thumbnail::loadProcParams()
{
    MyMutex::MyLock lock(mutex);

    pparamsValid = false;
    pparams->setDefaults();

    if (options.paramsLoadLocation == PLL_Input) {
        // try to load it from params file next to the image file
        const int ppres = pparams->load(fname + paramFileExtension);
        pparamsValid = !ppres && pparams->ppVersion >= 220;

        // if no success, try to load the cached version of the procparams
        if (!pparamsValid) {
            pparamsValid = !pparams->load(getCacheFileName("profiles", paramFileExtension));
        }
    } else {
        // try to load it from cache
        pparamsValid = !pparams->load(getCacheFileName("profiles", paramFileExtension));

        // if no success, try to load it from params file next to the image file
        if (!pparamsValid) {
            const int ppres = pparams->load(fname + paramFileExtension);
            pparamsValid = !ppres && pparams->ppVersion >= 220;
        }
    }
}

void Thumbnail::clearProcParams (int whoClearedIt)
{

    /*  Clarification on current "clear profile" functionality:
        a. if rank/colorlabel/inTrash are NOT set,
        the "clear profile" will delete the pp3 file (as before).

        b. if any of the rank/colorlabel/inTrash ARE set,
        the "clear profile" will lead to execution of ProcParams::setDefaults
        (the CPB is NOT called) to set the params values and will preserve
        rank/colorlabel/inTrash in the param file. */

    {
        MyMutex::MyLock lock(mutex);

        // preserve rank, colorlabel and inTrash across clear
        int rank = getRank();
        int colorlabel = getColorLabel();
        int inTrash = getStage();


        cfs.recentlySaved = false;
        pparamsValid = false;

        //TODO: run though customprofilebuilder?
        // probably not as this is the only option to set param values to default

        // reset the params to defaults
        pparams->setDefaults();

        // and restore rank and inTrash
        setRank(rank);
        pparamsValid = cfs.rating != rank;
        setColorLabel(colorlabel);
        setStage(inTrash);

        // params could get validated by rank/inTrash values restored above
        if (pparamsValid) {
            updateCache();
        } else {
            // remove param file from cache
            Glib::ustring fname_ = getCacheFileName ("profiles", paramFileExtension);
            g_remove (fname_.c_str ());

            // remove param file located next to the file
            fname_ = fname + paramFileExtension;
            g_remove (fname_.c_str ());

            fname_ = removeExtension(fname) + paramFileExtension;
            g_remove (fname_.c_str ());

            if (cfs.format == FT_Raw && options.internalThumbIfUntouched && cfs.thumbImgType != CacheImageData::QUICK_THUMBNAIL) {
                // regenerate thumbnail, ie load the quick thumb again. For the rare formats not supporting quick thumbs this will
                // be a bit slow as a new full thumbnail will be generated unnecessarily, but currently there is no way to pre-check
                // if the format supports quick thumbs.
                initial_ = true;
                _generateThumbnailImage();
                initial_ = false;
            }
        }

    } // end of mutex lock

    for (size_t i = 0; i < listeners.size(); i++) {
        listeners[i]->procParamsChanged (this, whoClearedIt);
    }
}

bool Thumbnail::hasProcParams () const
{

    return pparamsValid;
}

void Thumbnail::setProcParams (const ProcParams& pp, ParamsEdited* pe, int whoChangedIt, bool updateCacheNow, bool resetToDefault)
{
    const bool needsReprocessing =
           resetToDefault
        || pparams->toneCurve != pp.toneCurve
        || pparams->locallab != pp.locallab
        || pparams->labCurve != pp.labCurve
        || pparams->localContrast != pp.localContrast
        || pparams->rgbCurves != pp.rgbCurves
        || pparams->colorToning != pp.colorToning
        || pparams->vibrance != pp.vibrance
        || pparams->wb != pp.wb
        || pparams->colorappearance != pp.colorappearance
        || pparams->epd != pp.epd
        || pparams->fattal != pp.fattal
        || pparams->sh != pp.sh
        || pparams->toneEqualizer != pp.toneEqualizer
        || pparams->crop != pp.crop
        || pparams->coarse != pp.coarse
        || pparams->commonTrans != pp.commonTrans
        || pparams->rotate != pp.rotate
        || pparams->distortion != pp.distortion
        || pparams->lensProf != pp.lensProf
        || pparams->perspective != pp.perspective
        || pparams->gradient != pp.gradient
        || pparams->pcvignette != pp.pcvignette
        || pparams->cacorrection != pp.cacorrection
        || pparams->vignetting != pp.vignetting
        || pparams->chmixer != pp.chmixer
        || pparams->blackwhite != pp.blackwhite
        || pparams->icm != pp.icm
        || pparams->hsvequalizer != pp.hsvequalizer
        || pparams->filmSimulation != pp.filmSimulation
        || pparams->softlight != pp.softlight
        || pparams->dehaze != pp.dehaze
        || pparams->filmNegative != pp.filmNegative
        || whoChangedIt == FILEBROWSER
        || whoChangedIt == BATCHEDITOR;

    {
        MyMutex::MyLock lock(mutex);

        if (*pparams != pp) {
            cfs.recentlySaved = false;
        } else if (pparamsValid && !updateCacheNow) {
            // nothing to do
            return;
        }

        // do not update rank, colorlabel and inTrash
        const int rank = getRank();
        const int colorlabel = getColorLabel();
        const int inTrash = getStage();

        if (pe) {
            pe->combine(*pparams, pp, true);
        } else {
            *pparams = pp;
        }

        pparamsValid = true;

        setRank(rank);
        setColorLabel(colorlabel);
        setStage(inTrash);

        if (updateCacheNow) {
            updateCache();
        }
    } // end of mutex lock

    if (needsReprocessing) {
        for (size_t i = 0; i < listeners.size(); i++) {
            listeners[i]->procParamsChanged (this, whoChangedIt);
        }
    }
}

bool Thumbnail::isRecentlySaved () const
{

    return cfs.recentlySaved;
}

void Thumbnail::imageDeveloped ()
{

    cfs.recentlySaved = true;
    cfs.save (getCacheFileName ("data", ".txt"));

    if (options.saveParamsCache) {
        pparams->save (getCacheFileName ("profiles", paramFileExtension));
    }
}

void Thumbnail::imageEnqueued ()
{

    enqueueNumber++;
}

void Thumbnail::imageRemovedFromQueue ()
{

    enqueueNumber--;
}

bool Thumbnail::isEnqueued () const
{

    return enqueueNumber > 0;
}

bool Thumbnail::isPixelShift () const
{
    return cfs.isPixelShift;
}
bool Thumbnail::isHDR () const
{
    return cfs.isHDR;
}

void Thumbnail::increaseRef ()
{
    MyMutex::MyLock lock(mutex);
    ++ref;
}

void Thumbnail::decreaseRef ()
{
    {
        MyMutex::MyLock lock(mutex);

        if ( ref == 0 ) {
            return;
        }

        if ( --ref != 0 ) {
            return;
        }
    }
    cachemgr->closeThumbnail (this);
}

void Thumbnail::getThumbnailSize(int &w, int &h, const rtengine::procparams::ProcParams *pparams)
{
    MyMutex::MyLock lock(mutex);

    int tw_ = tw;
    int th_ = th;

    float imgRatio_ = imgRatio;

    if (pparams) {
        int ppCoarse = pparams->coarse.rotate;

        if (ppCoarse >= 180) {
            ppCoarse -= 180;
        }

        int thisCoarse = this->pparams->coarse.rotate;

        if (thisCoarse >= 180) {
            thisCoarse -= 180;
        }

        if (thisCoarse != ppCoarse) {
            // different orientation -> swapping width & height
            std::swap(th_, tw_);
            if (imgRatio_ >= 0.0001f) {
                imgRatio_ = 1.f / imgRatio_;
            }
        }
    }

    if (imgRatio_ > 0.) {
        w = imgRatio_ * static_cast<float>(h);
    } else {
        w = tw_ * h / th_;
    }

    if (w > options.maxThumbnailWidth) {
        const float s = static_cast<float>(options.maxThumbnailWidth) / w;
        w = options.maxThumbnailWidth;
        h = std::max<int>(h * s, 1);
    }
}

void Thumbnail::getFinalSize (const rtengine::procparams::ProcParams& pparams, int& w, int& h)
{
    MyMutex::MyLock lock(mutex);

    // WARNING: When downscaled, the ratio have loosed a lot of precision, so we can't get back the exact initial dimensions
    double fw = lastW * lastScale;
    double fh = lastH * lastScale;

    if (pparams.coarse.rotate == 90 || pparams.coarse.rotate == 270) {
        fh = lastW * lastScale;
        fw = lastH * lastScale;
    }

    if (!pparams.resize.enabled) {
        w = fw;
        h = fh;
    } else {
        w = (int)(fw + 0.5);
        h = (int)(fh + 0.5);
    }
}

void Thumbnail::getOriginalSize (int& w, int& h) const
{
    w = tw;
    h = th;
}

rtengine::IImage8* Thumbnail::processThumbImage (const rtengine::procparams::ProcParams& pparams, int h, double& scale)
{

    MyMutex::MyLock lock(mutex);

    if ( tpp == nullptr ) {
        _loadThumbnail();

        if ( tpp == nullptr ) {
            return nullptr;
        }
    }

    rtengine::IImage8* image = nullptr;

    if ( cfs.thumbImgType == CacheImageData::QUICK_THUMBNAIL ) {
        // RAW internal thumbnail, no profile yet: just do some rotation etc.
        image = tpp->quickProcessImage (pparams, h, rtengine::TI_Nearest);
    } else {
        // Full thumbnail: apply profile
        // image = tpp->processImage (pparams, h, rtengine::TI_Bilinear, cfs.getCamera(), cfs.focalLen, cfs.focalLen35mm, cfs.focusDist, cfs.shutter, cfs.fnumber, cfs.iso, cfs.expcomp, scale );
        image = tpp->processImage (pparams, static_cast<rtengine::eSensorType>(cfs.sensortype), h, rtengine::TI_Bilinear, &cfs, scale );
    }

    tpp->getDimensions(lastW, lastH, lastScale);

    delete tpp;
    tpp = nullptr;
    return image;
}

rtengine::IImage8* Thumbnail::upgradeThumbImage (const rtengine::procparams::ProcParams& pparams, int h, double& scale)
{

    MyMutex::MyLock lock(mutex);

    if ( cfs.thumbImgType != CacheImageData::QUICK_THUMBNAIL ) {
        return nullptr;
    }

    _generateThumbnailImage();

    if ( tpp == nullptr ) {
        return nullptr;
    }

    // rtengine::IImage8* image = tpp->processImage (pparams, h, rtengine::TI_Bilinear, cfs.getCamera(), cfs.focalLen, cfs.focalLen35mm, cfs.focusDist, cfs.shutter, cfs.fnumber, cfs.iso, cfs.expcomp,  scale );
    rtengine::IImage8* image = tpp->processImage (pparams, static_cast<rtengine::eSensorType>(cfs.sensortype), h, rtengine::TI_Bilinear, &cfs, scale );
    tpp->getDimensions(lastW, lastH, lastScale);

    delete tpp;
    tpp = nullptr;
    return image;
}

void Thumbnail::generateExifDateTimeStrings ()
{
    if (cfs.timeValid) {
        std::string dateFormat = options.dateFormat;
        std::ostringstream ostr;
        bool spec = false;

        for (size_t i = 0; i < dateFormat.size(); i++)
            if (spec && dateFormat[i] == 'y') {
                ostr << cfs.year;
                spec = false;
            } else if (spec && dateFormat[i] == 'm') {
                ostr << (int)cfs.month;
                spec = false;
            } else if (spec && dateFormat[i] == 'd') {
                ostr << (int)cfs.day;
                spec = false;
            } else if (dateFormat[i] == '%') {
                spec = true;
            } else {
                ostr << (char)dateFormat[i];
                spec = false;
            }

        ostr << " " << (int)cfs.hour;
        ostr << ":" << std::setw(2) << std::setfill('0') << (int)cfs.min;
        ostr << ":" << std::setw(2) << std::setfill('0') << (int)cfs.sec;

        dateTimeString = ostr.str ();
        dateTime = Glib::DateTime::create_local(cfs.year, cfs.month, cfs.day,
                                                cfs.hour, cfs.min, cfs.sec);
    }

    if (!dateTime.gobj() || !cfs.timeValid) {
        dateTimeString = "";
        dateTime = Glib::DateTime::create_now_utc(0);
    }

    if (!cfs.exifValid) {
        exifString = "";
        return;
    }

    exifString = Glib::ustring::compose ("f/%1 %2s %3%4 %5mm", Glib::ustring(rtengine::FramesData::apertureToString(cfs.fnumber)), Glib::ustring(rtengine::FramesData::shutterToString(cfs.shutter)), M("QINFO_ISO"), cfs.iso, Glib::ustring::format(std::setw(3), std::fixed, std::setprecision(2), cfs.focalLen));

    if (options.fbShowExpComp && cfs.expcomp != "0.00" && !cfs.expcomp.empty()) { // don't show exposure compensation if it is 0.00EV;old cache files do not have ExpComp, so value will not be displayed.
        exifString = Glib::ustring::compose ("%1 %2EV", exifString, cfs.expcomp);    // append exposure compensation to exifString
    }
}

const Glib::ustring& Thumbnail::getExifString () const
{

    return exifString;
}

const Glib::ustring& Thumbnail::getDateTimeString () const
{

    return dateTimeString;
}

const Glib::DateTime& Thumbnail::getDateTime () const
{

    return dateTime;
}

void Thumbnail::getAutoWB (double& temp, double& green, double equal, rtengine::StandardObserver observer, double tempBias)
{
    if (cfs.redAWBMul != -1.0) {
        rtengine::ColorTemp ct(cfs.redAWBMul, cfs.greenAWBMul, cfs.blueAWBMul, equal, observer);
        temp = ct.getTemp();
        green = ct.getGreen();
    } else {
        temp = green = -1.0;
    }
}


ThFileType Thumbnail::getType () const
{

    return (ThFileType) cfs.format;
}

int Thumbnail::infoFromImage (const Glib::ustring& fname)
{
    return infoFromImage(fname, cfs);
}

int Thumbnail::infoFromImage(const Glib::ustring &fname, CacheImageData &cfs)
{
    std::unique_ptr<rtengine::FramesMetaData> idata(rtengine::FramesMetaData::fromFile (fname));

    if (!idata) {
        return 0;
    }

    int deg = 0;
    cfs.timeValid = false;
    cfs.exifValid = false;

    if (idata->getDateTimeAsTS() > 0) {
        cfs.year         = 1900 + idata->getDateTime().tm_year;
        cfs.month        = idata->getDateTime().tm_mon + 1;
        cfs.day          = idata->getDateTime().tm_mday;
        cfs.hour         = idata->getDateTime().tm_hour;
        cfs.min          = idata->getDateTime().tm_min;
        cfs.sec          = idata->getDateTime().tm_sec;
        cfs.timeValid    = true;
    }

    if (idata->hasExif()) {
        cfs.shutter      = idata->getShutterSpeed ();
        cfs.fnumber      = idata->getFNumber ();
        cfs.focalLen     = idata->getFocalLen ();
        cfs.focalLen35mm = idata->getFocalLen35mm ();
        cfs.focusDist    = idata->getFocusDist ();
        cfs.iso          = idata->getISOSpeed ();
        cfs.expcomp      = idata->expcompToString (idata->getExpComp(), false); // do not mask Zero expcomp
        cfs.isHDR        = idata->getHDR ();
        cfs.isPixelShift = idata->getPixelShift ();
        cfs.frameCount   = idata->getFrameCount ();
        cfs.sampleFormat = idata->getSampleFormat ();
        cfs.lens         = idata->getLens();
        cfs.camMake      = idata->getMake();
        cfs.camModel     = idata->getModel();
        cfs.rating       = idata->getRating();
        cfs.exifValid    = true;

        if (idata->getOrientation() == "Rotate 90 CW") {
            deg = 90;
        } else if (idata->getOrientation() == "Rotate 180") {
            deg = 180;
        } else if (idata->getOrientation() == "Rotate 270 CW") {
            deg = 270;
        }
    } else {
        cfs.lens     = "Unknown";
        cfs.camMake  = "Unknown";
        cfs.camModel = "Unknown";
    }

    // get image filetype
    std::string::size_type idx;
    idx = fname.rfind('.');

    if(idx != std::string::npos) {
        cfs.filetype = fname.substr(idx + 1);
    } else {
        cfs.filetype = "";
    }

    idata->getDimensions(cfs.width, cfs.height);

    return deg;
}

/*
 * Read all thumbnail's data from the cache; build and save them if doesn't exist - NON PROTECTED
 * This includes:
 *  - image's bitmap (*.rtti)
 *  - auto exposure's histogram (full thumbnail only)
 *  - embedded profile (full thumbnail only)
 *  - LiveThumbData section of the data file
 */
void Thumbnail::_loadThumbnail(bool firstTrial)
{

    tw = -1;
    th = options.maxThumbnailHeight;
    delete tpp;
    tpp = new rtengine::Thumbnail ();
    tpp->isRaw = (cfs.format == (int) FT_Raw);

    // load supplementary data
    bool succ = tpp->readData (getCacheFileName ("data", ".txt"));

    if (succ) {
        tpp->getAutoWBMultipliers(cfs.redAWBMul, cfs.greenAWBMul, cfs.blueAWBMul);
    }

    // thumbnail image
    succ = succ && tpp->readImage (getCacheFileName ("images", ""));

    if (!succ && firstTrial) {
        _generateThumbnailImage ();

        if (cfs.supported) {
            _loadThumbnail (false);
        }

        if (tpp == nullptr) {
            return;
        }
    } else if (!succ) {
        delete tpp;
        tpp = nullptr;
        return;
    }

    if ( cfs.thumbImgType == CacheImageData::FULL_THUMBNAIL ) {
        // load embedded profile
        tpp->readEmbProfile (getCacheFileName ("embprofiles", ".icc"));

        tpp->init ();
    }

    if (!initial_) {
        tw = tpp->getImageWidth (getProcParamsU(), th, imgRatio);    // this might return 0 if image was just building
    }
}

/*
 * Save thumbnail's data to the cache - NON PROTECTED
 * This includes:
 *  - image's bitmap (*.rtti)
 *  - auto exposure's histogram (full thumbnail only)
 *  - embedded profile (full thumbnail only)
 *  - LiveThumbData section of the data file
 */
void Thumbnail::_saveThumbnail ()
{

    if (!tpp) {
        return;
    }

    g_remove (getCacheFileName ("images", ".rtti").c_str ());

    // save thumbnail image
    tpp->writeImage (getCacheFileName ("images", ""));

    // save embedded profile
    tpp->writeEmbProfile (getCacheFileName ("embprofiles", ".icc"));

    // save supplementary data
    tpp->writeData (getCacheFileName ("data", ".txt"));
}

/*
 * Save thumbnail's data to the cache - MUTEX PROTECTED
 * This includes:
 *  - image's bitmap (*.rtti)
 *  - auto exposure's histogram (full thumbnail only)
 *  - embedded profile (full thumbnail only)
 *  - LiveThumbData section of the data file
 */
void Thumbnail::saveThumbnail ()
{
    MyMutex::MyLock lock(mutex);
    _saveThumbnail();
}

/*
 * Update the cached files
 *  - updatePParams==true (default)        : write the procparams file (sidecar or cache, depending on the options)
 *  - updateCacheImageData==true (default) : write the CacheImageData values in the cache folder,
 *                                           i.e. some General, DateTime, ExifInfo, File info and ExtraRawInfo,
 */
void Thumbnail::updateCache (bool updatePParams, bool updateCacheImageData)
{

    if (updatePParams && pparamsValid) {
        pparams->save (
            options.saveParamsFile  ? fname + paramFileExtension : "",
            options.saveParamsCache ? getCacheFileName ("profiles", paramFileExtension) : "",
            true
        );
    }

    if (updateCacheImageData) {
        cfs.save (getCacheFileName ("data", ".txt"));
    }

    if (updatePParams && pparamsValid) {
        saveMetadata();
    }
}

Thumbnail::~Thumbnail ()
{
    mutex.lock();

    delete [] lastImg;
    delete tpp;
    mutex.unlock();
}

Glib::ustring Thumbnail::getCacheFileName (const Glib::ustring& subdir, const Glib::ustring& fext) const
{
    return cachemgr->getCacheFileName (subdir, fname, fext, cfs.md5);
}

void Thumbnail::setFileName (const Glib::ustring &fn)
{

    fname = fn;
    cfs.md5 = ::getMD5 (fname);
}

int Thumbnail::getRank  () const
{
    // prefer the user-set rank over the embedded Rating
    // pparams->rank == -1 means that there is no saved rank yet, so we should
    // next look for the embedded Rating metadata.
    if (pparams->rank != -1) {
        return pparams->rank;
    } else {
        return cfs.rating;
    }
}

void Thumbnail::setRank  (int rank)
{
    pparams->rank = rank;
    pparamsValid = true;
}

int Thumbnail::getColorLabel  () const
{
    return pparams->colorlabel;
}

void Thumbnail::setColorLabel  (int colorlabel)
{
    if (pparams->colorlabel != colorlabel) {
        pparams->colorlabel = colorlabel;
        pparamsValid = true;
    }
}

int Thumbnail::getStage () const
{
    return pparams->inTrash;
}

void Thumbnail::setStage (bool stage)
{
    if (pparams->inTrash != stage) {
        pparams->inTrash = stage;
        pparamsValid = true;
    }
}

void Thumbnail::addThumbnailListener (ThumbnailListener* tnl)
{

    increaseRef();
    listeners.push_back (tnl);
}

void Thumbnail::removeThumbnailListener (ThumbnailListener* tnl)
{

    std::vector<ThumbnailListener*>::iterator f = std::find (listeners.begin(), listeners.end(), tnl);

    if (f != listeners.end()) {
        listeners.erase (f);
        decreaseRef();
    }
}

// Calculates the standard filename for the automatically named batch result
// and opens it in OS default viewer
// destination: 1=Batch conf. file; 2=batch out dir; 3=RAW dir
// Return: Success?
bool Thumbnail::openDefaultViewer(int destination)
{

#ifdef _WIN32
    Glib::ustring openFName;

    if (destination == 1) {
        openFName = Glib::ustring::compose ("%1.%2", BatchQueue::calcAutoFileNameBase(fname), options.saveFormatBatch.format);

        if (Glib::file_test (openFName, Glib::FILE_TEST_EXISTS)) {
            wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (openFName.c_str(), -1, NULL, NULL, NULL);
            ShellExecuteW(NULL, L"open", wfilename, NULL, NULL, SW_SHOWMAXIMIZED );
            g_free(wfilename);
        } else {
            printf("%s not found\n", openFName.data());
            return false;
        }
    } else {
        openFName = destination == 3 ? fname
                    : Glib::ustring::compose ("%1.%2", BatchQueue::calcAutoFileNameBase(fname), options.saveFormatBatch.format);

        printf("Opening %s\n", openFName.c_str());

        if (Glib::file_test (openFName, Glib::FILE_TEST_EXISTS)) {
            // Output file exists, so open explorer and select output file
            wchar_t* org = (wchar_t*)g_utf8_to_utf16 (Glib::ustring::compose("/select,\"%1\"", openFName).c_str(), -1, NULL, NULL, NULL);
            wchar_t* par = new wchar_t[wcslen(org) + 1];
            wcscpy(par, org);

            // In this case the / disturbs
            wchar_t* p = par + 1; // skip the first backslash

            while (*p != 0) {
                if (*p == L'/') {
                    *p = L'\\';
                }

                p++;
            }

            ShellExecuteW(NULL, L"open", L"explorer.exe", par, NULL, SW_SHOWNORMAL );

            delete[] par;
            g_free(org);
        } else if (Glib::file_test (Glib::path_get_dirname(openFName), Glib::FILE_TEST_EXISTS)) {
            // Out file does not exist, but directory
            wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (Glib::path_get_dirname(openFName).c_str(), -1, NULL, NULL, NULL);
            ShellExecuteW(NULL, L"explore", wfilename, NULL, NULL, SW_SHOWNORMAL );
            g_free(wfilename);
        } else {
            printf("File and dir not found\n");
            return false;
        }
    }

    return true;

#else
    // TODO: Add more OSes here
    printf("Automatic opening not supported on this OS\n");
    return false;
#endif

}

bool Thumbnail::imageLoad(bool loading)
{
    MyMutex::MyLock lock(mutex);
    bool previous = imageLoading;

    if( loading && !previous ) {
        imageLoading = true;
        return true;
    } else if( !loading ) {
        imageLoading = false;
    }

    return false;
}

void Thumbnail::getCamWB(double& temp, double& green, rtengine::StandardObserver observer) const
{
    if (tpp) {
        tpp->getCamWB  (temp, green, observer);
    } else {
        temp = green = -1.0;
    }
}

void Thumbnail::saveMetadata()
{
    if (options.rtSettings.metadata_xmp_sync != rtengine::Settings::MetadataXmpSync::READ_WRITE) {
        return;
    }

    if (pparams->metadata.exif.empty() && pparams->metadata.iptc.empty()) {
        return;
    }

    auto fn = rtengine::Exiv2Metadata::xmpSidecarPath(fname);
    try {
        auto xmp = rtengine::Exiv2Metadata::getXmpSidecar(fname);
        rtengine::Exiv2Metadata meta;
        meta.xmpData() = std::move(xmp);
        meta.setExif(pparams->metadata.exif);
        meta.setIptc(pparams->metadata.iptc);
        meta.saveToXmp(fn);
        if (options.rtSettings.verbose) {
            std::cout << "saved edited metadata for " << fname << " to "
                      << fn << std::endl;
        }
    } catch (std::exception &exc) {
        std::cerr << "ERROR saving metadata for " << fname << " to " << fn
                  << ": " << exc.what() << std::endl;
    }
}
void Thumbnail::getSpotWB(int x, int y, int rect, double& temp, double& green)
{
    if (tpp) {
        tpp->getSpotWB (getProcParams(), x, y, rect, temp, green);
    } else {
        temp = green = -1.0;
    }
}

void Thumbnail::applyAutoExp (rtengine::procparams::ProcParams& pparams)
{
    if (tpp) {
        tpp->applyAutoExp (pparams);
    }
}

const CacheImageData* Thumbnail::getCacheImageData() const
{
    return &cfs;
}

std::string Thumbnail::getMD5() const
{
    return cfs.md5;
}

bool Thumbnail::isQuick() const
{
    return cfs.thumbImgType == CacheImageData::QUICK_THUMBNAIL;
}

bool Thumbnail::isPParamsValid() const
{
    return pparamsValid;
}
