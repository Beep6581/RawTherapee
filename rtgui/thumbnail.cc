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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "multilangmgr.h"
#include "thumbnail.h"
#include <sstream>
#include <iomanip>
#include "options.h"
#include "../rtengine/mytime.h"
#include <cstdio>
#include <cstdlib>
#include <glibmm.h>
#include "../rtengine/imagedata.h"
#include <glib/gstdio.h>
#include "guiutils.h"
#include "profilestore.h"
#include "batchqueue.h"
#include "../rtengine/safegtk.h"
#include "version.h"

using namespace rtengine::procparams;

Thumbnail::Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf) 
    : fname(fname), cfs(*cf), cachemgr(cm), ref(1), enqueueNumber(0), tpp(NULL),
      pparamsValid(false), needsReProcessing(true),imageLoading(false), lastImg(NULL),
		initial_(false), lastW(0), lastH(0), lastScale(0),pparamsId(-1),idata(NULL) {

    cfs.load (getCacheFileName ("data")+".txt");
    loadInfoFromImage ();
    //loadProcParams ();
    _loadThumbnail ();
    generateExifDateTimeStrings ();

 	delete tpp;
 	tpp = 0;
}

Thumbnail::Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5)
    : fname(fname), cachemgr(cm), ref(1), enqueueNumber(0), tpp(NULL), pparamsValid(false),
      needsReProcessing(true),imageLoading(false), lastImg(NULL),
		initial_(true),pparamsId(-1),idata(NULL) {


    cfs.md5 = md5;
    _generateThumbnailImage ();
    //loadProcParams ();
    cfs.recentlySaved = false;

	initial_ = false;
 
 	delete tpp;
 	tpp = 0;
}

void Thumbnail::_generateThumbnailImage () {

	//  delete everything loaded into memory
	delete tpp;
	tpp = NULL;
	delete [] lastImg;
	lastImg = NULL;
	tw = -1;
	th = options.maxThumbnailHeight;
	imgRatio = -1.;

	// generate thumbnail image
	Glib::ustring ext = getExtension (fname);
	if (ext=="") 
		return;
	cfs.supported = false;
	cfs.exifValid = false;
	cfs.timeValid = false;

	if (ext.lowercase()=="jpg") {
			tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, 1, loadInfoFromImage () );
			if (tpp)
					cfs.format = FT_Jpeg;
	}
	else if (ext.lowercase()=="png") {
			tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, 1, loadInfoFromImage ());
			if (tpp)
					cfs.format = FT_Png;
	}
	else if (ext.lowercase()=="tif" || ext.lowercase()=="tiff") {
			tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, 1, loadInfoFromImage ());
			if (tpp)
					cfs.format = FT_Tiff;
	}
    else {
 		// RAW works like this:
		//  1. if we are here it's because we aren't in the cache so load the JPG
		//     image out of the RAW. Mark as "quick".
		//  2. if we don't find that then just grab the real image.
		bool quick = false;
		rtengine::RawMetaDataLocation ri;
		if ( initial_ && options.internalThumbIfUntouched)
		{
			quick = true;
			tpp = rtengine::Thumbnail::loadQuickFromRaw (fname, ri, tw, th, 1, TRUE);
		}
		if ( tpp == NULL )
		{
			quick = false;
			tpp = rtengine::Thumbnail::loadFromRaw (fname, ri, tw, th, 1, TRUE);
		}
		if (tpp) {
			cfs.format = FT_Raw;
			cfs.thumbImgType = quick ? CacheImageData::QUICK_THUMBNAIL : CacheImageData::FULL_THUMBNAIL;
			loadInfoFromImage ();
		}
	}

    if (tpp)
	{
        _saveThumbnail ();
        cfs.supported = true;
        needsReProcessing = true;

        cfs.save (getCacheFileName ("data")+".txt");
    	if( idata && !idata->isEmbedded() )
    		idata->saveXMP();
        generateExifDateTimeStrings ();
    }
}

bool Thumbnail::isSupported () {
    return cfs.supported;
}

const ProcParams& Thumbnail::getProcParams () {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(mutex);
	#endif

    if (pparamsValid)
        return pparams;
    else {
        pparams = *(profileStore.getDefaultProcParams (getType()==FT_Raw));
        if (pparams.wb.method=="Camera") {
            double ct;
            getCamWB (ct, pparams.wb.green);
            pparams.wb.temperature = ct;
        }
        else if (pparams.wb.method=="Auto") {
            double ct;
            getAutoWB (ct, pparams.wb.green);
            pparams.wb.temperature = ct;
        }
    }
    return pparams; // there is no valid pp to return, but we have to return something
}

/*
 *  Create default params on demand and returns a new updatable object
 *  The loaded profile may be partial, but it return a complete ProcParams (i.e. without ParamsEdited)
 */
rtengine::procparams::ProcParams* Thumbnail::createProcParamsForUpdate(bool returnParams, bool forceCPB) {
    // try to load the last saved parameters from the cache or from the paramfile file
    ProcParams* ldprof = NULL;

    Glib::ustring defProf = getType()==FT_Raw ? options.defProfRaw : options.defProfImg;

    const CacheImageData* cfs=getCacheImageData();
    if (!options.customProfileBuilder.empty() && (!hasProcParams() || forceCPB) && cfs && cfs->exifValid) {
        // For the filename etc. do NOT use streams, since they are not UTF8 safe
        Glib::ustring cmdLine=Glib::ustring("\"") + options.customProfileBuilder + Glib::ustring("\" \"") + fname + Glib::ustring("\" \"")
        + options.rtdir + Glib::ustring("/") + options.profilePath + Glib::ustring("/") + defProf + Glib::ustring(".pp3") + Glib::ustring("\" ");

        // ustring doesn't know int etc formatting, so take these via (unsafe) stream
        std::ostringstream strm;
        strm << cfs->fnumber << Glib::ustring(" ") << cfs->shutter << Glib::ustring(" ");
        strm << cfs->focalLen << Glib::ustring(" ") << cfs->iso << Glib::ustring(" \"");
        strm << cfs->lens << Glib::ustring("\" \"") << cfs->camera << Glib::ustring("\"");
 
        bool success = safe_spawn_command_line_sync (cmdLine + strm.str());

        // Now they SHOULD be there (and potentially "partial"), so try to load them and store it as a full procparam
        if (success) loadInfoFromImage();
    }

    if (returnParams && hasProcParams()) {
        ldprof = new ProcParams ();
        *ldprof = getProcParams ();
    }

    return ldprof;
}

void Thumbnail::notifylisterners_procParamsChanged(int whoChangedIt){
	for (size_t i=0; i<listeners.size(); i++)
		listeners[i]->procParamsChanged (this, whoChangedIt);
}

void Thumbnail::clearProcParams (int whoClearedIt) {
    
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(mutex);
	#endif

    cfs.recentlySaved = false;
    pparamsValid = false;
    needsReProcessing = true;

    //TODO: run though customprofilebuilder?
    // probably not as this is the only option to set param values to default

    // reset the params to defaults
    pparams.setDefaults();

    idata->deleteAllSnapshots();

    for (size_t i=0; i<listeners.size(); i++)
        listeners[i]->procParamsChanged (this, whoClearedIt);
}

bool Thumbnail::hasProcParams () {
    
    return pparamsValid;
}

void Thumbnail::setProcParams (const ProcParams& pp, ParamsEdited* pe, int whoChangedIt, bool updateCacheNow) {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(mutex);
	#endif
    
    if (pparams!=pp) 
        cfs.recentlySaved = false;


    if (pe) {
    	// coarse.rotate works in ADD mode only, so we set it to 0 first
    	if (pe->coarse.rotate)
    		pparams.coarse.rotate = 0;
    	pe->combine(pparams, pp, true);
    }
    else pparams = pp;
    pparamsValid = true;
    needsReProcessing = true;


    if (updateCacheNow)
        updateCache ();

    for (size_t i=0; i<listeners.size(); i++)
        listeners[i]->procParamsChanged (this, whoChangedIt);
}

bool Thumbnail::isRecentlySaved () {
    
    return cfs.recentlySaved;
}

void Thumbnail::imageDeveloped () {
        
    cfs.recentlySaved = true;
    cfs.save (getCacheFileName ("data")+".txt");
}

void Thumbnail::imageEnqueued () {

    enqueueNumber++;
}

void Thumbnail::imageRemovedFromQueue () {

    enqueueNumber--;
}

bool Thumbnail::isEnqueued () {
    
    return enqueueNumber > 0;
}

void Thumbnail::increaseRef ()
{
	Glib::Mutex::Lock lock(mutex);
   	++ref;
}

void Thumbnail::decreaseRef () 
{
	{
		Glib::Mutex::Lock lock(mutex);
		if ( ref == 0 )
		{
			return;
		}
		if ( --ref != 0 )
		{
			return;
		}
	}
	cachemgr->closeThumbnail (this); 
}

void Thumbnail::getThumbnailSize (int &w, int &h) {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(mutex);
	#endif

	w=0;
	if (!initial_ && tpp) w = tpp->getImageWidth (getProcParams(), h, imgRatio);  // this might return 0 if image was just building
	if (w==0) {
		if (imgRatio > 0.)
			w = (int)(imgRatio * (float)h);
		else
			w = tw * h / th;
	}
}

void Thumbnail::getFinalSize (const rtengine::procparams::ProcParams& pparams, int& w, int& h) {
    // TODO: Check for Linux
    #ifdef WIN32
    Glib::Mutex::Lock lock(mutex);
    #endif

    // WARNING: When downscaled, the ratio have loosed a lot of precision, so we can't get back the exact initial dimensions
    double fw = lastW*lastScale;
    double fh = lastH*lastScale;

    if (pparams.coarse.rotate==90 || pparams.coarse.rotate==270) {
        fh = lastW*lastScale;
        fw = lastH*lastScale;
    }
    if (!pparams.resize.enabled) {
        w = fw;
        h = fh;
    }
    else {
        w = (int)(fw+0.5);
        h = (int)(fh+0.5);
    }
}


rtengine::IImage8* Thumbnail::processThumbImage (const rtengine::procparams::ProcParams& pparams, int h, double& scale) {

	Glib::Mutex::Lock lock(mutex);

     if ( tpp == 0 )
 	{
 		_loadThumbnail();
 		if ( tpp == 0 )
 		{
 			return 0;
 		}
 	}
 
 	rtengine::IImage8* image = 0;

	if ( cfs.thumbImgType == CacheImageData::QUICK_THUMBNAIL )
	{
		// RAW internal thumbnail, no profile yet: just do some rotation etc.
 		image = tpp->quickProcessImage (pparams, h, rtengine::TI_Nearest, scale);
	}
	else
	{
		// Full thumbnail: apply profile
 		image = tpp->processImage (pparams, h, rtengine::TI_Bilinear, cfs.camera, scale );
	}
 
    tpp->getDimensions(lastW,lastH,lastScale);

 	delete tpp;
 	tpp = 0;
 	return image;
}

rtengine::IImage8* Thumbnail::upgradeThumbImage (const rtengine::procparams::ProcParams& pparams, int h, double& scale) {

	Glib::Mutex::Lock lock(mutex);

	if ( cfs.thumbImgType != CacheImageData::QUICK_THUMBNAIL )
	{
		return 0;
	}

	_generateThumbnailImage();
 	if ( tpp == 0 )
 	{
 		return 0;
 	}
 
 	rtengine::IImage8* image = tpp->processImage (pparams, h, rtengine::TI_Bilinear, cfs.camera, scale );
    tpp->getDimensions(lastW,lastH,lastScale);
 
 	delete tpp;
 	tpp = 0;
 	return image;
}

void Thumbnail::generateExifDateTimeStrings () {

    exifString = "";
    dateTimeString = "";


    if( idata )
    exifString = Glib::ustring::compose ("f/%1 %2s %3%4 %5mm",
    		Glib::ustring(rtengine::ImageMetaData::apertureToString(idata->getFNumber())),
    		Glib::ustring(rtengine::ImageMetaData::shutterToString(idata->getShutterSpeed())),
    		M("QINFO_ISO"),
    		idata->getISOSpeed(),
    		idata->getFocalLen());

    if (options.fbShowExpComp && idata->getExpComp() ) // don't show exposure compensation if it is 0.00EV;
    	exifString = Glib::ustring::compose ("%1 %2EV", exifString, idata->expcompToString (idata->getExpComp(), false)); // append exposure compensation to exifString

    std::string dateFormat = options.dateFormat;
    std::ostringstream ostr;
    bool spec = false;
    for (size_t i=0; i<dateFormat.size(); i++)
        if (spec && dateFormat[i]=='y') {
            ostr << cfs.year;
            spec = false;
        }
        else if (spec && dateFormat[i]=='m') {
            ostr << (int)cfs.month;
            spec = false;
        }
        else if (spec && dateFormat[i]=='d') {
            ostr << (int)cfs.day;
            spec = false;
        }
        else if (dateFormat[i]=='%') 
            spec = true;
        else {
            ostr << (char)dateFormat[i];
            spec = false;
        }

    ostr << " " << (int)cfs.hour;
    ostr << ":" << std::setw(2) << std::setfill('0') << (int)cfs.min;
    ostr << ":" << std::setw(2) << std::setfill('0') << (int)cfs.sec;

    dateTimeString = ostr.str ();
}

const Glib::ustring& Thumbnail::getExifString () {

    return exifString;
}

const Glib::ustring& Thumbnail::getDateTimeString () {

    return dateTimeString;
}

ThFileType Thumbnail::getType () {

    return (ThFileType) cfs.format;
}

int Thumbnail::loadInfoFromImage ( ) {

	Glib::ustring xmpfname = removeExtension(fname) + paramFileExtension;
	Glib::ustring xmpfname2 = options.saveParamsCache? getCacheFileName ("profiles")+paramFileExtension :"";

	if (idata)
		delete idata;

	if( getExtension(fname).lowercase().compare("dng")==0 && options.embedXmpIntoDNG )
		idata = rtengine::ImageMetaData::fromFile (fname,"",xmpfname2,true);
	else if( getExtension(fname).lowercase().compare("png")==0 && options.embedXmpIntoPNG )
		idata = rtengine::ImageMetaData::fromFile (fname,"",xmpfname2,true);
	else if( getExtension(fname).lowercase().compare("jpg")==0 && options.embedXmpIntoJPG )
		idata = rtengine::ImageMetaData::fromFile (fname,"",xmpfname2,true);
	else if( (getExtension(fname).lowercase().compare("tif")==0 || getExtension(fname).lowercase().compare("tiff")==0) && options.embedXmpIntoTIFF )
		idata = rtengine::ImageMetaData::fromFile (fname,"",xmpfname2,true);
	else
       idata = rtengine::ImageMetaData::fromFile (fname,xmpfname,xmpfname2);
    if (!idata)
        return 0;

    int currentId = idata->getSnapshotId( );
    if( currentId >=0 ){
    	rtengine::SnapshotInfo si = idata->getSnapshot( currentId );
		pparams = si.params;
		pparamsId = currentId;
		pparamsValid=true;
    }

    int deg = 0;

    cfs.shutter = idata->getShutterSpeed();
	cfs.fnumber = idata->getFNumber();
	cfs.focalLen = idata->getFocalLen();
	cfs.iso = idata->getISOSpeed();
    cfs.expcomp  = idata->expcompToString (idata->getExpComp(), false); // do not mask Zero expcomp
	cfs.year = 1900 + idata->getDateTime().tm_year;
	cfs.month = idata->getDateTime().tm_mon + 1;
	cfs.day = idata->getDateTime().tm_mday;
	cfs.hour = idata->getDateTime().tm_hour;
	cfs.min = idata->getDateTime().tm_min;
	cfs.sec = idata->getDateTime().tm_sec;
	cfs.timeValid = true;
	cfs.exifValid = true;
	cfs.lens = idata->getLens();
	cfs.camera = idata->getCamera();

	if (idata->getOrientation()	== rtengine::ImageMetaData::ePhotoOrientationRotate90) {
		deg = 90;
	} else if (idata->getOrientation() == rtengine::ImageMetaData::ePhotoOrientationRotate180) {
		deg = 180;
	} else if (idata->getOrientation() == rtengine::ImageMetaData::ePhotoOrientationRotate270) {
		deg = 270;
	}

		// get image filetype
    std::string::size_type idx;
	idx = fname.rfind('.');
	if (idx != std::string::npos) {
		cfs.filetype = fname.substr(idx + 1);
	} else {
		cfs.filetype = "";
	}
	
    return deg;
}

void Thumbnail::_loadThumbnail(bool firstTrial) {

    needsReProcessing = true;
 	tw = -1;
 	th = options.maxThumbnailHeight;
    delete tpp;
    tpp = new rtengine::Thumbnail ();
    tpp->isRaw = (cfs.format == (int) FT_Raw);

    // load supplementary data
    bool succ = tpp->readData (getCacheFileName ("data")+".txt");
    
    // thumbnail image
    succ = succ && tpp->readImage (getCacheFileName ("images"));

    if (!succ && firstTrial) {
        _generateThumbnailImage ();
        if (cfs.supported && firstTrial)
            _loadThumbnail (false);

        if (tpp==NULL) return;
    }
    else if (!succ) {
        delete tpp;
        tpp = NULL;
		return;
    }
 
    if ( cfs.thumbImgType == CacheImageData::FULL_THUMBNAIL ) {
        // load aehistogram
        tpp->readAEHistogram (getCacheFileName ("aehistograms"));

        // load embedded profile
        tpp->readEmbProfile (getCacheFileName ("embprofiles")+".icc");

        tpp->init ();
        
    }
 
 	getThumbnailSize(tw,th);
}

void Thumbnail::loadThumbnail (bool firstTrial) {
	Glib::Mutex::Lock lock(mutex);
	_loadThumbnail(firstTrial);
}

void Thumbnail::_saveThumbnail () {

    if (!tpp)
        return;

    safe_g_remove (getCacheFileName ("images")+".cust");
    safe_g_remove (getCacheFileName ("images")+".cust16");
    safe_g_remove (getCacheFileName ("images")+".jpg");
    
    // save thumbnail image
    if (options.thumbnailFormat == FT_Custom) 
        tpp->writeImage (getCacheFileName ("images")+".cust", 1);
    else if (options.thumbnailFormat == FT_Custom16) 
        tpp->writeImage (getCacheFileName ("images")+".cust16", 2);
    else if (options.thumbnailFormat == FT_Jpeg) 
        tpp->writeImage (getCacheFileName ("images")+".jpg", 3);

    // save aehistogram
    tpp->writeAEHistogram (getCacheFileName ("aehistograms"));

    // save embedded profile
    tpp->writeEmbProfile (getCacheFileName ("embprofiles")+".icc");

    // save supplementary data
    tpp->writeData (getCacheFileName ("data")+".txt");
}

void Thumbnail::saveThumbnail () 
{
   	Glib::Mutex::Lock lock(mutex);
	_saveThumbnail(); 
}

void Thumbnail::updateCache (bool updatePParams, bool updateCacheImageData) {

    if (updatePParams && pparamsValid) {
        if( idata ){
			int idSnapshot=idata->getSnapshotId();
			if( idSnapshot<0 )
				idata->newSnapshot( rtengine::SnapshotInfo::kCurrentSnapshotName , pparams);
			else
				idata->updateSnapshot(idSnapshot, pparams );
        }
    }
    if (updateCacheImageData)
    cfs.save (getCacheFileName ("data")+".txt");
}

Thumbnail::~Thumbnail () {
	// TODO: Check for Linux
	#ifdef WIN32
	Glib::Mutex::Lock lock(mutex);
	#endif

    delete [] lastImg;
    delete tpp;
}

// Glib::ustring Thumbnail::getCacheFileName (Glib::ustring subdir) {
Glib::ustring Thumbnail::getCacheFileName (const Glib::ustring &subdir) const
{

    return cachemgr->getCacheFileName (subdir, fname, cfs.md5);
}

void Thumbnail::setFileName (const Glib::ustring fn) { 
    
    fname = fn; 
    cfs.md5 = cachemgr->getMD5 (fname); 
}

void Thumbnail::addThumbnailListener (ThumbnailListener* tnl) {

    increaseRef();
    listeners.push_back (tnl);
}

void Thumbnail::removeThumbnailListener (ThumbnailListener* tnl) {

    std::vector<ThumbnailListener*>::iterator f = std::find (listeners.begin(), listeners.end(), tnl);
    if (f!=listeners.end()) {
        listeners.erase (f);
        decreaseRef();
	}
}

// Calculates the standard filename for the automatically named batch result 
// and opens it in OS default viewer
// destination: 1=Batch conf. file; 2=batch out dir; 3=RAW dir
// Return: Success?
bool Thumbnail::openDefaultViewer(int destination) {

#ifdef WIN32 
    Glib::ustring openFName;
  
    if (destination==1) {
            openFName = Glib::ustring::compose ("%1.%2", BatchQueue::calcAutoFileNameBase(fname /* Merge uncertainty: was "this" */), options.saveFormatBatch.format);
            if (safe_file_test (openFName, Glib::FILE_TEST_EXISTS)) {
              wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (openFName.c_str(), -1, NULL, NULL, NULL);
              ShellExecuteW(NULL, L"open", wfilename, NULL, NULL, SW_SHOWMAXIMIZED );
              g_free(wfilename);
            } else {
                printf("File not found\n");
                return false;
            }
    } else {
        openFName = destination == 3 ? fname
            : Glib::ustring::compose ("%1.%2", BatchQueue::calcAutoFileNameBase(fname /* Merge uncertainty: was "this" */), options.saveFormatBatch.format);

        printf("Opening %s\n", openFName.c_str());

        if (safe_file_test (openFName, Glib::FILE_TEST_EXISTS)) {
            // Output file exists, so open explorer and select output file
            wchar_t* org=(wchar_t*)g_utf8_to_utf16 (Glib::ustring::compose("/select,\"%1\"", openFName).c_str(), -1, NULL, NULL, NULL);
            wchar_t* par=new wchar_t[wcslen(org)+1];
            wcscpy(par, org);

            // In this case the / disturbs
            wchar_t* p = par+1;  // skip the first backslash
            while (*p!=0) {
                if (*p==L'/') *p=L'\\';
                p++;
            }

            ShellExecuteW(NULL, L"open", L"explorer.exe", par, NULL, SW_SHOWNORMAL );

            delete[] par;
            g_free(org);
        } else if (safe_file_test (Glib::path_get_dirname(openFName), Glib::FILE_TEST_EXISTS)) {
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
	Glib::Mutex::Lock lock(mutex);
	bool previous = imageLoading;
	if( loading && !previous ){
		imageLoading = true;
		return true;
	}else if( !loading )
		imageLoading = false;
    return false;
}


/* Create a new snapshot with the given name and save inside it the procparams
 * Append the processing parameters to the end of the list then save the xmp file
 * Returns the id of the new snapshot
 */
int Thumbnail::newSnapshot(const Glib::ustring &newName, const rtengine::procparams::ProcParams& params, bool queued)
{
	int id =idata->newSnapshot( newName,params,queued);

    for (int i=0; i<listeners.size(); i++)
        listeners[i]->snapshotChanged(this,id);
    return id;
}

/* Delete the snapshot with the given id
 * then save the xmp file
 */
bool Thumbnail::deleteSnapshot( int id )
{
	return idata->deleteSnapshot(id );
}

/* Rename the snapshot with the given id
 * then save the xmp file
 */
bool Thumbnail::renameSnapshot(int id, const Glib::ustring &newname )
{
   bool b= idata->renameSnapshot(id, newname);
   for (int i=0; i<listeners.size(); i++)
       listeners[i]->snapshotChanged(this,id);
   return b;
}

/*
 *  Read all the snapshots from the xmpData
 */
rtengine::snapshotsList_t Thumbnail::getSnapshotsList()
{
    return idata->getSnapshotsList();
}

int Thumbnail::getNumQueued()
{
    return idata->getNumQueued();
}

rtengine::SnapshotInfo Thumbnail::getSnapshot( int id )
{
	return idata->getSnapshot( id );
}

bool  Thumbnail::setQueued( int id, bool inqueue )
{
    bool b= idata->setQueuedSnapshot( id, inqueue );
    for (int i=0; i<listeners.size(); i++)
        listeners[i]->snapshotChanged(this,id);
    return b;
}

int Thumbnail::getNumSaved()
{
    return idata->getSavedSnapshots();
}

bool Thumbnail::setSaved( int id, bool saved, const Glib::ustring &filename )
{
    bool b= idata->setSavedSnapshot( id, saved, filename );
    for (int i=0; i<listeners.size(); i++)
        listeners[i]->snapshotChanged(this,id);
    return b;
}

void Thumbnail::setRank  (int rank)
{
	if( idata ) idata->setRank(rank);
}
int Thumbnail::getRank  ()
{
    if( idata ) return idata->getRank();
    return 0;
}

void Thumbnail::setLabel( const Glib::ustring &colorlabel )
{
	if( idata ) idata->setLabel( colorlabel );
}

Glib::ustring Thumbnail::getLabel()
{
   if( idata ) return idata->getLabel();
   return "";
}
