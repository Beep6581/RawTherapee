/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include <multilangmgr.h>
#include <thumbnail.h>
#include <sstream>
#include <options.h>
#include <mytime.h>
#include <stdio.h>
#include <stdlib.h>
#include <glibmm.h>
#include <imagedata.h>
#include <glib/gstdio.h>
#include <guiutils.h>
#include <profilestore.h>

using namespace rtengine::procparams;

Thumbnail::Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf) 
    : cachemgr(cm), cfs(*cf), pparamsValid(false), fname(fname), 
        lastImg(NULL), ref(1), enqueueNumber(0), tpp(NULL), needsReProcessing(true) {

    mutex = new Glib::Mutex ();
    cfs.load (getCacheFileName ("data")+".txt");
    loadProcParams ();
    loadThumbnail ();
    generateExifDateTimeStrings ();
}

Thumbnail::Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5)
    : cachemgr(cm), pparamsValid(false), lastImg(NULL), fname(fname),
        ref(1), enqueueNumber(0), tpp(NULL), needsReProcessing(true) {

    mutex = new Glib::Mutex ();

    cfs.md5 = md5;
    generateThumbnailImage ();
    loadProcParams ();
    cfs.recentlySaved = false;
}

void Thumbnail::generateThumbnailImage (bool internal) {

    if (!internal)
        mutex->lock ();

//  delete everything loaded into memory
    delete tpp;
    tpp = NULL;
    delete [] lastImg;
    lastImg = NULL;
    tw = -1;
    th = options.maxThumbnailHeight;

// generate thumbnail image

    Glib::ustring ext = getExtension (fname);
    if (ext=="") 
        return;
    cfs.supported = false;
    cfs.exifValid = false;
    cfs.timeValid = false;
    
    delete tpp;
    tpp = NULL;
    if (ext.lowercase()=="jpg" || ext.lowercase()=="png" || ext.lowercase()=="tif" || ext.lowercase()=="tiff")
        tpp = rtengine::Thumbnail::loadFromImage (fname, tw, th, 1);
    if (tpp) {
        if (ext.lowercase()=="jpg") {
            cfs.format = FT_Jpeg;
            infoFromImage (fname);
        }
        else if (ext.lowercase()=="png")
            cfs.format = FT_Png;
        else if (ext.lowercase()=="tif" || ext.lowercase()=="tiff") {
            cfs.format = FT_Tiff;
            infoFromImage (fname);
        }
    }
    else {
        rtengine::RawMetaDataLocation ri;
        tpp = rtengine::Thumbnail::loadFromRaw (fname, ri, tw, th, 1);
        if (tpp) {
            cfs.format = FT_Raw;
            infoFromImage (fname, &ri);
        }
    }
    if (tpp) {
        // save thumbnail image to cache
        saveThumbnail ();
        cfs.supported = true;
    }
    needsReProcessing = true;

    cfs.save (getCacheFileName ("data")+".txt");

    generateExifDateTimeStrings ();
    
    if (!internal)
        mutex->unlock ();
}

bool Thumbnail::isSupported () {

    return cfs.supported;
}

const ProcParams& Thumbnail::getProcParams () {

    if (pparamsValid)
        return pparams;
    else {
        rtengine::procparams::ProcParams* pp = profileStore.getDefaultProcParams (getType()==FT_Raw);
        if (pp->wb.method=="Camera") {
            double ct;
            getCamWB (ct, pp->wb.green);
            pp->wb.temperature = ct;
        }
        else if (pp->wb.method=="Auto") {
            double ct;
            getAutoWB (ct, pp->wb.green);
            pp->wb.temperature = ct;
        }
        if (pp)
            return *pp;
    }
    return pparams; // there is no valid pp to return, but we have to return something
}

void Thumbnail::loadProcParams () {

    pparamsValid = false;
    if (options.paramsLoadLocation==PLL_Input) {
        // try to load it from pp2 file next to the image file    
        int ppres = pparams.load (fname + ".pp2");
        pparamsValid = !ppres && pparams.version>=220;
        if (!pparamsValid) 
                pparamsValid = !pparams.load (getCacheFileName ("profiles")+".pp2");
    }
    else {
        // try to load it from cache
        pparamsValid = !pparams.load (getCacheFileName ("profiles")+".pp2");
        // if no success, load it from pp2 file next to the image file    
        if (!pparamsValid) {
            int ppres = pparams.load (fname + ".pp2");
            pparamsValid = !ppres && pparams.version>=220;
        }
    }
}

void Thumbnail::clearProcParams (int whoClearedIt) {

    cfs.recentlySaved = false;
    pparamsValid = false;
    needsReProcessing = true;
    // remove pp2 file from cache
    Glib::ustring fname_ = getCacheFileName ("profiles")+".pp2";
    if (Glib::file_test (fname_, Glib::FILE_TEST_EXISTS))
        ::g_remove (fname_.c_str());
    // remove pp2 file located next to the file
//    fname_ = removeExtension(fname) + ".pp2";
    fname_ = fname + ".pp2";
    if (Glib::file_test (fname_, Glib::FILE_TEST_EXISTS))
        ::g_remove (fname_.c_str());
    fname_ = removeExtension(fname) + ".pp2";
    if (Glib::file_test (fname_, Glib::FILE_TEST_EXISTS))
        ::g_remove (fname_.c_str());

    for (int i=0; i<listeners.size(); i++)
        listeners[i]->procParamsChanged (this, whoClearedIt);
}

bool Thumbnail::hasProcParams () {
    
    return pparamsValid;
}

void Thumbnail::setProcParams (const ProcParams& pp, int whoChangedIt, bool updateCacheNow) {
    
    if (pparams!=pp) 
        cfs.recentlySaved = false;

    pparams = pp;
    pparamsValid = true;
    needsReProcessing = true;
    if (updateCacheNow)
        updateCache ();

    for (int i=0; i<listeners.size(); i++)
        listeners[i]->procParamsChanged (this, whoChangedIt);
}

bool Thumbnail::isRecentlySaved () {
    
    return cfs.recentlySaved;
}

void Thumbnail::imageDeveloped () {
        
    cfs.recentlySaved = true;
    cfs.save (getCacheFileName ("data")+".txt");
    pparams.save (getCacheFileName ("profiles")+".pp2");
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

void Thumbnail::increaseRef () { ref++; }
void Thumbnail::decreaseRef () { ref--; if (!ref) cachemgr->closeThumbnail (this); }

void Thumbnail::getThumbnailSize (int &w, int &h) {

    if (tpp) 
        w = tpp->getImageWidth (getProcParams(), h);
    else
        w = tw * h / th;
}

rtengine::IImage8* Thumbnail::processThumbImage (const rtengine::procparams::ProcParams& pparams, int h, double& scale) {

    mutex->lock ();

    if (!tpp)
        return NULL;

    rtengine::IImage8* res = tpp->processImage (pparams, h, rtengine::TI_Bilinear, scale);

    mutex->unlock ();
    return res;
}

void Thumbnail::generateExifDateTimeStrings () {

    exifString = "";
    dateTimeString = "";

    if (!cfs.exifValid)
        return;

    exifString = Glib::ustring::compose ("f/%1 %2s %3%4", Glib::ustring(rtengine::ImageData::apertureToString(cfs.fnumber)), Glib::ustring(rtengine::ImageData::shutterToString(cfs.shutter)), M("QINFO_ISO"), cfs.iso);

    std::string dateFormat = options.dateFormat;
    std::ostringstream ostr;
    bool spec = false;
    for (int i=0; i<dateFormat.size(); i++)
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
    ostr << " " << (int)cfs.hour << ":" << (int)cfs.min << ":" << (int)cfs.sec;

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

void Thumbnail::infoFromImage (const Glib::ustring& fname, rtengine::RawMetaDataLocation* rml) {

    rtengine::ImageMetaData* idata = rtengine::ImageMetaData::fromFile (fname, rml);
    if (!idata)
        return;
    cfs.timeValid = false;
    cfs.exifValid = false;
    if (idata->hasExif()) {
        cfs.shutter  = idata->getShutterSpeed ();
        cfs.fnumber  = idata->getFNumber ();
        cfs.focalLen = idata->getFocalLen ();
        cfs.iso      = idata->getISOSpeed ();
        cfs.year     = 1900 + idata->getDateTime().tm_year;
        cfs.month    = idata->getDateTime().tm_mon + 1;
        cfs.day      = idata->getDateTime().tm_mday;
        cfs.hour     = idata->getDateTime().tm_hour;
        cfs.min      = idata->getDateTime().tm_min;
        cfs.sec      = idata->getDateTime().tm_sec;
        cfs.timeValid = true;
        cfs.exifValid = true;
        cfs.lens     = idata->getLens();
        cfs.camera   = idata->getMake() + " " + idata->getModel();
    }
	else {
        cfs.lens     = "Unknown";
        cfs.camera   = "Unknown";
	}
    delete idata;
}

void Thumbnail::loadThumbnail (bool internal, bool firstTrial) {

    if (!internal)
        mutex->lock ();

    needsReProcessing = true;
    delete tpp;
    tpp = new rtengine::Thumbnail ();
    tpp->isRaw = (cfs.format == (int) FT_Raw);

    // load supplementary data
    bool succ = tpp->readData (getCacheFileName ("data")+".txt");
    
    // thumbnail image
    succ = succ && tpp->readImage (getCacheFileName ("images"));

    if (!succ && firstTrial) {
        generateThumbnailImage (true);
        if (cfs.supported && firstTrial)
            loadThumbnail (true, false);
    }
    else if (!succ) {
        delete tpp;
        tpp = NULL;
    }
    else {
        // load aehistogram
        tpp->readAEHistogram (getCacheFileName ("aehistograms"));

        // load embedded profile
        tpp->readEmbProfile (getCacheFileName ("embprofiles")+".icc");

        tpp->init ();
        
    }
    if (!internal)
        mutex->unlock ();
}

void Thumbnail::saveThumbnail () {

    if (!tpp)
        return;

    ::g_remove ((getCacheFileName ("images")+".cust").c_str());
    ::g_remove ((getCacheFileName ("images")+".cust16").c_str());
    ::g_remove ((getCacheFileName ("images")+".jpg").c_str());
    
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

void Thumbnail::updateCache () {

    if (pparamsValid) {
        if (options.saveParamsCache)
            pparams.save (getCacheFileName ("profiles")+".pp2");
        if (options.saveParamsFile)
//            pparams.save (removeExtension(fname) + ".pp2");
            pparams.save (fname + ".pp2");
    }
    cfs.save (getCacheFileName ("data")+".txt");
}

Thumbnail::~Thumbnail () {

    delete [] lastImg;
    delete tpp;
}

Glib::ustring Thumbnail::getCacheFileName (Glib::ustring subdir) {

    return cachemgr->getCacheFileName (subdir, fname, cfs.md5);
}

void Thumbnail::setFileName (const Glib::ustring fn) { 
    
    fname = fn; 
    cfs.md5 = cachemgr->getMD5 (fname); 
}

void Thumbnail::addThumbnailListener (ThumbnailListener* tnl) {

    listeners.push_back (tnl);
}

void Thumbnail::removeThumbnailListener (ThumbnailListener* tnl) {

    std::vector<ThumbnailListener*>::iterator f = std::find (listeners.begin(), listeners.end(), tnl);
    if (f!=listeners.end())
        listeners.erase (f);
}

