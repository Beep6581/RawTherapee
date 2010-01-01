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
#include <thumbnail.h>
#include <sstream>
#include <options.h>
#include <mytime.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <tiff.h>
#include <tiffio.h>
#include <glibmm.h>
#include <imagedata.h>
#include <glib/gstdio.h>
extern "C" {
#include <jpeglib.h>
extern jmp_buf jpeg_jmp_buf;
extern GLOBAL(struct jpeg_error_mgr *)
my_jpeg_std_error (struct jpeg_error_mgr * err);
extern GLOBAL(void)
my_jpeg_stdio_src (j_decompress_ptr cinfo, FILE * infile);
}
#include <guiutils.h>
#include <profilestore.h>

using namespace rtengine::procparams;
using namespace rtengine;

Thumbnail::Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheFileStruct* cf, int fpos) 
    : cachemgr(cm), cfs(cf), tImgData(NULL), tw(-1), th(-1), pparamsValid(false), filePos(fpos), fname(fname), 
        lastImg(NULL), ref(1), enqueueNumber(0), tpp(NULL), needsReProcessing(true) {

    mutex = new Glib::Mutex ();

    loadProcParams ();
    if (options.liveThumbnails!=cfs->thumbProcessed) {
        cfs->thumbProcessed = options.liveThumbnails;
        generateThumbnailImage (false);
    }
    if (!tImgData || !tpp)
        loadThumbnail (options.thumbCacheMemPolicy==MP_Memory);
}

Thumbnail::Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5, CacheFileStruct* cf, int fpos)
    : cachemgr(cm), cfs(cf), pparamsValid(false), filePos(fpos), lastImg(NULL), tImgData(NULL), fname(fname),
        ref(1), enqueueNumber(0), tpp(NULL), needsReProcessing(true) {

    mutex = new Glib::Mutex ();

    memcpy (cfs->MD5, md5.c_str(), 32);
    cfs->ppformat = 0;
    cfs->preinterp = 0;
    cfs->accessed = -1;
    cfs->rank = 0;
    cfs->stage = 0;
    cfs->thumbnail = FT_Invalid;
    cfs->recentlySaved = false;
    cfs->thumbProcessed = options.liveThumbnails;

    generateThumbnailImage (false);
}

void Thumbnail::generateThumbnailImage (bool useLock) {

    mutex->lock ();

    if (cfs->thumbProcessed)
        generateProcessedThumbnailImage (useLock);
    else
        generatePreviewThumbnailImage (useLock);

    mutex->unlock ();
}

void Thumbnail::generatePreviewThumbnailImage (bool useLock) {

//  delete everything loaded into memory
    delete tpp;
    tpp = NULL;
    delete [] lastImg;
    lastImg = NULL;
    tw = -1;
    th = options.maxThumbnailHeight;
    delete [] tImgData;
    tImgData = NULL;

// generate thumbnail image
    int lastdot = fname.find_last_of ('.');
    if (lastdot==Glib::ustring::npos) 
        return;

    cfs->thumbnail = (int) options.thumbnailFormat;
    Glib::ustring ext = fname.substr (lastdot);
    cfs->supported = false;
    cfs->exifValid = 0;
    cfs->timeValid = 0;
    if (ext.lowercase()==".jpg") {
        // try to load it as jpeg
        int thumbW, thumbH;
        int res = loadJPEGThumbnail (0, tImgData, tw, th);
        if (!res) {
            cfs->format = (int) FT_Jpeg;
            infoFromImage (fname);
            cfs->supported = true;
        }
        else {
            // try to let it load by dcraw
            rtengine::RawMetaDataLocation ri;
            res = getRawFileBasicInfo (fname, ri, cfs->rotate, thumbW, thumbH, cfs->thumbOffset, cfs->thumbImgType);
            if (!res) {
                cfs->format = (int) FT_Raw;
                infoFromImage (fname, &ri);
                cfs->supported = true;
                if (cfs->thumbImgType == 1)
                    loadJPEGThumbnail (cfs->thumbOffset, tImgData, tw, th, cfs->rotate);
                else if (cfs->thumbImgType == 2)
                    loadPPMThumbnail (cfs->thumbOffset, thumbW, thumbH, tImgData, tw, th, cfs->rotate);
                else {
                    cfs->thumbnail = (int) FT_None;
                    tw = th * thumbW / thumbH;
                }

            }
        }
    }
    else if (ext.lowercase()==".tif" || ext.lowercase()==".tiff") {
        // try to detect it as raw file first
        int thumbW, thumbH;
        rtengine::RawMetaDataLocation ri;
        if (!getRawFileBasicInfo (fname, ri, cfs->rotate, thumbW, thumbH, cfs->thumbOffset, cfs->thumbImgType)) {
            cfs->format = (int) FT_Raw;
            infoFromImage (fname, &ri);
            cfs->supported = true;
            if (cfs->thumbImgType == 1)
                loadJPEGThumbnail (cfs->thumbOffset, tImgData, tw, th, cfs->rotate);
            else if (cfs->thumbImgType == 2)
                loadPPMThumbnail (cfs->thumbOffset, thumbW, thumbH, tImgData, tw, th, cfs->rotate);           
            else {
                cfs->thumbnail = (int) FT_None;
                tw = th * thumbW / thumbH;
            }
        }
        else if (!loadTIFFThumbnail (tImgData, tw, th)) {
            cfs->format = (int) FT_Tiff;
            infoFromImage (fname);
            cfs->supported = true;
        }
    }
    else if (ext.lowercase()==".png") {
        int res = loadPNGThumbnail (tImgData, tw, th);
        if (!res) {
            cfs->format = (int) FT_Png;
            cfs->supported = true;
        }       
    }
    else {
        int thumbW, thumbH;
        rtengine::RawMetaDataLocation ri;
        if (!getRawFileBasicInfo (fname, ri, cfs->rotate, thumbW, thumbH, cfs->thumbOffset, cfs->thumbImgType)) {
            cfs->format = (int) FT_Raw;
            infoFromImage (fname, &ri);
            cfs->supported = true;
            if (cfs->thumbImgType == 1)
                loadJPEGThumbnail (cfs->thumbOffset, tImgData, tw, th, cfs->rotate);
            else if (cfs->thumbImgType == 2)
                loadPPMThumbnail (cfs->thumbOffset, thumbW, thumbH, tImgData, tw, th, cfs->rotate);
            else {
                cfs->thumbnail = (int) FT_None;
                tw = th * thumbW / thumbH;
            }
        }
    }
        
    cfs->thumbProcessed = 0;
// save thumbnail image to cache
    if (cfs->supported) {
        saveThumbnail ();
        if (options.thumbCacheMemPolicy==MP_Memory) {
            delete [] tImgData;
            tImgData = NULL;
        }
    }
    updateCache (true, true, useLock);
}

void Thumbnail::generateProcessedThumbnailImage (bool useLock) {

//  delete everything loaded into memory
    delete tpp;
    tpp = NULL;
    delete [] lastImg;
    lastImg = NULL;
    tw = -1;
    th = options.maxThumbnailHeight;
    delete [] tImgData;
    tImgData = NULL;

// generate thumbnail image
    int lastdot = fname.find_last_of ('.');
    if (lastdot==Glib::ustring::npos) 
        return;

    cfs->thumbnail = (int) options.thumbnailFormat;
    Glib::ustring ext = fname.substr (lastdot);
    cfs->supported = false;
    cfs->exifValid = 0;
    cfs->timeValid = 0;
    
    delete tpp;
    tpp = NULL;
    if (ext.lowercase()==".jpg" || ext.lowercase()==".tif" || ext.lowercase()==".tiff" || ext.lowercase()==".png")
        tpp = ThumbnailProcessingParameters::loadFromImage (fname, tw, th, 1);

    if (tpp) {
        if (ext.lowercase()==".jpg") {
            cfs->format = (int) FT_Jpeg;
            infoFromImage (fname);
        }
        else if (ext.lowercase()==".tif" || ext.lowercase()==".tiff") {
            cfs->format = (int) FT_Tiff;
            infoFromImage (fname);
        }
        else if (ext.lowercase()==".png")
            cfs->format = (int) FT_Png;
    }
    else {
        rtengine::RawMetaDataLocation ri;
        tpp = ThumbnailProcessingParameters::loadFromRaw (fname, ri, tw, th, 1);
        if (tpp) {
            cfs->format = (int) FT_Raw;
            infoFromImage (fname, &ri);
        }
    }
    if (tpp) {
        cfs->supported = true;
        cfs->camwbRed = tpp->camwbRed;
        cfs->camwbGreen = tpp->camwbGreen;
        cfs->camwbBlue = tpp->camwbBlue;
        cfs->autowbTemp = tpp->autowbTemp;
        cfs->autowbGreen = tpp->autowbGreen;
        cfs->aeHistCompression = tpp->aeHistCompression;
        cfs->embProfileLength = tpp->embProfileLength;
        cfs->redMultiplier = tpp->redMultiplier;
        cfs->greenMultiplier = tpp->greenMultiplier;
        cfs->blueMultiplier = tpp->blueMultiplier;
        cfs->scale = tpp->scale;
        cfs->defGain = tpp->defGain;
        cfs->scaleForSave = tpp->scaleForSave;
        cfs->gammaCorrected = tpp->gammaCorrected;
        memcpy (cfs->colorMatrix, tpp->colorMatrix, sizeof (cfs->colorMatrix));
/*        IImage8* res = tpp->processImage (getProcParams ());
        tw = res->getWidth ();
        th = res->getHeight ();
        tImgData = new unsigned char [tw*th*3];
        memcpy (tImgData, res->getData (), tw*th*3);
        res->free ();*/
    }
        
    cfs->thumbProcessed = 1;
// save thumbnail image to cache
    if (cfs->supported) {
        saveThumbnail ();
        if (options.thumbCacheMemPolicy==MP_Memory) {
            delete tpp;
            tpp = NULL;
//            delete [] tImgData;
//            tImgData = NULL;
        }
    }

    updateCache (true, true, useLock);
    needsReProcessing = true;
}

bool Thumbnail::isSupported () {

    return cfs->supported;
}

const ProcParams& Thumbnail::getProcParams () {

    if (pparamsValid)
        return pparams;
    else {
        rtengine::procparams::ProcParams* pp = profileStore.getDefaultProcParams (getType()==FT_Raw);
        if (pp)
            return *pp;
    }
    return pparams; // there is no valid pp to return, but we have to return something
}

void Thumbnail::loadProcParams () {

    pparamsValid = false;
    if (options.paramsLoadLocation==PLL_Input) {
        // try to load it from pp2 file next to the image file    
        int ppres = pparams.load (removeExtension(fname) + ".pp2");
        pparamsValid = !ppres && pparams.version>=220;
        // if no success, load it from the cache
        if (!pparamsValid)
            pparamsValid = !pparams.load (getCacheFileName ("profiles")+".pp2");
    }
    else {
        // try to load it from cache
        pparamsValid = !pparams.load (getCacheFileName ("profiles")+".pp2");
        // if no success, load it from pp2 file next to the image file    
        if (!pparamsValid) {
            int ppres = pparams.load (removeExtension(fname) + ".pp2");
            pparamsValid = !ppres && pparams.version>=220;
        }
    }
}

void Thumbnail::clearProcParams () {

    cfs->recentlySaved = false;
    cfs->ppformat = 0;
    pparamsValid = false;
    needsReProcessing = true;
    // remove pp2 file from cache
    Glib::ustring fname_ = getCacheFileName ("profiles")+".pp2";
    if (Glib::file_test (fname_, Glib::FILE_TEST_EXISTS))
        ::g_remove (fname_.c_str());
    // remove pp2 file located next to the file
    fname_ = removeExtension(fname) + ".pp2";
    if (Glib::file_test (fname_, Glib::FILE_TEST_EXISTS))
        ::g_remove (fname_.c_str());
}

bool Thumbnail::hasProcParams () {
    
    return pparamsValid;
}

void Thumbnail::setProcParams (const ProcParams& pp, bool updateCacheNow) {
    
    if (pparams!=pp) 
        cfs->recentlySaved = false;

    pparams = pp;
    pparamsValid = true;
    needsReProcessing = true;
    if (updateCacheNow)
        pparams.save (getCacheFileName ("profiles")+".pp2");
}

bool Thumbnail::isRecentlySaved () {
    
    return cfs->recentlySaved;
}

void Thumbnail::imageDeveloped () {
        
    cfs->recentlySaved = true;
    updateCache (false, false);
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

void Thumbnail::getThumbnailSize (int &w, int &h, int fixwh) { // fixwh = 0: fix w and calculate h, =1: fix h and calculate w

    mutex->lock ();

    if (fixwh==1) 
        w = tw * h / th;
    else if (fixwh==0) 
        h = th * w / tw;

    mutex->unlock ();
}

unsigned char* Thumbnail::getThumbnailImage (int &w, int &h, int fixwh) { // fixwh = 0: fix w and calculate h, =1: fix h and calculate w

    getThumbnailSize (w, h, fixwh);

    mutex->lock ();

    // if the result of last query fits the needs, we just have to return it
    if (lastImg && lastW==w && lastH==h && !(cfs->thumbProcessed && needsReProcessing)) {
        mutex->unlock ();
        return lastImg;
    }
    else 
        delete [] lastImg;
       
    if (!cfs->thumbProcessed) {
        // if thumbnail not loaed yet, load it
        if (!tImgData) {
            loadThumbnail ();
            if (!tImgData) {
                if (fixwh==1) 
                    w = tw * h / th;
                else if (fixwh==0) 
                    h = th * w / tw;
                mutex->unlock ();
                return NULL;
            }
        }
        lastW = w;
        lastH = h;
        lastImg = new unsigned char [lastW*lastH*3];
        thumbInterp (tImgData, tw, th, lastImg, w, h);
    }
    else {
        if (!tpp)
            loadThumbnail ();
        IImage8* res = tpp->processImage (getProcParams (), h, TI_Bilinear);
        lastW = w = res->getWidth ();
        lastH = h = res->getHeight ();           
        lastImg = new unsigned char [lastW*lastH*3];
        memcpy (lastImg, res->getData(), lastW*lastH*3);
        res->free ();
        needsReProcessing = false;
    }
        
        // if we have to spare with memory, delete full res thumb image from memory
    if (options.thumbCacheMemPolicy==MP_Memory) {
        delete [] tImgData;
        delete tpp;
        tpp = NULL;
        tImgData = NULL;
    }

    mutex->unlock ();
    return lastImg;
}

Glib::ustring Thumbnail::getExifString () {

    if (!cfs->exifValid)
        return "";

    std::ostringstream s;
    s << "f/" << ImageData::apertureToString(cfs->fnumber) << " ";
    s << ImageData::shutterToString(cfs->shutter) << "s ";
    s << "ISO" << cfs->iso;

    return s.str();
}

Glib::ustring Thumbnail::getDateTimeString () {

    if (!cfs->timeValid)
        return "";
    std::string dateFormat = options.dateFormat;
    std::ostringstream ostr;
    bool spec = false;
    for (int i=0; i<dateFormat.size(); i++)
    if (spec && dateFormat[i]=='y') {
        ostr << cfs->year;
        spec = false;
    }
    else if (spec && dateFormat[i]=='m') {
        ostr << (int)cfs->month;
        spec = false;
    }
    else if (spec && dateFormat[i]=='d') {
        ostr << (int)cfs->day;
        spec = false;
    }
    else if (dateFormat[i]=='%') 
        spec = true;
    else {
        ostr << (char)dateFormat[i];
        spec = false;
    }
    ostr << " " << (int)cfs->hour << ":" << (int)cfs->min << ":" << (int)cfs->sec;
    return ostr.str ();
}

ThFileType Thumbnail::getType () {

    return (ThFileType) cfs->format;
}

void Thumbnail::infoFromImage (const Glib::ustring& fname, rtengine::RawMetaDataLocation* rml) {

    ImageMetaData* idata = ImageMetaData::fromFile (fname, rml);
    if (!idata)
        return;
    cfs->timeValid = false;
    cfs->exifValid = false;
    if (idata->hasExif()) {
        cfs->shutter  = idata->getShutterSpeed ();
        cfs->fnumber  = idata->getFNumber ();
        cfs->focalLen = idata->getFocalLen ();
        cfs->iso      = idata->getISOSpeed ();
        cfs->year     = 1900 + idata->getDateTime().tm_year;
        cfs->month    = idata->getDateTime().tm_mon + 1;
        cfs->day      = idata->getDateTime().tm_mday;
        cfs->hour     = idata->getDateTime().tm_hour;
        cfs->min      = idata->getDateTime().tm_min;
        cfs->sec      = idata->getDateTime().tm_sec;
        cfs->timeValid = true;
        cfs->exifValid = true;
        cfs->lens     = idata->getLens();
        cfs->camera   = idata->getMake() + " " + idata->getModel();
    }
    delete idata;
}

void Thumbnail::loadThumbnail (bool sizeOnly) {

    mutex->lock ();

    delete tImgData;
    delete tpp;
    tImgData = NULL;
    tpp = NULL;
    needsReProcessing = true;
    
    if (cfs->thumbnail == FT_None)
        return;

    unsigned char* data = NULL;
    Glib::ustring fname = getCacheFileName ("images");
    int imgType = 0;
    if (Glib::file_test (fname+".cust", Glib::FILE_TEST_EXISTS))
        imgType = 1;
    else if (Glib::file_test (fname+".jpg", Glib::FILE_TEST_EXISTS))
        imgType = 2;


    FILE* f = g_fopen (fname.c_str(), "rb");
    if (!f || !imgType) {
        if (f)
            fclose(f);
        mutex->unlock ();
        generateThumbnailImage ();
        if (cfs->supported)
            loadThumbnail (sizeOnly);
        else
            return;
    }
    else if (imgType==1) {
        fread (&tw, 1, sizeof (int), f);
        fread (&th, 1, sizeof (int), f);
        if (!sizeOnly) {
            data = new unsigned char [tw*th*3];
            fread (data, tw*th, 3, f);
        }
        fclose (f);
    }
    else if (imgType==2) {
        struct jpeg_decompress_struct cinfo;
        struct jpeg_error_mgr jerr;
        if (!setjmp(jpeg_jmp_buf)) {
            cinfo.err = my_jpeg_std_error (&jerr);
            jpeg_create_decompress (&cinfo);
            my_jpeg_stdio_src (&cinfo,f);
            jpeg_read_header (&cinfo, TRUE);
            tw = cinfo.image_width;
            th = cinfo.image_height;
            if (!sizeOnly) {
                cinfo.dct_method = JDCT_FASTEST;
                cinfo.do_fancy_upsampling = 1;
                jpeg_start_decompress(&cinfo);
                data = new unsigned char [tw*th*3];
                while (cinfo.output_scanline < cinfo.output_height) {
                    unsigned char* row = data + cinfo.output_scanline*tw*3;
                    jpeg_read_scanlines (&cinfo, &row, 1);
                }
                jpeg_finish_decompress (&cinfo);
            }
            jpeg_destroy_decompress (&cinfo);
            fclose (f);
        }
        else {
            fclose (f);
            mutex->unlock ();
            return;
        }
    }
    if (!cfs->thumbProcessed) {
        tImgData = data;
    }
    else if (!sizeOnly) {
        tpp = new ThumbnailProcessingParameters ();
        tpp->data = data;
        tpp->width = tw;
        tpp->height = th;
        tpp->camwbRed = cfs->camwbRed;
        tpp->camwbGreen = cfs->camwbGreen;
        tpp->camwbBlue = cfs->camwbBlue;
        tpp->autowbTemp = cfs->autowbTemp;
        tpp->autowbGreen = cfs->autowbGreen;
        tpp->aeHistCompression = cfs->aeHistCompression;
        FILE* f = fopen (getCacheFileName ("aehistograms").c_str(), "rb");
        if (!f)
            tpp->aeHistogram = NULL;
        else {
            tpp->aeHistogram = new int[65536>>tpp->aeHistCompression];
            fread (tpp->aeHistogram, 1, (65536>>tpp->aeHistCompression)*sizeof(int), f);
            fclose (f);
        }
        tpp->embProfileLength = cfs->embProfileLength;
        Glib::ustring pfName = getCacheFileName ("embprofiles")+".icc";
        f = fopen (pfName.c_str(), "rb");
        if (!f) {
            tpp->embProfileData = NULL;
            tpp->embProfile = NULL;
        }
        else {
            tpp->embProfileData = new unsigned char[tpp->embProfileLength];
            fread (tpp->embProfileData, 1, tpp->embProfileLength, f);
            fclose (f);
            tpp->embProfile = cmsOpenProfileFromMem (tpp->embProfileData, tpp->embProfileLength);
        }
        tpp->redMultiplier = cfs->redMultiplier;
        tpp->greenMultiplier = cfs->greenMultiplier;
        tpp->blueMultiplier = cfs->blueMultiplier;
        tpp->scale = cfs->scale;
        tpp->defGain = cfs->defGain;
        tpp->scaleForSave = cfs->scaleForSave;
        tpp->gammaCorrected = cfs->gammaCorrected;
        tpp->isRaw = (cfs->format == (int) FT_Raw);
        memcpy (tpp->colorMatrix, cfs->colorMatrix, sizeof (tpp->colorMatrix));
        tpp->init ();
    }
    mutex->unlock ();
}

void Thumbnail::saveThumbnail () {

    if (!tImgData && !tpp)
        return;

    unsigned char* data = tImgData;
    int w = tw, h = th;
    
    if (cfs->thumbProcessed && tpp) {
        data = tpp->data;
        w = tpp->width;
        h = tpp->height;
    }

    if (data) {
        if (cfs->thumbnail == FT_Custom) {
            Glib::ustring fname = getCacheFileName ("images")+".cust";
            FILE* f = g_fopen (fname.c_str(), "wb");
            if (!f) 
                return;
            fwrite (&w, 1, sizeof (int), f);
            fwrite (&h, 1, sizeof (int), f);
            fwrite (data, w*h, 3, f);
            fclose (f);
        }
        else if (cfs->thumbnail == FT_Jpeg) {
            Glib::ustring fname = getCacheFileName ("images")+".jpg";
            FILE* f = g_fopen (fname.c_str(), "wb");
            if (!f) 
                return;
        	jpeg_compress_struct cinfo;
        	jpeg_error_mgr jerr;
	        cinfo.err = jpeg_std_error (&jerr);
	        jpeg_create_compress (&cinfo);
        	jpeg_stdio_dest (&cinfo, f);
            cinfo.image_width  = w;
	        cinfo.image_height = h;
	        cinfo.in_color_space = JCS_RGB;
	        cinfo.input_components = 3;
	        jpeg_set_defaults (&cinfo);
            cinfo.write_JFIF_header = FALSE;
            jpeg_set_quality (&cinfo, 85, true);
        	jpeg_start_compress(&cinfo, TRUE);
            int rowlen = w*3;
        	while (cinfo.next_scanline < cinfo.image_height) {
                unsigned char* row = data + cinfo.next_scanline*w*3;
        		if (jpeg_write_scanlines (&cinfo, &row, 1) < 1) {
                    jpeg_finish_compress (&cinfo);
	                jpeg_destroy_compress (&cinfo);
        	        fclose (f);
                    return;
                }
            }
        	jpeg_finish_compress (&cinfo);
        	jpeg_destroy_compress (&cinfo);
        	fclose (f);
        }
    }
    
    if (tpp) {
        // save aehistogram
        if (tpp->aeHistogram) {
            FILE* f = fopen (getCacheFileName ("aehistograms").c_str(), "wb");
            if (f) {
                fwrite (tpp->aeHistogram, 1, (65536>>tpp->aeHistCompression)*sizeof(int), f);
                fclose (f);
            }
        }
        // save embedded profile
        if (tpp->embProfileLength) {
            Glib::ustring pfName = getCacheFileName ("embprofiles")+".icc";
            FILE* f = fopen (pfName.c_str(), "wb");
            if (f) {
                fwrite (tpp->embProfileData, 1, tpp->embProfileLength, f);
                fclose (f);
            }
        }
    }
}

void Thumbnail::updateCache (bool savePparams, bool saveThumbImg, bool useLock) {

    if (filePos>=0) {

        if (pparamsValid) 
            cfs->ppformat = 1;
        else
            cfs->ppformat = 0;

        if (pparamsValid && savePparams) {
            if (options.saveParamsCache)
                pparams.save (getCacheFileName ("profiles"));
            if (options.saveParamsFile)
                pparams.save (removeExtension(fname) + ".pp2");
        }

        cachemgr->saveThumbnail (this, useLock);
    }
}

Thumbnail::~Thumbnail () {

    delete [] lastImg;
    delete [] tImgData;
}

Glib::ustring Thumbnail::getCacheFileName (Glib::ustring subdir) {

    Glib::ustring cfn = cachemgr->getBaseDir() + "/" + subdir + "/";
    for (int i=0; i<32; i++)
        cfn += cfs->MD5[i];
    return cfn;
}

int Thumbnail::loadJPEGThumbnail (int offset, unsigned char* &data, int &width, int &height, int degree, int fixwh) {

    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    data = NULL;

    FILE* f = g_fopen (fname.c_str(), "rb");

    if (!f) 
        return 2;
    fseek (f, offset, SEEK_SET);

    if (!setjmp(jpeg_jmp_buf)) {

        cinfo.err = my_jpeg_std_error (&jerr);
        jpeg_create_decompress (&cinfo);
        my_jpeg_stdio_src (&cinfo,f);
        jpeg_read_header (&cinfo, TRUE);
        cinfo.dct_method = JDCT_FASTEST;
        cinfo.scale_num = 1;

        if (fixwh==1) {
            if (degree==0 || degree==180) 
                width = height * cinfo.image_width / cinfo.image_height;
            else
                width = height * cinfo.image_height / cinfo.image_width;
        }
        else {
            if (degree==0 || degree==180) 
                height = width * cinfo.image_height / cinfo.image_width;
            else
                height = width * cinfo.image_width / cinfo.image_height;
        }

        if (cinfo.image_width > 8*width && cinfo.image_height > 8*height ) 
            cinfo.scale_denom = 8;
        else if (cinfo.image_width > 4*width && cinfo.image_height > 4*height ) 
            cinfo.scale_denom = 4;
        else if (cinfo.image_width > 2*width && cinfo.image_height > 2*height ) 
            cinfo.scale_denom = 2;
        else 
            cinfo.scale_denom = 1;

        cinfo.do_fancy_upsampling = 1;
        jpeg_start_decompress(&cinfo);
        unsigned char* fullData = new unsigned char [3*cinfo.output_height*cinfo.output_width];
        unsigned char* row = fullData;
        for (int i=0; i<cinfo.output_height; i++, row += 3*cinfo.output_width)
            jpeg_read_scanlines (&cinfo, &row, 1);

        unsigned char* tdata = new unsigned char [3*width*height];       
        if (degree==0)
            thumbInterp (fullData, cinfo.output_width, cinfo.output_height, tdata, width, height);
        else if (degree==90 || degree==270) {
            thumbInterp (fullData, cinfo.output_width, cinfo.output_height, tdata, height, width);
            rotate (tdata, height, width, degree);
        }
        else if (degree==180) {
            thumbInterp (fullData, cinfo.output_width, cinfo.output_height, tdata, width, height);
            rotate (tdata, width, height, degree);
        }
        delete [] fullData;
           
        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);
        fclose(f);

        data = tdata;
        return 0;
    }
    else {
        fclose(f);
        return 1;
    }

}

int Thumbnail::loadPPMThumbnail (int offset, int iw, int ih, unsigned char* &data, int &width, int &height, int degree, int fixwh) {

    data = NULL;

    FILE* f = g_fopen (fname.c_str(), "rb");
    if (!f) 
        return 2;

    fseek (f, offset, SEEK_SET);

    if (fixwh==1) {
        if (degree==0 || degree==180) 
            width = height * iw / ih;
        else
            width = height * ih / iw;
    }
    else {
        if (degree==0 || degree==180) 
            height = width * ih / iw;
        else
            height = width * iw / ih;
    }

    unsigned char* fullData = new unsigned char [3*iw*ih];
    fread (fullData, 1, iw*ih*3, f);
    unsigned char* tdata = new unsigned char [3*width*height];       
    if (degree==0)
        thumbInterp (fullData, iw, ih, tdata, width, height);
    else if (degree==90 || degree==270) {
        thumbInterp (fullData, iw, ih, tdata, height, width);
        rotate (tdata, height, width, degree);
    }
    else if (degree==180) {
        thumbInterp (fullData, iw, iw, tdata, width, height);
        rotate (tdata, width, height, degree);
    }
    delete [] fullData;
    fclose(f);
    data = tdata;
    return 0;
}


void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length) {
   /* fread() returns 0 on error, so it is OK to store this in a png_size_t
    * instead of an int, which is what fread() actually returns.
    */
   png_size_t check = (png_size_t)fread(data, (png_size_t)1, length, (FILE *)png_ptr->io_ptr);
   if (check != length)
      png_error(png_ptr, "Read Error");
}

int Thumbnail::loadPNGThumbnail (unsigned char* &data, int &width, int &height, int fixwh) {


    data = NULL;

    FILE *file=g_fopen(fname.c_str(),"rb");
    if (!file) return 2;

	//reading PNG header
	unsigned char header[8];
	fread(header,1,8,file);
	if (!png_check_sig(header,8)) {
		fclose(file);
		return 3;
	}
	//initializing main structures
	png_structp png= png_create_read_struct(PNG_LIBPNG_VER_STRING,0,0,0);
	if (!png) {
		fclose(file);
		return 4;
	}
	png_infop info= png_create_info_struct(png);
	png_infop end_info = png_create_info_struct(png);
	if (!end_info || !info) {
		png_destroy_read_struct(&png,&info,&end_info);
		fclose(file);
		return 5;
	}

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_read_struct(&png,&info,&end_info);
		fclose(file);
		return 3;
        }

	//set up png read
    png_set_read_fn(png, file, png_read_data);
	png_set_sig_bytes(png,8);

	png_read_info(png,info);

	//retrieving image information
    unsigned long fwidth, fheight;
	int bit_depth,color_type,interlace_type,compression_type,filter_method;
	png_get_IHDR(png,info,&fwidth,&fheight,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

	//converting to 32bpp format
	if (color_type==PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
	
	if (color_type==PNG_COLOR_TYPE_GRAY || color_type==PNG_COLOR_TYPE_GRAY_ALPHA)
          png_set_gray_to_rgb(png);

	if (png_get_valid(png,info,PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);
	
	if (interlace_type!=PNG_INTERLACE_NONE) {
          fclose (file);
          return 6;
        }

    if (color_type & PNG_COLOR_MASK_ALPHA)
        png_set_strip_alpha(png);

    if (bit_depth==16)
        png_set_strip_16(png);

	//setting gamma
	double gamma;
	if (png_get_gAMA(png,info,&gamma))
		png_set_gamma(png, 2.0, gamma);
	else
		png_set_gamma(png, 1.0, 1.0);

	//updating png info struct
	png_read_update_info(png,info);
	png_get_IHDR(png,info,&fwidth,&fheight,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

    if (color_type & PNG_COLOR_MASK_ALPHA)
        png_set_strip_alpha(png);
    if (bit_depth==16)
        png_set_strip_16(png);

	png_read_update_info(png,info);
	png_get_IHDR(png,info,&fwidth,&fheight,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

    png_timep time;
    cfs->timeValid = png_get_tIME (png, info, &time);
    cfs->exifValid = 0;
    cfs->year = time->year;
    cfs->month = time->month;
    cfs->day = time->day;
    cfs->hour = time->hour;
    cfs->min = time->minute;
    cfs->sec = time->second; 

    if (fixwh==1) 
         width = height * fwidth / fheight;
    else 
         height = width * fheight / fwidth;

    unsigned char* fullData = new unsigned char [3*fheight*fwidth];
    unsigned char* row = fullData;
    for (int i=0; i<fheight; i++, row += 3*fwidth)
        png_read_row (png, (png_byte*)row, NULL);

    unsigned char* tdata = new unsigned char [3*width*height];       
    thumbInterp (fullData, fwidth, fheight, tdata, width, height);
    delete [] fullData;

	png_read_end(png,0);
	png_destroy_read_struct(&png,&info,&end_info);
	
	fclose(file);
	data = tdata;
    return 0;
}

int Thumbnail::loadTIFFThumbnail (unsigned char* &data, int &width, int &height, int fixwh) {

    data = NULL;
#ifdef WIN32
    wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
    TIFF* in = TIFFOpenW (wfilename, "r");
    g_free (wfilename);
#else
    TIFF* in = TIFFOpen(fname.c_str(), "r");
#endif
	if (in == NULL)
		return 1;
    int fwidth, fheight;
	TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &fwidth);
	TIFFGetField(in, TIFFTAG_IMAGELENGTH, &fheight);
    uint16 bitspersample, samplesperpixel;
	TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
    uint16 photometric;
	if (!TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photometric) || photometric != PHOTOMETRIC_RGB || samplesperpixel < 3) {
        TIFFClose(in);
		return 2;
	}
    uint16 config;
	TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);
	if (config != PLANARCONFIG_CONTIG) {
        TIFFClose(in);
		return 3;
	}

    if (fixwh==1) 
        width = height * fwidth / fheight;
    else 
        height = width * fheight / fwidth;

    unsigned char* fullData = new unsigned char [3*fheight*fwidth];
    unsigned char* row = fullData;
    if (bitspersample == 16) {
        unsigned short* linebuffer = new unsigned short[TIFFScanlineSize(in)];
        int ix = 0;
        for (int i=0; i<fheight; i++) {
            TIFFReadScanline(in, linebuffer, i, 0);
            for (int j=0; j<3*fwidth; j++)
                fullData[ix++] = linebuffer[j] >> 8;
        }
        TIFFClose(in);
        delete [] linebuffer;
    }
    else if (bitspersample == 8) {
        for (int i=0; i<fheight; i++, row += 3*fwidth) 
            TIFFReadScanline(in, row, i, 0);
    }
    else {
      	TIFFClose(in);
        return 4;
    }
    unsigned char* tdata = new unsigned char [3*width*height];       
    thumbInterp (fullData, fwidth, fheight, tdata, width, height);
    delete [] fullData;
    data = tdata;
    return 0;
}

