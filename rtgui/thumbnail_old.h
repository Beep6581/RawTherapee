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
#ifndef _THUMBNAIL_
#define _THUMBNAIL_

#include <string>
#include <glibmm.h>
#include <cachemanager.h>
#include <options.h>
#include <rtengine.h>
#include <cachefilestruct.h>

class CacheManager;
class Thumbnail {

        Glib::Mutex*    mutex;

        Glib::ustring   fname;              // file name corresponding to the thumbnail
        CacheImageData* cfs;                // cache entry corresponding to the thumbnail
        CacheManager*   cachemgr;           // parent
        int             ref;                // variable for reference counting
        int             enqueueNumber;      // the number of instances in the batch queue corresponding to this thumbnail

        // if the thumbnail is in non-processed mode, these fields hold the thumbnail image:
        unsigned char*  tImgData;
        int             tw, th;     // dimensions of timgdata (it stores tpp->width and tpp->height in processed mode for simplicity)
        
        // if the thumbnail is in processed mode, this class holds its data:
        rtengine::ThumbnailProcessingParameters* tpp;
        
        rtengine::procparams::ProcParams      pparams;
        bool            pparamsValid;
        bool            pparamsSet;
        bool            needsReProcessing;
        
        // these are the data of the result image of the last getthumbnailimage  call (for caching purposes)
        unsigned char*  lastImg;
        int             lastW;
        int             lastH;
        
        void            infoFromImage (const Glib::ustring& fname, rtengine::RawMetaDataLocation* rml=NULL);
        void            loadThumbnail (bool sizeOnly=false);
        void            saveThumbnail ();
        void            generateProcessedThumbnailImage (bool useLock=true);
        void            generatePreviewThumbnailImage (bool useLock=true);


        Glib::ustring    getCacheFileName (Glib::ustring subdir);
        
    public:
        Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf);
        Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5);
        ~Thumbnail ();
        
        bool              hasProcParams ();
        const rtengine::procparams::ProcParams& getProcParams ();
        void              setProcParams (const rtengine::procparams::ProcParams& pp, bool updateCacheNow=true);
        void              clearProcParams ();
        void              loadProcParams ();

        bool              isRecentlySaved ();
        void              imageDeveloped ();
        void              imageEnqueued ();
        void              imageRemovedFromQueue ();
        bool              isEnqueued ();

        unsigned char*  getThumbnailImage (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
        void            getThumbnailSize  (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
        void            generateThumbnailImage (bool useLock=true);

        Glib::ustring   getExifString ();
        Glib::ustring   getDateTimeString ();
        
        ThFileType      getType ();
        Glib::ustring   getFileName () { return fname; }
        void            setFileName (const Glib::ustring fn) { fname = fn; }

        bool            isSupported ();

        void            saveToCache (FILE* f, bool full);

        void            setFilePos (int fp) { filePos = fp; }
        int             getFilePos ()       { return filePos; }
        
        int             getRank  () { return cfs->rank; }
        void            setRank  (int rank) { cfs->rank = rank;  }
        int             getStage () { return cfs->stage; }
        void            setStage (int stage) { cfs->stage = stage;  }

        void            increaseRef ();
        void            decreaseRef ();
        
        void            updateCache (bool pparams=true, bool thumbImg=true, bool useLock=true);

        const CacheFileStruct* getCacheFileStruct () const { return cfs; }
        
        int  loadJPEGThumbnail   (int offset, unsigned char* &data, int &width, int &height, int degree=0, int fixwh=1);
        int  loadPPMThumbnail    (int offset, int iw, int ih, unsigned char* &data, int &width, int &height, int degree=0, int fixwh=1);
        int  loadPNGThumbnail    (unsigned char* &data, int &width, int &height, int fixwh=1);
        int  loadTIFFThumbnail   (unsigned char* &data, int &width, int &height, int fixwh=1);
};


#endif

