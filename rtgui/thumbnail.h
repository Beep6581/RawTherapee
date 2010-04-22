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
#include <rtthumbnail.h>
#include <cacheimagedata.h>
#include <thumbnaillistener.h>

class CacheManager;
class Thumbnail {

        Glib::Mutex*    mutex;

        Glib::ustring   fname;              // file name corresponding to the thumbnail
        CacheImageData  cfs;                // cache entry corresponding to the thumbnail
        CacheManager*   cachemgr;           // parent
        int             ref;                // variable for reference counting
        int             enqueueNumber;      // the number of instances in the batch queue corresponding to this thumbnail
       
        // if the thumbnail is in processed mode, this class holds its data:
        rtengine::Thumbnail* tpp;
        int             tw, th;     // dimensions of timgdata (it stores tpp->width and tpp->height in processed mode for simplicity)
//        double          scale;      // portion of the sizes of the processed thumbnail image and the full scale image
        
        rtengine::procparams::ProcParams      pparams;
        bool            pparamsValid;
        bool            pparamsSet;
        bool            needsReProcessing;
        
        // these are the data of the result image of the last getthumbnailimage  call (for caching purposes)
        unsigned char*  lastImg;
        int             lastW;
        int             lastH;

        // exif & date/time strings
        Glib::ustring   exifString;
        Glib::ustring   dateTimeString;
        
        // vector of listeners
        std::vector<ThumbnailListener*> listeners;
        
        void            infoFromImage (const Glib::ustring& fname, rtengine::RawMetaDataLocation* rml=NULL);
        void            loadThumbnail (bool internal=false, bool firstTrial=true);
        void            saveThumbnail ();
        void            generateExifDateTimeStrings ();

        Glib::ustring    getCacheFileName (Glib::ustring subdir);
        
    public:
        Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf);
        Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5);
        ~Thumbnail ();
        
        bool              hasProcParams ();
        const rtengine::procparams::ProcParams& getProcParams ();
        void              setProcParams (const rtengine::procparams::ProcParams& pp, int whoChangedIt=-1, bool updateCacheNow=true);
        void              clearProcParams (int whoClearedIt=-1);
        void              loadProcParams ();

        bool              isRecentlySaved ();
        void              imageDeveloped ();
        void              imageEnqueued ();
        void              imageRemovedFromQueue ();
        bool              isEnqueued ();

//        unsigned char*  getThumbnailImage (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
        rtengine::IImage8* processThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale);
        void            processThumbImage2      (const rtengine::procparams::ProcParams& pparams, int h, rtengine::IImage8*& img, double& scale) { img = processThumbImage(pparams, h, scale); }
        void            getThumbnailSize        (int &w, int &h);
        void            getFinalSize            (const rtengine::procparams::ProcParams& pparams, int& w, int& h) { if (tpp) tpp->getFinalSize (pparams, w, h); }

        void            generateThumbnailImage (bool internal=false);

        const Glib::ustring&  getExifString ();
        const Glib::ustring&  getDateTimeString ();
        void                  getCamWB (double& temp, double& green) { if (tpp) tpp->getCamWB (temp, green); }
        void                  getAutoWB (double& temp, double& green) { if (tpp) tpp->getAutoWB (temp, green); }
        void                  getSpotWB (int x, int y, int rect, double& temp, double& green) { if (tpp) tpp->getSpotWB (getProcParams(), x, y, rect, temp, green); }
        void                  applyAutoExp (rtengine::procparams::ProcParams& pparams) { if (tpp) tpp->applyAutoExp (pparams); }
        
        ThFileType      getType ();
        Glib::ustring   getFileName () { return fname; }
        void            setFileName (const Glib::ustring fn);

        bool            isSupported ();

        const CacheImageData* getCacheImageData() { return &cfs; }
        std::string     getMD5   () { return cfs.md5; }

        int             getRank  () { return cfs.rank; }
        void            setRank  (int rank) { cfs.rank = rank;  }

        int             getStage () { return cfs.inTrash; }
        void            setStage (int stage) { cfs.inTrash = stage;  }

        void            addThumbnailListener (ThumbnailListener* tnl);
        void            removeThumbnailListener (ThumbnailListener* tnl);

        void            increaseRef ();
        void            decreaseRef ();
        
        void            updateCache ();
        void            reSaveThumbnail () { mutex->lock (); saveThumbnail(); mutex->unlock(); }
};


#endif

