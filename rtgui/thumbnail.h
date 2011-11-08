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
#include <rtXmp.h>


class SnapshotListener {
	public:
		virtual int  newSnapshot(const Glib::ustring &name, const rtengine::procparams::ProcParams& params, bool queued=false ){};
		virtual bool deleteSnapshot( int id ){};
		virtual bool renameSnapshot(int id, const Glib::ustring &newname ){};
};

class CacheManager;
class Thumbnail :public SnapshotListener{

        Glib::Mutex    mutex;

        Glib::ustring   fname;              // file name corresponding to the thumbnail
        CacheImageData  cfs;                // cache entry corresponding to the thumbnail
        CacheManager*   cachemgr;           // parent
        int             ref;                // variable for reference counting
        int             enqueueNumber;      // the number of instances in the batch queue corresponding to this thumbnail
       
        // if the thumbnail is in processed mode, this class holds its data:
        rtengine::Thumbnail* tpp;
        int             tw, th;             // dimensions of timgdata (it stores tpp->width and tpp->height in processed mode for simplicity)
        float           imgRatio;           // hack to avoid rounding error
//        double          scale;            // portion of the sizes of the processed thumbnail image and the full scale image

        rtengine::ImageMetaData  *idata;    // Metadata

        rtengine::procparams::ProcParams      pparams; // Current parameters for developing the shot
        int             pparamsId;	        // For identifing current pparams among the snapshots
        bool            pparamsValid;
        bool            pparamsSet;
        bool            needsReProcessing;
        bool            imageLoading;
        
        // these are the data of the result image of the last getthumbnailimage  call (for caching purposes)
        unsigned char*  lastImg;
        int             lastW;
        int             lastH;
        double          lastScale;

        // exif & date/time strings
        Glib::ustring   exifString;
        Glib::ustring   dateTimeString;
        
		bool            initial_;

        // vector of listeners
        std::vector<ThumbnailListener*> listeners;
        
        void            _loadThumbnail (bool firstTrial=true);
        void            _saveThumbnail ();
        void            _generateThumbnailImage ();

        void            loadThumbnail (bool firstTrial=true);
        void            generateExifDateTimeStrings ();
        int             loadInfoFromImage ( );

        int              initRTXMP();
    public:
        Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf);
        Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5);
        ~Thumbnail ();
        
        int loadXMP( );
        int saveXMP( ) const;

        bool              hasProcParams ();
        const rtengine::procparams::ProcParams& getProcParams ();
        Glib::ustring    getCacheFileName (const Glib::ustring &subdir) const;

        // Use this to create params on demand for update
        rtengine::procparams::ProcParams* createProcParamsForUpdate (bool returnParams, bool forceCPB);

        void              setProcParams (const rtengine::procparams::ProcParams& pp, int whoChangedIt=-1, bool updateCacheNow=true);
        void              clearProcParams (int whoClearedIt=-1);

        void              notifylisterners_procParamsChanged(int whoChangedIt);

		bool              isQuick() { return cfs.thumbImgType == CacheImageData::QUICK_THUMBNAIL; }
		bool              isPParamsValid() { return pparamsValid; }
        bool              isRecentlySaved ();
        void              imageDeveloped ();
        bool              isEnqueued ();

//        unsigned char*  getThumbnailImage (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
        rtengine::IImage8* processThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale);
        rtengine::IImage8* upgradeThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale);
        void            getThumbnailSize        (int &w, int &h);
        void            getFinalSize            (const rtengine::procparams::ProcParams& pparams, int& w, int& h);

        rtengine::ImageMetaData  *getMetadata() { return idata; }
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

        int             getRank  ();
        void            setRank  (int rank);

        Glib::ustring   getLabel  ();
        void            setLabel  (const Glib::ustring &colorlabel);

//        int             getStage ();              // stage => Rank=-1
//        void            setStage (int stage);

        void            addThumbnailListener (ThumbnailListener* tnl);
        void            removeThumbnailListener (ThumbnailListener* tnl);

        void            increaseRef ();
        void            decreaseRef ();
        
        void            updateCache (bool updatePParams = true, bool updateCacheImageData = true);
        void            saveThumbnail ();

        bool            openDefaultViewer(int destination);
        bool            imageLoad(bool loading);
		int 			newSnapshot(const Glib::ustring &name, const rtengine::procparams::ProcParams& params, bool queued=false );
		bool 			deleteSnapshot( int id );
		bool 			renameSnapshot(int id, const Glib::ustring &newname );
		rtengine::snapshotsList_t getSnapshotsList();
		rtengine::SnapshotInfo    getSnapshot( int id );
		int				getNumQueued();
		bool			setQueued( int id, bool inqueue=true );
		bool			setSaved( int id, bool saved=true, const Glib::ustring &filename="" );
};


#endif

