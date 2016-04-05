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
#include "cachemanager.h"
#include "options.h"
#include "../rtengine/rtengine.h"
#include "../rtengine/rtthumbnail.h"
#include "cacheimagedata.h"
#include "thumbnaillistener.h"
#include "threadutils.h"

class CacheManager;
class Thumbnail
{

    MyMutex         mutex;

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

    rtengine::procparams::ProcParams      pparams;
    bool            pparamsValid;
    bool            pparamsSet;
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

    void            _loadThumbnail ();
    void            _saveThumbnail ();
    void            _generateThumbnailImage ();
    int             infoFromImage (const Glib::ustring& fname, rtengine::RawMetaDataLocation* rml = NULL);
    void            loadThumbnail ();
    void            generateExifDateTimeStrings ();

    Glib::ustring    getCacheFileName (const Glib::ustring& subdir, const Glib::ustring& fext) const;

public:
    Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf);
    Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5);
    ~Thumbnail ();

    bool              hasProcParams ();
    const rtengine::procparams::ProcParams& getProcParams ();
    const rtengine::procparams::ProcParams& getProcParamsU ();  // Unprotected version

    // Use this to create params on demand for update ; if flaggingMode=true, the procparams is created for a file being flagged (inTrash, rank, colorLabel)
    rtengine::procparams::ProcParams* createProcParamsForUpdate (bool returnParams, bool forceCPB, bool flaggingMode = false);

    void              setProcParams (const rtengine::procparams::ProcParams& pp, ParamsEdited* pe = NULL, int whoChangedIt = -1, bool updateCacheNow = true);
    void              clearProcParams (int whoClearedIt = -1);
    void              loadProcParams ();

    void              notifylisterners_procParamsChanged(int whoChangedIt);

    bool              isQuick()
    {
        return cfs.thumbImgType == CacheImageData::QUICK_THUMBNAIL;
    }
    bool              isPParamsValid()
    {
        return pparamsValid;
    }
    bool              isRecentlySaved ();
    void              imageDeveloped ();
    void              imageEnqueued ();
    void              imageRemovedFromQueue ();
    bool              isEnqueued ();

//        unsigned char*  getThumbnailImage (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
    rtengine::IImage8* processThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale);
    rtengine::IImage8* upgradeThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale);
    void            getThumbnailSize        (int &w, int &h, const rtengine::procparams::ProcParams *pparams = NULL);
    void            getFinalSize            (const rtengine::procparams::ProcParams& pparams, int& w, int& h);

    const Glib::ustring&  getExifString ();
    const Glib::ustring&  getDateTimeString ();
    void                  getCamWB  (double& temp, double& green)
    {
        if (tpp) {
            tpp->getCamWB  (temp, green);
        } else {
            temp = green = -1.0;
        }
    }
    void                  getAutoWB (double& temp, double& green, double equal);
    void                  getSpotWB (int x, int y, int rect, double& temp, double& green)
    {
        if (tpp) {
            tpp->getSpotWB (getProcParams(), x, y, rect, temp, green);
        } else {
            temp = green = -1.0;
        }
    }
    void                  applyAutoExp (rtengine::procparams::ProcParams& pparams)
    {
        if (tpp) {
            tpp->applyAutoExp (pparams);
        }
    }

    ThFileType      getType ();
    Glib::ustring   getFileName ()
    {
        return fname;
    }
    void            setFileName (const Glib::ustring fn);

    bool            isSupported () const;

    const CacheImageData* getCacheImageData()
    {
        return &cfs;
    }
    std::string     getMD5   ()
    {
        return cfs.md5;
    }

    int             getRank  ()
    {
        return pparams.rank;
    }
    void            setRank  (int rank)
    {
        if (pparams.rank != rank) {
            pparams.rank = rank;
            pparamsValid = true;
        }
    }

    int             getColorLabel  ()
    {
        return pparams.colorlabel;
    }
    void            setColorLabel  (int colorlabel)
    {
        if (pparams.colorlabel != colorlabel) {
            pparams.colorlabel = colorlabel;
            pparamsValid = true;
        }
    }

    int             getStage ()
    {
        return pparams.inTrash;
    }
    void            setStage (int stage)
    {
        if (pparams.inTrash != stage) {
            pparams.inTrash = stage;
            pparamsValid = true;
        }
    }

    void            addThumbnailListener (ThumbnailListener* tnl);
    void            removeThumbnailListener (ThumbnailListener* tnl);

    void            increaseRef ();
    void            decreaseRef ();

    void            updateCache (bool updatePParams = true, bool updateCacheImageData = true);
    void            saveThumbnail ();

    bool            openDefaultViewer(int destination);
    bool            imageLoad(bool loading);
};


#endif

