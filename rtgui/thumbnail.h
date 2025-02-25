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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <memory>
#include <string>

#include <glibmm/ustring.h>
#include <glibmm/datetime.h>

#include "cacheimagedata.h"
#include "threadutils.h"
#include "thumbnaillistener.h"

namespace rtengine
{
class Thumbnail;

namespace procparams
{

class ProcParams;

}

}
class CacheManager;

struct ParamsEdited;

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
//  double          scale;              // portion of the sizes of the processed thumbnail image and the full scale image

    const std::unique_ptr<rtengine::procparams::ProcParams>      pparams;
    bool            pparamsValid;
    bool            imageLoading;

    // these are the data of the result image of the last getthumbnailimage  call (for caching purposes)
    unsigned char*  lastImg;            // pointer to the processed and scaled base ImageIO image
    int             lastW;              // non rotated width of the cached ImageIO image
    int             lastH;              // non rotated height of the cached ImageIO image
    double          lastScale;          // scale of the cached ImageIO image

    // exif & date/time strings
    Glib::ustring   exifString;
    Glib::ustring   dateTimeString;
    Glib::DateTime  dateTime;

    bool            initial_;

    // Properties holds values and edited states for rank, color and trashed
    struct Properties {
        template <class T> struct Property {
            T value;
            bool edited;
            Property(T v): value(v), edited(false) {}
            Property& operator=(T v)
            {
                value = v;
                edited = true;
                return *this;
            }
            operator T() const { return value; }
        };
        Property<int> rank;
        Property<int> color;
        Property<bool> trashed;

        explicit Properties(int r=0, int c=0, bool t=false):
            rank(r), color(c), trashed(t) {}
        bool edited() const { return rank.edited || color.edited || trashed.edited; }
    };
    Properties properties;

    // vector of listeners
    std::vector<ThumbnailListener*> listeners;

    void            _loadThumbnail (bool firstTrial = true);
    void            _saveThumbnail ();
    void            _generateThumbnailImage ();
    int             infoFromImage (const Glib::ustring& fname);
    void            generateExifDateTimeStrings ();

    Glib::ustring    getCacheFileName (const Glib::ustring& subdir, const Glib::ustring& fext) const;

    void saveMetadata();
    void loadProperties();
    void updateProcParamsProperties(bool forceUpdate = false);
    void saveXMPSidecarProperties();

public:
    Thumbnail (CacheManager* cm, const Glib::ustring& fname, CacheImageData* cf);
    Thumbnail (CacheManager* cm, const Glib::ustring& fname, const std::string& md5, const std::string &xmpSidecarMd5);
    ~Thumbnail ();

    static int infoFromImage(const Glib::ustring &fname, CacheImageData &cfs);
    static Glib::ustring xmpSidecarPath(const Glib::ustring &imagePath);

    bool              hasProcParams () const;
    const rtengine::procparams::ProcParams& getProcParams ();
    const rtengine::procparams::ProcParams& getProcParamsU ();  // Unprotected version

    // Use this to create params on demand for update ; if flaggingMode=true, the procparams is created for a file being flagged (inTrash, rank, colorLabel)
    rtengine::procparams::ProcParams* createProcParamsForUpdate (bool returnParams, bool force, bool flaggingMode = false);

    void              setProcParams (const rtengine::procparams::ProcParams& pp, ParamsEdited* pe = nullptr, int whoChangedIt = -1, bool updateCacheNow = true, bool resetToDefault = false);
    void              clearProcParams (int whoClearedIt = -1);
    void              loadProcParams ();

    void              notifylisterners_procParamsChanged(int whoChangedIt);

    bool              isQuick() const;
    bool              isPParamsValid() const;
    bool              isRecentlySaved () const;
    void              imageDeveloped ();
    void              imageEnqueued ();
    void              imageRemovedFromQueue ();
    bool              isEnqueued () const;
    bool              isPixelShift () const;
    bool              isHDR () const;

//        unsigned char*  getThumbnailImage (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
    rtengine::IImage8* processThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale);
    rtengine::IImage8* upgradeThumbImage    (const rtengine::procparams::ProcParams& pparams, int h, double& scale, bool forceUpgrade);
    void            getThumbnailSize        (int &w, int &h, const rtengine::procparams::ProcParams *pparams = nullptr);
    void            getFinalSize            (const rtengine::procparams::ProcParams& pparams, int& w, int& h);
    void            getOriginalSize         (int& w, int& h) const;

    const Glib::ustring&  getExifString () const;
    const Glib::ustring&  getDateTimeString () const;
    const Glib::DateTime& getDateTime () const;
    void                  getCamWB  (double& temp, double& green, rtengine::StandardObserver observer) const;
    void                  getAutoWB (double& temp, double& green, double equal, rtengine::StandardObserver observer, double tempBias);
    void                  getSpotWB (int x, int y, int rect, double& temp, double& green);
    void                  applyAutoExp (rtengine::procparams::ProcParams& pparams);

    ThFileType      getType () const;
    Glib::ustring   getFileName () const
    {
        return fname;
    }
    void            setFileName (const Glib::ustring &fn);

    bool            isSupported () const;

    const CacheImageData* getCacheImageData() const;
    std::string     getMD5   () const;

    int             getRank  () const;
    void            setRank  (int rank);

    int             getColorLabel  () const;
    void            setColorLabel  (int colorlabel);

    bool            getTrashed () const;
    void            setTrashed (bool trashed);

    void            addThumbnailListener (ThumbnailListener* tnl);
    void            removeThumbnailListener (ThumbnailListener* tnl);

    void            increaseRef ();
    void            decreaseRef ();

    void            updateCache (bool updatePParams = true, bool updateCacheImageData = true);
    void            saveThumbnail ();

    bool            openDefaultViewer(int destination);
    bool            imageLoad(bool loading);
};
