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
#ifndef _IMPROCCOORDINATOR_H_
#define _IMPROCCOORDINATOR_H_

#include <rtengine.h>
#include <improcfun.h>
#include <image8.h>
#include <image16.h>
#include <imagesource.h>
#include <procevents.h>
#include <dcrop.h>

namespace rtengine {

using namespace procparams;

class Crop;

class ImProcCoordinator : public StagedImageProcessor {

    friend class Crop;

    protected:
        Image16 *orig_prev;
        Image16 *oprevi;
        LabImage *oprevl;    
        LabImage *nprevl;    
        Image8 *previmg;
        ImageSource* imgsrc;
        
        int** buffer;
        
        SHMap* shmap;
        
        ColorTemp currWB;
        ColorTemp autoWB;

        bool awbComputed;

        ImProcFunctions ipf;

        int scale;
        bool allocated;
        
        void freeAll ();

        int tonecurve [65536];
        int lumacurve [65536];
        
        int vhist16[65536];
        int lhist16[65536];

        unsigned int rhist[256], ghist[256], bhist[256], Lhist[256];
        
        int fw, fh, tr, fullw, fullh;
        int pW, pH;

        ProgressListener* plistener;
        PreviewImageListener* imageListener;
        AutoExpListener* aeListener;
        HistogramListener* hListener;
        std::vector<SizeListener*> sizeListeners;
        
        std::vector<Crop*> crops;
        
        bool resultValid;
        
        Glib::Mutex minit;

        void progress (Glib::ustring str, int pr);
        void reallocAll ();
        void updateHistograms (int x1, int y1, int x2, int y2);
        void setScale (int prevscale, bool internal=false);
        void updatePreviewImage (int todo);

        Glib::Mutex mProcessing;
        ProcParams params;

        // members of the updater:
        Glib::Thread* thread;
        Glib::Mutex updaterThreadStart;
        Glib::Mutex paramsUpdateMutex;
        int  changeSinceLast;
        bool updaterRunning;
        ProcParams nextParams;
        bool destroying;

        void startProcessing ();
        void process ();
    
    public:

        ImProcCoordinator ();
        ~ImProcCoordinator ();
        void assign     (ImageSource* imgsrc);

        void        getParams (procparams::ProcParams* dst) { *dst = params; }

        ProcParams* getParamsForUpdate (ProcEvent change);
        void        paramsUpdateReady ();
        void        stopProcessing ();


        void setPreviewScale    (int scale) { setScale (scale); }
        int  getPreviewScale    () { return scale; }

        void fullUpdatePreviewImage  ();
        void fullUpdateDetailedCrops ();

        int getFullWidth ()     { return fullw; }
        int getFullHeight ()    { return fullh; }

        int getPreviewWidth ()     { return pW; }
        int getPreviewHeight ()    { return pH; }

        DetailedCrop* createCrop  ();

        void getAutoWB   (double& temp, double& green);
        void getCamWB    (double& temp, double& green);
        void getSpotWB   (int x, int y, int rectSize, double& temp, double& green);
        void getAutoCrop (double ratio, int &x, int &y, int &w, int &h);

        void setProgressListener (ProgressListener* pl)  { plistener = pl; }
        void setPreviewImageListener    (PreviewImageListener* il)     {imageListener = il; }
        void setSizeListener     (SizeListener* il)      {sizeListeners.push_back (il); }
        void delSizeListener     (SizeListener* il)      {std::vector<SizeListener*>::iterator it = std::find (sizeListeners.begin(), sizeListeners.end(), il); if (it!=sizeListeners.end()) sizeListeners.erase (it); }
        void setAutoExpListener  (AutoExpListener* ael)  {aeListener = ael; }
        void setHistogramListener(HistogramListener *h)  {hListener = h; }

        void saveInputICCReference (const Glib::ustring& fname);
        
        InitialImage*  getInitialImage () { return imgsrc; }
};
}
#endif
