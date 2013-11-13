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
#ifndef __CROPHANDLER__
#define __CROPHANDLER__

#include "../rtengine/rtengine.h"
#include "threadutils.h"
#include <gtkmm.h>

class CropHandlerListener {

    public:
        virtual void cropImageUpdated    () {}
        virtual void cropWindowChanged   () {}
        virtual void initialImageArrived () {}
};

class CropHandler;
struct CropHandlerIdleHelper {
    CropHandler* cropHandler;
    bool destroyed;
    int pending;
};

class CropHandler : public rtengine::DetailedCropListener, public rtengine::SizeListener {

    friend int createpixbufs (void* data);

    protected:
        int zoom;
        int ww, wh;             // size of the crop view on the screen
        int cx, cy, cw, ch;     // position and size of the requested crop
        int cropX, cropY, cropW, cropH; // position and size of the crop corresponding to cropPixbuf
        bool enabled;
        unsigned char* cropimg;
        unsigned char* cropimgtrue;
        int cropimg_width, cropimg_height, cix, ciy, ciw, cih, cis;
        bool initial;
        bool isLowUpdatePriority;

        rtengine::StagedImageProcessor* ipc;
        rtengine::DetailedCrop* crop;

        CropHandlerListener* listener;
        CropHandlerIdleHelper* chi;

        void    compDim ();

    public:

        void    update  ();


        rtengine::procparams::CropParams cropParams;
        rtengine::procparams::ColorManagementParams colorParams;
        Glib::RefPtr<Gdk::Pixbuf> cropPixbuf;
        Glib::RefPtr<Gdk::Pixbuf> cropPixbuftrue;

        MyMutex cimg;

        CropHandler ();
        ~CropHandler ();

        void    setCropHandlerListener (CropHandlerListener* l) { listener = l; }

        void    newImage    (rtengine::StagedImageProcessor* ipc_);
        void    setZoom     (int z, int centerx=-1, int centery=-1);
        double  getFitZoom  ();
        void    setWSize    (int w, int h);
        void    getWSize    (int& w, int &h);
        void    setPosition (int x, int y, bool update=true);
        void    getPosition (int& x, int& y);
        void    getSize     (int& w, int& h);
        void    getFullImageSize (int& w, int& h);

        void    setEnabled (bool e);
        bool    getEnabled ();

        // DetailedCropListener interface
        void    setDetailedCrop (rtengine::IImage8* im, rtengine::IImage8* imworking,rtengine::procparams::ColorManagementParams cmp,
                                 rtengine::procparams::CropParams cp, int cx, int cy, int cw, int ch, int skip);
        bool    getWindow (int& cwx, int& cwy, int& cww, int& cwh, int& cskip);
        // SizeListener interface
        void    sizeChanged  (int w, int h, int ow, int oh);

        void    cutRectToImgBounds (int& x, int& y, int& w, int& h);
};

#endif
