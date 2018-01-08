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
#include "edit.h"
#include "lockablecolorpicker.h"
#include <gtkmm.h>

class CropDisplayHandler
{

public:
    virtual ~CropDisplayHandler() {}
    virtual void cropImageUpdated    () {}
    virtual void cropWindowChanged   () {}
    virtual void initialImageArrived () {}
    virtual void setDisplayPosition  (int x, int y) {}
};

/**
 *  This class handle the displayed part of the image, ask for the initial data and process it so it can display it.
 *  Its position on the preview is handled not set by this class but by the CropHandlerListener (i.e. CropWindow) with which it works closely.
 */
class CropHandler final :
    public rtengine::DetailedCropListener,
    public rtengine::SizeListener
{
public:
    CropHandler ();
    ~CropHandler ();

    void    setDisplayHandler (CropDisplayHandler* l)
    {
        displayHandler = l;
    }
    void    setEditSubscriber      (rtedit::EditSubscriber* newSubscriber);

    void    newImage      (rtengine::StagedImageProcessor* ipc_, bool isDetailWindow);
    void    setZoom       (int z, int centerx = -1, int centery = -1);
    float   getZoomFactor ();
    double  getFitZoom    ();
    double  getFitCropZoom();
    bool    isFullDisplay ();
    void    setWSize      (int w, int h);
    void    getWSize      (int& w, int &h);
    void    getAnchorPosition (int& x, int& y);
    void    setAnchorPosition (int x, int y, bool update = true);
    void    moveAnchor    (int deltaX, int deltaY, bool update = true);
    void    centerAnchor  (bool update = true);
    void    getPosition   (int& x, int& y);
    void    getSize       (int& w, int& h);
    void    getFullImageSize (int& w, int& h);

    void    setEnabled (bool e);
    bool    getEnabled ();

    void    colorPick (const rtengine::Coord &pickerPos, float &r, float &g, float &b, float &rpreview, float &gpreview, float &bpreview, LockableColorPicker::Size size);

    rtengine::DetailedCrop* getCrop()
    {
        return crop;
    }

    // DetailedCropListener interface
    void    setDetailedCrop (rtengine::IImage8* im, rtengine::IImage8* imworking, rtengine::procparams::ColorManagementParams cmp,
                             rtengine::procparams::CropParams cp, int cx, int cy, int cw, int ch, int skip);
    bool    getWindow (int& cwx, int& cwy, int& cww, int& cwh, int& cskip);
    // SizeListener interface
    void    sizeChanged  (int w, int h, int ow, int oh);

    void    update  ();


    rtengine::procparams::CropParams cropParams;
    rtengine::procparams::ColorManagementParams colorParams;
    Glib::RefPtr<Gdk::Pixbuf> cropPixbuf;
    Glib::RefPtr<Gdk::Pixbuf> cropPixbuftrue;

    MyMutex cimg;

private:
    struct IdleHelper {
        CropHandler* cropHandler;
        bool destroyed;
        int pending;
    };

    void    compDim ();

    int zoom;               // scale factor (e.g. 5 if 1:5 scale) ; if 1:1 scale and bigger, factor is multiplied by 1000  (i.e. 1000 for 1:1 scale, 2000 for 2:1, etc...)
    int ww, wh;             // size of the crop's canvas on the screen ; might be bigger than the displayed image, but not smaller
    int cax, cay;           // clamped crop anchor's coordinate, i.e. point of the image that coincide to the center of the display area, expressed in image coordinates; cannot be outside the image's bounds; but if cax==cay==-1, designate the center of the image
    int cx, cy, cw, ch;     // position and size of the requested crop ; position expressed in image coordinates, so cx and cy might be negative and cw and ch higher than the image's 1:1 size
    int cropX, cropY, cropW, cropH; // cropPixbuf's displayed area (position and size), i.e. coordinates in 1:1 scale, i.e. cx, cy, cw & ch trimmed to the image's bounds
    bool enabled;
    unsigned char* cropimg;
    unsigned char* cropimgtrue;
    int cropimg_width, cropimg_height, cix, ciy, ciw, cih, cis;
    bool initial;
    bool isLowUpdatePriority;

    rtengine::StagedImageProcessor* ipc;
    rtengine::DetailedCrop* crop;

    CropDisplayHandler* displayHandler;
    IdleHelper* idle_helper;

    IdleRegister idle_register;
};

#endif
