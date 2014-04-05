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
#include <iomanip>

#include "cropwindow.h"
#include "options.h"
#include "guiutils.h"
#include "threadutils.h"
#include "../rtengine/mytime.h"
#include "imagearea.h"
#include "cursormanager.h"
#include "../rtengine/safegtk.h"
#include "../rtengine/rt_math.h"
#include "../rtengine/dcrop.h"

using namespace rtengine;

struct ZoomStep {   
    Glib::ustring label;
    double zoom;
    int czoom;
};

ZoomStep zoomSteps[] = {
		                {"  1%",  0.01,    100},
		                {"  2%",  0.02,    50},
		                {"  5%",  0.05,    20},
		                {"6.7%",  1.0/15.0,15},
		                {"  8%",  1.0/12.0,12},
		                {" 10%",  0.1,     10},
                        {"12.5%", 0.125,   8},
                        {" 14%",  1.0/7.0, 7},
                        {"16.6%", 1.0/6.0, 6},
                        {" 20%",  0.2,     5},
                        {" 25%",  0.25,    4},
                        {" 33%",  1.0/3.0, 3},
                        {" 50%",  0.5,     2},
                        {"100%",  1.0,     1000},
                        {"200%",  2.0,     2000},
                        {"300%",  3.0,     3000},
                        {"400%",  4.0,     4000},
                        {"500%",  5.0,     5000},
                        {"600%",  6.0,     6000},
                        {"700%",  7.0,     7000},
                        {"800%",  8.0,     8000}};
#define MAXZOOMSTEPS 20
#define ZOOM11INDEX  13

CropWindow::CropWindow (ImageArea* parent, rtengine::StagedImageProcessor* ipc_, bool isLowUpdatePriority_) 
    : onResizeArea(false), deleted(false), fitZoomEnabled(true), fitZoom(false), isLowUpdatePriority(isLowUpdatePriority_),
    backColor(options.bgcolor), decorated(true), titleHeight(30),
    sideBorderWidth(3), lowerBorderWidth(3), upperBorderWidth(1), sepWidth(2),
    xpos(30), ypos(30), imgX(0), imgY(0), imgW(1), imgH(1), iarea(parent),
    cropZoom(0), cropgl(NULL), pmlistener(NULL), observedCropWin(NULL), isFlawnOver(false) {

    Glib::RefPtr<Pango::Context> context = parent->get_pango_context () ;
    Pango::FontDescription fontd = context->get_font_description ();       
    fontd.set_weight (Pango::WEIGHT_BOLD);
    fontd.set_size(8*Pango::SCALE);
    context->set_font_description (fontd);
    cropLabel = "100%";
    Glib::RefPtr<Pango::Layout> cllayout = parent->create_pango_layout("1000%");
    exposeVersion = zoomVersion = 0;
 
    int iw, ih;
    cllayout->get_pixel_size (iw, ih);
    
    titleHeight = ih;
       
    resizeSurface = safe_create_from_png ("resize.png");
    bZoomIn  = new LWButton (safe_create_from_png ("gtk-zoom-in.png"),  0, NULL, LWButton::Left, LWButton::Center, "Zoom In");
    bZoomOut = new LWButton (safe_create_from_png ("gtk-zoom-out.png"), 1, NULL, LWButton::Left, LWButton::Center, "Zoom Out");
    bZoom100 = new LWButton (safe_create_from_png ("gtk-zoom-100.png"), 2, NULL, LWButton::Left, LWButton::Center, "Zoom 100/%");
    //bZoomFit = new LWButton (safe_create_from_png ("gtk-zoom-fit.png"), 3, NULL, LWButton::Left, LWButton::Center, "Zoom Fit");
    bClose   = new LWButton (safe_create_from_png ("gtk-close.png"),    4, NULL, LWButton::Right, LWButton::Center, "Close");
    
    buttonSet.add (bZoomIn);
    buttonSet.add (bZoomOut);
    buttonSet.add (bZoom100);
    buttonSet.add (bClose);
    
    buttonSet.setColors (Gdk::Color("black"), Gdk::Color("white"));
    buttonSet.setButtonListener (this);

    int bsw, bsh;
    buttonSet.getMinimalDimensions (bsw, bsh);

    if (bsh>titleHeight) 
        titleHeight = bsh;

    minWidth = bsw + iw + 2*sideBorderWidth;

    cropHandler.setCropHandlerListener (this);
    cropHandler.newImage (ipc_);
    cropHandler.setEnabled (true);

    state = SNormal;
}

void CropWindow::setPosition (int x, int y) {
    
    if (y<0)
        y = 0;
    xpos = x;
    ypos = y;
    if (decorated)
        buttonSet.arrangeButtons (xpos + sideBorderWidth, ypos + upperBorderWidth, width - 2*sideBorderWidth, titleHeight);
}

void CropWindow::getPosition (int& x, int& y) {
    
    x = xpos;
    y = ypos;
}

void CropWindow::getCropPosition (int& x, int& y) {
    
    int cropX, cropY;
    cropHandler.getPosition (cropX, cropY);
    if (state!=SCropImgMove) {
        x = cropX;
        y = cropY;
    }
    else {
        x = cropX + action_x;
        y = cropY + action_y;
    }
}

void CropWindow::getCropRectangle (int& x, int& y, int& w, int& h) {
    
    int cropX, cropY, cropW, cropH;
    cropHandler.getPosition (cropX, cropY);
    cropHandler.getSize (cropW, cropH);
    if (state!=SCropImgMove) {
        x = cropX;
        y = cropY;
    }
    else {
        x = cropX + action_x;
        y = cropY + action_y;
    }
    if (state!=SCropWinResize) {
        w = cropW;
        h = cropH;
    }
    else {
        w = imgAreaW;
        h = imgAreaH;
    }
    cropHandler.cutRectToImgBounds (x, y, w, h);
}

void CropWindow::setCropPosition (int x, int y) {
    
    cropHandler.setPosition (x, y);
    for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
        (*i)->cropPositionChanged (this);
}

void CropWindow::setSize (int w, int h, bool norefresh) {

    width = w;
    height = h;    
    
    fitZoom = false;
    
    if (width<minWidth)
        width = minWidth;
    if (height<64)
        height = 64;

    if (decorated) {
        imgAreaX = sideBorderWidth;
        imgAreaY = upperBorderWidth + titleHeight + sepWidth;
        imgAreaW = width - 2*sideBorderWidth;
        imgAreaH = height - lowerBorderWidth - titleHeight - sepWidth - upperBorderWidth;
        buttonSet.arrangeButtons (xpos + sideBorderWidth, ypos + upperBorderWidth, width - 2*sideBorderWidth, titleHeight);
    }
    else {
        imgAreaX = imgAreaY = 0;
        imgAreaW = width;
        imgAreaH = height;
    }
   
    if (!norefresh)
        cropHandler.setWSize (imgAreaW, imgAreaH);

    //iarea->redraw ();
}

void CropWindow::getSize (int& w, int& h) {
    
    w = width;
    h = height;
}

void CropWindow::getCropSize (int& w, int& h) {
    
    w = imgAreaW;
    h = imgAreaH;
}

bool CropWindow::isInside (int x, int y) {
    
    return x>=xpos && x<xpos+width && y>=ypos && y<ypos+height;
}

void CropWindow::leaveNotify (GdkEventCrossing* event) {
    EditSubscriber* subscriber = iarea->getCurrSubscriber();
    if (state==SNormal && subscriber && subscriber->getEditingType()==ET_PIPETTE) {
        iarea->pipetteVal[0] = iarea->pipetteVal[1] = iarea->pipetteVal[2] = -1.f;
        if (subscriber->mouseOver(0)) {
            iarea->redraw();
        }
    }
}

void CropWindow::flawnOver (bool isFlawnOver) {
    this->isFlawnOver = isFlawnOver;
}

void CropWindow::buttonPress (int button, int type, int bstate, int x, int y) {

    iarea->grabFocus (this);
    if (button==1 && type==GDK_2BUTTON_PRESS && onArea (CropImage, x, y) && (state==SNormal || state==SCropImgMove)) {
        if (fitZoomEnabled) {
            if (fitZoom) {
                state = SNormal;
                zoomVersion = exposeVersion;
                screenCoordToImage (x, y, action_x, action_y);
                changeZoom (ZOOM11INDEX, true, action_x, action_y);
                fitZoom = false;
            }
            else
                zoomFit ();
        }
        else 
            zoom11 ();
        state = SNormal;
    }
  //below code is no longer working/needed after adding buttons for each of the backColor values
  /*else if (button==1 && type==GDK_2BUTTON_PRESS && onArea (CropBorder, x, y)) {
        backColor = (backColor+1) % 3;
        options.bgcolor = backColor;
    }*/
    else if (button==1 && type==GDK_BUTTON_PRESS && state==SNormal && onArea (CropToolBar, x, y)) {
        if (!decorated || !buttonSet.pressNotify (x, y)) {
            state = SCropWinMove;
            action_x = x;
            action_y = y;
            press_x = xpos;
            press_y = ypos;
        }
    }
    else if (button==1 && type==GDK_BUTTON_PRESS && state==SNormal && onArea (CropResize, x, y)) {
        state = SCropWinResize;
        action_x = x;
        action_y = y;
        press_x = width;
        press_y = height;
    }
    else if (button==1 && type==GDK_BUTTON_PRESS && state==SNormal && onArea (CropImage, x, y)) {
        if (onArea (CropTopLeft, x, y)) {
            state = SResizeTL;
            press_x = x;
            action_x = cropHandler.cropParams.x;
            press_y = y;
            action_y = cropHandler.cropParams.y;
        }
        else if (onArea (CropTopRight, x, y)) {
            state = SResizeTR;
            press_x = x;
            action_x = cropHandler.cropParams.w;
            press_y = y;
            action_y = cropHandler.cropParams.y;
        }
        else if (onArea (CropBottomLeft, x, y)) {
            state = SResizeBL;
            press_x = x;
            action_x = cropHandler.cropParams.x;
            press_y = y;
            action_y = cropHandler.cropParams.h;
        }
        else if (onArea (CropBottomRight, x, y)) {
            state = SResizeBR;
            press_x = x;
            action_x = cropHandler.cropParams.w;
            press_y = y;
            action_y = cropHandler.cropParams.h;
        }
        else if (onArea (CropTop, x, y)) {
            state = SResizeH1;
            press_y = y;
            action_y = cropHandler.cropParams.y;
        }
        else if (onArea (CropBottom, x, y)) {
            state = SResizeH2;
            press_y = y;
            action_y = cropHandler.cropParams.h;
        }
        else if (onArea (CropLeft, x, y)) {
            state = SResizeW1;
            press_x = x;
            action_x = cropHandler.cropParams.x;
        }
        else if (onArea (CropRight, x, y)) {
            state = SResizeW2;
            press_x = x;
            action_x = cropHandler.cropParams.w;
        }
        else if ((bstate & GDK_SHIFT_MASK) && onArea (CropInside, x, y)) {
            state = SCropMove;
            press_x = x;
            press_y = y;
            action_x = cropHandler.cropParams.x;
            action_y = cropHandler.cropParams.y;
        }
        else if (iarea->getToolMode () == TMHand) {
            EditSubscriber *editSubscriber = iarea->getCurrSubscriber();

            if      (button==1 && editSubscriber && cropgl && cropgl->inImageArea(iarea->posImage.x, iarea->posImage.y) && (editSubscriber->getEditingType()==ET_OBJECTS && iarea->object>-1) ) {
                editSubscriber->button1Pressed(bstate);
                state=SEditDrag;
                press_x = x;
                press_y = y;
                action_x = 0;
                action_y = 0;
            }
            else if (onArea (CropObserved, x, y)) {
                state = SObservedMove;
                press_x = x;
                press_y = y;
                action_x = 0;
                action_y = 0;
            }
            else if (button==1 && editSubscriber && cropgl && cropgl->inImageArea(iarea->posImage.x, iarea->posImage.y) && (editSubscriber->getEditingType()==ET_PIPETTE && (bstate & GDK_CONTROL_MASK)) ) {
                editSubscriber->button1Pressed(bstate);
                state=SEditDrag;
                press_x = x;
                press_y = y;
                action_x = 0;
                action_y = 0;
            }
            else if(zoomSteps[cropZoom].zoom > cropHandler.getFitZoom()) { // only allow move when image is only partial visible
                state = SCropImgMove;
                press_x = x;
                press_y = y;
                action_x = 0;
                action_y = 0;
            }
        }
        else if (onArea (CropObserved, x, y)) {
            state = SObservedMove;
            press_x = x;
            press_y = y;
        }
        else if (iarea->getToolMode () == TMStraighten) {
            state = SRotateSelecting;
            press_x = x;
            press_y = y;
            action_x = x;
            action_y = y;
            rot_deg = 0;
        }
        else if (iarea->getToolMode () == TMSpotWB) {
            screenCoordToImage (x, y, action_x, action_y);
            iarea->spotWBSelected (action_x, action_y);
        }
        else if (iarea->getToolMode () == TMCropSelect && cropgl) {
            state = SCropSelecting;
            screenCoordToImage (x, y, press_x, press_y);
            cropHandler.cropParams.enabled = true;
            cropHandler.cropParams.x = press_x;
            cropHandler.cropParams.y = press_y;
            cropHandler.cropParams.w = cropHandler.cropParams.h = 1;
            cropgl->cropInit (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        }
    }
    if (button==3) {
        iarea->pipetteVal[0] = iarea->pipetteVal[1] = iarea->pipetteVal[2] = -1.f;
        EditSubscriber *editSubscriber = iarea->getCurrSubscriber();
        if (editSubscriber && editSubscriber->getEditingType()==ET_PIPETTE)
            editSubscriber->mouseOver(0);
        state = SNormal;
        iarea->setToolHand ();
        if (pmhlistener) {
          pmhlistener->toggleFreeze();
        }
    }
    iarea->redraw ();
    updateCursor (x, y);
}

void CropWindow::buttonRelease (int button, int num, int bstate, int x, int y) {

    EditSubscriber *editSubscriber = iarea->getCurrSubscriber();

    bool needRedraw = false;
    if (state==SCropWinResize) {
        setSize (press_x + x - action_x, press_y + y - action_y);
        state = SNormal;
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropWindowSizeChanged (this);
        needRedraw = true;
    }
    else if (state==SCropImgMove) {
        int cropX, cropY;
        cropHandler.getPosition (cropX, cropY);
        cropHandler.setPosition (cropX + action_x, cropY + action_y);
        cropHandler.getPosition (cropX, cropY);
        state = SNormal;
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropPositionChanged (this);
        needRedraw = true;
    }
    else if (state==SRotateSelecting) {
        iarea->straightenReady (rot_deg);
        iarea->setToolHand ();
        needRedraw = true;
    }
    else if (state==SObservedMove) {
        observedCropWin->remoteMoveReady ();
        state = SNormal;
        needRedraw = true;
    }
    else if (state==SEditDrag) {
        editSubscriber->button1Released();
        iarea->object = 0;
        iarea->deltaImage.set(0, 0);
        iarea->deltaScreen.set(0, 0);
        iarea->deltaPrevImage.set(0, 0);
        iarea->deltaPrevScreen.set(0, 0);
        state = SNormal;
        needRedraw = true;
    }
    if (cropgl && (state==SCropSelecting || state==SResizeH1 || state==SResizeH2 || state==SResizeW1 || state==SResizeW2 || state==SResizeTL || state==SResizeTR || state==SResizeBL || state==SResizeBR || state==SCropMove))
    {
        cropgl->cropManipReady ();
        iarea->setToolHand ();
        needRedraw = true;
    }

    if (decorated)
        buttonSet.releaseNotify (x, y);

    if (deleted) {
        iarea->flawnOverWindow = NULL;
        delete this;
        return;
    }

    state = SNormal;
    iarea->grabFocus (NULL);
    if (needRedraw)
        iarea->redraw ();
    updateCursor (x, y);
}

void CropWindow::pointerMoved (int bstate, int x, int y) {

    EditSubscriber *editSubscriber = iarea->getCurrSubscriber();

    if (state==SCropWinMove) {
        setPosition (press_x + x - action_x, press_y + y - action_y);
        iarea->redraw ();
    }
    else if (state==SCropWinResize) {
        setSize (press_x + x - action_x, press_y + y - action_y, true);
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropWindowSizeChanged (this);
        iarea->redraw ();
    }
    else if (state==SCropImgMove) {
        // multiplier is the amplification factor ; disabled if the user selected "1" (no amplification)
        double factor = options.panAccelFactor == 1 ? 1.0 : options.panAccelFactor * zoomSteps[cropZoom].zoom;
        // never move the preview slower than the cursor
        if (factor < 1.0)
            factor = 1.0;
        action_x =  (press_x - x) / zoomSteps[cropZoom].zoom * factor;
        action_y =  (press_y - y) / zoomSteps[cropZoom].zoom * factor;
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropPositionChanged (this);
        iarea->redraw ();
    }
    else if (state==SRotateSelecting) {
        action_x = x;
        action_y = y;
        iarea->redraw ();
    }
    else if (state==SNormal && iarea->getToolMode () == TMSpotWB) {
        action_x = x;
        action_y = y;
        iarea->redraw ();
    }
    else if (state==SResizeH1 && cropgl) { 
        int oy = cropHandler.cropParams.y;
        cropHandler.cropParams.y = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.h += oy - cropHandler.cropParams.y;
        cropgl->cropHeight1Resized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeH2 && cropgl) { 
        cropHandler.cropParams.h = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropgl->cropHeight2Resized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeW1 && cropgl) { 
        int ox = cropHandler.cropParams.x;
        cropHandler.cropParams.x = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.w += ox - cropHandler.cropParams.x;
        cropgl->cropWidth1Resized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeW2 && cropgl) { 
        cropHandler.cropParams.w = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropgl->cropWidth2Resized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeTL && cropgl) { 
        int ox = cropHandler.cropParams.x;
        cropHandler.cropParams.x = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.w += ox - cropHandler.cropParams.x;
        int oy = cropHandler.cropParams.y;
        cropHandler.cropParams.y = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.h += oy - cropHandler.cropParams.y;
        cropgl->cropTopLeftResized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeTR && cropgl) { 
        cropHandler.cropParams.w = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        int oy = cropHandler.cropParams.y;
        cropHandler.cropParams.y = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.h += oy - cropHandler.cropParams.y;
        cropgl->cropTopRightResized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeBL && cropgl) { 
        int ox = cropHandler.cropParams.x;
        cropHandler.cropParams.x = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.w += ox - cropHandler.cropParams.x;
        cropHandler.cropParams.h = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropgl->cropBottomLeftResized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SResizeBR && cropgl) { 
        cropHandler.cropParams.w = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.h = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropgl->cropBottomRightResized (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SCropMove && cropgl) { 
        cropHandler.cropParams.x = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.y = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropgl->cropMoved (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SCropSelecting && cropgl) { 
        screenCoordToImage (x, y, action_x, action_y);
        int cx1 = press_x, cy1 = press_y;
        int cx2 = action_x, cy2 = action_y;
        cropgl->cropResized (cx1, cy1, cx2, cy2);
        if (cx2 > cx1) {
            cropHandler.cropParams.x = cx1;
            cropHandler.cropParams.w = cx2 - cx1 + 1;
        }
        else {
            cropHandler.cropParams.x = cx2;
            cropHandler.cropParams.w = cx1 - cx2 + 1;
        }
        if (cy2 > cy1) {
            cropHandler.cropParams.y = cy1;
            cropHandler.cropParams.h = cy2 - cy1 + 1;
        }
        else {
            cropHandler.cropParams.y = cy2;
            cropHandler.cropParams.h = cy1 - cy2 + 1;
        }
        iarea->redraw ();
    }
    else if (state==SObservedMove) {
        observedCropWin->remoteMove ((x - press_x)/zoomSteps[cropZoom].zoom, (y - press_y)/zoomSteps[cropZoom].zoom);
        iarea->redraw ();
    }
    else if (editSubscriber) {
        rtengine::Crop* crop = static_cast<rtengine::Crop*>(cropHandler.getCrop());
        if (state==SNormal) {
            Coord imgPos;
            action_x = x;
            action_y = y;
            screenCoordToImage (x, y, imgPos.x, imgPos.y);

            iarea->posImage.set(imgPos.x, imgPos.y);
            iarea->posScreen.set(x, y);

            Coord cropPos;
            screenCoordToCropBuffer(x, y, cropPos.x, cropPos.y);
            if (editSubscriber->getEditingType()==ET_PIPETTE) {
                iarea->object = onArea (CropImage, x, y) && !onArea (CropObserved, x, y) ? 1 : 0;
                //iarea->object = cropgl && cropgl->inImageArea(iarea->posImage.x, iarea->posImage.y) ? 1 : 0;
                if (iarea->object) {
                    crop->getPipetteData(iarea->pipetteVal, cropPos.x, cropPos.y, iarea->getPipetteRectSize());
                    //printf("PipetteData:  %.3f  %.3f  %.3f\n", iarea->pipetteVal[0], iarea->pipetteVal[1], iarea->pipetteVal[2]);
                }
                else {
                    iarea->pipetteVal[0] = iarea->pipetteVal[1] = iarea->pipetteVal[2] = -1.f;
                }
            }
            else if (editSubscriber->getEditingType()==ET_OBJECTS) {
                if (onArea (CropImage, x, y))
                    iarea->object = crop->getObjectID(cropPos);
                else
                    iarea->object = -1;
            }

            if (editSubscriber->mouseOver(bstate))
                iarea->redraw ();
        }
        else if (state==SEditDrag) {
            Coord currPos;
            action_x = x;
            action_y = y;
            Coord oldPosImage = iarea->posImage+iarea->deltaImage;
            //printf(">>> IMG / ImgPrev(%d x %d) = (%d x %d) + (%d x %d)\n", oldPosImage.x, oldPosImage.y, iarea->posImage.x, iarea->posImage.y, iarea->deltaImage.x, iarea->deltaImage.y);
            screenCoordToImage (x, y, currPos.x, currPos.y);
            iarea->deltaImage     = currPos - iarea->posImage;
            iarea->deltaPrevImage = currPos - oldPosImage;
            //printf("          action_ & xy (%d x %d) -> (%d x %d) = (%d x %d) + (%d x %d) / deltaPrev(%d x %d)\n", action_x, action_y, currPos.x, currPos.y, iarea->posImage.x, iarea->posImage.y, iarea->deltaImage.x, iarea->deltaImage.y, iarea->deltaPrevImage.x, iarea->deltaPrevImage.y);

            Coord oldPosScreen = iarea->posScreen+iarea->deltaScreen;
            //printf(">>> SCR / ScrPrev(%d x %d) = (%d x %d) + (%d x %d)\n", oldPosScreen.x, oldPosScreen.y, iarea->posScreen.x, iarea->posScreen.y, iarea->deltaScreen.x, iarea->deltaScreen.y);
            currPos.set(x, y);
            iarea->deltaScreen     = currPos - iarea->posScreen;
            iarea->deltaPrevScreen = currPos - oldPosScreen;
            //printf("          action_ & xy (%d x %d) -> (%d x %d) = (%d x %d) + (%d x %d) / deltaPrev(%d x %d)\n", action_x, action_y, currPos.x, currPos.y, iarea->posScreen.x, iarea->posScreen.y, iarea->deltaScreen.x, iarea->deltaScreen.y, iarea->deltaPrevScreen.x, iarea->deltaPrevScreen.y);

            if (editSubscriber->drag(bstate))
                iarea->redraw ();
        }
    }
    updateCursor (x, y);
    
    bool oRA = onArea (CropResize, x, y);
    if (oRA!=onResizeArea) {
        onResizeArea = oRA;
        iarea->redraw ();
    }

    if (decorated)
        buttonSet.motionNotify (x, y);

    if (pmlistener) {
        int mx, my;
        screenCoordToImage (x, y, mx, my);
        if (!onArea (CropImage, x, y) || !cropHandler.cropPixbuf) {
            cropHandler.getFullImageSize(mx,my);
            pmlistener->pointerMoved (false, cropHandler.colorParams.working, mx, my, -1, -1, -1);
            if (pmhlistener) pmhlistener->pointerMoved (false, cropHandler.colorParams.working, mx, my, -1, -1, -1);
        }
        else {
            /*MyMutex::MyLock lock(cropHandler.cimg);

            int vx = x - xpos - imgX;
            int vy = y - ypos - imgY;
            guint8* pix = cropHandler.cropPixbuf->get_pixels() + vy*cropHandler.cropPixbuf->get_rowstride() + vx*3;
            if (vx < cropHandler.cropPixbuf->get_width() && vy < cropHandler.cropPixbuf->get_height())
            pmlistener->pointerMoved (true, mx, my, pix[0], pix[1], pix[2]);

            */

            cropHandler.cimg.lock ();
            int vx = x - xpos - imgX;
            int vy = y - ypos - imgY;
//          guint8* pix = cropHandler.cropPixbuf->get_pixels() + vy*cropHandler.cropPixbuf->get_rowstride() + vx*3;
//          if (vx < cropHandler.cropPixbuf->get_width() && vy < cropHandler.cropPixbuf->get_height())
//              pmlistener->pointerMoved (true, mx, my, pix[0], pix[1], pix[2]);
            int imwidth = cropHandler.cropPixbuf->get_width();
            int imheight = cropHandler.cropPixbuf->get_height();
            guint8* pix = cropHandler.cropPixbuftrue->get_pixels() + vy*cropHandler.cropPixbuf->get_rowstride() + vx*3;
            if (vx < imwidth && vy < imheight) {
                pmlistener->pointerMoved (true, cropHandler.colorParams.working, mx, my, pix[0], pix[1], pix[2]);
                if (pmhlistener)
                    pmhlistener->pointerMoved (true, cropHandler.colorParams.working, mx, my, pix[0], pix[1], pix[2]);
            }
            cropHandler.cimg.unlock ();
        }
    }
}

bool CropWindow::onArea (CursorArea a, int x, int y) {

    int CROPRESIZEBORDER = 6 / zoomSteps[cropZoom].zoom;
    int x1, y1, w, h;
    switch (a) {
        case CropWinButtons:
            return decorated && buttonSet.inside (x, y);
        case CropToolBar:
            return x>xpos && y>ypos && x<xpos+width-1 && y<ypos+imgAreaY;
        case CropImage:
            return x>=xpos+imgX && y>=ypos+imgY && x<xpos+imgX+imgW && y<ypos+imgY+imgH;
        case CropBorder:
            return 
                (x>=xpos+imgAreaX && y>=ypos+imgAreaY && x<xpos+imgAreaX+imgAreaW && y<ypos+imgAreaY+imgAreaH) &&
                !(x>=xpos+imgX && y>=ypos+imgY && x<xpos+imgX+imgW && y<ypos+imgY+imgH);
        case CropTopLeft:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>=cropHandler.cropParams.y-CROPRESIZEBORDER && 
                y1<=cropHandler.cropParams.y+CROPRESIZEBORDER &&
                y>=ypos+imgY &&
                x1>=cropHandler.cropParams.x-CROPRESIZEBORDER && 
                x1<=cropHandler.cropParams.x+CROPRESIZEBORDER &&
                x>=xpos+imgX;
        case CropTopRight:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>=cropHandler.cropParams.y-CROPRESIZEBORDER && 
                y1<=cropHandler.cropParams.y+CROPRESIZEBORDER &&
                y>=ypos+imgY &&
                x1>=cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                x1<=cropHandler.cropParams.x+cropHandler.cropParams.w-1+CROPRESIZEBORDER &&
                x<xpos+imgX+imgW;
        case CropBottomLeft:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>=cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                y1<=cropHandler.cropParams.y+cropHandler.cropParams.h-1+CROPRESIZEBORDER &&
                y<ypos+imgY+imgH &&
                x1>=cropHandler.cropParams.x-CROPRESIZEBORDER && 
                x1<=cropHandler.cropParams.x+CROPRESIZEBORDER &&
                x>=xpos+imgX;
        case CropBottomRight:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>=cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                y1<=cropHandler.cropParams.y+cropHandler.cropParams.h-1+CROPRESIZEBORDER &&
                y<ypos+imgY+imgH &&
                x1>=cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                x1<=cropHandler.cropParams.x+cropHandler.cropParams.w-1+CROPRESIZEBORDER &&
                x<xpos+imgX+imgW;
        case CropTop:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                x1>cropHandler.cropParams.x+CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                y1>cropHandler.cropParams.y-CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+CROPRESIZEBORDER &&
                y>=ypos+imgY;
        case CropBottom:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                x1>cropHandler.cropParams.x+CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                y1>cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1+CROPRESIZEBORDER &&
                y<ypos+imgY+imgH;
        case CropLeft:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>cropHandler.cropParams.y+CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                x1>cropHandler.cropParams.x-CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+CROPRESIZEBORDER &&
                x>=xpos+imgX;
        case CropRight:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>cropHandler.cropParams.y+CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                x1>cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1+CROPRESIZEBORDER &&
                x<xpos+imgX+imgW;
        case CropInside:
            screenCoordToImage (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>cropHandler.cropParams.y && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1 && 
                x1>cropHandler.cropParams.x && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1;
        case CropResize:
            return decorated && x>=xpos+width-16 && y>=ypos+height-16 && x<xpos+width && y<ypos+height;
        case CropObserved:
            if (!observedCropWin)
                return false;
            getObservedFrameArea (x1, y1, w, h);
            return x>x1-6 && y>y1-6 && x<x1+w-1+6 && y<y1+h-1+6 &&
                !(x>x1+2 && y>y1+2 && x<x1+w-1-2 && y<y1+h-1-2);
    }
    return false;
}

void CropWindow::updateCursor (int x, int y) {

    EditSubscriber *editSubscriber = iarea->getCurrSubscriber();
    ToolMode tm = iarea->getToolMode ();

    if (state==SNormal) {
        if (onArea (CropWinButtons, x, y)) 
            cursorManager.setCursor (iarea->get_window(), CSArrow);
        else if (onArea (CropToolBar, x, y))
            cursorManager.setCursor (iarea->get_window(), CSMove);
        else if (onArea (CropResize, x, y)) 
            cursorManager.setCursor (iarea->get_window(), CSResizeDiagonal);
        else if (tm==TMHand && (onArea (CropTopLeft, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeTopLeft);
        else if (tm==TMHand && (onArea (CropTopRight, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeTopRight);
        else if (tm==TMHand && (onArea (CropBottomLeft, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeBottomLeft);
        else if (tm==TMHand && (onArea (CropBottomRight, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeBottomRight);
        else if (tm==TMHand && (onArea (CropTop, x, y) || onArea (CropBottom, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeHeight);
        else if (tm==TMHand && (onArea (CropLeft, x, y) || onArea (CropRight, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeWidth);
        else if (onArea (CropImage, x, y)) {
            int objectID = -1;
            if (editSubscriber) {
                Coord cropPos;
                screenCoordToCropBuffer(iarea->posScreen.x, iarea->posScreen.y, cropPos.x, cropPos.y);
                objectID = static_cast<rtengine::Crop*>(cropHandler.getCrop())->getObjectID(cropPos);
            }
            if (objectID > -1) {
                cursorManager.setCursor (iarea->get_window(), editSubscriber->getCursor(objectID));
            }
            else if (tm==TMHand) {
                if (onArea (CropObserved, x, y))
                    cursorManager.setCursor (iarea->get_window(), CSMove);
                else
                    cursorManager.setCursor (iarea->get_window(), CSOpenHand);
            }
            else if (tm==TMSpotWB)
                cursorManager.setCursor (iarea->get_window(), CSSpotWB);
            else if (tm==TMCropSelect)
                cursorManager.setCursor (iarea->get_window(), CSCropSelect);
            else if (tm==TMStraighten)
                cursorManager.setCursor (iarea->get_window(), CSStraighten);
        }
        else
            cursorManager.setCursor (iarea->get_window(), CSArrow);
    }
    else if (state==SCropSelecting)
        cursorManager.setCursor (iarea->get_window(), CSCropSelect);
    else if (state==SRotateSelecting) 
        cursorManager.setCursor (iarea->get_window(), CSStraighten);
    else if (state==SCropMove || state==SCropWinMove || state==SObservedMove)
        cursorManager.setCursor (iarea->get_window(), CSMove);
    else if (state==SHandMove || state==SCropImgMove)
        cursorManager.setCursor (iarea->get_window(), CSClosedHand);
    else if (state==SResizeW1 || state==SResizeW2)
        cursorManager.setCursor (iarea->get_window(), CSResizeWidth);
    else if (state==SResizeH1 || state==SResizeH2)
        cursorManager.setCursor (iarea->get_window(), CSResizeHeight);
    else if (state==SResizeTL)
        cursorManager.setCursor (iarea->get_window(), CSResizeTopLeft);
    else if (state==SResizeTR)
        cursorManager.setCursor (iarea->get_window(), CSResizeTopRight);
    else if (state==SResizeBL)
        cursorManager.setCursor (iarea->get_window(), CSResizeBottomLeft);
    else if (state==SResizeBR)
        cursorManager.setCursor (iarea->get_window(), CSResizeBottomRight);
    else if (state==SCropWinResize)
        cursorManager.setCursor (iarea->get_window(), CSResizeDiagonal);
}

void CropWindow::expose (Cairo::RefPtr<Cairo::Context> cr) {
    MyMutex::MyLock lock(cropHandler.cimg);

    if (decorated)
        drawDecoration (cr);

    int x = xpos, y = ypos, h = height, w = width;
 
    // draw the background
    backColor = iarea->previewModePanel->GetbackColor();
    options.bgcolor = backColor;
    if (backColor==0) {
        Gdk::Color cback = iarea->get_style()->get_bg(Gtk::STATE_NORMAL);
        cr->set_source_rgb (cback.get_red_p(), cback.get_green_p(), cback.get_blue_p());
    }
    else if (backColor==1)
        cr->set_source_rgb (0,0,0);
    else if (backColor==2)
        cr->set_source_rgb (1,1,1);

    cr->set_line_width (0.);
    cr->rectangle (x+imgAreaX, y+imgAreaY, imgAreaW, imgAreaH);
    cr->stroke_preserve ();
    cr->fill ();


    // draw image
    if (state==SCropImgMove || state==SCropWinResize) {
        // draw a rough image
        int cropX, cropY;
        cropHandler.getPosition (cropX, cropY);
        if (state==SCropImgMove) {
            cropX += action_x;
            cropY += action_y;
        }
        Glib::RefPtr<Gdk::Pixbuf> rough = iarea->getPreviewHandler()->getRoughImage (cropX, cropY, imgAreaW, imgAreaH, zoomSteps[cropZoom].zoom);
        if (rough) {
            iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), rough, 0, 0, x+imgAreaX+(imgAreaW-rough->get_width())/2, y+imgAreaY+(imgAreaH-rough->get_height())/2, -1, -1, Gdk::RGB_DITHER_NORMAL, 0, 0);  
//            if (cropHandler.cropParams.enabled)
//                drawCrop (cr, x+imgX, y+imgY, imgW, imgH, cropX, cropY, zoomSteps[cropZoom].zoom, cropHandler.cropParams);
        }
        if (observedCropWin)
            drawObservedFrame (cr);
    }
    else {
        if (cropHandler.cropPixbuf) {
            imgW = cropHandler.cropPixbuf->get_width ();
            imgH = cropHandler.cropPixbuf->get_height ();
            imgX = imgAreaX + (imgAreaW-imgW)/2;
            imgY = imgAreaY + (imgAreaH-imgH)/2;
            exposeVersion++;

            bool showcs = iarea->indClippedPanel->showClippedShadows();
            bool showch = iarea->indClippedPanel->showClippedHighlights();
            const bool showR  = iarea->previewModePanel->showR(); // will show clipping if R channel is clipped
            const bool showG  = iarea->previewModePanel->showG(); // will show clipping if G channel is clipped
            const bool showB  = iarea->previewModePanel->showB(); // will show clipping if B channel is clipped
            const bool showL  = iarea->previewModePanel->showL(); // will show clipping if L value   is clipped
            const bool showFocusMask  = iarea->previewModePanel->showFocusMask();

            // While the Right-side ALT is pressed, auto-enable highlight and shadow clipping indicators
            // TODO: Add linux/MacOS specific functions for alternative
            #ifdef WIN32
            if (GetKeyState(VK_RMENU)<0) {
                showcs=true; showch=true;
            }
            #endif

            if (showcs || showch || showR || showG || showB || showL || showFocusMask) {
                Glib::RefPtr<Gdk::Pixbuf> tmp = cropHandler.cropPixbuf->copy ();
                guint8* pix = tmp->get_pixels();
                guint8* pixWrkSpace = cropHandler.cropPixbuftrue->get_pixels();

                const int pixRowStride = tmp->get_rowstride ();
                const int pixWSRowStride = cropHandler.cropPixbuftrue->get_rowstride ();

				const int bHeight = tmp->get_height();
				const int bWidth = tmp->get_width();

                if (showFocusMask) { // modulate preview to display focus mask
					const int blur_radius2 = 1;								// radius of small kernel. 1 => 3x3 kernel
					const int blur_dim2 = 2*blur_radius2+1;					// dimension of small kernel
					const int blur_radius = (blur_dim2*blur_dim2)/2; 		// radius of big kernel
					const float kernel_size = SQR(2.f*blur_radius+1.f); 	// count of pixels in the big blur kernel
					const float rkernel_size = 1.0f/kernel_size;			// reciprocal of kernel_size to avoid divisions
					const float kernel_size2 = SQR(2.f*blur_radius2+1.f);   // count of pixels in the small blur kernel
					const float rkernel_size2 = 1.0f/kernel_size2;			// reciprocal of kernel_size to avoid divisions

					// aloocate buffer for precalculated Luminance		
					float* tmpL = (float*)malloc(bHeight * bWidth * sizeof(float) );
					// aloocate buffers for sums and sums of squares of small kernel
					float* tmpLsum = (float*)malloc((bHeight) * (bWidth) * sizeof(float) );
					float* tmpLsumSq = (float*)malloc((bHeight) * (bWidth) * sizeof(float) );
					float* tmpstdDev2 = (float*)malloc((bHeight) * (bWidth) * sizeof(float) );
					float maxstdDev_L2 = 0.f;

#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
#pragma omp for
#endif
					// precalculate Luminance
					for(int i=0;i<bHeight;i++) {
						guint8* currWS = pixWrkSpace + i*pixWSRowStride;
						float*  currL = tmpL + i*bWidth;
						for(int j=0;j<bWidth;j++) {
							*currL = 0.299f*(currWS)[0]+0.587f*(currWS)[1]+0.114f*(currWS)[2];
							currL++;
							currWS+=3;
						}
					}

					float maxthrstdDev_L2 = 0.f;
#ifdef _OPENMP
#pragma omp for nowait
#endif
					// precalculate sum and sum of squares of small kernel
					for(int i=blur_radius2;i<bHeight-blur_radius2;i++) {
						for(int j=blur_radius2;j<bWidth-blur_radius2;j++) {
							float sumL = 0.f;
							float sumLSqu = 0.f;
							for(int kh=-blur_radius2;kh<=blur_radius2;kh++) {
								for(int kw=-blur_radius2;kw<=blur_radius2;kw++) {
									float curL = tmpL[(i+kh)*bWidth + j + kw];
									sumL += curL;
									sumLSqu += SQR(curL);
								}
							}
							tmpLsum[i*bWidth+j] = sumL;
							tmpLsumSq[i*bWidth+j] = sumLSqu;
							float stdDev_L2 = rkernel_size2 * sqrtf(sumLSqu * kernel_size2 - sumL * sumL);
							if(stdDev_L2 > maxthrstdDev_L2)
								maxthrstdDev_L2 = stdDev_L2;
							tmpstdDev2[i*bWidth+j] = stdDev_L2;
						}
					}
#pragma omp critical
{
					if(maxthrstdDev_L2 > maxstdDev_L2)
						maxstdDev_L2 = maxthrstdDev_L2;
}
}

					const float focus_thresh = 80.f;
					maxstdDev_L2 = std::min(maxstdDev_L2,focus_thresh);
					const float focus_threshby10 = focus_thresh / 10.f;
					#ifdef _OPENMP
					#pragma omp parallel for schedule(dynamic,16)
					#endif
					for (int i=blur_radius+1; i<bHeight-blur_radius; i++) {
						guint8* curr = pix + i*pixRowStride + 3*(blur_radius+1);
						guint8* currWs = pixWrkSpace + i*pixWSRowStride + 3*(blur_radius+1);
						for (int j=blur_radius+1; j<bWidth-blur_radius; j++) {

							//*************
							// Copyright (c) 2011 Michael Ezra michael@michaelezra.com
							// determine if pixel is in the sharp area of the image using
							// standard deviation analysis on two different scales
							//float focus_thresh2;
							//float opacity = 0.9;//TODO: implement opacity
							//TODO: evaluate effects of altering sampling frequency


							//TODO: dynamically determine appropriate values based on image analysis

							// calculate average in +-blur_radius pixels area around the current pixel
							// speed up: calculate sum of squares in the same loops

							float sum_L = 0.f;
							float sumsq_L = 0.f;
							// use precalculated values of small kernel to reduce number of iterations
							for (int kh=-blur_radius+blur_radius2; kh<=blur_radius-blur_radius2;kh+=blur_dim2){
								float* currLsum = &tmpLsum[(i+kh)*bWidth+j-blur_radius+1];
								float* currLsumSqu = &tmpLsumSq[(i+kh)*bWidth+j-blur_radius+1];
								for (int k=-blur_radius+blur_radius2; k<=blur_radius-blur_radius2;k+=blur_dim2,currLsum+=blur_dim2,currLsumSqu+=blur_dim2){
									sum_L += *currLsum;
									sumsq_L += *currLsumSqu;
								}
							}
							float sum_L2 = tmpLsum[i*bWidth+j];
							float sumsq_L2 = tmpLsumSq[i*bWidth+j];
							//*************
							// averages
							// Optimized formulas to avoid divisions
							float stdDev_L = rkernel_size * sqrtf(sumsq_L * kernel_size - sum_L * sum_L);
							float stdDev_L2 = tmpstdDev2[i*bWidth+j];
//							float stdDev_L2 = rkernel_size2 * sqrtf(sumsq_L2 * kernel_size2 - sum_L2 * sum_L2);

									//TODO: try to normalize by average L of the entire (preview) image

									//detection method 1: detect focus in features
									//there is no strict condition between stdDev_L and stdDev_L2 themselves
	/*                                if (stdDev_L2>focus_thresh2
									&& (stdDev_L <focus_thresh)){ // this excludes false positives due to high contrast edges

										curr[1]=255;
										curr[0]=0;
										curr[2]=0;

									}*/

									//detection method 2: detect focus in texture
									// key point is std deviation on lower scale is higher than for the larger scale
									// plus some boundary conditions
							if (focus_thresh >= stdDev_L2 //TODO: could vary this to bypass noise better
													&& stdDev_L2 > stdDev_L //this is the key to select fine detail within lower contrast on larger scale
																&& stdDev_L > focus_threshby10 //options.highlightThreshold
								){
									// transpareny depends on sdtDev_L2 and maxstdDev_L2
									float transparency = 1.f - std::min(stdDev_L2 / maxstdDev_L2, 1.0f) ;
									// first row of circle
									guint8* currtmp = &curr[0] + (-3*pixRowStride);
									guint8* currtmpWS = &currWs[0] + (-3*pixWSRowStride);
									for(int jj=-3;jj<=3;jj+=3) {
										guint8* currtmpl=currtmp+jj;
										guint8* currtmpWSl=currtmpWS+jj;
										//transparent green
										currtmpl[0] = transparency * currtmpWSl[0];
										currtmpl[1] = transparency * currtmpWSl[1] + (1.f-transparency)*255.f;
										currtmpl[2] = transparency * currtmpWSl[2];
									}
									// second row of circle
									currtmp = &curr[0] + (-2*pixRowStride);
									currtmpWS = &currWs[0] + (-2*pixWSRowStride);
									for(int jj=-6;jj<=6;jj+=3) {
										guint8* currtmpl=currtmp+jj;
										guint8* currtmpWSl=currtmpWS+jj;
										//transparent green
										currtmpl[0] = transparency * currtmpWSl[0];
										currtmpl[1] = transparency * currtmpWSl[1] + (1.f-transparency)*255.f;
										currtmpl[2] = transparency * currtmpWSl[2];
									}
									// three middle row of circle
									for(int ii=-1;ii<=1;ii++) {
										currtmp = &curr[0] + (ii*pixRowStride);
										currtmpWS = &currWs[0] + (ii*pixWSRowStride);
										for(int jj=-9;jj<=9;jj+=3) {
											guint8* currtmpl=currtmp+jj;
											guint8* currtmpWSl=currtmpWS+jj;
											//transparent green
											currtmpl[0] = transparency * currtmpWSl[0];
											currtmpl[1] = transparency * currtmpWSl[1] + (1.f-transparency)*255.f;
											currtmpl[2] = transparency * currtmpWSl[2];
										}
									}
									// second last row of circle
									currtmp = &curr[0] + (2*pixRowStride);
									currtmpWS = &currWs[0] + (2*pixWSRowStride);
									for(int jj=-6;jj<=6;jj+=3) {
										guint8* currtmpl=currtmp+jj;
										guint8* currtmpWSl=currtmpWS+jj;
										//transparent green
										currtmpl[0] = transparency * currtmpWSl[0];
										currtmpl[1] = transparency * currtmpWSl[1] + (1.f-transparency)*255.f;
										currtmpl[2] = transparency * currtmpWSl[2];
									}
									// last row of circle
									currtmp = &curr[0] + (3*pixRowStride);
									currtmpWS = &currWs[0] + (3*pixWSRowStride);
									for(int jj=-3;jj<=3;jj+=3) {
										guint8* currtmpl=currtmp+jj;
										guint8* currtmpWSl=currtmpWS+jj;
										//transparent green
										currtmpl[0] = transparency * currtmpWSl[0];
										currtmpl[1] = transparency * currtmpWSl[1] + (1.f-transparency)*255.f;
										currtmpl[2] = transparency * currtmpWSl[2];
									}
								 }
							curr+=3; currWs+=3;
						}
					}
					free(tmpL);
					free(tmpLsum);
					free(tmpLsumSq);
					free(tmpstdDev2);
					
                } else { // !showFocusMask
				
					const int hlThreshold = options.highlightThreshold;
					const int shThreshold = options.shadowThreshold;
	                const float ShawdowFac = 64 / (options.shadowThreshold+1);
					const float HighlightFac = 64 / (256-options.highlightThreshold);
					const bool showclippedAny = (!showR && !showG && !showB && !showL); // will show clipping if any of RGB chanels is clipped

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,16)
#endif
					for (int i=0; i<bHeight; i++) {
						guint8* curr = pix + i*pixRowStride;
						guint8* currWS = pixWrkSpace + i*pixWSRowStride;

						for (int j=0; j<bWidth; j++) {
                            // we must compare clippings in working space, since the cropPixbuf is in sRGB, with mon profile

                            bool changedHL = false;
                            bool changedSH = false;
							int delta = 0;
                            // for efficiency, pre-calculate currWS_L as it may be needed in both
                            // if (showch) and if (showcs) branches
                            int currWS_L;
                            if (showL && (showch||showcs))
								currWS_L = (int)(0.299f*currWS[0]+0.587f*currWS[1]+0.114f*currWS[2]);
    
                            if (showch) {
                                if ((showclippedAny || showR) && currWS[0]>=hlThreshold ) { delta += 255-currWS[0]; changedHL=true; }
                                if ((showclippedAny || showG) && currWS[1]>=hlThreshold ) { delta += 255-currWS[1]; changedHL=true; }
                                if ((showclippedAny || showB) && currWS[2]>=hlThreshold ) { delta += 255-currWS[2]; changedHL=true; }
                                if (showL && currWS_L>= hlThreshold )                     { delta += 255-currWS_L ; changedHL=true; }

                                if (changedHL) { 
                                    delta *= HighlightFac; 
                                    if (showclippedAny) curr[0]=curr[1]=curr[2]=delta; // indicate clipped highlights in gray
                                    else {curr[0]=255; curr[1]=curr[2]=delta;}         // indicate clipped highlights in red
                                }
                            }
                            if (showcs) {
                                if ((showclippedAny || showR) && currWS[0]<=shThreshold ) { delta += currWS[0]; changedSH=true; }
                                if ((showclippedAny || showG) && currWS[1]<=shThreshold ) { delta += currWS[1]; changedSH=true; }
                                if ((showclippedAny || showB) && currWS[2]<=shThreshold ) { delta += currWS[2]; changedSH=true; }
                                if (showL && currWS_L <=shThreshold )                     { delta += currWS_L ; changedSH=true; }

                                if (changedSH) {
                                    if (showclippedAny) {
                                        delta = 255 - (delta * ShawdowFac);
                                        curr[0]=curr[1]=curr[2]=delta;      // indicate clipped shadows in gray
                                    }
                                    else {
                                        delta *= ShawdowFac;
                                        curr[2]=255; curr[0]=curr[1]=delta; // indicate clipped shadows in blue
                                    }
                                }
                            } //if (showcs)
                            
                            // modulate the preview of channels & L;
                            if (!changedHL && !changedSH && !showclippedAny){          //This condition allows clipping indicators for RGB channels to remain in color
                                if (showR) curr[1]=curr[2]=curr[0]; //Red   channel in grayscale
                                if (showG) curr[0]=curr[2]=curr[1]; //Green channel in grayscale
                                if (showB) curr[0]=curr[1]=curr[2]; //Blue  channel in grayscale
                                if (showL) {                        //Luminosity
                                    // see http://en.wikipedia.org/wiki/HSL_and_HSV#Lightness for more info
                                    //int L = (int)(0.212671*curr[0]+0.715160*curr[1]+0.072169*curr[2]);
                                    int L = (int)(0.299*curr[0]+0.587*curr[1]+0.114*curr[2]); //Lightness - this matches Luminosity mode in Photoshop CS5
                                    curr[0]=curr[1]=curr[2]=L; 
                                }
                            }

                        /*
                            if (showch && (currWS[0]>=options.highlightThreshold || currWS[1]>=options.highlightThreshold || currWS[2]>=options.highlightThreshold))
                                curr[0] = curr[1] = curr[2] = 0;
                            else if (showcs && (currWS[0]<=options.shadowThreshold || currWS[1]<=options.shadowThreshold || currWS[2]<=options.shadowThreshold))
                                curr[0] = curr[1] = curr[2] = 255;
                            //if (showch && ((0.299*curr[0]+0.587*curr[1]+0.114*curr[2])>=options.highlightThreshold))
                            //    curr[0] = curr[1] = curr[2] = 0;
                            //else if (showcs && ((0.299*curr[0]+0.587*curr[1]+0.114*curr[2])<=options.shadowThreshold))
                            //    curr[0] = curr[1] = curr[2] = 255;
                        */

							curr+=3; currWS+=3;
						}
					}
                }
//printf("zoomSteps[cropZoom].zoom=%d\n",zoomSteps[cropZoom].zoom);
                iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), tmp, 0, 0, x+imgX, y+imgY, -1, -1, Gdk::RGB_DITHER_NONE, 0, 0);
            }
            else
                iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), cropHandler.cropPixbuf, 0, 0, x+imgX, y+imgY, -1, -1, Gdk::RGB_DITHER_NONE, 0, 0);

            if (cropHandler.cropParams.enabled) {
                int cropX, cropY;
                cropHandler.getPosition (cropX, cropY);
                drawCrop (cr, x+imgX, y+imgY, imgW, imgH, cropX, cropY, zoomSteps[cropZoom].zoom, cropHandler.cropParams);
            }
            if (observedCropWin)
                drawObservedFrame (cr);

            EditSubscriber *editSubscriber = iarea->getCurrSubscriber();
            rtengine::Crop* crop = static_cast<rtengine::Crop*>(cropHandler.getCrop());
            if (editSubscriber && crop->bufferCreated()) {

                // clip the region
                if (this != iarea->mainCropWindow) {
                    cr->set_line_width (0.);
                    cr->rectangle (x+imgX, y+imgY, imgW, imgH);
                    cr->clip();
                }

                // drawing Subscriber's visible geometry
                const std::vector<Geometry*> visibleGeom = editSubscriber->getVisibleGeometry();
                cr->set_antialias(Cairo::ANTIALIAS_DEFAULT); // ANTIALIAS_SUBPIXEL ?
                cr->set_line_cap(Cairo::LINE_CAP_SQUARE);
                cr->set_line_join(Cairo::LINE_JOIN_ROUND);

                // drawing outer lines
                for (std::vector<Geometry*>::const_iterator i = visibleGeom.begin(); i != visibleGeom.end(); ++i)
                    (*i)->drawOuterGeometry(cr, crop, *this);

                // drawing inner lines
                for (std::vector<Geometry*>::const_iterator i = visibleGeom.begin(); i != visibleGeom.end(); ++i)
                    (*i)->drawInnerGeometry(cr, crop, *this);

                if (this != iarea->mainCropWindow)
                    cr->reset_clip();

                // drawing to the "mouse over" channel
                if (editSubscriber->getEditingType() == ET_OBJECTS) {
                    const std::vector<Geometry*> mouseOverGeom = editSubscriber->getMouseOverGeometry();
                    if (mouseOverGeom.size()) {
                        //printf("ObjectMap (%d x %d)\n", crop->getObjectMap()->get_width(), crop->getObjectMap()->get_height());
                        Cairo::RefPtr<Cairo::Context> crMO = Cairo::Context::create(crop->getObjectMap());
                        crMO->set_antialias(Cairo::ANTIALIAS_NONE);
                        crMO->set_line_cap(Cairo::LINE_CAP_SQUARE);
                        crMO->set_line_join(Cairo::LINE_JOIN_ROUND);
                        crMO->set_operator(Cairo::OPERATOR_SOURCE);

                        // clear the bitmap
                        crMO->set_source_rgba(0., 0., 0., 0.);
                        crMO->rectangle(0., 0., crop->getObjectMap()->get_width(), crop->getObjectMap()->get_height());
                        crMO->set_line_width(0.);
                        crMO->fill();

                        Cairo::RefPtr<Cairo::Context> crMO2;
                        if (crop->getObjectMode() > OM_255) {
                            crMO2 = Cairo::Context::create(crop->getObjectMap2());
                            crMO2->set_antialias(Cairo::ANTIALIAS_NONE);
                            crMO2->set_line_cap(Cairo::LINE_CAP_SQUARE);
                            crMO2->set_line_join(Cairo::LINE_JOIN_ROUND);
                            crMO2->set_operator(Cairo::OPERATOR_SOURCE);

                            // clear the bitmap
                            crMO2->set_source_rgba(0., 0., 0., 0.);
                            crMO2->rectangle(0., 0., crop->getObjectMap2()->get_width(), crop->getObjectMap2()->get_height());
                            crMO2->set_line_width(0.);
                            crMO2->fill();
                        }
                        std::vector<Geometry*>::const_iterator i;
                        int a;
                        for (a=0, i=mouseOverGeom.begin(); i != mouseOverGeom.end(); ++i, ++a) {
                            (*i)->drawToMOChannel(crMO, crMO2, a, crop, *this);
                        }

                        // Debug code: save the "mouse over" image to a new file at each occurrence
                        #if 0
                        {
                        static unsigned int count=0;
                        int w = crop->getObjectMap()->get_width();
                        int h = crop->getObjectMap()->get_height();
                        Glib::RefPtr<Gdk::Pixbuf> img = Gdk::Pixbuf::create(Gdk::COLORSPACE_RGB, false, 8, w, h);
                        guint8 *dst = img->get_pixels();
                        unsigned char *src1 = crop->getObjectMap()->get_data();
                        unsigned char *src2 = crop->getObjectMode() > OM_255 ? crop->getObjectMap()->get_data() : NULL;
                        memcpy(dst, src1, w*h);
                        for (int n=0, n3=0; n<w*h;) {
                            dst[n3++] = src1[n];
                            if (src2)
                                dst[n3++] = src2[n];
                            else
                                dst[n3++] = 0;
                            dst[n3++] = 0;
                            ++n;
                        }
                        img->save(Glib::ustring::compose("mouseOverImage-%1.png", count++), "png");
                        }
                        #endif
                    }
                }
            }
        }
        else {
            // cropHandler.cropPixbuf is null
            int cropX, cropY;
            cropHandler.getPosition (cropX, cropY);
            Glib::RefPtr<Gdk::Pixbuf> rough = iarea->getPreviewHandler()->getRoughImage (cropX, cropY, imgAreaW, imgAreaH, zoomSteps[cropZoom].zoom);
            if (rough) {
                iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), rough, 0, 0, x+imgAreaX+(imgAreaW-rough->get_width())/2, y+imgAreaY+(imgAreaH-rough->get_height())/2, -1, -1, Gdk::RGB_DITHER_NORMAL, 0, 0);  
                if (cropHandler.cropParams.enabled) {
                    drawCrop (cr, x+imgAreaX+(imgAreaW-rough->get_width())/2, y+imgAreaY+(imgAreaH-rough->get_height())/2, rough->get_width(), rough->get_height(), cropX, cropY, zoomSteps[cropZoom].zoom, cropHandler.cropParams);
                }
                if (observedCropWin)
                    drawObservedFrame (cr, rough->get_width(), rough->get_height());
            }
        }
    }

    // if cursor stays above resize area, draw the icon
    if (decorated && (state==SCropWinResize || onResizeArea)) {
        int rw = resizeSurface->get_width ();
        int rh = resizeSurface->get_height ();
        cr->set_source_rgb (0.5,0.5,0.5);
        cr->rectangle (x+w-1.5-rw-1, y+h-1.5-rh-1, rw+1, rh+1);
        cr->stroke_preserve ();
        cr->fill ();
        cr->set_source (resizeSurface, x+w-1.5-rw, y+h-1.5-rh);
        cr->paint ();
        cr->set_source_rgb (0,0,0);
        cr->move_to (x+w-2.5-rw, y+h-1.5);
        cr->line_to (x+w-2.5-rw, y+h-2.5-rh);
        cr->line_to (x+w-1.5, y+h-2.5-rh);
        cr->stroke ();
    }
    if (state==SRotateSelecting) 
        drawStraightenGuide (cr);
    if (state==SNormal && isFlawnOver) {
        EditSubscriber *editSubscriber = iarea->getCurrSubscriber();
        if (iarea->getToolMode () == TMHand && editSubscriber && editSubscriber->getEditingType()==ET_PIPETTE && iarea->object)
            drawUnscaledSpotRectangle (cr, iarea->getPipetteRectSize ());
        else if (iarea->getToolMode () == TMSpotWB)
            drawScaledSpotRectangle (cr, iarea->getSpotWBRectSize ());
    }

    //t2.set ();
//    printf ("etime --> %d, %d\n", t2.etime (t1), t4.etime (t3));
}

// calculate the center of the zoomed in/out preview given a cursor position
void CropWindow::findCenter  (int deltaZoom, int& x, int& y) {
	int cursorX, cursorY;
	screenCoordToImage(x, y, cursorX, cursorY);

	int cropX, cropY, cropW, cropH, skip;
	cropHandler.getWindow (cropX, cropY, cropW, cropH, skip);

	int currCenterX = cropX + cropW/2;
	int currCenterY = cropY + cropH/2;

	int deltaX = currCenterX - cursorX;
	int deltaY = currCenterY - cursorY;

	double factor = zoomSteps[cropZoom].zoom / zoomSteps[cropZoom+deltaZoom].zoom;
	x = cursorX + (int)((double)(deltaX)*factor);
	y = cursorY + (int)((double)(deltaY)*factor);
}

// zoom* is called from the zoomPanel or the scroll wheel in the preview area
void CropWindow::zoomIn (bool toCursor, int cursorX, int cursorY) {

    int x = -1;
    int y = -1;

    if (toCursor) {
        x = cursorX;
        y = cursorY;
    } else {
        if (zoomSteps[cropZoom].zoom <= cropHandler.getFitZoom()) {
            if (cropHandler.cropParams.enabled) {
                x = cropHandler.cropParams.x + cropHandler.cropParams.w / 2;
                y = cropHandler.cropParams.y + cropHandler.cropParams.h / 2;
            } else {
                int fw, fh;
                cropHandler.getFullImageSize(fw, fh);
                x = fw / 2;
                y = fh / 2;
            }
            zoomVersion = exposeVersion;
        } else if (zoomVersion != exposeVersion) {
            screenCoordToImage(xpos+imgX+imgW/2, ypos+imgY+imgH/2, x, y);
            if (cropHandler.cropParams.enabled) {
                // add some gravity towards crop center
                int x1 = cropHandler.cropParams.x + cropHandler.cropParams.w / 2;
                int y1 = cropHandler.cropParams.y + cropHandler.cropParams.h / 2;
                double cropd = sqrt(cropHandler.cropParams.h*cropHandler.cropParams.h+cropHandler.cropParams.w*cropHandler.cropParams.w) * zoomSteps[cropZoom].zoom;
                double imd = sqrt(imgW*imgW+imgH+imgH);
                double d;
                // the more we can see of the crop, the more gravity towards crop center
                if (cropd > imd) {
                    d = 0.8;
                } else if (cropd < imd * 0.5) {
                    d = 0.0;
                } else {
                    d = 1.6 * (cropd - imd * 0.5) / imd;
                }
                x = d * x + (1.0 - d) * x1;
                y = d * y + (1.0 - d) * y1;
            }
            zoomVersion = exposeVersion;
        }
    }

    changeZoom (cropZoom+1, true, x, y);
    fitZoom = false;
}

void CropWindow::zoomOut (bool toCursor, int cursorX, int cursorY) {

    int x = -1;
    int y = -1;

    if (toCursor) {
        x = cursorX;
        y = cursorY;
    }

    zoomVersion = exposeVersion;
    changeZoom (cropZoom-1, true, x, y);
    fitZoom = false;
}

void CropWindow::zoom11 () {

    int x = -1;
    int y = -1;
    if (zoomSteps[cropZoom].zoom <= cropHandler.getFitZoom()) {
        if (cropHandler.cropParams.enabled) {
            x = cropHandler.cropParams.x + cropHandler.cropParams.w / 2;
            y = cropHandler.cropParams.y + cropHandler.cropParams.h / 2;
        } else {
            int fw, fh;
            cropHandler.getFullImageSize(fw, fh);
            x = fw/2;
            y = fh/2;
        }
        zoomVersion = exposeVersion;
    }
    changeZoom (ZOOM11INDEX, true, x, y);
    fitZoom = false;
}

double CropWindow::getZoom () {

    return zoomSteps[cropZoom].zoom;
}

bool CropWindow::isMinZoom () {
    return cropZoom <= 0;
}

bool CropWindow::isMaxZoom () {
    return cropZoom >= MAXZOOMSTEPS;
}

void CropWindow::setZoom (double zoom) {
    int cz = MAXZOOMSTEPS;
    if (zoom < zoomSteps[0].zoom)
        cz = 0;
    else 
        for (int i=0; i<MAXZOOMSTEPS; i++)
            if (zoomSteps[i].zoom <= zoom && zoomSteps[i+1].zoom > zoom) {
                cz = i;
                break;
            }
    changeZoom (cz, false);
}

void CropWindow::zoomFit () {

    double z = cropHandler.getFitZoom ();
    int cz = MAXZOOMSTEPS;
    if (z < zoomSteps[0].zoom)
        cz = 0;
    else 
        for (int i=0; i<MAXZOOMSTEPS; i++)
            if (zoomSteps[i].zoom <= z && zoomSteps[i+1].zoom > z) {
                cz = i;
                break;
            }
    zoomVersion = exposeVersion;
    changeZoom (cz);
    fitZoom = true;
}

void CropWindow::buttonPressed (LWButton* button, int actionCode, void* actionData) {

    if (button==bZoomIn) // zoom in
        zoomIn ();
    else if (button==bZoomOut) // zoom out
        zoomOut ();
    else if (button==bZoom100) // zoom 100
        zoom11 ();
    else if (button==bClose) {// close
        deleted = true;
        iarea->cropWindowClosed (this);
    }
}

void CropWindow::redrawNeeded (LWButton* button) {
    
    iarea->redraw ();
}

void CropWindow::changeZoom  (int zoom, bool notify, int centerx, int centery) {

    if (zoom<0)
        zoom = 0;
    else if (zoom>MAXZOOMSTEPS)
        zoom = MAXZOOMSTEPS;

    if (cropZoom == zoom) {
    	// We are already at the start/end of the zoom range, so we do nothing
    	return;
    }

    	cropZoom = zoom;

    cropLabel = zoomSteps[cropZoom].label;
    cropHandler.setZoom (zoomSteps[cropZoom].czoom, centerx, centery);
    if (notify)
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropZoomChanged (this);
    iarea->redraw ();
}

void CropWindow::screenCoordToCropBuffer (int phyx, int phyy, int& cropx, int& cropy) {

    rtengine::Crop* crop = static_cast<rtengine::Crop*>(cropHandler.getCrop());
    cropx = phyx - xpos - imgX;
    cropy = phyy - ypos - imgY;
    if (zoomSteps[cropZoom].zoom > 1.) {
        cropx = int(double(cropx)/zoomSteps[cropZoom].zoom);
        cropy = int(double(cropy)/zoomSteps[cropZoom].zoom);
    }
    cropx += crop->getLeftBorder();
    cropy += crop->getUpperBorder();
}

void CropWindow::screenCoordToImage (int phyx, int phyy, int& imgx, int& imgy) {

    int cropX, cropY;
    cropHandler.getPosition (cropX, cropY);
    imgx = cropX + (phyx - xpos - imgX)/zoomSteps[cropZoom].zoom;
    imgy = cropY + (phyy - ypos - imgY)/zoomSteps[cropZoom].zoom;
}

void CropWindow::screenCoordToPreview (int phyx, int phyy, int& prevx, int& prevy) {

    prevx = phyx - xpos - imgX;
    prevy = phyy - ypos - imgY;
}

void CropWindow::imageCoordToScreen (int imgx, int imgy, int& phyx, int& phyy) {

    int cropX, cropY;
    cropHandler.getPosition (cropX, cropY);
    phyx = (imgx - cropX)*zoomSteps[cropZoom].zoom + xpos + imgX;
    phyy = (imgy - cropY)*zoomSteps[cropZoom].zoom + ypos + imgY;
    //printf("imgx:%d  /  imgy:%d  /  cropX:%d  /  cropY:%d  /  xpos:%d  /  ypos:%d  /  imgX:%d  /  imgY:%d  /  leftBorder: %d  /  upperBorder:%d  /  phyx:%d  /  phyy:%d\n", imgx, imgy, cropX, cropY, xpos, ypos, imgX, imgY, crop->getLeftBorder(), crop->getUpperBorder(), phyx, phyy);
}

void CropWindow::imageCoordToCropBuffer (int imgx, int imgy, int& phyx, int& phyy) {
    int cropX, cropY;
    rtengine::Crop* crop = static_cast<rtengine::Crop*>(cropHandler.getCrop());
    cropHandler.getPosition (cropX, cropY);
    phyx = (imgx - cropX)*zoomSteps[cropZoom].zoom + /*xpos + imgX +*/ crop->getLeftBorder();
    phyy = (imgy - cropY)*zoomSteps[cropZoom].zoom + /*ypos + imgY +*/ crop->getUpperBorder();
    //printf("imgx:%d  /  imgy:%d  /  cropX:%d  /  cropY:%d  /  xpos:%d  /  ypos:%d  /  imgX:%d  /  imgY:%d  /  leftBorder: %d  /  upperBorder:%d  /  phyx:%d  /  phyy:%d\n", imgx, imgy, cropX, cropY, xpos, ypos, imgX, imgY, crop->getLeftBorder(), crop->getUpperBorder(), phyx, phyy);
}

int CropWindow::scaleValueToImage (int value) {
    return int(double(value)/zoomSteps[cropZoom].zoom);
}

float CropWindow::scaleValueToImage (float value) {
    return float(double(value)/zoomSteps[cropZoom].zoom);
}

double CropWindow::scaleValueToImage (double value) {
    return value/zoomSteps[cropZoom].zoom;
}

int CropWindow::scaleValueToScreen (int value) {
    return int(double(value)*zoomSteps[cropZoom].zoom);
}

float CropWindow::scaleValueToScreen (float value) {
    return float(double(value)*zoomSteps[cropZoom].zoom);
}

double CropWindow::scaleValueToScreen (double value) {
    return value*zoomSteps[cropZoom].zoom;
}

void CropWindow::drawDecoration (Cairo::RefPtr<Cairo::Context> cr) {

    int x = xpos, y = ypos;
    // prepare label
    Glib::RefPtr<Pango::Context> context = iarea->get_pango_context () ;
    Pango::FontDescription fontd = context->get_font_description ();       
    fontd.set_weight (Pango::WEIGHT_BOLD);
    fontd.set_size(8*Pango::SCALE);
    context->set_font_description (fontd);
    Glib::RefPtr<Pango::Layout> cllayout = iarea->create_pango_layout(cropLabel);
    int iw, ih;
    cllayout->get_pixel_size (iw, ih);  

    // draw decoration (border)
    int h = height, w = width;
    cr->set_source_rgb (0,0,0);
    cr->set_line_width (1.0);
    cr->move_to (x+0.5, y+h-0.5);
    cr->line_to (x+0.5, y+0.5);
    cr->line_to (x+w-0.5, y+0.5);
    cr->stroke ();
    cr->set_source_rgb (1,1,1);
    cr->move_to (x+w-0.5, y+0.5);
    cr->line_to (x+w-0.5, y+h-0.5);
    cr->line_to (x+0.5, y+h-0.5);
    cr->stroke ();
    cr->set_source_rgb (0.5,0.5,0.5);
    cr->rectangle (x+1.5, y+1.5+titleHeight, w-3, h-titleHeight-3);
    cr->stroke ();
    cr->set_source_rgb (1,1,1);
    cr->move_to (x+2.5, y+h-2.5);
    cr->line_to (x+2.5, y+titleHeight+2.5);
    cr->line_to (x+w-2.5, y+titleHeight+2.5);
    cr->stroke ();
    cr->set_source_rgb (0,0,0);
    cr->move_to (x+w-2.5, y+titleHeight+2.5);
    cr->line_to (x+w-2.5, y+h-2.5);
    cr->line_to (x+2.5, y+h-2.5);
    cr->stroke ();
    cr->set_source_rgb (0.5,0.5,0.5);
    cr->rectangle (x+1.5, y+1.5, w-3, titleHeight);
    cr->stroke_preserve ();
    cr->fill ();
    
    // draw label
    cr->set_source_rgb (1,1,1);
    cr->move_to (x+6+sideBorderWidth+bZoomIn->getIcon()->get_width()+bZoomOut->getIcon()->get_width()+bZoom100->getIcon()->get_width(), y+upperBorderWidth+(titleHeight-ih)/2);
    cllayout->add_to_cairo_context (cr);
    cr->fill ();

    buttonSet.redraw (cr);
}

void CropWindow::drawStraightenGuide (Cairo::RefPtr<Cairo::Context> cr) {

    if (action_x!=press_x || action_y!=press_y) {
        double arg = (press_x-action_x) / sqrt(double((press_x-action_x)*(press_x-action_x)+(press_y-action_y)*(press_y-action_y)));
        double sol1, sol2;
        double pi = M_PI;
        if (press_y>action_y) {
            sol1 = acos(arg)*180/pi;
            sol2 = -acos(-arg)*180/pi;
        }
        else {
            sol1 = acos(-arg)*180/pi;
            sol2 = -acos(arg)*180/pi;
        }
        if (fabs(sol1)<fabs(sol2))
            rot_deg = sol1;
        else
           rot_deg = sol2;

        if (rot_deg<-45)
           rot_deg = 90.0 + rot_deg;
        else if (rot_deg>45)
           rot_deg = - 90.0 + rot_deg;
    }
    else
        rot_deg = 0;

    Glib::RefPtr<Pango::Context> context = iarea->get_pango_context () ;
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_BOLD);
    fontd.set_size (8*Pango::SCALE);
    context->set_font_description (fontd);
    Glib::RefPtr<Pango::Layout> deglayout = iarea->create_pango_layout(Glib::ustring::compose ("%1 deg", Glib::ustring::format(std::setprecision(2), rot_deg)));

    int x1 = press_x;
    int y1 = press_y;
    int y2 = action_y;
    int x2 = action_x;
/*    if (x1<0) x1 = 0;
    if (y1<0) y1 = 0;
    if (x2<0) x2 = 0;
    if (y2<0) y2 = 0;
    if (x2>=image->getWidth()) x2 = image->getWidth()-1;
    if (y2>=image->getHeight()) y2 = image->getHeight()-1;
    if (x1>=image->getWidth()) x1 = image->getWidth()-1;
    if (y1>=image->getHeight()) y1 = image->getHeight()-1;
*/

    cr->set_line_width (1);
    cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
    cr->move_to (x1+0.5, y1+0.5);
    cr->line_to (x2+0.5, y2+0.5);
    cr->stroke ();
    cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (x1+0.5, y1+0.5);
    cr->line_to (x2+0.5, y2+0.5);
    cr->stroke ();

    if (press_x!=action_x && press_y!=action_y) {
        cr->set_source_rgb (0.0, 0.0, 0.0);
        cr->move_to ((x1+x2)/2+1, (y1+y2)/2+1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1+x2)/2+1, (y1+y2)/2-1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1+x2)/2-1, (y1+y2)/2+1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1+x2)/2+1, (y1+y2)/2+1);
        deglayout->add_to_cairo_context (cr);
        cr->fill ();
        cr->set_source_rgb (1.0, 1.0, 1.0);
        cr->move_to ((x1+x2)/2, (y1+y2)/2);
        deglayout->add_to_cairo_context (cr);
        cr->fill ();
    }
}

void CropWindow::drawScaledSpotRectangle (Cairo::RefPtr<Cairo::Context> cr, int rectSize) {

    int x1 = action_x/zoomSteps[cropZoom].zoom - rectSize;
    int y1 = action_y/zoomSteps[cropZoom].zoom - rectSize;
    int y2 = action_y/zoomSteps[cropZoom].zoom + rectSize;
    int x2 = action_x/zoomSteps[cropZoom].zoom + rectSize;

    cr->set_line_width (1.0);
    cr->rectangle (xpos+imgX-0.5, ypos+imgY-0.5, imgW, imgH);
    cr->clip ();

    cr->set_source_rgb (1.0, 1.0, 1.0);
    cr->rectangle (x1*zoomSteps[cropZoom].zoom-1.5, y1*zoomSteps[cropZoom].zoom-1.5, x2*zoomSteps[cropZoom].zoom-x1*zoomSteps[cropZoom].zoom+2, y2*zoomSteps[cropZoom].zoom-y1*zoomSteps[cropZoom].zoom+2);
    cr->stroke ();
    cr->set_source_rgb (0.0, 0.0, 0.0);
    cr->rectangle (x1*zoomSteps[cropZoom].zoom-0.5, y1*zoomSteps[cropZoom].zoom-0.5, x2*zoomSteps[cropZoom].zoom-x1*zoomSteps[cropZoom].zoom, y2*zoomSteps[cropZoom].zoom-y1*zoomSteps[cropZoom].zoom);
    cr->stroke ();
    
    cr->reset_clip ();
}

void CropWindow::drawUnscaledSpotRectangle (Cairo::RefPtr<Cairo::Context> cr, int rectSize) {

    int x1 = action_x - rectSize;
    int y1 = action_y - rectSize;
    int y2 = action_y + rectSize;
    int x2 = action_x + rectSize;

    cr->set_line_width (1.0);
    cr->rectangle (xpos+imgX-0.5, ypos+imgY-0.5, imgW, imgH);
    cr->clip ();

    cr->set_source_rgb (1.0, 1.0, 1.0);
    cr->rectangle (x1-1.5, y1-1.5, x2-x1+2, y2-y1+2);
    cr->stroke ();
    cr->set_source_rgb (0.0, 0.0, 0.0);
    cr->rectangle (x1-0.5, y1-0.5, x2-x1, y2-y1);
    cr->stroke ();

    cr->reset_clip ();
}

void CropWindow::getObservedFrameArea (int& x, int& y, int& w, int& h, int rw, int rh) {

    int cropX, cropY, cropW, cropH;
    observedCropWin->getCropRectangle (cropX, cropY, cropW, cropH);
    int myCropX, myCropY, myCropW, myCropH;
    getCropRectangle (myCropX, myCropY, myCropW, myCropH);
    
    // translate it to screen coordinates
    if (rw) {
		x = xpos + imgAreaX+(imgAreaW-rw)/2 + (cropX-myCropX)*zoomSteps[cropZoom].zoom;
		y = ypos + imgAreaY+(imgAreaH-rh)/2 + (cropY-myCropY)*zoomSteps[cropZoom].zoom;
    }
    else {
		x = xpos + imgX + (cropX-myCropX)*zoomSteps[cropZoom].zoom;
		y = ypos + imgY + (cropY-myCropY)*zoomSteps[cropZoom].zoom;
    }
    w = cropW * zoomSteps[cropZoom].zoom;
    h = cropH * zoomSteps[cropZoom].zoom;
}

void CropWindow::drawObservedFrame (Cairo::RefPtr<Cairo::Context> cr, int rw, int rh) {

    int x, y, w, h;
    getObservedFrameArea (x, y, w, h, rw, rh);

    cr->set_source_rgb (1.0, 1.0, 1.0);
    cr->set_line_width (4);
    cr->rectangle (x-2, y-2, w+4, h+4);
    cr->stroke ();
    cr->set_source_rgb (1.0, 0.0, 0.0);
    cr->set_line_width (2);
    cr->rectangle (x-2, y-2, w+4, h+4);
    cr->stroke ();
}

void CropWindow::cropImageUpdated () {

    iarea->redraw ();
}

void CropWindow::cropWindowChanged () {

    if (!decorated)
        iarea->syncBeforeAfterViews ();
    iarea->redraw ();
}

void CropWindow::initialImageArrived () {

    for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
        (*i)->initialImageArrived (this);
}


void CropWindow::remoteMove (int deltaX, int deltaY) {

    state = SCropImgMove;
    action_x =  deltaX;
    action_y =  deltaY;
    for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
        (*i)->cropPositionChanged (this);
}

void CropWindow::remoteMoveReady () {

    int cropX, cropY;
    cropHandler.getPosition (cropX, cropY);
    cropHandler.setPosition (cropX + action_x, cropY + action_y);
    cropHandler.getPosition (cropX, cropY);
    state = SNormal;
    for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
        (*i)->cropPositionChanged (this);
}

void CropWindow::delCropWindowListener (CropWindowListener* l) {

    std::list<CropWindowListener*>::iterator i=listeners.begin();
    while (i!=listeners.end()) 
        if (*i==l)
            i = listeners.erase (i);
        else
            i++;
}

EditDataProvider* CropWindow::getImageArea() {
    return iarea;
}
