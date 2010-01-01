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
#include <cropwindow.h>
#include <options.h>
#include <iomanip>
#include <guiutils.h>
#include <mytime.h>
#include <imagearea.h>
#include <cursormanager.h>

struct ZoomStep {   
    Glib::ustring label;
    double zoom;
    int czoom;
};

ZoomStep zoomSteps[] = {"10%",   0.1,     10,
                        "12.5%", 0.125,   8,
                        "16.6%", 1.0/6.0, 6,
                        "20%",   0.2,     5,
                        "25%",   0.25,    4,
                        "33%",   1.0/3.0, 3,
                        "50%",   0.5,     2,
                        "100%",  1.0,     1000,
                        "200%",  2.0,     2000,
                        "300%",  3.0,     3000,
                        "400%",  4.0,     4000,
                        "500%",  5.0,     5000,
                        "600%",  6.0,     6000,
                        "700%",  7.0,     7000,
                        "800%",  8.0,     8000};
#define MAXZOOMSTEPS 14
#define ZOOM11INDEX  7

CropWindow::CropWindow (ImageArea* parent, rtengine::StagedImageProcessor* ipc_) 
    : iarea(parent), cropgl(NULL), xpos(30), ypos(30), imgX(0), imgY(0), imgW(1), imgH(1),
    titleHeight(30), sideBorderWidth(3), upperBorderWidth(1), lowerBorderWidth(3), sepWidth(2),
    cropZoom(ZOOM11INDEX), deleted(false), onResizeArea(false), fitZoom(false),
    fitZoomEnabled(true), decorated(true), backColor(0), observedCropWin(NULL),
	pmlistener(NULL)	{

    Glib::RefPtr<Pango::Context> context = parent->get_pango_context () ;
    Pango::FontDescription fontd = context->get_font_description ();       
    fontd.set_weight (Pango::WEIGHT_BOLD);
    fontd.set_size(8*Pango::SCALE);
    context->set_font_description (fontd);
    cropLabel = "100%";
    Glib::RefPtr<Pango::Layout> cllayout = parent->create_pango_layout("1000%");
 
    int iw, ih;
    cllayout->get_pixel_size (iw, ih);
    
    titleHeight = ih;
       
    resizeSurface = Cairo::ImageSurface::create_from_png (argv0+"/images/resize.png");
    bZoomIn  = new LWButton (Cairo::ImageSurface::create_from_png (argv0+"/images/gtk-zoom-in.png"),  0, NULL, LWButton::Left, LWButton::Center, "Zoom In");
    bZoomOut = new LWButton (Cairo::ImageSurface::create_from_png (argv0+"/images/gtk-zoom-out.png"), 1, NULL, LWButton::Left, LWButton::Center, "Zoom Out");
    bZoom100 = new LWButton (Cairo::ImageSurface::create_from_png (argv0+"/images/gtk-zoom-100.png"), 2, NULL, LWButton::Left, LWButton::Center, "Zoom 100/%");
    bZoomFit = new LWButton (Cairo::ImageSurface::create_from_png (argv0+"/images/gtk-zoom-fit.png"), 3, NULL, LWButton::Left, LWButton::Center, "Zoom Fit");
    bClose   = new LWButton (Cairo::ImageSurface::create_from_png (argv0+"/images/gtk-close.png"),    4, NULL, LWButton::Right, LWButton::Center, "Close");
    
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
    
    setSize (300, 300);
    cropHandler.newImage (ipc_);
    cropHandler.setPosition (0,0);
    cropHandler.setEnabled (true);
    cropHandler.setCropHandlerListener (this);
    state = SNormal;
}

CropWindow::~CropWindow () {

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
    iarea->redraw ();
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
        
void CropWindow::buttonPress (int button, int type, int bstate, int x, int y) {

    iarea->grabFocus (this);
    if (button==1 && type==GDK_2BUTTON_PRESS && onArea (CropImage, x, y) && (state==SNormal || state==SCropImgMove)) {
        if (fitZoomEnabled) {
            if (fitZoom) {
                translateCoord (x, y, action_x, action_y);
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
    else if (button==1 && type==GDK_2BUTTON_PRESS && onArea (CropBorder, x, y)) {
        backColor = (backColor+1) % 3;
    }
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
        if (onArea (CropTop, x, y)) {
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
        else if (onArea (CropObserved, x, y)) {
            state = SObservedMove;
            press_x = x;
            press_y = y;
        }
        else if ((bstate & GDK_SHIFT_MASK) && onArea (CropInside, x, y)) {
            state = SCropMove;
            press_x = x;
            press_y = y;
            action_x = cropHandler.cropParams.x;
            action_y = cropHandler.cropParams.y;
        }
        else if (iarea->getToolMode () == TMHand) {
            state = SCropImgMove;
            action_x = 0;
            action_y = 0;
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
            translateCoord (x, y, action_x, action_y);
            iarea->spotWBSelected (action_x, action_y);
        }
        else if (iarea->getToolMode () == TMCropSelect && cropgl) {
            state = SCropSelecting;
            translateCoord (x, y, press_x, press_y);
            cropHandler.cropParams.x = press_x;
            cropHandler.cropParams.y = press_y;
            cropHandler.cropParams.w = cropHandler.cropParams.h = 1;
            cropgl->cropInit (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        }
    }
    if (button==3) {
        state = SNormal;
        iarea->setToolHand ();
    }
    iarea->redraw ();
    updateCursor (x, y);
}

void CropWindow::buttonRelease (int button, int num, int bstate, int x, int y) {

    if (state==SCropWinResize) {
        setSize (press_x + x - action_x, press_y + y - action_y);
        state = SNormal;
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropWindowSizeChanged (this);
    }
    else if (state==SCropImgMove) {
        int cropX, cropY;
        cropHandler.getPosition (cropX, cropY);
        cropHandler.setPosition (cropX + action_x, cropY + action_y);
        cropHandler.getPosition (cropX, cropY);
        state = SNormal;
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropPositionChanged (this);
    }
    else if (state==SRotateSelecting) {
        iarea->straightenReady (rot_deg);
        iarea->setToolHand ();
    }
    else if (state==SObservedMove) {
        observedCropWin->remoteMoveReady ();
        state = SNormal;
    }
    if (cropgl && (state==SCropSelecting || state==SResizeH1 || state==SResizeH2 || state==SResizeW1 || state==SResizeW2 || state==SCropMove)) {
        cropgl->cropManipReady ();
        iarea->setToolHand ();
    }

    if (decorated)
        buttonSet.releaseNotify (x, y);
    if (deleted)
        return;

    state = SNormal;
    iarea->grabFocus (NULL);
    iarea->redraw ();
    updateCursor (x, y);
}

void CropWindow::pointerMoved (int x, int y) {

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
        action_x =  (press_x - x) / zoomSteps[cropZoom].zoom;
        action_y =  (press_y - y) / zoomSteps[cropZoom].zoom;
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
    else if (state==SCropMove && cropgl) { 
        cropHandler.cropParams.x = action_x + (x-press_x) / zoomSteps[cropZoom].zoom;
        cropHandler.cropParams.y = action_y + (y-press_y) / zoomSteps[cropZoom].zoom;
        cropgl->cropMoved (cropHandler.cropParams.x, cropHandler.cropParams.y, cropHandler.cropParams.w, cropHandler.cropParams.h);
        iarea->redraw ();
    }
    else if (state==SCropSelecting && cropgl) { 
        translateCoord (x, y, action_x, action_y);
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
		translateCoord (x, y, mx, my);
		if (!onArea (CropImage, x, y) || !cropHandler.cropPixbuf) 
			pmlistener->pointerMoved (false, mx, my, -1, -1, -1);
		else {
			cropHandler.cimg.lock ();
			int vx = x - xpos - imgX;
			int vy = y - ypos - imgY;
			guint8* pix = cropHandler.cropPixbuf->get_pixels() + vy*cropHandler.cropPixbuf->get_rowstride() + vx*3;
			if (vx < cropHandler.cropPixbuf->get_width() && vy < cropHandler.cropPixbuf->get_height())
				pmlistener->pointerMoved (true, mx, my, pix[0], pix[1], pix[2]);
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
        case CropTop:
            translateCoord (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                x1>cropHandler.cropParams.x+CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                y1>cropHandler.cropParams.y-CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+CROPRESIZEBORDER;
        case CropBottom:
            translateCoord (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                x1>cropHandler.cropParams.x+CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                y1>cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1+CROPRESIZEBORDER;
        case CropLeft:
            translateCoord (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>cropHandler.cropParams.y+CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                x1>cropHandler.cropParams.x-CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+CROPRESIZEBORDER;
        case CropRight:
            translateCoord (x, y, x1, y1);
            return cropHandler.cropParams.enabled && 
                y1>cropHandler.cropParams.y+CROPRESIZEBORDER && 
                y1<cropHandler.cropParams.y+cropHandler.cropParams.h-1-CROPRESIZEBORDER && 
                x1>cropHandler.cropParams.x+cropHandler.cropParams.w-1-CROPRESIZEBORDER && 
                x1<cropHandler.cropParams.x+cropHandler.cropParams.w-1+CROPRESIZEBORDER;
        case CropInside:
            translateCoord (x, y, x1, y1);
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
    
    ToolMode tm = iarea->getToolMode ();
    
    if (state==SNormal) {
        if (onArea (CropWinButtons, x, y)) 
            cursorManager.setCursor (iarea->get_window(), CSArrow);
        else if (onArea (CropToolBar, x, y)) 
            cursorManager.setCursor (iarea->get_window(), CSMove);
        else if (onArea (CropResize, x, y)) 
            cursorManager.setCursor (iarea->get_window(), CSResizeDiagonal);
        else if (tm==TMHand && (onArea (CropTop, x, y) || onArea (CropBottom, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeHeight);
        else if (tm==TMHand && (onArea (CropLeft, x, y) || onArea (CropRight, x, y))) 
            cursorManager.setCursor (iarea->get_window(), CSResizeWidth);
        else if (onArea (CropImage, x, y)) { 
            if (tm==TMHand) {
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
    else if (state==SCropWinResize)
        cursorManager.setCursor (iarea->get_window(), CSResizeDiagonal);
}

void CropWindow::expose (Cairo::RefPtr<Cairo::Context> cr) {

    MyTime t1, t2, t3, t4;
    
    t1.set ();

    if (decorated)
        drawDecoration (cr);

    int x = xpos, y = ypos, h = height, w = width;
 
    // draw border
    if (backColor==0) {
        Gdk::Color cback = iarea->get_style()->get_bg(Gtk::STATE_NORMAL);
        cr->set_source_rgb (cback.get_red_p(), cback.get_green_p(), cback.get_blue_p());
    }
    else if (backColor==1)
        cr->set_source_rgb (0,0,0);
    else if (backColor==2)
        cr->set_source_rgb (1,1,1);
        
    cr->rectangle (x+imgAreaX+0.5, y+imgAreaY+0.5, imgAreaW, imgAreaH);
    cr->stroke_preserve ();
    cr->fill ();

    cropHandler.cimg.lock ();
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
    }
    else {         
        if (cropHandler.cropPixbuf) {
            imgW = cropHandler.cropPixbuf->get_width ();
            imgH = cropHandler.cropPixbuf->get_height ();
            imgX = imgAreaX + (imgAreaW-imgW)/2;
            imgY = imgAreaY + (imgAreaH-imgH)/2;
// PERFORMANCE BOTTLENECK STARTS HERE
    t3.set ();
			bool showcs = iarea->indClippedPanel->showClippedShadows();
			bool showch = iarea->indClippedPanel->showClippedHighlights();
			if (showcs || showch) {
				Glib::RefPtr<Gdk::Pixbuf> tmp = cropHandler.cropPixbuf->copy ();
				guint8* pix = tmp->get_pixels();
				for (int i=0; i<tmp->get_height(); i++)
					for (int j=0; j<tmp->get_width(); j++) {
						guint8* curr = pix + i*tmp->get_rowstride () + j*3;
						if (showch && (curr[0]>=options.highlightThreshold || curr[1]>=options.highlightThreshold || curr[2]>=options.highlightThreshold))
							curr[0] = curr[1] = curr[2] = 0;
						else if (showcs && (curr[0]<=options.shadowThreshold || curr[1]<=options.shadowThreshold || curr[2]<=options.shadowThreshold))
							curr[0] = curr[1] = curr[2] = 255;
					}
				iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), tmp, 0, 0, x+imgX, y+imgY, -1, -1, Gdk::RGB_DITHER_NONE, 0, 0);  
			}
			else
				iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), cropHandler.cropPixbuf, 0, 0, x+imgX, y+imgY, -1, -1, Gdk::RGB_DITHER_NONE, 0, 0);  
    t4.set ();
// END OF BOTTLENECK
            if (cropHandler.cropParams.enabled) {
                int cropX, cropY;
                cropHandler.getPosition (cropX, cropY);
                drawCrop (cr, x+imgX, y+imgY, imgW, imgH, cropX, cropY, zoomSteps[cropZoom].zoom, cropHandler.cropParams);
            }
        }
        else {
            int cropX, cropY;
            cropHandler.getPosition (cropX, cropY);
            Glib::RefPtr<Gdk::Pixbuf> rough = iarea->getPreviewHandler()->getRoughImage (cropX, cropY, imgAreaW, imgAreaH, zoomSteps[cropZoom].zoom);
            if (rough) {
                iarea->get_window()->draw_pixbuf (iarea->get_style()->get_base_gc(Gtk::STATE_NORMAL), rough, 0, 0, x+imgAreaX+(imgAreaW-rough->get_width())/2, y+imgAreaY+(imgAreaH-rough->get_height())/2, -1, -1, Gdk::RGB_DITHER_NORMAL, 0, 0);  
                if (cropHandler.cropParams.enabled) {
                    int cropX, cropY;
                    cropHandler.getPosition (cropX, cropY);
                    drawCrop (cr, x+imgX, y+imgY, imgW, imgH, cropX, cropY, zoomSteps[cropZoom].zoom, cropHandler.cropParams);
                }
            }
        }
    }

    if (observedCropWin)
        drawObservedFrame (cr);

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
    if (state==SNormal && iarea->getToolMode () == TMSpotWB)
        drawSpotWBRectangle (cr);

    t2.set ();
    cropHandler.cimg.unlock ();
//    printf ("etime --> %d, %d\n", t2.etime (t1), t4.etime (t3));
}

void CropWindow::zoomIn () {

    changeZoom (cropZoom+1);
    fitZoom = false;
}

void CropWindow::zoomOut () {

    changeZoom (cropZoom-1);
    fitZoom = false;
}

void CropWindow::zoom11 () {

    changeZoom (ZOOM11INDEX);
    fitZoom = false;
}

double CropWindow::getZoom () {

    return zoomSteps[cropZoom].zoom;
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

    cropZoom = zoom;
    if (cropZoom<0)
        cropZoom = 0;
    else if (cropZoom>MAXZOOMSTEPS)
        cropZoom = MAXZOOMSTEPS;

    cropLabel = zoomSteps[cropZoom].label;
    cropHandler.setZoom (zoomSteps[cropZoom].czoom, centerx, centery);
    if (notify)
        for (std::list<CropWindowListener*>::iterator i=listeners.begin(); i!=listeners.end(); i++)
            (*i)->cropZoomChanged (this);
    iarea->redraw ();
}

void CropWindow::translateCoord (int phyx, int phyy, int& imgx, int& imgy) {

    int cropX, cropY;
    cropHandler.getPosition (cropX, cropY);
    imgx = cropX + (phyx - xpos - imgX)/zoomSteps[cropZoom].zoom;
    imgy = cropY + (phyy - ypos - imgY)/zoomSteps[cropZoom].zoom;
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
        double arg = (press_x-action_x) / sqrt((press_x-action_x)*(press_x-action_x)+(press_y-action_y)*(press_y-action_y));
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

    cr->set_line_width (1.5);
    cr->set_source_rgb (1.0, 1.0, 1.0);
    cr->move_to (x1, y1);
    cr->line_to (x2, y2);
    cr->stroke ();
    cr->set_source_rgb (0.0, 0.0, 0.0);
    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (x1, y1);
    cr->line_to (x2, y2);
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

void CropWindow::drawSpotWBRectangle (Cairo::RefPtr<Cairo::Context> cr) {

    int rectsize = iarea->getSpotWBRectSize ();
    int x1 = action_x/zoomSteps[cropZoom].zoom - rectsize;
    int y1 = action_y/zoomSteps[cropZoom].zoom - rectsize;
    int y2 = action_y/zoomSteps[cropZoom].zoom + rectsize;
    int x2 = action_x/zoomSteps[cropZoom].zoom + rectsize;

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

void CropWindow::getObservedFrameArea (int& x, int& y, int& w, int& h) {
    
    int cropX, cropY, cropW, cropH;
    observedCropWin->getCropRectangle (cropX, cropY, cropW, cropH);
    int myCropX, myCropY, myCropW, myCropH;
    getCropRectangle (myCropX, myCropY, myCropW, myCropH);
    
    // translate it to screen coordinates
    x = xpos + imgX + (cropX-myCropX)*zoomSteps[cropZoom].zoom;
    y = ypos + imgY + (cropY-myCropY)*zoomSteps[cropZoom].zoom;
    w = cropW * zoomSteps[cropZoom].zoom;
    h = cropH * zoomSteps[cropZoom].zoom;
}

void CropWindow::drawObservedFrame (Cairo::RefPtr<Cairo::Context> cr) {
    
    int x, y, w, h;
    getObservedFrameArea (x, y, w, h);

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
        iarea->updateScrollbars ();
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
