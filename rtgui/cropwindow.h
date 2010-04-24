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
#ifndef _CROPWINDOW_
#define _CROPWINDOW_

#include <rtengine.h>
#include <gtkmm.h>
#include <lwbutton.h>
#include <lwbuttonset.h>
#include <editenums.h>
#include <crophandler.h>
#include <list>
#include <cropguilistener.h>
#include <pointermotionlistener.h>

class CropWindow;
class CropWindowListener {
    
    public:
        virtual void cropPositionChanged   (CropWindow*) {}
        virtual void cropWindowSizeChanged (CropWindow*) {}
        virtual void cropZoomChanged       (CropWindow*) {}
        virtual void initialImageArrived   (CropWindow*) {}
};

class ImageArea;
class CropWindow : public LWButtonListener, public CropHandlerListener {

        // state management
        ImgEditState state;					// current state of user (see enum State)
        int action_x, action_y, press_x, press_y;
        double rot_deg;
        bool onResizeArea;
        bool deleted;
        bool fitZoomEnabled;
        bool fitZoom;

        // decoration
        Cairo::RefPtr<Cairo::ImageSurface> resizeSurface;
        LWButton *bZoomIn, *bZoomOut, *bZoom100, *bZoomFit, *bClose;
        LWButtonSet buttonSet;
        Glib::ustring cropLabel;
        int backColor;
        bool decorated;
        
        // sizes, positions
        int titleHeight, sideBorderWidth, lowerBorderWidth, upperBorderWidth, sepWidth, minWidth;
        int imgAreaX, imgAreaY, imgAreaW, imgAreaH;
        int imgX, imgY, imgW, imgH;
        int xpos, ypos, width, height;
        
        // image handling
        CropHandler cropHandler;
        ImageArea* iarea;
        int cropZoom; // *1000
        
        // crop gui listener
        CropGUIListener* cropgl;
		PointerMotionListener* pmlistener;
        std::list<CropWindowListener*> listeners;
        
        CropWindow* observedCropWin;

        bool        onArea              (CursorArea a, int x, int y);
        void        updateCursor        (int x, int y);
        void        drawDecoration      (Cairo::RefPtr<Cairo::Context> cr);
        void        drawStraightenGuide (Cairo::RefPtr<Cairo::Context> cr);
        void        drawSpotWBRectangle (Cairo::RefPtr<Cairo::Context> cr);
        void        drawObservedFrame   (Cairo::RefPtr<Cairo::Context> cr);
        void        translateCoord      (int phyx, int phyy, int& imgx, int& imgy);
        void        changeZoom          (int zoom, bool notify=true, int centerx=-1, int centery=-1); 
        void        getObservedFrameArea(int& x, int& y, int& w, int& h);
    
    public:
        CropWindow (ImageArea* parent, rtengine::StagedImageProcessor* ipc_);
        ~CropWindow ();
        
        void setDecorated       (bool decorated)    { this->decorated = decorated; }
        void setFitZoomEnabled  (bool fze)          { fitZoomEnabled = fze; }
        void setObservedCropWin (CropWindow* cw)    { observedCropWin = cw; }
               
        void setPosition (int x, int y);
        void getPosition (int& x, int& y);
        void setSize     (int w, int h, bool norefresh=false);
        void getSize     (int& w, int& h);
        
        // zoomlistener interface
        void zoomIn      ();
        void zoomOut     ();
        void zoom11      ();
        void zoomFit     ();
        double getZoom   ();
        void   setZoom   (double zoom);

        bool isInside    (int x, int y);

        void buttonPress   (int button, int num, int state, int x, int y);
        void buttonRelease (int button, int num, int state, int x, int y);
        void pointerMoved  (int x, int y);

        void expose        (Cairo::RefPtr<Cairo::Context> cr);

        // interface lwbuttonlistener
        void buttonPressed (LWButton* button, int actionCode, void* actionData);
        void redrawNeeded  (LWButton* button);
        
        // crop handling
        void getCropRectangle (int& x, int& y, int& w, int& h);
        void getCropPosition  (int& x, int& y);
        void setCropPosition  (int x, int y);
        void getCropSize      (int& w, int& h);
        
        // listeners
        void setCropGUIListener 	  (CropGUIListener* cgl) { cropgl = cgl; }
		void setPointerMotionListener (PointerMotionListener* pml) { pmlistener = pml; }
        
        // crop window listeners
        void addCropWindowListener (CropWindowListener* l) { listeners.push_back (l); }
        void delCropWindowListener (CropWindowListener* l);
        
        // crophandlerlistener interface
        void cropImageUpdated ();
        void cropWindowChanged ();
        void initialImageArrived ();
        
        void remoteMove      (int deltaX, int deltaY);
        void remoteMoveReady ();
};

#endif
