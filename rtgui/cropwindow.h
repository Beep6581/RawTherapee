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

#include "../rtengine/rtengine.h"
#include <gtkmm.h>
#include "lwbutton.h"
#include "lwbuttonset.h"
#include "editenums.h"
#include "crophandler.h"
#include <list>
#include "cropguilistener.h"
#include "pointermotionlistener.h"
#include "edit.h"

class CropWindow;
class CropWindowListener
{

public:
    virtual ~CropWindowListener() {}
    virtual void cropPositionChanged   (CropWindow*) {}
    virtual void cropWindowSizeChanged (CropWindow*) {}
    virtual void cropZoomChanged       (CropWindow*) {}
    virtual void initialImageArrived   (CropWindow*) {}
};

class ImageArea;
class CropWindow : public LWButtonListener, public CropDisplayHandler, public EditCoordSystem, public ObjectMOBuffer
{

    // state management
    ImgEditState state;                  // current state of user (see enum State)
    int press_x, press_y;               // position of the cursor in the GUI space on button press
    int action_x, action_y;             // parameter that will evolve during a pan or drag action
    int pickedObject;
    int pickModifierKey;
    double rot_deg;
    bool onResizeArea;
    bool deleted;
    bool fitZoomEnabled;
    bool fitZoom;
    bool isLowUpdatePriority;

    // color pickers
    std::vector<LockableColorPicker*> colorPickers;
    LockableColorPicker* hoveredPicker;

    // decoration
    LWButton *bZoomIn, *bZoomOut, *bZoom100, /**bZoomFit,*/ *bClose;
    LWButtonSet buttonSet;
    Glib::ustring cropLabel;
    int backColor;
    bool decorated;
    bool isFlawnOver;

    // crop frame description
    int titleHeight, sideBorderWidth, lowerBorderWidth, upperBorderWidth, sepWidth, minWidth;
    // size & position of the crop relative to the top left corner
    // of the main preview area
    int xpos, ypos, width, height;
    // size & pos of the drawable area relative to the top left corner of the crop
    int imgAreaX, imgAreaY, imgAreaW, imgAreaH;
    // size & pos of the piece of preview image relative to the top left corner of the crop
    int imgX, imgY, imgW, imgH;

    // image handling

    ImageArea* iarea;
    int cropZoom; // *1000
    unsigned int zoomVersion, exposeVersion;

    // crop gui listener
    CropGUIListener* cropgl;
    PointerMotionListener* pmlistener;
    PointerMotionListener* pmhlistener;
    std::list<CropWindowListener*> listeners;

    CropWindow* observedCropWin;  // Pointer to the currently active detail CropWindow

    bool onArea                    (CursorArea a, int x, int y);
    void updateCursor              (int x, int y);
    void drawDecoration            (Cairo::RefPtr<Cairo::Context> cr);
    void drawStraightenGuide       (Cairo::RefPtr<Cairo::Context> cr);
    void drawScaledSpotRectangle   (Cairo::RefPtr<Cairo::Context> cr, int rectSize);
    void drawUnscaledSpotRectangle (Cairo::RefPtr<Cairo::Context> cr, int rectSize);
    void drawObservedFrame         (Cairo::RefPtr<Cairo::Context> cr, int rw = 0, int rh = 0);
    void changeZoom                (int zoom, bool notify = true, int centerx = -1, int centery = -1);
    void updateHoveredPicker       (rtengine::Coord &imgPos);
    void cycleRGB                  ();
    void cycleHSV                  ();

    LockableColorPicker::Validity checkValidity (LockableColorPicker*  picker, const rtengine::Coord &pos);

    // Used by the mainCropWindow only
    void getObservedFrameArea      (int& x, int& y, int& w, int& h, int rw = 0, int rh = 0);

public:
    CropHandler cropHandler;
    CropWindow (ImageArea* parent, bool isLowUpdatePriority_, bool isDetailWindow);
    ~CropWindow ();

    void setDecorated       (bool decorated)
    {
        this->decorated = decorated;
    }
    void setFitZoomEnabled  (bool fze)
    {
        fitZoomEnabled = fze;
    }
    void setObservedCropWin (CropWindow* cw)
    {
        observedCropWin = cw;
    }
    void deleteColorPickers ();

    void screenCoordToCropBuffer (int phyx, int phyy, int& cropx, int& cropy);
    void screenCoordToImage (int phyx, int phyy, int& imgx, int& imgy);
    void screenCoordToCropCanvas (int phyx, int phyy, int& prevx, int& prevy);
    void imageCoordToCropCanvas (int imgx, int imgy, int& phyx, int& phyy);
    void imageCoordToScreen (int imgx, int imgy, int& phyx, int& phyy);
    void imageCoordToCropBuffer (int imgx, int imgy, int& phyx, int& phyy);
    void imageCoordToCropImage (int imgx, int imgy, int& phyx, int& phyy);
    int scaleValueToImage (int value);
    float scaleValueToImage (float value);
    double scaleValueToImage (double value);
    int scaleValueToCanvas (int value);
    float scaleValueToCanvas (float value);
    double scaleValueToCanvas (double value);
    double getZoomFitVal ();
    void setPosition (int x, int y);
    void getPosition (int& x, int& y);
    void setSize     (int w, int h, bool norefresh = false);
    void getSize     (int& w, int& h);
    void enable      ();

    void leaveNotify (GdkEventCrossing* event);
    void flawnOver   (bool isFlawnOver);

    // zoomlistener interface
    void zoomIn      (bool toCursor = false, int cursorX = -1, int cursorY = -1);
    void zoomOut     (bool toCursor = false, int cursorX = -1, int cursorY = -1);
    void zoom11      ();
    void zoomFit     ();
    void zoomFitCrop ();
    double getZoom   ();
    bool isMinZoom   ();
    bool isMaxZoom   ();
    void setZoom     (double zoom);

    bool isInside    (int x, int y);


    void scroll        (int state, GdkScrollDirection direction, int x, int y);
    void buttonPress   (int button, int num, int state, int x, int y);
    void buttonRelease (int button, int num, int state, int x, int y);
    void pointerMoved  (int bstate, int x, int y);

    void expose        (Cairo::RefPtr<Cairo::Context> cr);

    void setEditSubscriber (EditSubscriber* newSubscriber);

    // interface lwbuttonlistener
    void buttonPressed (LWButton* button, int actionCode, void* actionData);
    void redrawNeeded  (LWButton* button);

    // crop handling
    void getCropRectangle      (int& x, int& y, int& w, int& h);
    void getCropPosition       (int& x, int& y);
    void setCropPosition       (int x, int y, bool update = true);
    void centerCrop            (bool update = true);
    void getCropSize           (int& w, int& h);
    void getCropAnchorPosition (int& w, int& h);
    void setCropAnchorPosition (int& w, int& h);

    // listeners
    void setCropGUIListener       (CropGUIListener* cgl);
    void setPointerMotionListener (PointerMotionListener* pml);
    PointerMotionListener* getPointerMotionListener ();
    void setPointerMotionHListener (PointerMotionListener* pml);

    // crop window listeners
    void addCropWindowListener (CropWindowListener* l);
    void delCropWindowListener (CropWindowListener* l);

    // crophandlerlistener interface
    void cropImageUpdated ();
    void cropWindowChanged ();
    void initialImageArrived ();
    void setDisplayPosition (int x, int y);

    void remoteMove      (int deltaX, int deltaY);
    void remoteMoveReady ();

    ImageArea* getImageArea();
};

#endif
