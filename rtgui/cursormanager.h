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
#ifndef _CURSORMANAGER_
#define _CURSORMANAGER_

#include <gtkmm.h>

enum CursorShape {
    CSArrow, CSOpenHand, CSClosedHand, CSMove, CSMoveLeft,
    CSMoveRight, CSResizeWidth, CSResizeHeight, CSResizeDiagonal,
    CSResizeTopLeft, CSResizeTopRight, CSResizeBottomLeft, CSResizeBottomRight,
    CSMove2D, CSMove1DH, CSMove1DV, CSMoveRotate,
    CSSpotWB, CSCropSelect, CSStraighten, CSPlus, CSWait, CSEmpty
};

class CursorManager {

    protected:
        Gdk::Cursor* cResizeWidth;
        Gdk::Cursor* cResizeHeight;
        Gdk::Cursor* cResizeDiag;
        Gdk::Cursor* cResizeTopLeft;
        Gdk::Cursor* cResizeTopRight;
        Gdk::Cursor* cResizeBottomLeft;
        Gdk::Cursor* cResizeBottomRight;
        Gdk::Cursor* cCropMove;
        Gdk::Cursor* cCropMoving;
        Gdk::Cursor* cLeftTanMove;
        Gdk::Cursor* cRightTanMove;
        Gdk::Cursor* cNormal;
        Gdk::Cursor* cCropSelection;
        Gdk::Cursor* cAdd;
        Gdk::Cursor* cWait;
        Gdk::Cursor* cHand;
        Gdk::Cursor* cClosedHand;
        Gdk::Cursor* cWB;
        Gdk::Cursor* cHidden;
        Gdk::Cursor* cMove2D;
        Gdk::Cursor* cMove1DH;
        Gdk::Cursor* cMove1DV;
        Gdk::Cursor* cMoveRotate;
        Glib::RefPtr<Gdk::Window> mainWindow;

    public:
        void init        (Glib::RefPtr<Gdk::Window> mainWin);
        void setCursor   (Glib::RefPtr<Gdk::Window> window, CursorShape shape);
        void setCursor   (CursorShape shape);
};

extern CursorManager cursorManager;

#endif

