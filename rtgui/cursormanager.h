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

enum CursorShape {CSArrow, CSOpenHand, CSClosedHand, CSMove, CSResizeWidth, CSResizeHeight, CSResizeDiagonal, CSSpotWB, CSCropSelect, CSStraighten};

class CursorManager {

    protected:
	    Gdk::Cursor* cResizeWidth;
        Gdk::Cursor* cResizeHeight;
        Gdk::Cursor* cResizeDiag;
        Gdk::Cursor* cCropMove;
        Gdk::Cursor* cCropMoving;
        Gdk::Cursor* cNormal;
        Gdk::Cursor* cCropSelection;
        Gdk::Cursor* cHand;
        Gdk::Cursor* cClosedHand;
        Gdk::Cursor* cWB;
        Glib::RefPtr<Gdk::Window> mainWindow;

    public:
        void init        (Glib::RefPtr<Gdk::Window> mainWin);
        void setCursor   (Glib::RefPtr<Gdk::Window> window, CursorShape shape);
};

extern CursorManager cursorManager;

#endif

