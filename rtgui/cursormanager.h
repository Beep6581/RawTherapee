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
    CSUndefined, CSArrow, CSOpenHand, CSClosedHand, CSMove, CSMoveLeft,
    CSMoveRight, CSResizeWidth, CSResizeHeight, CSResizeDiagonal,
    CSResizeTopLeft, CSResizeTopRight, CSResizeBottomLeft, CSResizeBottomRight,
    CSMove2D, CSMove1DH, CSMove1DV, CSMoveRotate,
    CSSpotWB, CSCropSelect, CSStraighten, CSPlus, CSWait, CSEmpty
};

class CursorManager
{

private:
    Glib::RefPtr<Gdk::Cursor> cResizeWidth;
    Glib::RefPtr<Gdk::Cursor> cResizeHeight;
    Glib::RefPtr<Gdk::Cursor> cResizeDiag;
    Glib::RefPtr<Gdk::Cursor> cResizeTopLeft;
    Glib::RefPtr<Gdk::Cursor> cResizeTopRight;
    Glib::RefPtr<Gdk::Cursor> cResizeBottomLeft;
    Glib::RefPtr<Gdk::Cursor> cResizeBottomRight;
    Glib::RefPtr<Gdk::Cursor> cCropMove;
    Glib::RefPtr<Gdk::Cursor> cCropMoving;
    Glib::RefPtr<Gdk::Cursor> cLeftTanMove;
    Glib::RefPtr<Gdk::Cursor> cRightTanMove;
    Glib::RefPtr<Gdk::Cursor> cNormal;
    Glib::RefPtr<Gdk::Cursor> cCropSelection;
    Glib::RefPtr<Gdk::Cursor> cAdd;
    Glib::RefPtr<Gdk::Cursor> cWait;
    Glib::RefPtr<Gdk::Cursor> cHand;
    Glib::RefPtr<Gdk::Cursor> cClosedHand;
    Glib::RefPtr<Gdk::Cursor> cWB;
    Glib::RefPtr<Gdk::Cursor> cHidden;
    Glib::RefPtr<Gdk::Cursor> cMove2D;
    Glib::RefPtr<Gdk::Cursor> cMove1DH;
    Glib::RefPtr<Gdk::Cursor> cMove1DV;
    Glib::RefPtr<Gdk::Cursor> cMoveRotate;

    Glib::RefPtr<Gdk::Display> display;
    Glib::RefPtr<Gdk::Window> window;

    void setCursor (CursorShape shape);
    void setCursor (Glib::RefPtr<Gdk::Window> window, CursorShape shape);

public:
    void init                         (Glib::RefPtr<Gdk::Window> mainWindow);
    static void setWidgetCursor       (Glib::RefPtr<Gdk::Window> window, CursorShape shape);
    static void setCursorOfMainWindow (Glib::RefPtr<Gdk::Window> window, CursorShape shape);
};

extern CursorManager mainWindowCursorManager;
extern CursorManager editWindowCursorManager;

#endif

