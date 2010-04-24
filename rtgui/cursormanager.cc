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
#include <cursormanager.h>
#include <options.h>
#include <safegtk.h>

CursorManager cursorManager;

void CursorManager::init (Glib::RefPtr<Gdk::Window> mainWin) {

    cResizeWidth = new Gdk::Cursor (Gdk::SB_H_DOUBLE_ARROW);
    cResizeHeight = new Gdk::Cursor (Gdk::SB_V_DOUBLE_ARROW);
    cResizeDiag = new Gdk::Cursor (Gdk::BOTTOM_RIGHT_CORNER);
    cCropMove = new Gdk::Cursor (Gdk::FLEUR);
    cCropMoving = new Gdk::Cursor (Gdk::HAND2);
    cCropSelection = new Gdk::Cursor (Gdk::CROSSHAIR);
#ifdef _WIN32
    cNormal = new Gdk::Cursor (Gdk::LAST_CURSOR);
#else
    cNormal = new Gdk::Cursor (Gdk::ARROW);
#endif
		Glib::RefPtr<Gdk::Pixbuf> hand = safe_create_from_file(argv0+"/images/openhand22.png");
		Glib::RefPtr<Gdk::Pixbuf> close_hand = safe_create_from_file(argv0+"/images/closedhand22.png");
		Glib::RefPtr<Gdk::Pixbuf> wbpick = safe_create_from_file(argv0+"/images/wbpicker16.png");
		cHand = hand ? new Gdk::Cursor (cNormal->get_display(), hand, 10, 10) : new Gdk::Cursor (Gdk::HAND2);
    cClosedHand = close_hand ? new Gdk::Cursor (cNormal->get_display(), close_hand, 10, 10) : new Gdk::Cursor (Gdk::HAND2);
    cWB = wbpick ? new Gdk::Cursor (cNormal->get_display(), wbpick, 1, 12) : new Gdk::Cursor (Gdk::ARROW);

    mainWindow = mainWin;
}

void CursorManager::setCursor (Glib::RefPtr<Gdk::Window> window, CursorShape shape) {

    if (shape==CSArrow)
        window->set_cursor (*cNormal);
    else if (shape==CSOpenHand)
        window->set_cursor (*cHand);
    else if (shape==CSClosedHand)
        window->set_cursor (*cClosedHand);
    else if (shape==CSMove)
        window->set_cursor (*cCropMove);
    else if (shape==CSResizeWidth)
        window->set_cursor (*cResizeWidth);
    else if (shape==CSResizeHeight)
        window->set_cursor (*cResizeHeight);
    else if (shape==CSResizeDiagonal)
        window->set_cursor (*cResizeDiag);
    else if (shape==CSSpotWB)
        window->set_cursor (*cWB);
    else if (shape==CSCropSelect)
        window->set_cursor (*cCropSelection);
    else if (shape==CSStraighten)
        window->set_cursor (*cCropSelection);
}


