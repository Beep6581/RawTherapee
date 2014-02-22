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
#include "cursormanager.h"
#include "options.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"

CursorManager cursorManager;

void CursorManager::init (Glib::RefPtr<Gdk::Window> mainWin) {

    cResizeWidth = new Gdk::Cursor (Gdk::SB_H_DOUBLE_ARROW);
    cResizeHeight = new Gdk::Cursor (Gdk::SB_V_DOUBLE_ARROW);
    cResizeDiag = new Gdk::Cursor (Gdk::BOTTOM_RIGHT_CORNER);
    cResizeTopLeft = new Gdk::Cursor (Gdk::TOP_LEFT_CORNER);
    cResizeTopRight = new Gdk::Cursor (Gdk::TOP_RIGHT_CORNER);
    cResizeBottomLeft = new Gdk::Cursor (Gdk::BOTTOM_LEFT_CORNER);
    cResizeBottomRight = new Gdk::Cursor (Gdk::BOTTOM_RIGHT_CORNER);
    cCropMove = new Gdk::Cursor (Gdk::FLEUR);
    cCropMoving = new Gdk::Cursor (Gdk::HAND2);
    cCropSelection = new Gdk::Cursor (Gdk::CROSSHAIR);
    cLeftTanMove = new Gdk::Cursor (Gdk::SB_LEFT_ARROW);
    cRightTanMove = new Gdk::Cursor (Gdk::SB_RIGHT_ARROW);
    cAdd = new Gdk::Cursor (Gdk::PLUS);
    cWait = new Gdk::Cursor (Gdk::CLOCK);

    Glib::RefPtr<Gdk::Pixbuf> hand = safe_create_from_file("cross.png");
    Glib::RefPtr<Gdk::Pixbuf> close_hand = safe_create_from_file("closedhand.png");
    Glib::RefPtr<Gdk::Pixbuf> wbpick = safe_create_from_file("gtk-color-picker-small.png");
    Glib::RefPtr<Gdk::Pixbuf> empty = safe_create_from_file("empty.png");
    Glib::RefPtr<Gdk::Pixbuf> move2D = safe_create_from_file("move-2D.png");
    Glib::RefPtr<Gdk::Pixbuf> move1DH = safe_create_from_file("move-1D-h.png");
    Glib::RefPtr<Gdk::Pixbuf> move1DV = safe_create_from_file("move-1D-v.png");
    Glib::RefPtr<Gdk::Pixbuf> moveRotate = safe_create_from_file("move-rotate.png");

    cHand = hand ? new Gdk::Cursor (cAdd->get_display(), hand, 10, 10) : new Gdk::Cursor (Gdk::HAND2);
    cClosedHand = close_hand ? new Gdk::Cursor (cAdd->get_display(), close_hand, 10, 10) : new Gdk::Cursor (Gdk::HAND2);
    cWB = wbpick ? new Gdk::Cursor (cAdd->get_display(), wbpick, 1, 12) : new Gdk::Cursor (Gdk::ARROW);
    cHidden = empty ? new Gdk::Cursor (cAdd->get_display(), empty, 12, 12) : new Gdk::Cursor (Gdk::FLEUR);
    cMove2D = move2D ?  new Gdk::Cursor (cAdd->get_display(), move2D, 11, 11) : new Gdk::Cursor (Gdk::FLEUR);
    cMove1DH = move1DH ?  new Gdk::Cursor (cAdd->get_display(), move1DH, 11, 11) : new Gdk::Cursor (Gdk::FLEUR);
    cMove1DV = move1DV ?  new Gdk::Cursor (cAdd->get_display(), move1DV, 11, 11) : new Gdk::Cursor (Gdk::FLEUR);
    cMoveRotate = moveRotate ?  new Gdk::Cursor (cAdd->get_display(), moveRotate, 11, 11) : new Gdk::Cursor (Gdk::CIRCLE);

    mainWindow = mainWin;
}

/* Set the cursor of the given window */
void CursorManager::setCursor (Glib::RefPtr<Gdk::Window> window, CursorShape shape) {

    if (shape==CSArrow)
       // set_cursor without any arguments to select system default
        window->set_cursor ();
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
    else if (shape==CSResizeTopLeft)
        window->set_cursor (*cResizeTopLeft);
    else if (shape==CSResizeTopRight)
        window->set_cursor (*cResizeTopRight);
    else if (shape==CSResizeBottomLeft)
        window->set_cursor (*cResizeBottomLeft);
    else if (shape==CSResizeBottomRight)
        window->set_cursor (*cResizeBottomRight);
    else if (shape==CSMove2D)
        window->set_cursor (*cMove2D);
    else if (shape==CSMove1DH)
        window->set_cursor (*cMove1DH);
    else if (shape==CSMove1DV)
        window->set_cursor (*cMove1DV);
    else if (shape==CSMoveRotate)
        window->set_cursor (*cMoveRotate);
    else if (shape==CSSpotWB)
        window->set_cursor (*cWB);
    else if (shape==CSCropSelect)
        window->set_cursor (*cCropSelection);
    else if (shape==CSMoveLeft)
        window->set_cursor (*cLeftTanMove);
    else if (shape==CSMoveRight)
        window->set_cursor (*cRightTanMove);
    else if (shape==CSStraighten)
        window->set_cursor (*cCropSelection);
    else if (shape==CSWait)
        window->set_cursor (*cWait);
    else if (shape==CSPlus)
        window->set_cursor (*cAdd);
    else if (shape==CSEmpty)
        window->set_cursor (*cHidden);
}

/* Set the cursor of the main window */
void CursorManager::setCursor (CursorShape shape) {
	setCursor(mainWindow, shape);
}

