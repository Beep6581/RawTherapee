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
#include "rtimage.h"

CursorManager mainWindowCursorManager;
CursorManager editWindowCursorManager;

void CursorManager::init (Glib::RefPtr<Gdk::Window> mainWindow)
{

    display = Gdk::Display::get_default ();
#ifndef NDEBUG

    if (!display) {
        printf("Error: no default display!\n");
    }

#endif

    cResizeWidth = Gdk::Cursor::create (display, Gdk::SB_H_DOUBLE_ARROW);
    cResizeHeight = Gdk::Cursor::create (display, Gdk::SB_V_DOUBLE_ARROW);
    cResizeDiag = Gdk::Cursor::create (display, Gdk::BOTTOM_RIGHT_CORNER);
    cResizeTopLeft = Gdk::Cursor::create (display, Gdk::TOP_LEFT_CORNER);
    cResizeTopRight = Gdk::Cursor::create (display, Gdk::TOP_RIGHT_CORNER);
    cResizeBottomLeft = Gdk::Cursor::create (display, Gdk::BOTTOM_LEFT_CORNER);
    cResizeBottomRight = Gdk::Cursor::create (display, Gdk::BOTTOM_RIGHT_CORNER);
    cCropMove = Gdk::Cursor::create (display, Gdk::FLEUR);
    cCropMoving = Gdk::Cursor::create (display, Gdk::HAND2);
    cCropSelection = Gdk::Cursor::create (display, Gdk::CROSSHAIR);
    cLeftTanMove = Gdk::Cursor::create (display, Gdk::SB_LEFT_ARROW);
    cRightTanMove = Gdk::Cursor::create (display, Gdk::SB_RIGHT_ARROW);
    cAdd = Gdk::Cursor::create (display, Gdk::PLUS);
    cWait = Gdk::Cursor::create (display, Gdk::CLOCK);

    Glib::RefPtr<Gdk::Pixbuf> hand = RTImage::createFromFile ("cross.png");
    Glib::RefPtr<Gdk::Pixbuf> close_hand = RTImage::createFromFile ("closedhand.png");
    Glib::RefPtr<Gdk::Pixbuf> wbpick = RTImage::createFromFile ("gtk-color-picker-small.png");
    Glib::RefPtr<Gdk::Pixbuf> empty = RTImage::createFromFile ("empty.png");
    Glib::RefPtr<Gdk::Pixbuf> move2D = RTImage::createFromFile ("move-2D.png");
    Glib::RefPtr<Gdk::Pixbuf> move1DH = RTImage::createFromFile ("move-1D-h.png");
    Glib::RefPtr<Gdk::Pixbuf> move1DV = RTImage::createFromFile ("move-1D-v.png");
    Glib::RefPtr<Gdk::Pixbuf> moveRotate = RTImage::createFromFile ("move-rotate.png");

    cHand = hand ? Gdk::Cursor::create (cAdd->get_display(), hand, 10, 10) : Gdk::Cursor::create (cAdd->get_display(), Gdk::HAND2);
    cClosedHand = close_hand ? Gdk::Cursor::create (cAdd->get_display(), close_hand, 10, 10) : Gdk::Cursor::create (cAdd->get_display(), Gdk::HAND2);
    cWB = wbpick ? Gdk::Cursor::create (cAdd->get_display(), wbpick, 1, 12) : Gdk::Cursor::create (cAdd->get_display(), Gdk::ARROW);
    cHidden = empty ? Gdk::Cursor::create (cAdd->get_display(), empty, 12, 12) : Gdk::Cursor::create (cAdd->get_display(), Gdk::FLEUR);
    cMove2D = move2D ?  Gdk::Cursor::create (cAdd->get_display(), move2D, 11, 11) : Gdk::Cursor::create (cAdd->get_display(), Gdk::FLEUR);
    cMove1DH = move1DH ?  Gdk::Cursor::create (cAdd->get_display(), move1DH, 11, 11) : Gdk::Cursor::create (cAdd->get_display(), Gdk::FLEUR);
    cMove1DV = move1DV ?  Gdk::Cursor::create (cAdd->get_display(), move1DV, 11, 11) : Gdk::Cursor::create (cAdd->get_display(), Gdk::FLEUR);
    cMoveRotate = moveRotate ?  Gdk::Cursor::create (cAdd->get_display(), moveRotate, 11, 11) : Gdk::Cursor::create (cAdd->get_display(), Gdk::CIRCLE);

    window = mainWindow;
}

/* Set the cursor of the given window */
void CursorManager::setCursor (Glib::RefPtr<Gdk::Window> window, CursorShape shape)
{

    if (shape == CSArrow)
        // set_cursor without any arguments to select system default
    {
        window->set_cursor ();
    } else if (shape == CSOpenHand) {
        window->set_cursor (cHand);
    } else if (shape == CSClosedHand) {
        window->set_cursor (cClosedHand);
    } else if (shape == CSMove) {
        window->set_cursor (cCropMove);
    } else if (shape == CSResizeWidth) {
        window->set_cursor (cResizeWidth);
    } else if (shape == CSResizeHeight) {
        window->set_cursor (cResizeHeight);
    } else if (shape == CSResizeDiagonal) {
        window->set_cursor (cResizeDiag);
    } else if (shape == CSResizeTopLeft) {
        window->set_cursor (cResizeTopLeft);
    } else if (shape == CSResizeTopRight) {
        window->set_cursor (cResizeTopRight);
    } else if (shape == CSResizeBottomLeft) {
        window->set_cursor (cResizeBottomLeft);
    } else if (shape == CSResizeBottomRight) {
        window->set_cursor (cResizeBottomRight);
    } else if (shape == CSMove2D) {
        window->set_cursor (cMove2D);
    } else if (shape == CSMove1DH) {
        window->set_cursor (cMove1DH);
    } else if (shape == CSMove1DV) {
        window->set_cursor (cMove1DV);
    } else if (shape == CSMoveRotate) {
        window->set_cursor (cMoveRotate);
    } else if (shape == CSSpotWB) {
        window->set_cursor (cWB);
    } else if (shape == CSCropSelect) {
        window->set_cursor (cHand);
    } else if (shape == CSMoveLeft) {
        window->set_cursor (cLeftTanMove);
    } else if (shape == CSMoveRight) {
        window->set_cursor (cRightTanMove);
    } else if (shape == CSStraighten) {
        window->set_cursor (cHand);
    } else if (shape == CSWait) {
        window->set_cursor (cWait);
    } else if (shape == CSPlus) {
        window->set_cursor (cAdd);
    } else if (shape == CSEmpty) {
        window->set_cursor (cHidden);
    }
}

void CursorManager::setWidgetCursor (Glib::RefPtr<Gdk::Window> window, CursorShape shape)
{
    if (window->get_display() == mainWindowCursorManager.display) {
        mainWindowCursorManager.setCursor(window, shape);
    } else if (window->get_display() == editWindowCursorManager.display) {
        editWindowCursorManager.setCursor(window, shape);
    }

#ifndef NDEBUG
    else {
        printf("CursorManager::setWidgetCursor  /  Error: Display not found!\n");
    }

#endif
}

void CursorManager::setCursorOfMainWindow (Glib::RefPtr<Gdk::Window> window, CursorShape shape)
{
    if (window->get_display() == mainWindowCursorManager.display) {
        mainWindowCursorManager.setCursor(shape);
    } else if (window->get_display() == editWindowCursorManager.display) {
        editWindowCursorManager.setCursor(shape);
    }

#ifndef NDEBUG
    else {
        printf("CursorManager::setCursorOfMainWindow  /  Error: Display not found!\n");
    }

#endif
}

/* Set the cursor of the main window */
void CursorManager::setCursor (CursorShape shape)
{
    setCursor (window, shape);
}

