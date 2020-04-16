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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "cursormanager.h"

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

    Glib::RefPtr<Gdk::Pixbuf> add           = RTImage::createPixbufFromFile("crosshair-small.png");
    Glib::RefPtr<Gdk::Pixbuf> colPick       = RTImage::createPixbufFromFile("color-picker-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> colPickAdd    = RTImage::createPixbufFromFile("color-picker-add-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> cropDraw      = RTImage::createPixbufFromFile("crop-point-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> crosshair     = RTImage::createPixbufFromFile("crosshair-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> empty         = RTImage::createPixbufFromFile("empty.png");
    Glib::RefPtr<Gdk::Pixbuf> handClosed    = RTImage::createPixbufFromFile("hand-closed-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> handOpen      = RTImage::createPixbufFromFile("hand-open-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveBL        = RTImage::createPixbufFromFile("node-move-sw-ne-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveBR        = RTImage::createPixbufFromFile("node-move-nw-se-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveL         = RTImage::createPixbufFromFile("node-move-x-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveR         = RTImage::createPixbufFromFile("node-move-x-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveTL        = RTImage::createPixbufFromFile("node-move-nw-se-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveTR        = RTImage::createPixbufFromFile("node-move-sw-ne-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveX         = RTImage::createPixbufFromFile("node-move-x-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveXY        = RTImage::createPixbufFromFile("node-move-xy-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> moveY         = RTImage::createPixbufFromFile("node-move-y-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> rotate        = RTImage::createPixbufFromFile("rotate-aroundnode-hicontrast.png");
    Glib::RefPtr<Gdk::Pixbuf> wait          = RTImage::createPixbufFromFile("gears.png"); // Currently unused, create *-hicontrast once used.

    double s = RTScalable::getTweakedDPI() / RTScalable::baseDPI;  // RTScalable::getDPI() might be preferable, however it imply a lot of work to support this option

    cAdd = add                  ? Gdk::Cursor::create(display, add, (int)(8.*s), (int)(8.*s))           : Gdk::Cursor::create(display, Gdk::PLUS);
    cAddPicker = colPickAdd     ? Gdk::Cursor::create(display, colPickAdd, (int)(4.*s), (int)(21.*s))   : Gdk::Cursor::create(display, Gdk::PLUS);
    cCropDraw = cropDraw        ? Gdk::Cursor::create(display, cropDraw, (int)(3.*s), (int)(3.*s))      : Gdk::Cursor::create(display, Gdk::DIAMOND_CROSS);
    cCrosshair = crosshair      ? Gdk::Cursor::create(display, crosshair, (int)(12.*s), (int)(12.*s))   : Gdk::Cursor::create(display, Gdk::CROSSHAIR);
    cEmpty = empty              ? Gdk::Cursor::create(display, empty, 12, 12) /* PNG: do not scale */   : Gdk::Cursor::create(display, Gdk::BLANK_CURSOR);
    cHandClosed = handClosed    ? Gdk::Cursor::create(display, handClosed, (int)(12.*s), (int)(12.*s))  : Gdk::Cursor::create(display, Gdk::HAND1);
    cHandOpen = handOpen        ? Gdk::Cursor::create(display, handOpen, (int)(12.*s), (int)(12.*s))    : Gdk::Cursor::create(display, Gdk::HAND2);
    cMoveBL = moveBL            ? Gdk::Cursor::create(display, moveBL, (int)(12.*s), (int)(12.*s))      : Gdk::Cursor::create(display, Gdk::BOTTOM_LEFT_CORNER);
    cMoveBR = moveBR            ? Gdk::Cursor::create(display, moveBR, (int)(12.*s), (int)(12.*s))      : Gdk::Cursor::create(display, Gdk::BOTTOM_RIGHT_CORNER);
    cMoveL = moveL              ? Gdk::Cursor::create(display, moveL, (int)(12.*s), (int)(12.*s))       : Gdk::Cursor::create(display, Gdk::SB_LEFT_ARROW);
    cMoveR = moveR              ? Gdk::Cursor::create(display, moveR, (int)(12.*s), (int)(12.*s))       : Gdk::Cursor::create(display, Gdk::SB_RIGHT_ARROW);
    cMoveTL = moveTL            ? Gdk::Cursor::create(display, moveTL, (int)(12.*s), (int)(12.*s))      : Gdk::Cursor::create(display, Gdk::TOP_LEFT_CORNER);
    cMoveTR = moveTR            ? Gdk::Cursor::create(display, moveTR, (int)(12.*s), (int)(12.*s))      : Gdk::Cursor::create(display, Gdk::TOP_RIGHT_CORNER);
    cMoveX = moveX              ? Gdk::Cursor::create(display, moveX, (int)(12.*s), (int)(12.*s))       : Gdk::Cursor::create(display, Gdk::SB_H_DOUBLE_ARROW);
    cMoveXY = moveXY            ? Gdk::Cursor::create(display, moveXY, (int)(12.*s), (int)(12.*s))      : Gdk::Cursor::create(display, Gdk::FLEUR);
    cMoveY = moveY              ? Gdk::Cursor::create(display, moveY, (int)(12.*s), (int)(12.*s))       : Gdk::Cursor::create(display, Gdk::SB_V_DOUBLE_ARROW);
    cRotate = rotate            ? Gdk::Cursor::create(display, rotate, (int)(12.*s), (int)(12.*s))      : Gdk::Cursor::create(display, Gdk::EXCHANGE);
    cWB = colPick               ? Gdk::Cursor::create(display, colPick, (int)(4.*s), (int)(21.*s))      : Gdk::Cursor::create(display, Gdk::TARGET);
    cWait = wait                ? Gdk::Cursor::create(display, wait, (int)(12.*s), (int)(12.*s))        : Gdk::Cursor::create(display, Gdk::CLOCK);

    window = mainWindow;
}

void CursorManager::cleanup()
{
    cAdd.reset();
    cAddPicker.reset();
    cCropDraw.reset();
    cCrosshair.reset();
    cHandClosed.reset();
    cHandOpen.reset();
    cEmpty.reset();
    cMoveBL.reset();
    cMoveBR.reset();
    cMoveL.reset();
    cMoveR.reset();
    cMoveTL.reset();
    cMoveTR.reset();
    cMoveX.reset();
    cMoveY.reset();
    cMoveXY.reset();
    cRotate.reset();
    cWB.reset();
    cWait.reset();
}

/* Set the cursor of the given window */
void CursorManager::setCursor (Glib::RefPtr<Gdk::Window> window, CursorShape shape)
{
    switch (shape)
    {
        case CursorShape::CSAddColPicker:
            window->set_cursor(cAddPicker);
            break;
        case CursorShape::CSArrow:
            window->set_cursor(); // set_cursor without any arguments to select system default
            break;
        case CursorShape::CSCropSelect:
            window->set_cursor(cCropDraw);
            break;
        case CursorShape::CSCrosshair:
            window->set_cursor(cCrosshair);
            break;
        case CursorShape::CSEmpty:
            window->set_cursor(cEmpty);
            break;
        case CursorShape::CSHandClosed:
            window->set_cursor(cHandClosed);
            break;
        case CursorShape::CSHandOpen:
            window->set_cursor(cHandOpen);
            break;
        case CursorShape::CSMove:
            window->set_cursor(cHandClosed);
            break;
        case CursorShape::CSMove1DH:
            window->set_cursor(cMoveX);
            break;
        case CursorShape::CSMove1DV:
            window->set_cursor(cMoveY);
            break;
        case CursorShape::CSMove2D:
            window->set_cursor(cMoveXY);
            break;
        case CursorShape::CSMoveLeft:
            window->set_cursor(cMoveL);
            break;
        case CursorShape::CSMoveRight:
            window->set_cursor(cMoveR);
            break;
        case CursorShape::CSMoveRotate:
            window->set_cursor(cRotate);
            break;
        case CursorShape::CSPlus:
            window->set_cursor(cAdd);
            break;
        case CursorShape::CSResizeBottomLeft:
            window->set_cursor(cMoveBL);
            break;
        case CursorShape::CSResizeBottomRight:
            window->set_cursor(cMoveBR);
            break;
        case CursorShape::CSResizeDiagonal:
            window->set_cursor(cMoveXY);
            break;
        case CursorShape::CSResizeHeight:
            window->set_cursor(cMoveY);
            break;
        case CursorShape::CSResizeTopLeft:
            window->set_cursor(cMoveTL);
            break;
        case CursorShape::CSResizeTopRight:
            window->set_cursor(cMoveTR);
            break;
        case CursorShape::CSResizeWidth:
            window->set_cursor(cMoveX);
            break;
        case CursorShape::CSSpotWB:
            window->set_cursor(cWB);
            break;
        case CursorShape::CSStraighten:
            window->set_cursor(cRotate);
            break;
        case CursorShape::CSUndefined:
            break;
        case CursorShape::CSWait:
            window->set_cursor(cWait);
            break;
        default:
            window->set_cursor(cCrosshair);
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

