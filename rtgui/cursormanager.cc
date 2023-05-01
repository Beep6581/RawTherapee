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
#include "rtsurface.h"

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

    auto createCursor = [this] (const Glib::ustring &name, const Gdk::CursorType &fb_cursor,
        const double offX = 0., const double offY = 0.) -> Glib::RefPtr<Gdk::Cursor>
    {
        // Notes:
        // - Gdk Cursor Theme is not supported on some OS (ex : MacOS).
        // Cursor is retrieved from theme thanks to an RTSurface
        // - By default, cursor hotspot is located at middle of surface.
        // Use (offX, offY) between -1 and 0.99 to move cursor hotspot
        auto cursor_surf = RTSurface(name, Gtk::ICON_SIZE_MENU);
        double offXb = std::min(std::max(-1., offX), 0.99); // offX should belong to (-1; 0.99)
        double offYb = std::min(std::max(-1., offY), 0.99); // offY should belong to (-1; 0.99)
        auto cursor = Gdk::Cursor::create(this->display,
            cursor_surf.get(),
            cursor_surf.getWidth() / 2 * (1. + offXb),
            cursor_surf.getHeight() / 2 * (1. + offYb));

        if (!cursor) {
            cursor = Gdk::Cursor::create(this->display, fb_cursor);
        }

        return cursor;
    };

    cAdd        = createCursor("crosshair-small", Gdk::PLUS);
    cAddPicker  = createCursor("color-picker-add-hicontrast", Gdk::PLUS, -0.333, 0.75);
    cCropDraw   = createCursor("crop-point-hicontrast", Gdk::DIAMOND_CROSS);
    cCrosshair  = createCursor("crosshair-hicontrast", Gdk::CROSSHAIR);
    cEmpty      = createCursor("empty", Gdk::BLANK_CURSOR);
    cHandClosed = createCursor("hand-closed-hicontrast", Gdk::HAND1);
    cHandOpen   = createCursor("hand-open-hicontrast", Gdk::HAND2);
    cMoveBL     = createCursor("node-move-sw-ne-hicontrast", Gdk::BOTTOM_LEFT_CORNER);
    cMoveBR     = createCursor("node-move-nw-se-hicontrast", Gdk::BOTTOM_RIGHT_CORNER);
    cMoveL      = createCursor("node-move-x-hicontrast", Gdk::SB_LEFT_ARROW);
    cMoveR      = createCursor("node-move-x-hicontrast", Gdk::SB_RIGHT_ARROW);
    cMoveTL     = createCursor("node-move-nw-se-hicontrast", Gdk::TOP_LEFT_CORNER);
    cMoveTR     = createCursor("node-move-sw-ne-hicontrast", Gdk::TOP_RIGHT_CORNER);
    cMoveX      = createCursor("node-move-x-hicontrast", Gdk::SB_H_DOUBLE_ARROW);
    cMoveXY     = createCursor("node-move-xy-hicontrast", Gdk::FLEUR);
    cMoveY      = createCursor("node-move-y-hicontrast", Gdk::SB_V_DOUBLE_ARROW);
    cRotate     = createCursor("rotate-aroundnode-hicontrast", Gdk::EXCHANGE);
    cWB         = createCursor("color-picker-hicontrast", Gdk::TARGET, -0.333, 0.75);
    cWait       = createCursor("gears", Gdk::CLOCK);

    window = mainWindow;
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

