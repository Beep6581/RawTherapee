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
#include "filebrowserentry.h"

#include <cstring>
#include <iomanip>

#include "cropguilistener.h"
#include "cursormanager.h"
#include "guiutils.h"
#include "inspector.h"
#include "rtimage.h"
#include "threadutils.h"
#include "thumbbrowserbase.h"
#include "thumbnail.h"
#include "toolbar.h"

#include "../rtengine/procparams.h"

#define CROPRESIZEBORDER 4

//extern Glib::Threads::Thread* mainThread;

Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::editedIcon;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::recentlySavedIcon;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::enqueuedIcon;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::hdr;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::ps;

FileBrowserEntry::FileBrowserEntry (Thumbnail* thm, const Glib::ustring& fname)
    : ThumbBrowserEntryBase (fname, thm), wasInside(false), iatlistener(nullptr), press_x(0), press_y(0), action_x(0), action_y(0), rot_deg(0.0), landscape(true), cropParams(new rtengine::procparams::CropParams), cropgl(nullptr), state(SNormal), crop_custom_ratio(0.f)
{
    feih = new FileBrowserEntryIdleHelper;
    feih->fbentry = this;
    feih->destroyed = false;
    feih->pending = 0;

    italicstyle = thumbnail->getType() != FT_Raw;
    datetimeline = thumbnail->getDateTimeString ();
    exifline = thumbnail->getExifString ();

    scale = 1;

    thumbnail->addThumbnailListener (this);
}

FileBrowserEntry::~FileBrowserEntry ()
{
    idle_register.destroy();

    // so jobs arriving now do nothing
    if (feih->pending) {
        feih->destroyed = true;
    } else {
        delete feih;
        feih = nullptr;
    }

    thumbImageUpdater->removeJobs (this);

    if (thumbnail) {
        thumbnail->removeThumbnailListener (this);
        thumbnail->decreaseRef ();
    }
}

void FileBrowserEntry::init ()
{
    editedIcon = RTImage::createPixbufFromFile ("tick-small.png");
    recentlySavedIcon = RTImage::createPixbufFromFile ("save-small.png");
    enqueuedIcon = RTImage::createPixbufFromFile ("gears-small.png");
    hdr = RTImage::createPixbufFromFile ("filetype-hdr.png");
    ps = RTImage::createPixbufFromFile ("filetype-ps.png");
}

void FileBrowserEntry::refreshThumbnailImage ()
{

    if (!thumbnail) {
        return;
    }

    thumbImageUpdater->add (this, &updatepriority, false, this);
}

void FileBrowserEntry::refreshQuickThumbnailImage ()
{

    if (!thumbnail) {
        return;
    }

    // Only make a (slow) processed preview if the picture has been edited at all
    bool upgrade_to_processed = (!options.internalThumbIfUntouched || thumbnail->isPParamsValid());
    thumbImageUpdater->add(this, &updatepriority, upgrade_to_processed, this);
}

void FileBrowserEntry::calcThumbnailSize ()
{
    if (thumbnail) {
        int ow = prew, oh = preh;
        thumbnail->getThumbnailSize(prew, preh);
        if (ow != prew || oh != preh || preview.size() != static_cast<std::size_t>(prew * preh * 3)) {
            preview.clear();
        }
    }
}

std::vector<Glib::RefPtr<Gdk::Pixbuf>> FileBrowserEntry::getIconsOnImageArea ()
{
    if (!thumbnail) {
        return {};
    }

    std::vector<Glib::RefPtr<Gdk::Pixbuf>> ret;

    if (thumbnail->hasProcParams() && editedIcon) {
        ret.push_back(editedIcon);
    }

    if (thumbnail->isRecentlySaved() && recentlySavedIcon) {
        ret.push_back(recentlySavedIcon);
    }

    if (thumbnail->isEnqueued () && enqueuedIcon) {
        ret.push_back(enqueuedIcon);
    }

    return ret;
}

std::vector<Glib::RefPtr<Gdk::Pixbuf>> FileBrowserEntry::getSpecificityIconsOnImageArea ()
{
    if (!thumbnail) {
        return {};
    }

    std::vector<Glib::RefPtr<Gdk::Pixbuf>> ret;

    if (thumbnail->isHDR() && hdr) {
        ret.push_back (hdr);
    }

    if (thumbnail->isPixelShift() && ps) {
        ret.push_back (ps);
    }

    return ret;
}

void FileBrowserEntry::customBackBufferUpdate (Cairo::RefPtr<Cairo::Context> c)
{
    if(scale != 1.0 && cropParams->enabled) { // somewhere in pipeline customBackBufferUpdate is called when scale == 1.0, which is nonsense for a thumb
        if (state == SCropSelecting || state == SResizeH1 || state == SResizeH2 || state == SResizeW1 || state == SResizeW2 || state == SResizeTL || state == SResizeTR || state == SResizeBL || state == SResizeBR || state == SCropMove) {
            drawCrop (c, prex, prey, prew, preh, 0, 0, scale, *cropParams, true, false);
        } else {
            rtengine::procparams::CropParams cparams = thumbnail->getProcParams().crop;
            switch (options.cropGuides) {
            case Options::CROP_GUIDE_NONE:
                cparams.guide = rtengine::procparams::CropParams::Guide::NONE;
                break;
            case Options::CROP_GUIDE_FRAME:
                cparams.guide = rtengine::procparams::CropParams::Guide::FRAME;
                break;
            default:
                break;
            }

            if (cparams.enabled && !thumbnail->isQuick()) { // Quick thumb have arbitrary sizes, so don't apply the crop
                drawCrop (c, prex, prey, prew, preh, 0, 0, scale, cparams, true, false);
            }
        }
    }
}

void FileBrowserEntry::getIconSize (int& w, int& h) const
{

    w = editedIcon->get_width ();
    h = editedIcon->get_height ();
}

FileThumbnailButtonSet* FileBrowserEntry::getThumbButtonSet ()
{

    return (static_cast<FileThumbnailButtonSet*>(buttonSet));
}

void FileBrowserEntry::procParamsChanged (Thumbnail* thm, int whoChangedIt)
{

    if ( thumbnail->isQuick() ) {
        refreshQuickThumbnailImage ();
    } else {
        refreshThumbnailImage ();
    }
}

void FileBrowserEntry::updateImage(rtengine::IImage8* img, double scale, const rtengine::procparams::CropParams& cropParams)
{
    if (!feih) {
        return;
    }
    redrawRequests++;
    feih->pending++;

    idle_register.add(
        [this, img, scale, cropParams]() -> bool
        {
            if (feih->destroyed) {
                if (feih->pending == 1) {
                    delete feih;
                } else {
                    --feih->pending;
                }

                delete img;
                return false;
            }

            feih->fbentry->_updateImage(img, scale, cropParams);
            --feih->pending;

            return false;
        },
        G_PRIORITY_LOW
    );
}

void FileBrowserEntry::_updateImage(rtengine::IImage8* img, double s, const rtengine::procparams::CropParams& cropParams)
{
    MYWRITERLOCK(l, lockRW);

    redrawRequests--;
    scale = s;
    *this->cropParams = cropParams;

    bool newLandscape = img->getWidth() > img->getHeight();
    bool rotated = false;

    if (preh == img->getHeight()) {
        prew = img->getWidth ();

        // Check if image has been rotated since last time
        rotated = !preview.empty() && newLandscape != landscape;

        preview.resize(prew * preh * 3);
        std::copy(img->getData(), img->getData() + preview.size(), preview.begin());
        {
            GThreadLock lock;
            updateBackBuffer ();
        }
    }

    landscape = newLandscape;

    delete img;

    if (parent) {
        if (rotated) {
            parent->thumbRearrangementNeeded();
        } else if (redrawRequests == 0) {
            parent->redrawNeeded(this);
        }
    }
}

bool FileBrowserEntry::motionNotify (int x, int y)
{

    const bool b = ThumbBrowserEntryBase::motionNotify(x, y);

    const int ix = x - startx - ofsX;
    const int iy = y - starty - ofsY;

    Inspector* inspector = parent->getInspector();

    if (inspector && inspector->isActive() && (!parent->isInTabMode() || options.inspectorWindow)) {
        const rtengine::Coord2D coord(getPosInImgSpace(x, y));

        if (coord.x != -1.) {
            if (!wasInside) {
                inspector->switchImage(filename);
                wasInside = true;
            }
            inspector->mouseMove(coord, 0);
        } else {
            wasInside = false;
        }
    }

    if (inside(x, y)) {
        updateCursor(ix, iy);
    }

    if (state == SRotateSelecting) {
        action_x = x;
        action_y = y;
        parent->redrawNeeded (this);
    } else if (state == SResizeH1 && cropgl) {
        int oy = cropParams->y;
        cropParams->y = action_y + (y - press_y) / scale;
        cropParams->h += oy - cropParams->y;
        cropgl->cropHeight1Resized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeH2 && cropgl) {
        cropParams->h = action_y + (y - press_y) / scale;
        cropgl->cropHeight2Resized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeW1 && cropgl) {
        int ox = cropParams->x;
        cropParams->x = action_x + (x - press_x) / scale;
        cropParams->w += ox - cropParams->x;
        cropgl->cropWidth1Resized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeW2 && cropgl) {
        cropParams->w = action_x + (x - press_x) / scale;
        cropgl->cropWidth2Resized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeTL && cropgl) {
        int ox = cropParams->x;
        cropParams->x = action_x + (x - press_x) / scale;
        cropParams->w += ox - cropParams->x;
        int oy = cropParams->y;
        cropParams->y = action_y + (y - press_y) / scale;
        cropParams->h += oy - cropParams->y;
        cropgl->cropTopLeftResized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeTR && cropgl) {
        cropParams->w = action_x + (x - press_x) / scale;
        int oy = cropParams->y;
        cropParams->y = action_y + (y - press_y) / scale;
        cropParams->h += oy - cropParams->y;
        cropgl->cropTopRightResized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeBL && cropgl) {
        int ox = cropParams->x;
        cropParams->x = action_x + (x - press_x) / scale;
        cropParams->w += ox - cropParams->x;
        cropParams->h = action_y + (y - press_y) / scale;
        cropgl->cropBottomLeftResized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SResizeBR && cropgl) {
        cropParams->w = action_x + (x - press_x) / scale;
        cropParams->h = action_y + (y - press_y) / scale;
        cropgl->cropBottomRightResized (cropParams->x, cropParams->y, cropParams->w, cropParams->h, crop_custom_ratio);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SCropMove && cropgl) {
        cropParams->x = action_x + (x - press_x) / scale;
        cropParams->y = action_y + (y - press_y) / scale;
        cropgl->cropMoved (cropParams->x, cropParams->y, cropParams->w, cropParams->h);
        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    } else if (state == SCropSelecting && cropgl) {
        int cx1 = press_x, cy1 = press_y;
        int cx2 = (ix - prex) / scale, cy2 = (iy - prey) / scale;
        cropgl->cropResized (cx1, cy1, cx2, cy2);

        if (cx2 > cx1) {
            cropParams->x = cx1;
            cropParams->w = cx2 - cx1 + 1;
        } else {
            cropParams->x = cx2;
            cropParams->w = cx1 - cx2 + 1;
        }

        if (cy2 > cy1) {
            cropParams->y = cy1;
            cropParams->h = cy2 - cy1 + 1;
        } else {
            cropParams->y = cy2;
            cropParams->h = cy1 - cy2 + 1;
        }

        {
        MYREADERLOCK(l, lockRW);
        updateBackBuffer ();
        }
        parent->redrawNeeded (this);
    }

    return b;
}

bool FileBrowserEntry::pressNotify   (int button, int type, int bstate, int x, int y)
{

    bool b = ThumbBrowserEntryBase::pressNotify (button, type, bstate, x, y);

    if (!iatlistener || !iatlistener->getToolBar()) {
        return true;
    }

    ToolMode tm = iatlistener->getToolBar()->getTool ();
    int ix = x - startx - ofsX;
    int iy = y - starty - ofsY;

    if (tm == TMNone) {
        return b;
    }

    crop_custom_ratio = 0.f;

    if (!b && selected && inside (x, y)) {
        if (button == 1 && type == GDK_BUTTON_PRESS && state == SNormal) {
            if ((bstate & GDK_SHIFT_MASK) && cropParams->w > 0 && cropParams->h > 0) {
                crop_custom_ratio = float(cropParams->w) / float(cropParams->h);
            }
            if (onArea (CropTopLeft, ix, iy)) {
                state = SResizeTL;
                press_x = x;
                action_x = cropParams->x;
                press_y = y;
                action_y = cropParams->y;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropTopRight, ix, iy)) {
                state = SResizeTR;
                press_x = x;
                action_x = cropParams->w;
                press_y = y;
                action_y = cropParams->y;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropBottomLeft, ix, iy)) {
                state = SResizeBL;
                press_x = x;
                action_x = cropParams->x;
                press_y = y;
                action_y = cropParams->h;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropBottomRight, ix, iy)) {
                state = SResizeBR;
                press_x = x;
                action_x = cropParams->w;
                press_y = y;
                action_y = cropParams->h;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropTop, ix, iy)) {
                state = SResizeH1;
                press_y = y;
                action_y = cropParams->y;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropBottom, ix, iy)) {
                state = SResizeH2;
                press_y = y;
                action_y = cropParams->h;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropLeft, ix, iy)) {
                state = SResizeW1;
                press_x = x;
                action_x = cropParams->x;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropRight, ix, iy)) {
                state = SResizeW2;
                press_x = x;
                action_x = cropParams->w;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if ((bstate & GDK_SHIFT_MASK) && onArea (CropInside, ix, iy)) {
                state = SCropMove;
                press_x = x;
                press_y = y;
                action_x = cropParams->x;
                action_y = cropParams->y;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            } else if (onArea (CropImage, ix, iy)) {
                if (tm == TMStraighten) {
                    state = SRotateSelecting;
                    press_x = x;
                    press_y = y;
                    action_x = x;
                    action_y = y;
                    rot_deg = 0;
                    b = true;
                } else if (tm == TMSpotWB) {
                    iatlistener->spotWBselected ((ix - prex) / scale, (iy - prey) / scale, thumbnail);
                    b = true;
                } else if (tm == TMCropSelect) {
                    cropgl = iatlistener->startCropEditing (thumbnail);

                    if (cropgl) {
                        state = SCropSelecting;
                        press_x = cropParams->x = (ix - prex) / scale;
                        press_y = cropParams->y = (iy - prey) / scale;
                        cropParams->w = cropParams->h = 1;
                        cropgl->cropInit (cropParams->x, cropParams->y, cropParams->w, cropParams->h);
                        b = true;
                    }
                }
            }
        }

        updateCursor (ix, iy);
    }

    return b;
}

bool FileBrowserEntry::releaseNotify (int button, int type, int bstate, int x, int y)
{

    bool b = ThumbBrowserEntryBase::releaseNotify (button, type, bstate, x, y);

    int ix = x - startx - ofsX;
    int iy = y - starty - ofsY;

    if (!b) {
        if (state == SRotateSelecting) {
            iatlistener->rotateSelectionReady (rot_deg, thumbnail);

            if (iatlistener->getToolBar()) {
                iatlistener->getToolBar()->setTool (TMHand);
            }
        } else if (cropgl && (state == SCropSelecting || state == SResizeH1 || state == SResizeH2 || state == SResizeW1 || state == SResizeW2 || state == SResizeTL || state == SResizeTR || state == SResizeBL || state == SResizeBR || state == SCropMove)) {
            cropgl->cropManipReady ();
            cropgl = nullptr;
            iatlistener->cropSelectionReady ();

            if (iatlistener->getToolBar()) {
                iatlistener->getToolBar()->setTool (TMHand);
            }
        }

        state = SNormal;

        if (parent) {
            parent->redrawNeeded (this);
        }

        updateCursor (ix, iy);
    }

    return b;
}

bool FileBrowserEntry::onArea (CursorArea a, int x, int y)
{

    MYREADERLOCK(l, lockRW);
    if (!drawable || preview.empty()) {
        return false;
    }

    int x1 = (x - prex) / scale;
    int y1 = (y - prey) / scale;
    int cropResizeBorder = CROPRESIZEBORDER / scale;

    switch (a) {
    case CropImage:
        return x >= prex && x < prex + prew && y >= prey && y < prey + preh;

    case CropTopLeft:
        return cropParams->enabled &&
               y1 >= cropParams->y - cropResizeBorder &&
               y1 <= cropParams->y + cropResizeBorder &&
               x1 >= cropParams->x - cropResizeBorder &&
               x1 <= cropParams->x + cropResizeBorder;

    case CropTopRight:
        return cropParams->enabled &&
               y1 >= cropParams->y - cropResizeBorder &&
               y1 <= cropParams->y + cropResizeBorder &&
               x1 >= cropParams->x + cropParams->w - 1 - cropResizeBorder &&
               x1 <= cropParams->x + cropParams->w - 1 + cropResizeBorder;

    case CropBottomLeft:
        return cropParams->enabled &&
               y1 >= cropParams->y + cropParams->h - 1 - cropResizeBorder &&
               y1 <= cropParams->y + cropParams->h - 1 + cropResizeBorder &&
               x1 >= cropParams->x - cropResizeBorder &&
               x1 <= cropParams->x + cropResizeBorder;

    case CropBottomRight:
        return cropParams->enabled &&
               y1 >= cropParams->y + cropParams->h - 1 - cropResizeBorder &&
               y1 <= cropParams->y + cropParams->h - 1 + cropResizeBorder &&
               x1 >= cropParams->x + cropParams->w - 1 - cropResizeBorder &&
               x1 <= cropParams->x + cropParams->w - 1 + cropResizeBorder;

    case CropTop:
        return cropParams->enabled &&
               x1 > cropParams->x + cropResizeBorder &&
               x1 < cropParams->x + cropParams->w - 1 - cropResizeBorder &&
               y1 > cropParams->y - cropResizeBorder &&
               y1 < cropParams->y + cropResizeBorder;

    case CropBottom:
        return cropParams->enabled &&
               x1 > cropParams->x + cropResizeBorder &&
               x1 < cropParams->x + cropParams->w - 1 - cropResizeBorder &&
               y1 > cropParams->y + cropParams->h - 1 - cropResizeBorder &&
               y1 < cropParams->y + cropParams->h - 1 + cropResizeBorder;

    case CropLeft:
        return cropParams->enabled &&
               y1 > cropParams->y + cropResizeBorder &&
               y1 < cropParams->y + cropParams->h - 1 - cropResizeBorder &&
               x1 > cropParams->x - cropResizeBorder &&
               x1 < cropParams->x + cropResizeBorder;

    case CropRight:
        return cropParams->enabled &&
               y1 > cropParams->y + cropResizeBorder &&
               y1 < cropParams->y + cropParams->h - 1 - cropResizeBorder &&
               x1 > cropParams->x + cropParams->w - 1 - cropResizeBorder &&
               x1 < cropParams->x + cropParams->w - 1 + cropResizeBorder;

    case CropInside:
        return cropParams->enabled &&
               y1 > cropParams->y &&
               y1 < cropParams->y + cropParams->h - 1 &&
               x1 > cropParams->x &&
               x1 < cropParams->x + cropParams->w - 1;
    default: /* do nothing */ ;
    }

    return false;
}


void FileBrowserEntry::updateCursor (int x, int y)
{

    if (!iatlistener || !iatlistener->getToolBar()) {
        return;
    }

    CursorShape newCursor = CSUndefined;

    ToolMode tm = iatlistener->getToolBar()->getTool ();
    Glib::RefPtr<Gdk::Window> w = parent->getDrawingArea ()->get_window();

    if (!selected) {
        newCursor = CSArrow;
    } else if (state == SNormal) {
        if (tm == TMHand && (onArea (CropTop, x, y) || onArea (CropBottom, x, y))) {
            newCursor = CSResizeHeight;
        } else if (tm == TMHand && (onArea (CropLeft, x, y) || onArea (CropRight, x, y))) {
            newCursor = CSResizeWidth;
        } else if (tm == TMHand && (onArea (CropTopLeft, x, y))) {
            newCursor = CSResizeTopLeft;
        } else if (tm == TMHand && (onArea (CropTopRight, x, y))) {
            newCursor = CSResizeTopRight;
        } else if (tm == TMHand && (onArea (CropBottomLeft, x, y))) {
            newCursor = CSResizeBottomLeft;
        } else if (tm == TMHand && (onArea (CropBottomRight, x, y))) {
            newCursor = CSResizeBottomRight;
        } else if (onArea (CropImage, x, y)) {
            if (tm == TMHand) {
                newCursor = CSArrow;
            } else if (tm == TMSpotWB) {
                newCursor = CSSpotWB;
            } else if (tm == TMCropSelect) {
                newCursor = CSCropSelect;
            } else if (tm == TMStraighten) {
                newCursor = CSStraighten;
            }
        } else {
            newCursor = CSArrow;
        }
    } else if (state == SCropSelecting) {
        newCursor = CSCropSelect;
    } else if (state == SRotateSelecting) {
        newCursor = CSStraighten;
    } else if (state == SCropMove) {
        newCursor = CSMove;
    } else if (state == SResizeW1 || state == SResizeW2) {
        newCursor = CSResizeWidth;
    } else if (state == SResizeH1 || state == SResizeH2) {
        newCursor = CSResizeHeight;
    } else if (state == SResizeTL) {
        newCursor = CSResizeTopLeft;
    } else if (state == SResizeTR) {
        newCursor = CSResizeTopRight;
    } else if (state == SResizeBL) {
        newCursor = CSResizeBottomLeft;
    } else if (state == SResizeBR) {
        newCursor = CSResizeBottomRight;
    }

    if (newCursor != cursor_type) {
        cursor_type = newCursor;
        CursorManager::setCursorOfMainWindow (w, cursor_type);
    }
}

void FileBrowserEntry::draw (Cairo::RefPtr<Cairo::Context> cc)
{

    ThumbBrowserEntryBase::draw (cc);

    if (state == SRotateSelecting) {
        drawStraightenGuide (cc);
    }
}

void FileBrowserEntry::drawStraightenGuide (Cairo::RefPtr<Cairo::Context> cr)
{

    if (action_x != press_x || action_y != press_y) {
        double arg = (press_x - action_x) / sqrt(double((press_x - action_x) * (press_x - action_x) + (press_y - action_y) * (press_y - action_y)));
        double sol1, sol2;
        double pi = rtengine::RT_PI;

        if (press_y > action_y) {
            sol1 = acos(arg) * 180 / pi;
            sol2 = -acos(-arg) * 180 / pi;
        } else {
            sol1 = acos(-arg) * 180 / pi;
            sol2 = -acos(arg) * 180 / pi;
        }

        if (fabs(sol1) < fabs(sol2)) {
            rot_deg = sol1;
        } else {
            rot_deg = sol2;
        }

        if (rot_deg < -45) {
            rot_deg = 90.0 + rot_deg;
        } else if (rot_deg > 45) {
            rot_deg = - 90.0 + rot_deg;
        }
    } else {
        rot_deg = 0;
    }

    Glib::RefPtr<Pango::Context> context = parent->getDrawingArea()->get_pango_context () ;
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_BOLD);
    fontd.set_size (8 * Pango::SCALE);
    context->set_font_description (fontd);
    Glib::RefPtr<Pango::Layout> deglayout = parent->getDrawingArea()->create_pango_layout(Glib::ustring::compose ("%1 deg", Glib::ustring::format(std::setprecision(2), rot_deg)));

    int x1 = press_x;
    int y1 = press_y;
    int y2 = action_y;
    int x2 = action_x;

    {
    MYREADERLOCK(l, lockRW);
    if (x2 < prex + ofsX + startx) {
        y2 = y1 - (double)(y1 - y2) * (x1 - (prex + ofsX + startx)) / (x1 - x2);
        x2 = prex + ofsX + startx;
    } else if (x2 >= prew + prex + ofsX + startx) {
        y2 = y1 - (double)(y1 - y2) * (x1 - (prew + prex + ofsX + startx - 1)) / (x1 - x2);
        x2 = prew + prex + ofsX + startx - 1;
    }

    if (y2 < prey + ofsY + starty) {
        x2 = x1 - (double)(x1 - x2) * (y1 - (prey + ofsY + starty)) / (y1 - y2);
        y2 = prey + ofsY + starty;
    } else if (y2 >= preh + prey + ofsY + starty) {
        x2 = x1 - (double)(x1 - x2) * (y1 - (preh + prey + ofsY + starty - 1)) / (y1 - y2);
        y2 = preh + prey + ofsY + starty - 1;
    }
    }

    cr->set_line_width (1.5);
    cr->set_source_rgb (1.0, 1.0, 1.0);
    cr->move_to (x1, y1);
    cr->line_to (x2, y2);
    cr->stroke ();
    cr->set_source_rgb (0.0, 0.0, 0.0);
    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (x1, y1);
    cr->line_to (x2, y2);
    cr->stroke ();

    if (press_x != action_x && press_y != action_y) {
        cr->set_source_rgb (0.0, 0.0, 0.0);
        cr->move_to ((x1 + x2) / 2 + 1, (y1 + y2) / 2 + 1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1 + x2) / 2 + 1, (y1 + y2) / 2 - 1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1 + x2) / 2 - 1, (y1 + y2) / 2 + 1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1 + x2) / 2 + 1, (y1 + y2) / 2 + 1);
        deglayout->add_to_cairo_context (cr);
        cr->fill ();
        cr->set_source_rgb (1.0, 1.0, 1.0);
        cr->move_to ((x1 + x2) / 2, (y1 + y2) / 2);
        deglayout->add_to_cairo_context (cr);
        cr->fill ();
    }
}

