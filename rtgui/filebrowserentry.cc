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
#include <filebrowserentry.h>
#include <thumbbrowserbase.h>
#include <cursormanager.h>
#include <iomanip>
#include <guiutils.h>
#include <safegtk.h>

#define CROPRESIZEBORDER 4

bool FileBrowserEntry::iconsLoaded = false;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::editedIcon;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::recentlySavedIcon;
Glib::RefPtr<Gdk::Pixbuf> FileBrowserEntry::enqueuedIcon;

FileBrowserEntry::FileBrowserEntry (Thumbnail* thm, const Glib::ustring& fname) 
    : ThumbBrowserEntryBase (fname), thumbnail(thm), iatlistener(NULL), state(SNormal), cropgl(NULL) {

    feih = new FileBrowserEntryIdleHelper;
    feih->fbentry = this;
    feih->destroyed = false;
    feih->pending = 0;
    
    italicstyle = thumbnail->getType() != FT_Raw;
    datetimeline = thumbnail->getDateTimeString ();
    exifline = thumbnail->getExifString ();

    scale = 1;
    
    if (!iconsLoaded) {
        editedIcon = safe_create_from_file (argv0+"/images/edited.png");
        recentlySavedIcon = safe_create_from_file (argv0+"/images/saved.png");
        enqueuedIcon = safe_create_from_file (argv0+"/images/processing.png");
    }
    
    if (thm)
        thm->addThumbnailListener (this);
}

FileBrowserEntry::~FileBrowserEntry () {

    thumbImageUpdater.removeJobs (this);
    if (thumbnail)
        thumbnail->removeThumbnailListener (this);

    if (feih->pending)
        feih->destroyed = true;
    else
        delete feih;
}

void FileBrowserEntry::refreshThumbnailImage () {

    if (!thumbnail)
        return;

    thumbImageUpdater.add (thumbnail, thumbnail->getProcParams(), preh, &updatepriority, this);    
    thumbImageUpdater.process ();
}

void FileBrowserEntry::calcThumbnailSize () {

    if (thumbnail)
        thumbnail->getThumbnailSize (prew, preh);
}

std::vector<Glib::RefPtr<Gdk::Pixbuf> > FileBrowserEntry::getIconsOnImageArea () {

    std::vector<Glib::RefPtr<Gdk::Pixbuf> > ret;
    
    if (!thumbnail)
        return ret;

    if (thumbnail->hasProcParams() && editedIcon)
        ret.push_back (editedIcon);
    if (thumbnail->isRecentlySaved() && recentlySavedIcon)
        ret.push_back (recentlySavedIcon);
    if (thumbnail->isEnqueued () && enqueuedIcon)
        ret.push_back (enqueuedIcon);

   return ret;
}

void FileBrowserEntry::customBackBufferUpdate (Cairo::RefPtr<Cairo::Context> c) {
    
    if (state==SCropSelecting || state==SResizeH1 || state==SResizeH2 || state==SResizeW1 || state==SResizeW2 || state==SCropMove)
        drawCrop (c, prex, prey, prew, preh, 0, 0, scale, cropParams);
    else {
        rtengine::procparams::CropParams cparams = thumbnail->getProcParams().crop;
        if (cparams.enabled)
            drawCrop (c, prex, prey, prew, preh, 0, 0, scale, cparams);
    }
}

void FileBrowserEntry::getIconSize (int& w, int& h) {

    w = editedIcon->get_width ();
    h = editedIcon->get_height ();
}

FileThumbnailButtonSet* FileBrowserEntry::getThumbButtonSet () {

    return (FileThumbnailButtonSet*)buttonSet;
}

void FileBrowserEntry::procParamsChanged (Thumbnail* thm, int whoChangedIt) {

    refreshThumbnailImage ();
}

struct tiupdate {
    FileBrowserEntryIdleHelper* feih;
    rtengine::IImage8* img;
    double scale;
    rtengine::procparams::CropParams cropParams;
};

int fbeupdate (void* data) {
    
    gdk_threads_enter ();
    tiupdate* params = (tiupdate*)data;
    FileBrowserEntryIdleHelper* feih = params->feih;

    if (feih->destroyed) {
        if (feih->pending == 1)
            delete feih;
        else    
            feih->pending--;
        params->img->free ();
        delete params;
        gdk_threads_leave ();
        return 0;
    }
    
    feih->fbentry->_updateImage (params->img, params->scale, params->cropParams);
    feih->pending--;

    gdk_threads_leave ();
    delete params;
    
    return 0;
}

void FileBrowserEntry::updateImage (rtengine::IImage8* img, double scale, rtengine::procparams::CropParams cropParams) {

    redrawRequests++;
    feih->pending++;
    tiupdate* param = new tiupdate ();
    param->feih = feih;
    param->img = img;
    param->scale = scale;
    param->cropParams = cropParams;
    g_idle_add (fbeupdate, param);
}

void FileBrowserEntry::_updateImage (rtengine::IImage8* img, double s, rtengine::procparams::CropParams cropParams) {

    redrawRequests--; 
    scale = s;
    this->cropParams = cropParams;
    if (preh == img->getHeight ()) {
        prew = img->getWidth ();
        guint8* temp = preview;
        preview = NULL;
        delete [] temp;
        temp = new guint8 [prew*preh*3];
        memcpy (temp, img->getData(), prew*preh*3);
        preview = temp;
        updateBackBuffer ();
    }
    if (redrawRequests==0 && parent) 
        parent->redrawNeeded (this);
    img->free ();
}

bool FileBrowserEntry::motionNotify (int x, int y) {
    
    bool b = ThumbBrowserEntryBase::motionNotify (x, y);
    
    int ix = x - startx - ofsX;
    int iy = y - starty - ofsY;

    if (inside (x,y))
        updateCursor (ix, iy);

    if (state==SRotateSelecting) {
        action_x = x;
        action_y = y;
        parent->redrawNeeded (this);
    }
    else if (state==SResizeH1 && cropgl) { 
        int oy = cropParams.y;
        cropParams.y = action_y + (y-press_y) / scale;
        cropParams.h += oy - cropParams.y;
        cropgl->cropHeight1Resized (cropParams.x, cropParams.y, cropParams.w, cropParams.h);
        updateBackBuffer ();
        parent->redrawNeeded (this);
    }
    else if (state==SResizeH2 && cropgl) { 
        cropParams.h = action_y + (y-press_y) / scale;
        cropgl->cropHeight2Resized (cropParams.x, cropParams.y, cropParams.w, cropParams.h);
        updateBackBuffer ();
        parent->redrawNeeded (this);
    }
    else if (state==SResizeW1 && cropgl) { 
        int ox = cropParams.x;
        cropParams.x = action_x + (x-press_x) / scale;
        cropParams.w += ox - cropParams.x;
        cropgl->cropWidth1Resized (cropParams.x, cropParams.y, cropParams.w, cropParams.h);
        updateBackBuffer ();
        parent->redrawNeeded (this);
    }
    else if (state==SResizeW2 && cropgl) { 
        cropParams.w = action_x + (x-press_x) / scale;
        cropgl->cropWidth2Resized (cropParams.x, cropParams.y, cropParams.w, cropParams.h);
        updateBackBuffer ();
        parent->redrawNeeded (this);
    }
    else if (state==SCropMove && cropgl) { 
        cropParams.x = action_x + (x-press_x) / scale;
        cropParams.y = action_y + (y-press_y) / scale;
        cropgl->cropMoved (cropParams.x, cropParams.y, cropParams.w, cropParams.h);
        updateBackBuffer ();
        parent->redrawNeeded (this);
    }
    else if (state==SCropSelecting && cropgl) { 
        int cx1 = press_x, cy1 = press_y;
        int cx2 = (ix-prex) / scale, cy2 = (iy-prey) / scale;
        cropgl->cropResized (cx1, cy1, cx2, cy2);
        if (cx2 > cx1) {
            cropParams.x = cx1;
            cropParams.w = cx2 - cx1 + 1;
        }
        else {
            cropParams.x = cx2;
            cropParams.w = cx1 - cx2 + 1;
        }
        if (cy2 > cy1) {
            cropParams.y = cy1;
            cropParams.h = cy2 - cy1 + 1;
        }
        else {
            cropParams.y = cy2;
            cropParams.h = cy1 - cy2 + 1;
        }
        updateBackBuffer ();
        parent->redrawNeeded (this);
    }
       
    return b;
}

bool FileBrowserEntry::pressNotify   (int button, int type, int bstate, int x, int y) {

    bool b = ThumbBrowserEntryBase::pressNotify (button, type, bstate, x, y);

    ToolMode tm = iatlistener->getToolBar()->getTool ();
    int ix = x - startx - ofsX;
    int iy = y - starty - ofsY;
    if (!b && selected && inside (x,y)) {
        if (button==1 && type==GDK_BUTTON_PRESS && state==SNormal) {
            if (onArea (CropTop, ix, iy)) {
                state = SResizeH1;
                press_y = y;
                action_y = cropParams.y;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            }
            else if (onArea (CropBottom, ix, iy)) {
                state = SResizeH2;
                press_y = y;
                action_y = cropParams.h;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            }
            else if (onArea (CropLeft, ix, iy)) {
                state = SResizeW1;
                press_x = x;
                action_x = cropParams.x;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            }
            else if (onArea (CropRight, ix, iy)) {
                state = SResizeW2;
                press_x = x;
                action_x = cropParams.w;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            }
            else if ((bstate & GDK_SHIFT_MASK) && onArea (CropInside, ix, iy)) {
                state = SCropMove;
                press_x = x;
                press_y = y;
                action_x = cropParams.x;
                action_y = cropParams.y;
                cropgl = iatlistener->startCropEditing (thumbnail);
                b = true;
            }
            else if (onArea (CropImage, ix, iy)) {
                if (tm == TMStraighten) {
                    state = SRotateSelecting;
                    press_x = x;
                    press_y = y;
                    action_x = x;
                    action_y = y;
                    rot_deg = 0;
                    b = true;
                }
                else if (tm == TMSpotWB) {
                    iatlistener->spotWBselected ((ix-prex)/scale, (iy-prey)/scale, thumbnail);
                    b = true;
                }
                else if (tm == TMCropSelect) {
                    cropgl = iatlistener->startCropEditing (thumbnail);
                    if (cropgl) {
                        state = SCropSelecting;
                        press_x = cropParams.x = (ix-prex) / scale;
                        press_y = cropParams.y = (iy-prey) / scale;
                        cropParams.w = cropParams.h = 1;
                        cropgl->cropInit (cropParams.x, cropParams.y, cropParams.w, cropParams.h);
                        b = true;
                    }
                }
            }
        }
        updateCursor (ix, iy);
    }
    return b;
}

bool FileBrowserEntry::releaseNotify (int button, int type, int bstate, int x, int y) {

    bool b = ThumbBrowserEntryBase::releaseNotify (button, type, bstate, x, y);
 
    int ix = x - startx - ofsX;
    int iy = y - starty - ofsY;
    if (!b) {
        if (state==SRotateSelecting) {
            iatlistener->rotateSelectionReady (rot_deg, thumbnail);
            iatlistener->getToolBar()->setTool (TMHand);
        }
        else if (cropgl && (state==SCropSelecting || state==SResizeH1 || state==SResizeH2 || state==SResizeW1 || state==SResizeW2 || state==SCropMove)) {
            cropgl->cropManipReady ();
            cropgl = NULL;
            iatlistener->cropSelectionReady ();
            iatlistener->getToolBar()->setTool (TMHand);
        }
        state = SNormal;
        if (parent)
            parent->redrawNeeded (this);
        updateCursor (ix, iy);
    }   
    
    return b;
}

bool FileBrowserEntry::onArea (CursorArea a, int x, int y) {

    if (!drawable || !preview)
        return false;

    int x1 = (x-prex) / scale;
    int y1 = (y-prey) / scale;
    int cropResizeBorder = CROPRESIZEBORDER / scale;
    switch (a) {
        case CropImage:
            return x>=prex && x<prex+prew && y>=prey && y<prey+preh;
        case CropTop:
            return cropParams.enabled && 
                x1>cropParams.x+cropResizeBorder && 
                x1<cropParams.x+cropParams.w-1-cropResizeBorder && 
                y1>cropParams.y-cropResizeBorder && 
                y1<cropParams.y+cropResizeBorder;
        case CropBottom:
            return cropParams.enabled && 
                x1>cropParams.x+cropResizeBorder && 
                x1<cropParams.x+cropParams.w-1-cropResizeBorder && 
                y1>cropParams.y+cropParams.h-1-cropResizeBorder && 
                y1<cropParams.y+cropParams.h-1+cropResizeBorder;
        case CropLeft:
            return cropParams.enabled && 
                y1>cropParams.y+cropResizeBorder && 
                y1<cropParams.y+cropParams.h-1-cropResizeBorder && 
                x1>cropParams.x-cropResizeBorder && 
                x1<cropParams.x+cropResizeBorder;
        case CropRight:
            return cropParams.enabled && 
                y1>cropParams.y+cropResizeBorder && 
                y1<cropParams.y+cropParams.h-1-cropResizeBorder && 
                x1>cropParams.x+cropParams.w-1-cropResizeBorder && 
                x1<cropParams.x+cropParams.w-1+cropResizeBorder;
        case CropInside:
            return cropParams.enabled && 
                y1>cropParams.y && 
                y1<cropParams.y+cropParams.h-1 && 
                x1>cropParams.x && 
                x1<cropParams.x+cropParams.w-1;
    }
    return false;
}


void FileBrowserEntry::updateCursor (int x, int y) {
    
    if (!iatlistener)
        return;
        
    ToolMode tm = iatlistener->getToolBar()->getTool ();
    Glib::RefPtr<Gdk::Window> w = parent->getDrawingArea ()->get_window();
    
    if (!selected) {
        cursorManager.setCursor (w, CSArrow);
        return;
    }
    
    if (state==SNormal) {
        if (tm==TMHand && (onArea (CropTop, x, y) || onArea (CropBottom, x, y))) 
            cursorManager.setCursor (w, CSResizeHeight);
        else if (tm==TMHand && (onArea (CropLeft, x, y) || onArea (CropRight, x, y))) 
            cursorManager.setCursor (w, CSResizeWidth);
        else if (onArea (CropImage, x, y)) { 
            if (tm==TMHand)
                cursorManager.setCursor (w, CSArrow);
            else if (tm==TMSpotWB)
                cursorManager.setCursor (w, CSSpotWB);
            else if (tm==TMCropSelect)
                cursorManager.setCursor (w, CSCropSelect);
            else if (tm==TMStraighten)
                cursorManager.setCursor (w, CSStraighten);
        }
        else
            cursorManager.setCursor (w, CSArrow);
    }
    else if (state==SCropSelecting)
        cursorManager.setCursor (w, CSCropSelect);
    else if (state==SRotateSelecting) 
        cursorManager.setCursor (w, CSStraighten);
    else if (state==SCropMove)
        cursorManager.setCursor (w, CSMove);
    else if (state==SResizeW1 || state==SResizeW2)
        cursorManager.setCursor (w, CSResizeWidth);
    else if (state==SResizeH1 || state==SResizeH2)
        cursorManager.setCursor (w, CSResizeHeight);
}

void FileBrowserEntry::draw () {

    ThumbBrowserEntryBase::draw ();
    if (state==SRotateSelecting) {
        Cairo::RefPtr<Cairo::Context> cr = parent->getDrawingArea ()->get_window()->create_cairo_context();
        drawStraightenGuide (cr);
    }
}

void FileBrowserEntry::drawStraightenGuide (Cairo::RefPtr<Cairo::Context> cr) {

    if (action_x!=press_x || action_y!=press_y) {
        double arg = (press_x-action_x) / sqrt((press_x-action_x)*(press_x-action_x)+(press_y-action_y)*(press_y-action_y));
        double sol1, sol2;
        double pi = M_PI;
        if (press_y>action_y) {
            sol1 = acos(arg)*180/pi;
            sol2 = -acos(-arg)*180/pi;
        }
        else {
            sol1 = acos(-arg)*180/pi;
            sol2 = -acos(arg)*180/pi;
        }
        if (fabs(sol1)<fabs(sol2))
            rot_deg = sol1;
        else
           rot_deg = sol2;

        if (rot_deg<-45)
           rot_deg = 90.0 + rot_deg;
        else if (rot_deg>45)
           rot_deg = - 90.0 + rot_deg;
    }
    else
        rot_deg = 0;

    Glib::RefPtr<Pango::Context> context = parent->getDrawingArea()->get_pango_context () ;
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_BOLD);
    fontd.set_size (8*Pango::SCALE);
    context->set_font_description (fontd);
    Glib::RefPtr<Pango::Layout> deglayout = parent->getDrawingArea()->create_pango_layout(Glib::ustring::compose ("%1 deg", Glib::ustring::format(std::setprecision(2), rot_deg)));

    int x1 = press_x;
    int y1 = press_y;
    int y2 = action_y;
    int x2 = action_x;
    
    if (x2<prex+ofsX+startx) {
        y2 = y1 - (double)(y1-y2)*(x1 - (prex+ofsX+startx)) / (x1-x2);
        x2 = prex+ofsX+startx;
    }
    else if (x2>=prew+prex+ofsX+startx) {
        y2 = y1 - (double)(y1-y2)*(x1 - (prew+prex+ofsX+startx-1)) / (x1-x2);
        x2 = prew+prex+ofsX+startx-1;
    }
    if (y2<prey+ofsY+starty) {
        x2 = x1 - (double)(x1-x2)*(y1 - (prey+ofsY+starty)) / (y1-y2);
        y2 = prey+ofsY+starty;
    }
    else if (y2>=preh+prey+ofsY+starty) {
        x2 = x1 - (double)(x1-x2)*(y1 - (preh+prey+ofsY+starty-1)) / (y1-y2);
        y2 = preh+prey+ofsY+starty-1;
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

    if (press_x!=action_x && press_y!=action_y) {
        cr->set_source_rgb (0.0, 0.0, 0.0);
        cr->move_to ((x1+x2)/2+1, (y1+y2)/2+1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1+x2)/2+1, (y1+y2)/2-1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1+x2)/2-1, (y1+y2)/2+1);
        deglayout->add_to_cairo_context (cr);
        cr->move_to ((x1+x2)/2+1, (y1+y2)/2+1);
        deglayout->add_to_cairo_context (cr);
        cr->fill ();
        cr->set_source_rgb (1.0, 1.0, 1.0);
        cr->move_to ((x1+x2)/2, (y1+y2)/2);
        deglayout->add_to_cairo_context (cr);
        cr->fill ();
    }
}

