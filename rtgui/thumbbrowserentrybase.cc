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
#include "thumbbrowserentrybase.h"
#include "thumbbrowserbase.h"
#include "options.h"
#include "../rtengine/mytime.h"

ThumbBrowserEntryBase::ThumbBrowserEntryBase (const Glib::ustring& fname)
    : fnlabw(0), fnlabh(0), dtlabw(0), dtlabh(0), exlabw(0), exlabh(0), prew(0), preh(0),
      prex(0), prey(0), upperMargin(6), borderWidth(1), textGap(6), sideMargin(8), lowerMargin(8),
      preview(NULL), dispname(Glib::path_get_basename (fname)), buttonSet(NULL), width(0), height(0),
      exp_width(0), exp_height(0), startx(0), starty(0), ofsX(0), ofsY(0), redrawRequests(0),
      parent(NULL), original(NULL), bbSelected(false), bbFramed(false), bbPreview(NULL), cursor_type(CSUndefined),
      thumbnail(NULL), filename(fname), shortname(dispname), exifline(""), datetimeline(""),
      selected(false), drawable(false), filtered(false), framed(false), processing(false), italicstyle(false),
      edited(false), recentlysaved(false), updatepriority(false), withFilename(WFNAME_NONE) {}

ThumbBrowserEntryBase::~ThumbBrowserEntryBase ()
{

    if (preview) {
        delete [] preview;
    }

    delete buttonSet;
}

void ThumbBrowserEntryBase::addButtonSet (LWButtonSet* bs)
{

    buttonSet = bs;
}

void ThumbBrowserEntryBase::updateBackBuffer ()
{

    if (!parent) {
        return;
    }

    Gtk::Widget* w = parent->getDrawingArea ();

    Glib::RefPtr<Gdk::Window> win = w->get_window();

    if (!win)
        // Nothing to draw on, so we return
    {
        return;
    }

    if (backBuffer && (backBuffer->getWidth() != exp_width || backBuffer->getHeight() != exp_height )) {
        // deleting the existing BackBuffer
        backBuffer.reset();
    }
    if (!backBuffer) {
        backBuffer = Glib::RefPtr<BackBuffer> ( new BackBuffer (exp_width, exp_height, win) );
    }

    // If thumbnail is hidden by a filter, drawing to it will crash
    // if either with or height is zero then return early
    if (!backBuffer->getWidth() || !backBuffer->getHeight()) {
        return;
    }

    Cairo::RefPtr<Cairo::ImageSurface> surface = backBuffer->getSurface();

    bbSelected = selected;
    bbFramed = framed;
    bbPreview = preview;

    Cairo::RefPtr<Cairo::Context> cc = Cairo::Context::create(surface);

    Glib::RefPtr<Gtk::StyleContext> style = parent->getStyle();
    Gdk::RGBA textn = parent->getNormalTextColor();
    Gdk::RGBA texts = parent->getSelectedTextColor();
    Gdk::RGBA bgn = parent->getNormalBgColor();
    Gdk::RGBA bgs = parent->getSelectedBgColor();

    // clear area, draw frames and background
    style->render_background(cc, 0., 0., exp_width, exp_height);
    /*
    cc->set_line_width(0.);
    cc->set_line_cap(Cairo::LINE_CAP_BUTT);
    cc->set_antialias(Cairo::ANTIALIAS_NONE);
    cc->set_source_rgb(bgn.get_red(), bgn.get_green(), bgn.get_blue());
    cc->rectangle(0., 0., exp_width, exp_height);
    cc->fill();
    */

    cc->set_antialias(Cairo::ANTIALIAS_SUBPIXEL);

    drawFrame (cc, bgs, bgn);

    // calculate height of button set
    int bsHeight = 0;

    if (buttonSet) {
        int tmp;
        buttonSet->getAllocatedDimensions (tmp, bsHeight);
    }

    // draw preview frame
    //backBuffer->draw_rectangle (cc, false, (exp_width-prew)/2, upperMargin+bsHeight, prew+1, preh+1);
    // draw thumbnail image
    if (preview) {
        prex = borderWidth + (exp_width - prew) / 2;
        prey = upperMargin + bsHeight + borderWidth;
        backBuffer->copyRGBCharData(preview, 0, 0, prew, preh, prew * 3, prex, prey);
    }

    customBackBufferUpdate (cc);

    // draw icons onto the thumbnail area
    bbIcons = getIconsOnImageArea ();

    int infow, infoh;
    getTextSizes (infow, infoh);

    int iofs_x = 4, iofs_y = 4;
    int igap = 2;
    int istartx = prex;
    int istarty = prey;

    if ((parent->getLocation() != ThumbBrowserBase::THLOC_EDITOR && options.showFileNames && options.overlayedFileNames)
            || (parent->getLocation() == ThumbBrowserBase::THLOC_EDITOR && options.filmStripShowFileNames && options.filmStripOverlayedFileNames)) {
        cc->begin_new_path ();
        cc->rectangle (istartx, istarty, prew, fnlabh + dtlabh + exlabh + 2 * iofs_y);

        if ((texts.get_red() + texts.get_green() + texts.get_blue()) / 3 > 0.5) {
            cc->set_source_rgba (0, 0, 0, 0.5);
        } else {
            cc->set_source_rgba (1, 1, 1, 0.5);
        }

        cc->fill ();
    }

    istartx += iofs_x;
    istarty += iofs_y;

    if (!bbIcons.empty()) {
        int iwidth = 0;
        int iheight = 0;

        for (size_t i = 0; i < bbIcons.size(); i++) {
            iwidth += bbIcons[i]->get_width() + (i > 0 ? igap : 0);

            if (bbIcons[i]->get_height() > iheight) {
                iheight = bbIcons[i]->get_height();
            }
        }

        if ((parent->getLocation() != ThumbBrowserBase::THLOC_EDITOR && (!options.showFileNames || !options.overlayedFileNames))
                || (parent->getLocation() == ThumbBrowserBase::THLOC_EDITOR && (!options.filmStripShowFileNames || !options.filmStripOverlayedFileNames))) {
            // Draw the transparent black background around icons
            cc->begin_new_path ();
            cc->move_to(istartx - igap, istarty);
            cc->rel_line_to(igap, -igap);
            cc->rel_line_to(iwidth, 0);
            cc->rel_line_to(igap, igap);
            cc->rel_line_to(0, iheight);
            cc->rel_line_to(-igap, igap);
            cc->rel_line_to(-iwidth, 0);
            cc->rel_line_to(-igap, -igap);
            cc->rel_line_to(0, -iheight);
            cc->set_source_rgba (0, 0, 0, 0.6);
            cc->fill ();
        }

        for (size_t i = 0; i < bbIcons.size(); i++) {
            // Draw the image at 110, 90, except for the outermost 10 pixels.
            Gdk::Cairo::set_source_pixbuf(cc, bbIcons[i], istartx, istarty);
            cc->rectangle(istartx, istarty, bbIcons[i]->get_width(), bbIcons[i]->get_height());
            cc->fill();
            istartx += bbIcons[i]->get_width() + igap;
        }
    }

    if ( ( (parent->getLocation() != ThumbBrowserBase::THLOC_EDITOR && options.showFileNames)
            || (parent->getLocation() == ThumbBrowserBase::THLOC_EDITOR && options.filmStripShowFileNames))
            && withFilename > WFNAME_NONE) {
        int textposx_fn, textposx_ex, textposx_dt, textposy, textw;

        if (! ((parent->getLocation() != ThumbBrowserBase::THLOC_EDITOR && options.overlayedFileNames)
                || (parent->getLocation() == ThumbBrowserBase::THLOC_EDITOR && options.filmStripOverlayedFileNames)) ) {
            textposx_fn = exp_width / 2 - fnlabw / 2;

            if (textposx_fn < 0) {
                textposx_fn = 0;
            }

            textposx_ex = exp_width / 2 - exlabw / 2;

            if (textposx_ex < 0) {
                textposx_ex = 0;
            }

            textposx_dt = exp_width / 2 - dtlabw / 2;

            if (textposx_dt < 0) {
                textposx_dt = 0;
            }

            textposy = upperMargin + bsHeight + 2 * borderWidth + preh + borderWidth + textGap;
            textw = exp_width - 2 * textGap;

            if (selected) {
                cc->set_source_rgb(texts.get_red(), texts.get_green(), texts.get_blue());
            } else {
                cc->set_source_rgb(textn.get_red(), textn.get_green(), textn.get_blue());
            }
        } else {
            textposx_fn = istartx;
            textposx_ex = istartx;
            textposx_dt = istartx;
            textposy = istarty;
            textw = prew - (istartx - prex);
            cc->set_source_rgb(texts.get_red(), texts.get_green(), texts.get_blue());
        }

        // draw file name
        Glib::RefPtr<Pango::Context> context = w->get_pango_context () ;
        Pango::FontDescription fontd = context->get_font_description ();
        fontd.set_weight (Pango::WEIGHT_BOLD);

        if (italicstyle) {
            fontd.set_style (Pango::STYLE_ITALIC);
        } else {
            fontd.set_style (Pango::STYLE_NORMAL);
        }

        context->set_font_description (fontd);
        Glib::RefPtr<Pango::Layout> fn = w->create_pango_layout (dispname);
        fn->set_width (textw * Pango::SCALE);
        fn->set_ellipsize (Pango::ELLIPSIZE_MIDDLE);
        cc->move_to(textposx_fn, textposy);
        fn->add_to_cairo_context (cc);
        cc->fill();

        fontd.set_weight (Pango::WEIGHT_NORMAL);
        fontd.set_style (Pango::STYLE_NORMAL);
        context->set_font_description (fontd);

        if (withFilename == WFNAME_FULL) {
            // draw date/time label
            int tpos = fnlabh;

            if (options.fbShowDateTime && datetimeline != "") {
                fn = w->create_pango_layout (datetimeline);
                fn->set_width (textw * Pango::SCALE);
                fn->set_ellipsize (Pango::ELLIPSIZE_MIDDLE);
                cc->move_to(textposx_dt, textposy + tpos);
                fn->add_to_cairo_context (cc);
                cc->fill();
                tpos += dtlabh;
            }

            // draw basic exif info
            if (options.fbShowBasicExif && exifline != "") {
                fn = w->create_pango_layout (exifline);
                fn->set_width (textw * Pango::SCALE);
                fn->set_ellipsize (Pango::ELLIPSIZE_MIDDLE);
                cc->move_to(textposx_ex, textposy + tpos);
                fn->add_to_cairo_context (cc);
                cc->fill();
                tpos += exlabh;
            }
        }
    }

    backBuffer->setDirty(false);
}

void ThumbBrowserEntryBase::getTextSizes (int& infow, int& infoh)
{

    if (!parent) {
        return;
    }

    Gtk::Widget* w = parent->getDrawingArea ();

    // calculate dimensions of the text based fields
    dispname = shortname;

    Glib::RefPtr<Pango::Context> context = w->get_pango_context () ;
    context->set_font_description (w->get_style_context()->get_font());


    // filename:
    Pango::FontDescription fontd = context->get_font_description ();
    fontd.set_weight (Pango::WEIGHT_BOLD);
    context->set_font_description (fontd);
    Glib::RefPtr<Pango::Layout> fn = w->create_pango_layout(shortname);
    fn->get_pixel_size (fnlabw, fnlabh);

    // calculate cummulated height of all info fields
    infoh = fnlabh;
    infow = 0;

    if (withFilename == WFNAME_FULL) {
        // datetime
        fontd.set_weight (Pango::WEIGHT_NORMAL);
        context->set_font_description (fontd);
        fn = w->create_pango_layout (datetimeline);
        fn->get_pixel_size (dtlabw, dtlabh);

        // basic exif data
        fn = w->create_pango_layout (exifline);
        fn->get_pixel_size (exlabw, exlabh);

        // add date/tile size:
        if (options.fbShowDateTime) {
            infoh += dtlabh;

            if (dtlabw + 2 * sideMargin > infow) {
                infow = dtlabw + 2 * sideMargin;
            }
        } else {
            dtlabw = dtlabh = 0;
        }

        if (options.fbShowBasicExif) {
            infoh += exlabh;

            if (exlabw + 2 * sideMargin > infow) {
                infow = exlabw + 2 * sideMargin;
            }
        } else {
            exlabw = exlabh = 0;
        }
    } else {
        dtlabw = dtlabh = exlabw = exlabh = 0;
    }
}

void ThumbBrowserEntryBase::resize (int h)
{

#if PROTECT_VECTORS
    MYWRITERLOCK(l, lockRW);
#endif

    height = h;
    int old_preh = preh, old_width = width;

    // dimensions of the button set
    int bsw = 0, bsh = 0;

    if (buttonSet) {
        buttonSet->getMinimalDimensions (bsw, bsh);
    }

    if (parent->getLocation() == ThumbBrowserBase::THLOC_FILEBROWSER) {
        if (options.showFileNames) {
            withFilename = WFNAME_FULL;
        } else {
            withFilename = WFNAME_NONE;
        }
    } else if (parent->getLocation() == ThumbBrowserBase::THLOC_BATCHQUEUE) {
        withFilename = WFNAME_REDUCED;
    } else {
        if (options.filmStripShowFileNames) {
            withFilename = WFNAME_REDUCED;
        } else {
            withFilename = WFNAME_NONE;
        }
    }

    // calculate the height remaining for the thumbnail image
    preh = height - upperMargin - 2 * borderWidth - lowerMargin - bsh;
    int infow = 0;
    int infoh = 0;

    if (    (parent->getLocation() != ThumbBrowserBase::THLOC_EDITOR && options.showFileNames && !options.overlayedFileNames)
            || (parent->getLocation() == ThumbBrowserBase::THLOC_EDITOR && options.filmStripShowFileNames && !options.filmStripOverlayedFileNames)) {
        // dimensions of the info text
        getTextSizes (infow, infoh);
        infoh += textGap;
        //preh -= infoh;
        height += infoh;
    }

    // Minimum size for thumbs
    if (preh < 24) {
        preh = 24;
        height = preh + (upperMargin + 2 * borderWidth + lowerMargin) + bsh + infoh;
    }

    calcThumbnailSize ();  // recalculates prew

    width = prew + 2 * sideMargin + 2 * borderWidth;

    if (    (parent->getLocation() != ThumbBrowserBase::THLOC_EDITOR && options.showFileNames && !options.overlayedFileNames)
            || (parent->getLocation() == ThumbBrowserBase::THLOC_EDITOR && options.filmStripShowFileNames && !options.filmStripOverlayedFileNames)) {
        width = prew + 2 * sideMargin + 2 * borderWidth;

        if (width < infow + 2 * sideMargin + 2 * borderWidth) {
            width = infow + 2 * sideMargin + 2 * borderWidth;
        }
    }

    if (width < bsw + 2 * sideMargin + 2 * borderWidth) {
        width = bsw + 2 * sideMargin + 2 * borderWidth;
    }

    if ( preh != old_preh || width != old_width ) {
        delete [] preview;
        preview = NULL;
        refreshThumbnailImage ();
    } else if (backBuffer) {
        backBuffer->setDirty(true);    // This will force a backBuffer update on queue_draw
    }

    drawable = true;
}

void ThumbBrowserEntryBase::drawFrame (Cairo::RefPtr<Cairo::Context> cc, const Gdk::RGBA& bg, const Gdk::RGBA& fg)
{

    int radius = 4;

    if (selected || framed) {
        cc->move_to (radius, 0);
        cc->arc (exp_width - 1 - radius, radius, radius, -M_PI / 2, 0);
        cc->arc (exp_width - 1 - radius, exp_height - 1 - radius, radius, 0, M_PI / 2);
        cc->arc (radius, exp_height - 1 - radius, radius, M_PI / 2, M_PI);
        cc->arc (radius, radius, radius, M_PI, -M_PI / 2);
        cc->close_path ();

        if (selected) {
            cc->set_source_rgb (bg.get_red(), bg.get_green(), bg.get_blue());
            cc->fill_preserve ();
        }

        cc->set_source_rgb (bg.get_red() * 2 / 3, bg.get_green() * 2 / 3, bg.get_blue() * 2 / 3);
        cc->set_line_width (1.0);
        cc->stroke ();
    }

    if (framed) {
        cc->move_to (+2 + 0.5 + radius, +2 + 0.5);
        cc->arc (-2 + 0.5 + exp_width - 1 - radius, +2 + 0.5 + radius, radius, -M_PI / 2, 0);
        cc->arc (-2 + 0.5 + exp_width - 1 - radius, -2 + 0.5 + exp_height - 1 - radius, radius, 0, M_PI / 2);
        cc->arc (+2 + 0.5 + radius, -2 + exp_height - 1 - radius, radius, M_PI / 2, M_PI);
        cc->arc (+2 + 0.5 + radius, +2 + radius, radius, M_PI, -M_PI / 2);
        cc->close_path ();
        cc->set_source_rgb (fg.get_red(), fg.get_green(), fg.get_blue());
        cc->set_line_width (2.0);
        cc->stroke ();
    }
}

void ThumbBrowserEntryBase::draw (Cairo::RefPtr<Cairo::Context> cc)
{

    if (!drawable || !parent) {
        return;
    }

#if PROTECT_VECTORS
    MYREADERLOCK(l, lockRW);  // No resizes, position moves etc. inbetween
#endif

    int bbWidth, bbHeight;

    if (backBuffer) {
        bbWidth = backBuffer->getWidth();
        bbHeight = backBuffer->getHeight();
    }

    if (!backBuffer || selected != bbSelected || framed != bbFramed || preview != bbPreview
            || exp_width != bbWidth || exp_height != bbHeight || getIconsOnImageArea () != bbIcons || backBuffer->isDirty()) {
        updateBackBuffer ();
    }

    int w_ = startx + ofsX;
    int h_ = starty + ofsY;
    cc->set_source(backBuffer->getSurface(), w_, h_);
    cc->rectangle(w_, h_, backBuffer->getWidth(), backBuffer->getHeight());
    cc->fill();

    // check icon set changes!!!

//    drawProgressBar (window, cc, selected ? texts : textn, selected ? bgs : bgn, ofsX+startx, exp_width, ofsY+starty + upperMargin+bsHeight+borderWidth+preh+borderWidth+textGap+tpos, fnlabh);

    // redraw button set above the thumbnail
    if (buttonSet) {
        buttonSet->setColors (selected ? parent->getSelectedBgColor() : parent->getNormalBgColor(), selected ? parent->getNormalBgColor() : parent->getSelectedBgColor());
        buttonSet->redraw (cc);
    }
}

void ThumbBrowserEntryBase::setPosition (int x, int y, int w, int h)
{

#if PROTECT_VECTORS
    MYWRITERLOCK(l, lockRW);
#endif

    exp_width = w;
    exp_height = h;
    startx = x;
    starty = y;

    if (buttonSet) {
        buttonSet->arrangeButtons (ofsX + x + sideMargin, ofsY + y + upperMargin, w - 2 * sideMargin, -1);
    }
}

void ThumbBrowserEntryBase::setOffset (int x, int y)
{

#if PROTECT_VECTORS
    MYWRITERLOCK(l, lockRW);
#endif

    ofsX = -x;
    ofsY = -y;

    if (buttonSet) {
        buttonSet->move (ofsX + startx + sideMargin, ofsY + starty + upperMargin);
    }
}

bool ThumbBrowserEntryBase::inside (int x, int y)
{

    return x > ofsX + startx && x < ofsX + startx + exp_width && y > ofsY + starty && y < ofsY + starty + exp_height;
}

void ThumbBrowserEntryBase::getPosInImgSpace (int x, int y, rtengine::Coord2D &coord)
{

    coord.x = coord.y = -1.;

    if (preview) {
        x -= ofsX + startx;
        y -= ofsY + starty;

        if (x >= prex && x <= prex + prew && y >= prey && y <= prey + preh) {
            coord.x = double(x - prex) / double(prew);
            coord.y = double(y - prey) / double(preh);
        }
    }
}

bool ThumbBrowserEntryBase::insideWindow (int x, int y, int w, int h)
{

    return !(ofsX + startx > x + w || ofsX + startx + exp_width < x || ofsY + starty > y + h || ofsY + starty + exp_height < y);
}

std::vector<Glib::RefPtr<Gdk::Pixbuf> > ThumbBrowserEntryBase::getIconsOnImageArea()
{
    return std::vector<Glib::RefPtr<Gdk::Pixbuf> >();
}

void ThumbBrowserEntryBase::getIconSize(int& w, int& h)
{
    w = 0;
    h = 0;
}

bool ThumbBrowserEntryBase::motionNotify  (int x, int y)
{

    return buttonSet ? buttonSet->motionNotify (x, y) : false;
}

bool ThumbBrowserEntryBase::pressNotify   (int button, int type, int bstate, int x, int y)
{

    return buttonSet ? buttonSet->pressNotify (x, y) : false;
}

bool ThumbBrowserEntryBase::releaseNotify (int button, int type, int bstate, int x, int y)
{

    return buttonSet ? buttonSet->releaseNotify (x, y) : false;
}

Glib::ustring ThumbBrowserEntryBase::getToolTip (int x, int y)
{
    Glib::ustring tooltip = "";

    if (buttonSet) {
        tooltip = buttonSet->getToolTip (x, y);
    }

    // if the fileinfo is not shown anyway, make a tooltip with the info
    if (withFilename < WFNAME_FULL && inside(x, y) && tooltip.empty()) {
        tooltip = dispname;

        if (options.fbShowDateTime && datetimeline != "") {
            tooltip += Glib::ustring("\n") + datetimeline;
        }

        if (options.fbShowBasicExif && exifline != "") {
            tooltip += Glib::ustring("\n") + exifline;
        }
    }

    return tooltip;
}


