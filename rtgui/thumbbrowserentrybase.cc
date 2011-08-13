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
#include <thumbbrowserentrybase.h>
#include <thumbbrowserbase.h>
#include <options.h>
#include <mytime.h>

ThumbBrowserEntryBase::ThumbBrowserEntryBase (const Glib::ustring& fname) 
    : preh(0), preview(NULL), buttonSet(NULL), exp_width(0), exp_height(0), redrawRequests(0),
      parent(NULL), filename(fname), exifline(""), datetimeline(""), selected(false),
      drawable(false),framed(false), processing(false), italicstyle(false),
      updatepriority(false) {

    shortname = Glib::path_get_basename (fname);
    dispname = shortname;

    upperMargin = 6;
    borderWidth = 1;
    sideMargin = 8;
    lowerMargin = 8;
    textGap = 6;
    
    ofsX = ofsY = 0;
}

ThumbBrowserEntryBase::~ThumbBrowserEntryBase () {

    delete [] preview;
    delete buttonSet;
}

void ThumbBrowserEntryBase::addButtonSet (LWButtonSet* bs) {

    buttonSet = bs;
}

void ThumbBrowserEntryBase::updateBackBuffer () {

    if (!parent)
        return;

    Gtk::Widget* w = parent->getDrawingArea ();

    backBuffer = Gdk::Pixmap::create (w->get_window(), exp_width, exp_height);

    // If thumbnail is hidden by a filter drawing to it will crash
    int backbuffer_w=0, backbuffer_h=0;
    backBuffer->get_size(backbuffer_w, backbuffer_h);
    // if either with or height is zero then return early
    if (backbuffer_w * backbuffer_h == 0) return;

    bbSelected = selected;
    bbFramed = framed;
    bbPreview = preview;

    Glib::RefPtr<Gdk::GC> gc_ = Gdk::GC::create (backBuffer);

    Gdk::Color textn = w->get_style()->get_text(Gtk::STATE_NORMAL);
    Gdk::Color texts = w->get_style()->get_text(Gtk::STATE_SELECTED);
    Gdk::Color bgn = w->get_style()->get_bg(Gtk::STATE_NORMAL);
    Gdk::Color bgs = w->get_style()->get_bg(Gtk::STATE_SELECTED);   

    // clear area, draw frames and background
    gc_->set_foreground (bgn);
    gc_->set_background (bgn);
    backBuffer->draw_rectangle (gc_, true, 0, 0, exp_width, exp_height);
    Cairo::RefPtr<Cairo::Context> cr = backBuffer->create_cairo_context();
    drawFrame (cr, bgs, bgn);

    // calculate height of button set
    int bsHeight = 0;
    if (buttonSet) {
        int tmp;
        buttonSet->getAllocatedDimensions (tmp, bsHeight);
    }

    // draw preview frame
    backBuffer->draw_rectangle (gc_, false, (exp_width-prew)/2, upperMargin+bsHeight, prew+1, preh+1);
    // draw thumbnail image
    if (preview) {
        prex = borderWidth + (exp_width-prew)/2;
        prey = upperMargin+bsHeight+borderWidth;
        backBuffer->draw_rgb_image (gc_, prex, prey, prew, preh, Gdk::RGB_DITHER_NONE, preview, prew*3);
    }

    customBackBufferUpdate (cr);
    
    // draw icons onto the thumbnail area
    bbIcons = getIconsOnImageArea ();

    int infow, infoh;
    getTextSizes (infow, infoh);

    int iofs_x = 4, iofs_y = 4;
    int igap = 2;
    int istartx = prex;
    int istarty = prey;

    if (options.showFileNames && options.overlayedFileNames) {
        cr->begin_new_path ();
        cr->rectangle (istartx, istarty, prew, fnlabh+dtlabh+exlabh+2*iofs_y);
        if ((texts.get_red_p()+texts.get_green_p()+texts.get_blue_p())/3 > 0.5)
            cr->set_source_rgba (0, 0, 0, 0.5);
        else
            cr->set_source_rgba (1, 1, 1, 0.5);
        cr->fill ();
    }

    istartx += iofs_x;
    istarty += iofs_y;
    
    if (bbIcons.size()>0) {
        int iwidth = igap;
        int iheight = 0;
        for (int i=0; i<bbIcons.size(); i++) {
            iwidth += bbIcons[i]->get_width() + igap;
            if (bbIcons[i]->get_height() > iheight)
                iheight = bbIcons[i]->get_height();
        }
        iheight += 2*igap;
        if (!options.overlayedFileNames) {
            cr->begin_new_path ();
            cr->rectangle (istartx-igap, istarty-igap, iwidth, iheight);
            cr->set_source_rgba (0, 0, 0, 0.75);
            cr->fill ();
        }
        for (int i=0; i<bbIcons.size(); i++) {
            backBuffer->draw_pixbuf (gc_, bbIcons[i], 0, 0, istartx, istarty, bbIcons[i]->get_width(), bbIcons[i]->get_height(), Gdk::RGB_DITHER_NONE, 0, 0);
            istartx += bbIcons[i]->get_width() + igap;
        }
    }

    if( options.showFileNames ){
		int textposx_fn, textposx_ex, textposx_dt, textposy, textw;
		if (!options.overlayedFileNames) {
			textposx_fn = exp_width/2 - fnlabw/2;
			if (textposx_fn<0) textposx_fn = 0;
			textposx_ex = exp_width/2 - exlabw/2;
			if (textposx_ex<0) textposx_ex = 0;
			textposx_dt = exp_width/2 - dtlabw/2;
			if (textposx_dt<0) textposx_dt = 0;
			textposy = upperMargin + bsHeight + 2*borderWidth + preh + borderWidth + textGap;
			textw = exp_width - 2*textGap;
			gc_->set_foreground (selected ? texts : textn);
		}
		else {
			textposx_fn = istartx;
			textposx_ex = istartx;
			textposx_dt = istartx;
			textposy = istarty;
			textw = prew - (istartx - prex);
			gc_->set_foreground (texts);
		}

		// draw file name
		Glib::RefPtr<Pango::Context> context = w->get_pango_context () ;
		Pango::FontDescription fontd = context->get_font_description ();
		fontd.set_weight (Pango::WEIGHT_BOLD);
		if (italicstyle)
			fontd.set_style (Pango::STYLE_ITALIC);
		else
			fontd.set_style (Pango::STYLE_NORMAL);

		context->set_font_description (fontd);
		Glib::RefPtr<Pango::Layout> fn = w->create_pango_layout (dispname);
		fn->set_width (textw*Pango::SCALE);
		fn->set_ellipsize (Pango::ELLIPSIZE_MIDDLE);
		backBuffer->draw_layout(gc_, textposx_fn, textposy, fn);

		fontd.set_weight (Pango::WEIGHT_NORMAL);
		fontd.set_style (Pango::STYLE_NORMAL);
		context->set_font_description (fontd);

		// draw date/time label
		int tpos = fnlabh;
		if (options.fbShowDateTime && datetimeline!="") {
			fn = w->create_pango_layout (datetimeline);
			fn->set_width (textw*Pango::SCALE);
			fn->set_ellipsize (Pango::ELLIPSIZE_MIDDLE);
			backBuffer->draw_layout(gc_, textposx_dt, textposy + tpos, fn);
			tpos += dtlabh;
		}
		// draw basic exif info
		if (options.fbShowBasicExif && exifline!="") {
			fn = w->create_pango_layout (exifline);
			fn->set_width (textw*Pango::SCALE);
			fn->set_ellipsize (Pango::ELLIPSIZE_MIDDLE);
			backBuffer->draw_layout (gc_, textposx_ex, textposy + tpos, fn);
			tpos += exlabh;
		}
    }
}

void ThumbBrowserEntryBase::getTextSizes (int& infow, int& infoh) {

    if (!parent)
        return;
        
    Gtk::Widget* w = parent->getDrawingArea ();

    // calculate dimensions of the text based fields
    dispname = shortname;

    Glib::RefPtr<Pango::Context> context = w->get_pango_context () ;
    context->set_font_description (w->get_style()->get_font());
   
    // filename:
    Pango::FontDescription fontd = context->get_font_description ();     
    fontd.set_weight (Pango::WEIGHT_BOLD);
    context->set_font_description (fontd);
    Glib::RefPtr<Pango::Layout> fn = w->create_pango_layout(shortname);
    fn->get_pixel_size (fnlabw, fnlabh);

    // datetime
    fontd.set_weight (Pango::WEIGHT_NORMAL);
    context->set_font_description (fontd);
    fn = w->create_pango_layout (datetimeline);
    fn->get_pixel_size (dtlabw, dtlabh);
    
    // basic exif data
    fn = w->create_pango_layout (exifline);
    fn->get_pixel_size (exlabw, exlabh);

    // calculate cummulated height of all info fields
    infoh = fnlabh;
    infow = 0;
    // add date/tile size:
    if (options.fbShowDateTime) {
        infoh += dtlabh;
        if (dtlabw + 2*sideMargin > infow)
            infow = dtlabw + 2*sideMargin;
    }
    if (options.fbShowBasicExif) {
        infoh += exlabh;
        if (exlabw + 2*sideMargin > infow)
            infow = exlabw + 2*sideMargin;
    }
}

void ThumbBrowserEntryBase::resize (int h) {
    Glib::RWLock::WriterLock l(lockRW);

    height = h;
    int old_preh = preh, old_width = width;

    // dimensions of the button set
    int bsw=0, bsh=0;
    if (buttonSet)
        buttonSet->getMinimalDimensions (bsw, bsh);
    
    // calculate the height remaining for the thumbnail image
    preh = height - upperMargin - 2*borderWidth - lowerMargin - bsh;
    int infow=0;
    int infoh=0;
    if (options.showFileNames && !options.overlayedFileNames) {
        // dimensions of the info text
        getTextSizes (infow, infoh);
        infoh += textGap;
        preh -= infoh;
    }
    // Minimum size for thumbs
    if (preh<24){
    	preh=24;
    	height = preh + (upperMargin + 2*borderWidth + lowerMargin)+ bsh + infoh;
    }

    calcThumbnailSize ();  // recalculates prew

    width = prew + 2*sideMargin + 2*borderWidth;
    if (options.showFileNames && !options.overlayedFileNames) {
        width = prew + 2*sideMargin + 2*borderWidth;
        if (width<infow + 2*sideMargin + 2*borderWidth)
            width = infow + 2*sideMargin + 2*borderWidth;
    }
    if (width < bsw + 2*sideMargin + 2*borderWidth)
        width = bsw + 2*sideMargin + 2*borderWidth;

    if ( preh != old_preh || width != old_width )
    {
        delete [] preview;
        preview = NULL;
        refreshThumbnailImage ();
    } // causes skewed thumb sometimes: else updateBackBuffer();

    drawable = true;
}

void ThumbBrowserEntryBase::drawFrame (Cairo::RefPtr<Cairo::Context> cr, const Gdk::Color& bg, const Gdk::Color& fg) {

    int radius = 8;
    if (selected || framed) {
        cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        cr->move_to (radius, 0);
        cr->arc (exp_width-1-radius, radius, radius, -M_PI/2, 0);
        cr->arc (exp_width-1-radius, exp_height-1-radius, radius, 0, M_PI/2);
        cr->arc (radius, exp_height-1-radius, radius, M_PI/2, M_PI);
        cr->arc (radius, radius, radius, M_PI, -M_PI/2);
        cr->close_path ();
        if (selected) {
            cr->set_source_rgb (bg.get_red_p(), bg.get_green_p(), bg.get_blue_p());
            cr->fill_preserve ();
        }
        cr->set_source_rgb (bg.get_red_p()*2/3, bg.get_green_p()*2/3, bg.get_blue_p()*2/3);
        cr->set_line_width (1.0);
        cr->stroke ();
        
    }

    if (framed) {
        cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        cr->move_to (+2+0.5+radius, +2+0.5);
        cr->arc (-2+0.5+exp_width-1-radius, +2+0.5+radius, radius, -M_PI/2, 0);
        cr->arc (-2+0.5+exp_width-1-radius, -2+0.5+exp_height-1-radius, radius, 0, M_PI/2);
        cr->arc (+2+0.5+radius, -2+exp_height-1-radius, radius, M_PI/2, M_PI);
        cr->arc (+2+0.5+radius, +2+radius, radius, M_PI, -M_PI/2);
        cr->close_path ();
        cr->set_source_rgb (fg.get_red_p(), fg.get_green_p(), fg.get_blue_p());
        cr->set_line_width (2.0);
        cr->stroke ();
    }   
}

void ThumbBrowserEntryBase::draw () {

    if (!drawable || !parent)
        return;

    Glib::RWLock::ReaderLock l(lockRW);  // No resizes, position moves etc. inbetween

    int bbWidth, bbHeight;
    if (backBuffer)
        backBuffer->get_size (bbWidth, bbHeight);

    if (!backBuffer || selected!=bbSelected || framed!=bbFramed || preview!=bbPreview 
        || exp_width!=bbWidth || exp_height!=bbHeight || getIconsOnImageArea ()!=bbIcons)
        updateBackBuffer ();

    Gtk::Widget* w = parent->getDrawingArea ();

    Glib::RefPtr<Gdk::GC> gc_ = Gdk::GC::create (w->get_window());

 //   Gdk::Color textn = w->get_style()->get_text(Gtk::STATE_NORMAL);
  //  Gdk::Color texts = w->get_style()->get_text(Gtk::STATE_SELECTED);
    Gdk::Color bgn = w->get_style()->get_bg(Gtk::STATE_NORMAL);
    Gdk::Color bgs = w->get_style()->get_bg(Gtk::STATE_SELECTED);   

    w->get_window()->draw_drawable (gc_, backBuffer, 0, 0, startx + ofsX, starty + ofsY);

    // check icon set changes!!!

//    drawProgressBar (window, gc_, selected ? texts : textn, selected ? bgs : bgn, ofsX+startx, exp_width, ofsY+starty + upperMargin+bsHeight+borderWidth+preh+borderWidth+textGap+tpos, fnlabh);

    // redraw button set above the thumbnail   
    if (buttonSet) {
        buttonSet->setColors (selected ? bgs : bgn, selected ? bgn : bgs);
        buttonSet->redraw (w->get_window()->create_cairo_context());
    }
}

void ThumbBrowserEntryBase::setPosition (int x, int y, int w, int h) {
    Glib::RWLock::WriterLock l(lockRW);

    exp_width = w;
    exp_height = h;
    startx = x;
    starty = y;
    
    if (buttonSet) 
        buttonSet->arrangeButtons (ofsX+x+sideMargin, ofsY+y+upperMargin, w-2*sideMargin, -1);
}

void ThumbBrowserEntryBase::setOffset (int x, int y) {
    Glib::RWLock::WriterLock l(lockRW);

    ofsX = -x;
    ofsY = -y;
    
    if (buttonSet) 
        buttonSet->move (ofsX+startx+sideMargin, ofsY+starty+upperMargin);
}

bool ThumbBrowserEntryBase::inside (int x, int y) {

    return x>ofsX+startx && x<ofsX+startx+exp_width && y>ofsY+starty && y<ofsY+starty+exp_height; 
} 

bool ThumbBrowserEntryBase::insideWindow (int x, int y, int w, int h) {

    return !(ofsX+startx>x+w || ofsX+startx+exp_width<x || ofsY+starty>y+h || ofsY+starty+exp_height<y); 
} 

bool ThumbBrowserEntryBase::motionNotify  (int x, int y) {

    return buttonSet ? buttonSet->motionNotify (x, y) : false;
}

bool ThumbBrowserEntryBase::pressNotify   (int button, int type, int bstate, int x, int y) {

    return buttonSet ? buttonSet->pressNotify (x, y) : false;
}

bool ThumbBrowserEntryBase::releaseNotify (int button, int type, int bstate, int x, int y) {

    return buttonSet ? buttonSet->releaseNotify (x, y) : false;
}

Glib::ustring ThumbBrowserEntryBase::getToolTip (int x, int y) {
    Glib::ustring tooltip = "";
    
    if (buttonSet) tooltip = buttonSet->getToolTip (x, y);

    // if the fileinfo is not shown anyway, make a tooltip with the info
    if (!options.showFileNames && inside(x,y) && tooltip.empty()) {
        tooltip = dispname;
        if (options.fbShowDateTime && datetimeline!="") tooltip += Glib::ustring("\n") + datetimeline;
        if (options.fbShowBasicExif && exifline!="") tooltip += Glib::ustring("\n") + exifline;
    }

    return tooltip;
}


