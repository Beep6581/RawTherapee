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
#include "histogrampanel.h"
#include "multilangmgr.h"
#include "guiutils.h"
#include "options.h"
#include <cstring>
#include "../rtengine/LUT.h"
#include "rtimage.h"
#include "../rtengine/improccoordinator.h"
#include "../rtengine/color.h"

extern Glib::ustring argv0;
extern Options options;


//
//
// HistogramPanel
HistogramPanel::HistogramPanel () {

   histogramArea = Gtk::manage (new HistogramArea (this));
   histogramRGBArea = Gtk::manage (new HistogramRGBArea ());
   histogramRGBArea->show();

   gfxVBox = Gtk::manage (new Gtk::VBox (false, 2));
   histogramRGBArea->setParent(gfxVBox);
   gfxVBox->pack_start (*histogramArea, Gtk::PACK_EXPAND_WIDGET, 0);
   if (options.histogramBar)
      gfxVBox->pack_start (*histogramRGBArea, Gtk::PACK_SHRINK, 0);

   redImage   = new RTImage ("histRed.png");
   greenImage = new RTImage ("histGreen.png");
   blueImage  = new RTImage ("histBlue.png");
   valueImage = new RTImage ("histValue.png");
   chroImage  = new RTImage ("histChro.png");
   rawImage   = new RTImage ("histRaw.png");
   fullImage  = new RTImage ("histFull.png");
   barImage   = new RTImage ("histBar.png");

   redImage_g   = new RTImage ("histRedg.png");
   greenImage_g = new RTImage ("histGreeng.png");
   blueImage_g  = new RTImage ("histBlueg.png");
   valueImage_g = new RTImage ("histValueg.png");
   chroImage_g  = new RTImage ("histChrog.png");
   rawImage_g   = new RTImage ("histRawg.png");
   fullImage_g  = new RTImage ("histFullg.png");
   barImage_g   = new RTImage ("histBarg.png");

   showRed   = Gtk::manage (new Gtk::ToggleButton ());
   showGreen = Gtk::manage (new Gtk::ToggleButton ());
   showBlue  = Gtk::manage (new Gtk::ToggleButton ());
   showValue = Gtk::manage (new Gtk::ToggleButton ());
   showChro  = Gtk::manage (new Gtk::ToggleButton ());
   showRAW   = Gtk::manage (new Gtk::ToggleButton ());
   showFull  = Gtk::manage (new Gtk::ToggleButton ());
   showBAR   = Gtk::manage (new Gtk::ToggleButton ());

   showRed->set_name("histButton");   showRed->set_can_focus(false);
   showGreen->set_name("histButton"); showGreen->set_can_focus(false);
   showBlue->set_name("histButton");  showBlue->set_can_focus(false);
   showValue->set_name("histButton"); showValue->set_can_focus(false);
   showChro->set_name("histButton");  showChro->set_can_focus(false);
   showRAW->set_name("histButton");   showRAW->set_can_focus(false);
   showFull->set_name("fullButton");  showFull->set_can_focus(false);
   showBAR->set_name("histButton");   showBAR->set_can_focus(false);

   showRed->set_relief (Gtk::RELIEF_NONE);
   showGreen->set_relief (Gtk::RELIEF_NONE);
   showBlue->set_relief (Gtk::RELIEF_NONE);
   showValue->set_relief (Gtk::RELIEF_NONE);
   showChro->set_relief (Gtk::RELIEF_NONE);
   showRAW->set_relief (Gtk::RELIEF_NONE);
   showFull->set_relief (Gtk::RELIEF_NONE);
   showBAR->set_relief (Gtk::RELIEF_NONE);

   showRed->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_R"));
   showGreen->set_tooltip_text (M("HISTOGRAM_TOOLTIP_G"));
   showBlue->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_B"));
   showValue->set_tooltip_text (M("HISTOGRAM_TOOLTIP_L"));
   showChro->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_CHRO"));
   showRAW->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_RAW"));
   showFull->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_FULL"));
   showBAR->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_BAR"));

   buttonVBox = Gtk::manage (new Gtk::VBox (false, 2));
   showRed->set_active (true);
   showGreen->set_active (true);
   showBlue->set_active (true);
   showValue->set_active (true);
   showChro->set_active (false);//unactive by default

   showRAW->set_active (false);
   showFull->set_active (!options.histogramFullMode);
   showBAR->set_active (options.histogramBar);

   showRed->set_image   (showRed->get_active()   ? *redImage   : *redImage_g);
   showGreen->set_image (showGreen->get_active() ? *greenImage : *greenImage_g);
   showBlue->set_image  (showBlue->get_active()  ? *blueImage  : *blueImage_g);
   showValue->set_image (showValue->get_active() ? *valueImage : *valueImage_g);
   showChro->set_image  (showChro->get_active()   ? *chroImage : *chroImage_g);
   showRAW->set_image   (showRAW->get_active()   ? *rawImage   : *rawImage_g);
   showFull->set_image  (showFull->get_active()  ? *fullImage  : *fullImage_g);
   showBAR->set_image   (showBAR->get_active()   ? *barImage   : *barImage_g);

   showRed->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::red_toggled), showRed );
   showGreen->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::green_toggled), showGreen );
   showBlue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::blue_toggled), showBlue );
   showValue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::value_toggled), showValue );
   showChro->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::chro_toggled), showChro );
   showRAW->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::raw_toggled), showRAW );
   showFull->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::full_toggled), showFull );
   showBAR->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::bar_toggled), showBAR );

   buttonVBox->pack_start (*showRed, Gtk::PACK_SHRINK, 0);
   buttonVBox->pack_start (*showGreen, Gtk::PACK_SHRINK, 0);
   buttonVBox->pack_start (*showBlue, Gtk::PACK_SHRINK, 0);
   buttonVBox->pack_start (*showValue, Gtk::PACK_SHRINK, 0);
   buttonVBox->pack_start (*showChro, Gtk::PACK_SHRINK, 0);
	buttonVBox->pack_start (*showRAW, Gtk::PACK_SHRINK, 0);
   buttonVBox->pack_start (*showFull, Gtk::PACK_SHRINK, 0);
   buttonVBox->pack_start (*showBAR, Gtk::PACK_SHRINK, 0);

   // Put the button vbox next to the window's border to be less disturbing
   if (options.histogramPosition == 1) {
      pack_start (*buttonVBox, Gtk::PACK_SHRINK, 2);
      pack_start (*gfxVBox,Gtk::PACK_EXPAND_WIDGET, 2);
   }
   else {
      pack_start (*gfxVBox,Gtk::PACK_EXPAND_WIDGET, 2);
      pack_start (*buttonVBox, Gtk::PACK_SHRINK, 2);
   }

   show_all ();

   rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &HistogramPanel::resized) );
}

HistogramPanel::~HistogramPanel () {
  delete redImage;
  delete greenImage;
  delete blueImage;
  delete valueImage;
  delete chroImage;
  delete rawImage;
  delete fullImage;
  delete barImage;

  delete redImage_g;
  delete greenImage_g;
  delete blueImage_g;
  delete valueImage_g;
  delete chroImage_g;
  delete rawImage_g;
  delete fullImage_g;
  delete barImage_g;
  
}

void HistogramPanel::resized (Gtk::Allocation& req) {

    rconn.block (true);

    int gHeight = req.get_width()/2;
    if (gHeight > 150) gHeight = 150; else if (gHeight < 100) gHeight = 100;
    int bHeight = req.get_width()/30;
    if (bHeight >  10) bHeight =  10; else if (bHeight < 5  ) bHeight = 5;
    histogramArea->set_size_request (req.get_width(), gHeight);
    histogramRGBArea->set_size_request (req.get_width(), bHeight);

    rconn.block (false);

    histogramArea->renderHistogram ();
    histogramArea->queue_draw ();

    if (histogramRGBArea->getFreeze()==true) {
        histogramRGBArea->updateFreeze(false);
        // set histogramRGBArea invalid;
        histogramRGBArea->renderRGBMarks(-1, -1, -1);
        // re-set freeze to old state
        histogramRGBArea->updateFreeze(true);
        histogramRGBArea->queue_draw ();
    }
    else {
        // set histogramRGBArea invalid;
        histogramRGBArea->renderRGBMarks(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    }
}

void HistogramPanel::red_toggled () {
    showRed->set_image(showRed->get_active() ? *redImage : *redImage_g);
    rgbv_toggled();
}
void HistogramPanel::green_toggled () {
    showGreen->set_image(showGreen->get_active() ? *greenImage : *greenImage_g);
    rgbv_toggled();
}
void HistogramPanel::blue_toggled () {
    showBlue->set_image(showBlue->get_active() ? *blueImage : *blueImage_g);
    rgbv_toggled();
}
void HistogramPanel::value_toggled () {
	removeIfThere(showValue, valueImage, false);
	removeIfThere(showValue, valueImage_g, false);
    showValue->set_image(showValue->get_active() ? *valueImage : *valueImage_g);
    rgbv_toggled();
}
void HistogramPanel::chro_toggled () {
	removeIfThere(showChro, chroImage, false);
	removeIfThere(showChro, chroImage_g, false);
    showChro->set_image(showChro->get_active() ? *chroImage : *chroImage_g);
    rgbv_toggled();
}

void HistogramPanel::raw_toggled () {
	if (showRAW->get_active()) {
      showRAW->set_image(*rawImage);
      showValue->set_sensitive(false);
      showChro->set_sensitive(false);
	}
	else {
      showRAW->set_image(*rawImage_g);
      showValue->set_sensitive(true);
      showChro->set_sensitive(true);
	}
    rgbv_toggled();
}
void HistogramPanel::full_toggled () {
    options.histogramFullMode = !showFull->get_active();
    showFull->set_image(showFull->get_active() ? *fullImage : *fullImage_g);
    rgbv_toggled();
}
void HistogramPanel::bar_toggled () {
    showBAR->set_image(showBAR->get_active() ? *barImage : *barImage_g);
    rgbv_toggled();
}
void HistogramPanel::rgbv_toggled () {
  // Update Display
  histogramArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active(), showFull->get_active(), showChro->get_active());
  histogramArea->queue_draw ();

  histogramRGBArea->updateOptions (showRed->get_active(), showGreen->get_active(), showBlue->get_active(), showValue->get_active(), showRAW->get_active(), showBAR->get_active(), showChro->get_active());
  histogramRGBArea->renderRGBMarks (0,0,0);
  histogramArea->queue_draw ();
}

void HistogramPanel::setHistRGBInvalid () {
    // do something to un-show vertical bars
    histogramRGBArea->renderRGBMarks(-1, -1, -1);
    histogramRGBArea->queue_draw ();
}

// "Freeze" is not a button, but a RMB-click, so this is not in the RGBV-Toggle method 
void HistogramPanel::toggleFreeze () {
    if (histogramRGBArea->getFreeze()==true) {
        histogramRGBArea->updateFreeze(false);
    }
    else if (histogramRGBArea->getShow()==true) {
        histogramRGBArea->updateFreeze(true);
    }
    return;
}

void HistogramPanel::pointerMoved (bool validPos, Glib::ustring profile, int x, int y, int r, int g, int b) {

    if (!validPos) {
        // do something to un-show vertical bars
        histogramRGBArea->renderRGBMarks(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    }
    else {
        // do something to show vertical bars
        histogramRGBArea->renderRGBMarks(r, g, b, profile);
        histogramRGBArea->queue_draw ();
    }
}

/*
 * Move the vertical button bar to the right side
 * only allowed values for align are Gtk::ALIGN_LEFT and Gtk::ALIGN_RIGHT
 */
void HistogramPanel::reorder (Gtk::AlignmentEnum align) {
	if (align == Gtk::ALIGN_LEFT)
		reorder_child(*buttonVBox, 0);
	else
		reorder_child(*buttonVBox, 1);
}

// FullModeListener interface:
void HistogramPanel::toggle_button_full () {
    showFull->set_active (!showFull->get_active ());
    showFull->set_image(showFull->get_active() ? *fullImage : *fullImage_g);
}

//
//
//
// HistogramRGBArea
HistogramRGBArea::HistogramRGBArea () ://needChroma unactive by default
  frozen(false), valid(false), needRed(true), needGreen(true), needBlue(true), needLuma(true), rawMode(false), showMode(options.histogramBar), barDisplayed(options.histogramBar), needChroma(false){

  harih = new HistogramRGBAreaIdleHelper;
  harih->harea = this;
  harih->destroyed = false;
  harih->pending = 0;
}

HistogramRGBArea::~HistogramRGBArea () {

    if (harih->pending)
         harih->destroyed = true;
    else
         delete harih;
}

bool HistogramRGBArea::getFreeze() {
    return(frozen);
}

bool HistogramRGBArea::getShow() {
    return(showMode);
}

void HistogramRGBArea::updateFreeze (bool f) {
    frozen = f;
    return;
}

void HistogramRGBArea::renderRGBMarks (int r, int g, int b, Glib::ustring profile) {

  if (!is_realized ())
        return;
  
  if (frozen)  {
        return;
  }

  // Mostly not necessary, but should be in some case
  GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected

  Glib::RefPtr<Gdk::Window> window = get_window();
  int winx, winy, winw, winh, wind;
  window->get_geometry(winx, winy, winw, winh, wind);

  overlay = Gdk::Pixmap::create (window, winw, winh, -1);
  Glib::RefPtr<Gdk::GC> ovrl = Gdk::GC::create(overlay);

  Glib::RefPtr<Gtk::Style> style = get_style ();

  if (!showMode)  {
    ovrl->set_foreground (style->get_bg (Gtk::STATE_NORMAL));
    overlay->draw_rectangle (ovrl, true, 0, 0, winw, winh);
    if (rgbgc_ && overlay) {
      window->draw_drawable (rgbgc_, overlay, 0, 0, 0, 0, -1, -1); }
    return; }
  else {
    ovrl->set_foreground (mgray);
    overlay->draw_rectangle (ovrl, true, 0, 0, winw, winh);
    if (rgbgc_ && overlay) {
      window->draw_drawable (rgbgc_, overlay, 0, 0, 0, 0, -1, -1); }
  }

  Cairo::RefPtr<Cairo::Context> cr = overlay->create_cairo_context();
  cr->set_line_width (1.0);

  if ( r != -1 && g != -1 && b != -1 ) {
    if (needRed) {
      // Red
      cr->set_source_rgb(1.0, 0.0, 0.0);
      cr->move_to((int)(r*(winw/256.0)), 0);
      cr->line_to((int)(r*(winw/256.0)), winh-0);
      cr->stroke();
    }
    if (needGreen) {
      // Green
      cr->set_source_rgb(0.0, 1.0, 0.0);
      cr->move_to((int)(g*(winw/256.0)), 0);
      cr->line_to((int)(g*(winw/256.0)), winh-0);
      cr->stroke();
  }
    if (needBlue) {
      // Blue
      cr->set_source_rgb(0.0, 0.0, 1.0);
      cr->move_to((int)(b*(winw/256.0)), 0);
      cr->line_to((int)(b*(winw/256.0)), winh-0);
      cr->stroke();
  }
  if(needLuma || needChroma) {
	float Lab_L,Lab_a,Lab_b;
	rgb2lab( profile, r,g,b,Lab_L,Lab_a,Lab_b);
    if (needLuma) {
      // Luma
      cr->set_source_rgb(1.0, 1.0, 1.0);
      cr->move_to((int)((Lab_L)*(winw/100.0)), 0);
      cr->line_to((int)((Lab_L)*(winw/100.0)), winh-0);
      cr->stroke();
    }
    if (needChroma) {
      // Chroma
      float chromaval = sqrt(Lab_a*Lab_a + Lab_b*Lab_b)/1.8;
      cr->set_source_rgb(0.0, 0.0, 0.0);
      cr->move_to((int)(chromaval*(winw/100.0)), 0);
      cr->line_to((int)(chromaval*(winw/100.0)), winh-0);
      cr->stroke();
    }
  }
  }
}

void HistogramRGBArea::rgb2lab (Glib::ustring profile, int r, int g, int b, float &LAB_l, float &LAB_a, float &LAB_b) {
	double xyz_rgb[3][3];
	const double ep=216.0/24389.0;
	const double ka=24389.0/27.0;

	double var_R = r / 255.0;
	double var_G = g / 255.0;
	double var_B = b / 255.0;

	if (profile=="sRGB") {//apply sRGB inverse gamma

    // 
// if you want display = working space
// today as the gamma output can not be configured
// it is better that the user has the gamma of the output space
	if ( var_R > 0.04045 ) 
		var_R = pow ( ( ( var_R + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
	else    
		var_R = var_R / 12.92;
	if ( var_G > 0.04045 ) 
		var_G = pow ( ( ( var_G + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
	else                   
		var_G = var_G / 12.92;
	if ( var_B > 0.04045 ) 
		var_B = pow ( ( ( var_B + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
	else    
		var_B = var_B / 12.92;
	} 
// if you want display = output space
	else 
	if (profile=="ProPhoto") {// apply inverse gamma 1.8
		var_R = pow ( var_R, 1.8);
		var_G = pow ( var_G, 1.8);
		var_B = pow ( var_B, 1.8);
	}
	else {// apply inverse gamma 2.2
		var_R = pow ( var_R, 2.2);
		var_G = pow ( var_G, 2.2);
		var_B = pow ( var_B, 2.2);
	}

	/*for (int i=0; i<numprofiles; i++) {
		if (profile==wpnames[i]) {
			for (int m=0; m<3; m++) 
				for (int n=0; n<3; n++) {
					xyz_rgb[m][n] = wprofiles[i][m][n];
			}
			break;
		}
	}*/

    TMatrix wprof = rtengine::ICCStore::getInstance()->workingSpaceMatrix (profile);

	for (int m=0; m<3; m++) 
		for (int n=0; n<3; n++) {
			xyz_rgb[m][n] = wprof[m][n];
		}

	double varxx,varyy,varzz;
	double var_X = ( xyz_rgb[0][0]*var_R + xyz_rgb[0][1]*var_G + xyz_rgb[0][2]*var_B ) / rtengine::Color::D50x;
	double var_Y = ( xyz_rgb[1][0]*var_R + xyz_rgb[1][1]*var_G + xyz_rgb[1][2]*var_B ) ;
	double var_Z = ( xyz_rgb[2][0]*var_R + xyz_rgb[2][1]*var_G + xyz_rgb[2][2]*var_B ) / rtengine::Color::D50z;

	varxx = var_X>ep?cbrt(var_X):( ka * var_X  +  16.0) / 116.0 ;
	varyy = var_Y>ep?cbrt(var_Y):( ka * var_Y  +  16.0) / 116.0 ;
	varzz = var_Z>ep?cbrt(var_Z):( ka * var_Z  +  16.0) / 116.0 ;
	LAB_l = ( 116 * varyy ) - 16;
	LAB_a = 500 * ( varxx - varyy );
	LAB_b = 200 * ( varyy - varzz );

}


int histrgbupdate (void* data) {

    HistogramRGBAreaIdleHelper* harih = static_cast<HistogramRGBAreaIdleHelper*>(data);

    if (harih->destroyed) {
        if (harih->pending == 1)
            delete harih;
        else    
            harih->pending--;
        return 0;
    }
    
    harih->harea->renderRGBMarks(-1,-1,-1);
    harih->harea->queue_draw ();

    harih->pending--;

    return 0;
}

void HistogramRGBArea::update (int valh, int rh, int  gh, int bh) {
	
    if (valh) {
        val=valh;
        r=rh;
        g=gh;
        b=bh;
        valid = true;
    }
    else
        valid = false;

    harih->pending++;
    g_idle_add (histrgbupdate, harih);
}

void HistogramRGBArea::updateOptions (bool r, bool g, bool b, bool l, bool raw, bool bar, bool c) {

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needLuma  = l;
    rawMode   = raw;
    showMode  = bar;
	needChroma = c;

    // Histogram RGB BAR button logic goes here

    if (bar && !barDisplayed) {
        // Toggled on, add (show) the widget
    	parent->pack_start(*this, Gtk::PACK_SHRINK, 0);
        options.histogramBar = true;
        barDisplayed = true;
    }
    else if (!bar && barDisplayed){
        // Toggled off, remove (hide) the widget
    	removeIfThere(parent, this, false);
        options.histogramBar = false;
        barDisplayed = false;
        // unfreeze
        updateFreeze(false);
    }

    // Disable (but don't hide it) the bar button when RAW histogram is displayed
    if (rawMode) {
        showMode = false;
   }
}

void HistogramRGBArea::on_realize () {

  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  rgbgc_ = Gdk::GC::create(window);
  add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK);

  Glib::RefPtr<Gdk::Colormap> rgbcolormap = get_default_colormap();
  black = Gdk::Color ("black");
  red = Gdk::Color ("red");
  green = Gdk::Color ("green");
  blue = Gdk::Color ("blue");
  lgray = Gdk::Color ("gray75");
  mgray = Gdk::Color ("gray50");
  dgray = Gdk::Color ("gray25");
  rgbcolormap->alloc_color(black);
  rgbcolormap->alloc_color(white);
  rgbcolormap->alloc_color(red);
  rgbcolormap->alloc_color(green);
  rgbcolormap->alloc_color(blue);
  rgbcolormap->alloc_color(lgray);
  rgbcolormap->alloc_color(mgray);
  rgbcolormap->alloc_color(dgray);

}

bool HistogramRGBArea::on_expose_event(GdkEventExpose* event) {

  Glib::RefPtr<Gdk::Window> window = get_window();

  // on_realize & RenderRGBMarks have to be called before 
  if (rgbgc_ && overlay) {
      window->draw_drawable (rgbgc_, overlay, 0, 0, 0, 0, -1, -1);
  }

  return true;
}

bool HistogramRGBArea::on_button_press_event (GdkEventButton* event) {

  if (event->type==GDK_2BUTTON_PRESS && event->button==1) {
  // do something. Maybe un-freeze ? 
  }
  return true;
}

void HistogramRGBArea::on_style_changed (const Glib::RefPtr<Gtk::Style>& style) {

  white = get_style()->get_base(Gtk::STATE_NORMAL);
  queue_draw ();
}

//
//
//
// HistogramArea
HistogramArea::HistogramArea (FullModeListener *fml) : //needChroma unactive by default
    valid(false), fullMode(options.histogramFullMode), myFullModeListener(fml), oldwidth(-1), oldheight(-1), needLuma(true), needRed(true), needGreen(true), needBlue(true), rawMode(false), needChroma(false) {

    lhist(256);
    rhist(256);
    ghist(256);
    bhist(256);
    chist(256);

    haih = new HistogramAreaIdleHelper;
    haih->harea = this;
    haih->destroyed = false;
    haih->pending = 0;
}

HistogramArea::~HistogramArea () {

    if (haih->pending)
        haih->destroyed = true;
    else
        delete haih;

}

void HistogramArea::updateOptions (bool r, bool g, bool b, bool l, bool raw, bool full, bool c) {

    needRed   = r;
    needGreen = g;
    needBlue  = b;
    needLuma  = l;
    rawMode   = raw;
    fullMode  = !full;
	needChroma = c;

    renderHistogram ();
}

int histupdateUI (void* data) {

    HistogramAreaIdleHelper* haih = static_cast<HistogramAreaIdleHelper*>(data);

    if (haih->destroyed) {
        if (haih->pending == 1)
            delete haih;
        else    
            haih->pending--;

        return 0;
    }
    
    haih->harea->renderHistogram ();
    haih->harea->queue_draw ();

    haih->pending--;

    return 0;
}

void HistogramArea::update (LUTu &histRed, LUTu &histGreen, LUTu &histBlue, LUTu &histLuma, LUTu &histRedRaw, LUTu &histGreenRaw, LUTu &histBlueRaw, LUTu &histChroma) {
	
    if (histRed) {
        lhist=histLuma;chist=histChroma;
        rhist=histRed; ghist=histGreen; bhist=histBlue;
        rhistRaw=histRedRaw; ghistRaw=histGreenRaw; bhistRaw=histBlueRaw;

        valid = true;
    }
    else
        valid = false;

    haih->pending++;
    // Can be done outside of the GUI thread
    g_idle_add (histupdateUI, haih);
}

void HistogramArea::renderHistogram () {

  if (!is_realized ())
        return;

  Glib::RefPtr<Gdk::Window> window = get_window();
  int winx, winy, winw, winh, wind;
  window->get_geometry(winx, winy, winw, winh, wind);

  backBuffer = Gdk::Pixmap::create (window, winw, winh, -1);

  Glib::RefPtr<Gdk::GC> bgc = Gdk::GC::create(backBuffer);
  
  bgc->set_foreground (white);
  backBuffer->draw_rectangle (bgc, true, 0, 0, winw, winh);

  if (valid) {
    // For RAW mode use the other hists
    LUTu& rh = rawMode ? rhistRaw : rhist;
    LUTu& gh = rawMode ? ghistRaw : ghist;
    LUTu& bh = rawMode ? bhistRaw : bhist;

	// make double copies of LUT, one for faster access, another one to scale down the raw histos
	LUTu rhchanged(256),ghchanged(256),bhchanged(256);
	unsigned int lhisttemp[256],chisttemp[256],rhtemp[256],ghtemp[256],bhtemp[256];
	const int scale = (rawMode ? 8 : 1);
	for(int i=0;i<256;i++) {
		if(needLuma)
			lhisttemp[i] = lhist[i];
		if(needChroma)
			chisttemp[i] = chist[i];
		if(needRed)
			rhchanged[i] = rhtemp[i] = rh[i] / scale;
		if(needGreen)
			ghchanged[i] = ghtemp[i] = gh[i] / scale;
		if(needBlue)
			bhchanged[i] = bhtemp[i] = bh[i] / scale;
	}

    // compute height of the full histogram (realheight) and
    // does not take into account 0 and 255 values
    // them are handled separately

    int fullhistheight = 0;
    for (int i=1; i<255; i++) {
        if (needLuma && lhisttemp[i]>fullhistheight)
	        fullhistheight = lhisttemp[i];
        if (needChroma && chisttemp[i]>fullhistheight)
	        fullhistheight = chisttemp[i];
        if (needRed && rhtemp[i]>fullhistheight)
	        fullhistheight = rhtemp[i];
        if (needGreen && ghtemp[i]>fullhistheight)
	        fullhistheight = ghtemp[i];
        if (needBlue && bhtemp[i]>fullhistheight)
	        fullhistheight = bhtemp[i];
    }

    int realhistheight = fullhistheight;
    
    // though much faster than before, this still takes a lot of time especially for big files if rawMode is true
    if (!fullMode) {
        int area = 0;
        if(!rawMode)
			for (int i=0; i<fullhistheight; i++) {
				for (int j=0; j<256; j++)
					if ((needLuma && lhisttemp[j]>i) || (needChroma && chisttemp[j]>i) || (needRed && rhtemp[j]>i) || (needGreen && ghtemp[j]>i) || (needBlue && bhtemp[j]>i))
						area++;
				if ((double)area / (256*(i+1)) < 0.3) {
					realhistheight = i;
					break;
				}
			}
		else
			for (int i=0; i<fullhistheight; i++) {
				for (int j=0; j<256; j++)
					if ((needRed && rhtemp[j]>i) || (needGreen && ghtemp[j]>i) || (needBlue && bhtemp[j]>i))
						area++;
				if ((double)area / (256*(i+1)) < 0.3) {
					realhistheight = i;
					break;
				}
			}
			
    }

    if (realhistheight<winh-2)
        realhistheight = winh-2;

    Cairo::RefPtr<Cairo::Context> cr = backBuffer->create_cairo_context();
    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_width (1.0);

    int ui = 0, oi = 0;

    if (needLuma && !rawMode) {
        drawCurve(cr, lhist, realhistheight, winw, winh);
        cr->set_source_rgb (0.75, 0.75, 0.75);
		cr->fill_preserve ();
        cr->set_source_rgb (0.5, 0.5, 0.5);
		cr->stroke ();

        drawMarks(cr, lhist, realhistheight, winw, ui, oi);
    }
    if (needChroma && !rawMode) {
        drawCurve(cr, chist, realhistheight, winw, winh);
        cr->set_source_rgb (0.6, 0.6, 0.6);
	//	cr->fill_preserve ();
     //   cr->set_source_rgb (0.2, 0.2, 0.1);
		cr->stroke ();

        drawMarks(cr, chist, realhistheight, winw, ui, oi);
    }
    if (needRed) {
        drawCurve(cr, rhchanged, realhistheight, winw, winh);
        cr->set_source_rgb (1.0, 0.0, 0.0);
    	cr->stroke ();

        drawMarks(cr, rhchanged, realhistheight, winw, ui, oi);
    }
    if (needGreen) {
        drawCurve(cr, ghchanged, realhistheight, winw, winh);
        cr->set_source_rgb (0.0, 1.0, 0.0);
    	cr->stroke ();

        drawMarks(cr, ghchanged, realhistheight, winw, ui, oi);
    }
    if (needBlue) {
        drawCurve(cr, bhchanged, realhistheight, winw, winh);
        cr->set_source_rgb (0.0, 0.0, 1.0);
    	cr->stroke ();

        drawMarks(cr, bhchanged, realhistheight, winw, ui, oi);
    }
  }

  bgc->set_foreground (mgray);
  backBuffer->draw_rectangle (bgc, false, 0, 0, winw-1, winh-1);

  bgc->set_line_attributes (1, Gdk::LINE_ON_OFF_DASH, Gdk::CAP_NOT_LAST, Gdk::JOIN_MITER);

  backBuffer->draw_line (bgc, winw/4, 0, winw/4, winh);
  backBuffer->draw_line (bgc, 2*winw/4, 0, 2*winw/4, winh);
  backBuffer->draw_line (bgc, 3*winw/4, 0, 3*winw/4, winh);
  backBuffer->draw_line (bgc, 0, winh/4, winw, winh/4);
  backBuffer->draw_line (bgc, 0, 2*winh/4, winw, 2*winh/4);
  backBuffer->draw_line (bgc, 0, 3*winh/4, winw, 3*winh/4);

  bgc->set_line_attributes (1, Gdk::LINE_SOLID, Gdk::CAP_NOT_LAST, Gdk::JOIN_MITER);

  oldwidth = winw;
  oldheight = winh;
}

void HistogramArea::on_realize () {

  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  gc_ = Gdk::GC::create(window);
  add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK);
  Glib::RefPtr<Gdk::Colormap> colormap = get_default_colormap();

  black = Gdk::Color ("black");
  red = Gdk::Color ("red");
  green = Gdk::Color ("green");
  blue = Gdk::Color ("blue");
  lgray = Gdk::Color ("gray75");
  mgray = Gdk::Color ("gray50");
  dgray = Gdk::Color ("gray25");
  colormap->alloc_color(black);
  colormap->alloc_color(white);
  colormap->alloc_color(red);
  colormap->alloc_color(green);
  colormap->alloc_color(blue);
  colormap->alloc_color(lgray);
  colormap->alloc_color(mgray);
  colormap->alloc_color(dgray);

}

void HistogramArea::drawCurve(Cairo::RefPtr<Cairo::Context> &cr,
    LUTu & data, double scale, int hsize, int vsize)
{
    cr->move_to (0, vsize-1);
    for (int i = 0; i < 256; i++) {
        double val = data[i] * (double)(vsize-2) / scale;
        if (val > vsize - 1)
            val = vsize - 1;
        cr->line_to ((i/255.0)*(hsize - 1), vsize - 1 - val);
	}
    cr->line_to (hsize - 1, vsize - 1);
}

void HistogramArea::drawMarks(Cairo::RefPtr<Cairo::Context> &cr,
    LUTu & data, double scale, int hsize, int & ui, int & oi)
{
    int s = 8;
    
    if(data[0] > scale) {
            cr->rectangle(0, (ui++)*s, s, s);
    }
    if(data[255] > scale) {
            cr->rectangle(hsize - s, (oi++)*s, s, s);
    }
    cr->fill();
}

void HistogramArea::on_style_changed (const Glib::RefPtr<Gtk::Style>& style) {

    white = get_style()->get_base(Gtk::STATE_NORMAL);
    queue_draw ();
}

bool HistogramArea::on_expose_event(GdkEventExpose* event) {

    Glib::RefPtr<Gdk::Window> window = get_window();

    int winx, winy, winw, winh, wind;
    window->get_geometry(winx, winy, winw, winh, wind);

    if (winw!=oldwidth && winh!=oldheight)
    renderHistogram ();
    window->draw_drawable (gc_, backBuffer, 0, 0, 0, 0, -1, -1);

    return true;
}

bool HistogramArea::on_button_press_event (GdkEventButton* event) {

    if (event->type==GDK_2BUTTON_PRESS && event->button==1) {
        fullMode = !fullMode;
        options.histogramFullMode = fullMode;
        if (myFullModeListener)
            myFullModeListener->toggle_button_full ();
        renderHistogram ();
        queue_draw ();
    }
    return true;
}

