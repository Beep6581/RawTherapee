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
#include <crop.h>
#include <options.h>
#include <guiutils.h>
using namespace rtengine;
using namespace rtengine::procparams;

extern Glib::ustring argv0;
extern Options options;

class RefreshSpinHelper {

    public:
        Crop* crop;
        bool  notify;
        RefreshSpinHelper (Crop* _crop, bool _notify)
            : crop(_crop), notify(_notify) {}
};

Crop::Crop () {

  clistener = NULL;

  maxw = 3000;
  maxh = 2000;

  enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
  enabled->set_active (false);
  pack_start(*enabled);

  pack_start(*Gtk::manage (new  Gtk::HSeparator()));

  Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());

  hb1->pack_start (*Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("TP_CROP_X") +": ")));
  x = Gtk::manage (new Gtk::SpinButton ());
  x->set_size_request (60, -1);
  hb1->pack_start (*x);

  hb1->pack_start (*Gtk::manage (new Gtk::Label (Glib::ustring("   ") + M("TP_CROP_Y") + ": ")));
  y = Gtk::manage (new Gtk::SpinButton ());
  y->set_size_request (60, -1);
  hb1->pack_start (*y);

  pack_start (*hb1, Gtk::PACK_SHRINK, 2);

  Gtk::HBox* hb2 = Gtk::manage (new Gtk::HBox ());

  hb2->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_W") + ": ")));
  w = Gtk::manage (new Gtk::SpinButton ());
  w->set_size_request (60, -1);
  hb2->pack_start (*w);

  hb2->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_H") + ": ")));
  h = Gtk::manage (new Gtk::SpinButton ());
  h->set_size_request (60, -1);
  hb2->pack_start (*h);
  
  pack_start (*hb2, Gtk::PACK_SHRINK, 4);

  selectCrop = Gtk::manage (new Gtk::Button (M("TP_CROP_SELECTCROP")));
  selectCrop->set_image (*Gtk::manage (new Gtk::Image (argv0+"/images/crop16.png")));
  
  pack_start (*selectCrop, Gtk::PACK_SHRINK, 2);

  Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());

  fixr = Gtk::manage (new Gtk::CheckButton (M("TP_CROP_FIXRATIO")));
  fixr->set_active (1);

  hb3->pack_start (*fixr, Gtk::PACK_SHRINK, 4);

  ratio = Gtk::manage (new Gtk::ComboBoxText ());
  hb3->pack_start (*ratio, Gtk::PACK_SHRINK, 4);

  orientation = Gtk::manage (new Gtk::ComboBoxText ());
  hb3->pack_start (*orientation);

  pack_start (*hb3, Gtk::PACK_SHRINK, 4);

  Gtk::HBox* hb31 = Gtk::manage (new Gtk::HBox ());

  hb31->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_GUIDETYPE"))), Gtk::PACK_SHRINK, 4);
  guide = Gtk::manage (new Gtk::ComboBoxText ());
  hb31->pack_start (*guide);

  pack_start (*hb31, Gtk::PACK_SHRINK, 4);
  
  dpibox = Gtk::manage (new Gtk::VBox());
  dpibox->pack_start (*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_SHRINK, 2);

  Gtk::HBox* hb4 = Gtk::manage (new Gtk::HBox ());
  hb4->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_DPI"))));
  dpi = Gtk::manage (new Gtk::SpinButton ());
  dpi->set_size_request (60, -1);
  hb4->pack_start (*dpi);

  dpibox->pack_start (*hb4, Gtk::PACK_SHRINK, 2);

  sizecm = Gtk::manage (new Gtk::Label (M("GENERAL_NA") + " cm x " + M("GENERAL_NA") + " cm"));
  sizein = Gtk::manage (new Gtk::Label (M("GENERAL_NA") + " in x " + M("GENERAL_NA") + " in"));

  dpibox->pack_start (*sizecm, Gtk::PACK_SHRINK, 1);
  dpibox->pack_start (*sizein, Gtk::PACK_SHRINK, 1);
  
  pack_start (*dpibox, Gtk::PACK_SHRINK, 0);

  dpi->set_value (300);

  ratio->append_text ("3:2");
  ratio->append_text ("4:3");
  ratio->append_text ("16:9");
  ratio->append_text ("16:10");
  ratio->append_text ("5:4");
  ratio->append_text ("2:1");
  ratio->append_text ("1:1");
  ratio->append_text ("DIN");
  ratio->set_active (0);

  orientation->append_text (M("GENERAL_LANDSCAPE"));
  orientation->append_text (M("GENERAL_PORTRAIT"));
  orientation->set_active (0);

  guide->append_text (M("TP_CROP_GTNONE"));
  guide->append_text (M("TP_CROP_GTRULETHIRDS"));
  guide->append_text (M("TP_CROP_GTDIAGONALS"));
  guide->append_text (M("TP_CROP_GTHARMMEANS1"));
  guide->append_text (M("TP_CROP_GTHARMMEANS2"));
  guide->append_text (M("TP_CROP_GTHARMMEANS3"));
  guide->append_text (M("TP_CROP_GTHARMMEANS4"));
  guide->set_active (0);

  w->set_range (0, maxw);
  h->set_range (0, maxh);
  x->set_range (0, maxw);
  y->set_range (0, maxh);

  x->set_digits (0);
  x->set_increments (1,100);
  x->set_value (0);

  y->set_digits (0);
  y->set_increments (1,100);
  y->set_value (0);

  w->set_digits (0);
  w->set_increments (1,100);
  w->set_value (200);

  h->set_digits (0);
  h->set_increments (1,100);
  h->set_value (200);

  dpi->set_digits (0);
  dpi->set_increments (1,100);
  dpi->set_range (50, 12000);
  dpi->set_value (300);

  xconn = x->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::positionChanged), true);
  yconn = y->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::positionChanged), true);
  wconn = w->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::widthChanged), true);
  hconn = h->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::heightChanged), true);
  econn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Crop::enabledChanged) );
  fconn = fixr->signal_toggled().connect( sigc::mem_fun(*this, &Crop::ratioChanged) );
  rconn = ratio->signal_changed().connect( sigc::mem_fun(*this, &Crop::ratioChanged) );
  oconn = orientation->signal_changed().connect( sigc::mem_fun(*this, &Crop::ratioChanged) );
  gconn = guide->signal_changed().connect( sigc::mem_fun(*this, &Crop::notifyListener) );
  selectCrop->signal_pressed().connect( sigc::mem_fun(*this, &Crop::selectPressed) );
  dpi->signal_value_changed().connect( sigc::mem_fun(*this, &Crop::refreshSize) );

  nx = ny = nw = nh = 0;
  nsx = nsy = nsw = nsh = 0;
  lastScale = 1.0;
  lastRotationDeg = 0;
  show_all ();
}

void Crop::writeOptions () {

    options.cropDPI = (int)dpi->get_value ();
}

void Crop::readOptions () {

    disableListener ();

    dpi->set_value (options.cropDPI);

    enableListener ();
}

void Crop::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    
    xconn.block (true);
    yconn.block (true);
    wconn.block (true);
    hconn.block (true);
    rconn.block (true);
    fconn.block (true);
    oconn.block (true);
    gconn.block (true);
    enabled->set_active (pp->crop.enabled);

    // check if the new values are larger than the maximum
    double tmp, maxw, maxh;
    w->get_range (tmp, maxw);
    h->get_range (tmp, maxh);
    if (pp->crop.x + pp->crop.w > maxw || pp->crop.y + pp->crop.h > maxh)
        setDimensions (pp->crop.x + pp->crop.w, pp->crop.y + pp->crop.h);

    ratio->set_active_text (pp->crop.ratio);
    fixr->set_active (pp->crop.fixratio);

    if (pp->crop.orientation == "Landscape")
        orientation->set_active (0);
    else if (pp->crop.orientation == "Portrait")
        orientation->set_active (1);
    
    if (pp->crop.guide == "None")
        guide->set_active (0);
    else if (pp->crop.guide == "Rule of thirds")
        guide->set_active (1);
    else if (pp->crop.guide == "Rule of diagonals")
        guide->set_active (2);
    else if (pp->crop.guide == "Harmonic means 1")
        guide->set_active (3);
    else if (pp->crop.guide == "Harmonic means 2")
        guide->set_active (4);
    else if (pp->crop.guide == "Harmonic means 3")
        guide->set_active (5);
    else if (pp->crop.guide == "Harmonic means 4")
        guide->set_active (6);

    x->set_value (pp->crop.x);
    y->set_value (pp->crop.y);
    w->set_value (pp->crop.w);
    h->set_value (pp->crop.h);

    nx = pp->crop.x;
    ny = pp->crop.y;
    nw = pp->crop.w;
    nh = pp->crop.h;
    
    if (pp->resize.enabled)
        lastScale = pp->resize.scale;
    else
        lastScale = 1.0;
    nsx = nx / lastScale;
    nsy = ny / lastScale;
    nsw = nw / lastScale;
    nsh = nh / lastScale;
    lastRotationDeg = pp->coarse.rotate;

    wDirty = false;
    hDirty = false;
    xDirty = false;
    yDirty = false;

    if (pedited) {
        wDirty = pedited->crop.w;
        hDirty = pedited->crop.h;
        xDirty = pedited->crop.x;
        yDirty = pedited->crop.y;
        if (!pedited->crop.ratio)
            ratio->set_active (8);
        if (!pedited->crop.orientation) 
            orientation->set_active (2);
        if (!pedited->crop.guide) 
            guide->set_active (7);
        enabled->set_inconsistent (!pedited->crop.enabled);
        fixr->set_inconsistent (!pedited->crop.fixratio);
    }

    lastEnabled = pp->crop.enabled;
    lastAspect  = pp->crop.fixratio;

    xconn.block (false);
    yconn.block (false);
    wconn.block (false);
    hconn.block (false);
    rconn.block (false);
    fconn.block (false);
    oconn.block (false);
    gconn.block (false);

    enableListener ();
}

void Crop::write (ProcParams* pp, ParamsEdited* pedited) {

  pp->crop.enabled = enabled->get_active ();
  pp->crop.x = nx;
  pp->crop.y = ny;
  pp->crop.w = nw;
  pp->crop.h = nh;
  pp->crop.fixratio = fixr->get_active ();
  pp->crop.ratio = ratio->get_active_text ();

  if (orientation->get_active_row_number()==0)
      pp->crop.orientation = "Landscape";
  else if (orientation->get_active_row_number()==1)
      pp->crop.orientation = "Portrait";

  if (guide->get_active_row_number()==0)
    pp->crop.guide = "None";
  else if (guide->get_active_row_number()==1)
    pp->crop.guide = "Rule of thirds";
  else if (guide->get_active_row_number()==2)
    pp->crop.guide = "Rule of diagonals";
  else if (guide->get_active_row_number()==3)
    pp->crop.guide = "Harmonic means 1";
  else if (guide->get_active_row_number()==4)
    pp->crop.guide = "Harmonic means 2";
  else if (guide->get_active_row_number()==5)
    pp->crop.guide = "Harmonic means 3";
  else if (guide->get_active_row_number()==6)
    pp->crop.guide = "Harmonic means 4";

    if (pedited) {
        pedited->crop.enabled       = !enabled->get_inconsistent();
        pedited->crop.ratio         = ratio->get_active_row_number() != 8;
        pedited->crop.orientation   = orientation->get_active_row_number() != 2;
        pedited->crop.guide         = guide->get_active_row_number() != 7;
        pedited->crop.fixratio      = !fixr->get_inconsistent();
        pedited->crop.w             = wDirty;
        pedited->crop.h             = hDirty;
        pedited->crop.x             = xDirty;
        pedited->crop.y             = yDirty;
    }

}

void Crop::selectPressed () {

  if (clistener)
    clistener->cropSelectRequested ();
}

void Crop::notifyListener () {

    if (listener && enabled->get_active ()) 
        listener->panelChanged (EvCrop, Glib::ustring::compose ("%1=%2, %3=%4\n%5=%6, %7=%8", M("TP_CROP_X"), nx, M("TP_CROP_Y"), ny, M("TP_CROP_W"), nw, M("TP_CROP_H"), nh));
}

void Crop::enabledChanged () {

    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            econn.block (true);
            enabled->set_active (false);
            econn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvCrop, Glib::ustring::compose ("%1=%2, %3=%4\n%5=%6, %7=%8", M("TP_CROP_X"), nx, M("TP_CROP_Y"), ny, M("TP_CROP_W"), nw, M("TP_CROP_H"), nh));
        else
            listener->panelChanged (EvCrop, M("GENERAL_DISABLED"));
    }
}

int notifylistener (void* data) {

    gdk_threads_enter ();
    ((Crop*)data)->notifyListener ();
    gdk_threads_leave ();
    return 0;
}

int refreshspins (void* data) {

    gdk_threads_enter ();
    RefreshSpinHelper* rsh = (RefreshSpinHelper*) data;
    rsh->crop->refreshSpins (rsh->notify);
    delete rsh;
    gdk_threads_leave ();
    return 0;
}

void Crop::resizeScaleChanged (double rsc) {

    lastScale = rsc;

    nx = (int)round (nsx * rsc);
    ny = (int)round (nsy * rsc);
    nw = (int)round (nsw * rsc);
    nh = (int)round (nsh * rsc);

    if (nx+nw > maxw || ny+nh > maxh) 
        setDimensions (nx+nw, ny+nh);

    g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
}

void Crop::hFlipCrop () {

    nx = maxw - nx - nw;
    nsx = nx / lastScale;
    g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
}

void Crop::vFlipCrop () {

    ny = maxh - ny - nh;
    nsy = ny / lastScale;
    g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
}

void Crop::rotateCrop (int deg) {

    int tmp;
    switch ((360+deg-lastRotationDeg)%360) {
        case 90:
            tmp = nx;
            nx = maxh - ny - nh;
            ny = tmp;
            tmp = nw;
            nw = nh;
            nh = tmp;
            break;
        case 270:
            tmp = ny;
            ny = maxw - nx - nw;
            nx = tmp;
            tmp = nw;
            nw = nh;
            nh = tmp;
            break;
        case 180:
            nx = maxw - nx - nw;
            ny = maxh - ny - nh;
            break;
    }
    nsx = nx / lastScale;
    nsy = ny / lastScale;
    nsw = nw / lastScale;
    nsh = nh / lastScale;

    lastRotationDeg = deg;
    g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
}

void Crop::positionChanged () {

    xDirty = true;
    yDirty = true;

    int X = (int)x->get_value ();
    int Y = (int)y->get_value ();
    int W = nw;
    int H = nh;
    cropMoved (X, Y, W, H);
    g_idle_add (notifylistener, this);
}

void Crop::widthChanged () {

    wDirty = true;

    int X = nx;
    int Y = ny;
    int W = (int)w->get_value ();
    int H = nh;
    cropWidth2Resized (X, Y, W, H);
    g_idle_add (notifylistener, this);
}

void Crop::heightChanged () {

    hDirty = true;

    int X = nx;
    int Y = ny;
    int W = nw;
    int H = (int)h->get_value ();
    cropHeight2Resized (X, Y, W, H);
    g_idle_add (notifylistener, this);
}


void Crop::ratioChanged () {

    if (batchMode && lastAspect != fixr->get_active ()) {
        if (fixr->get_inconsistent()) {
            fixr->set_inconsistent (false);
            fconn.block (true);
            fixr->set_active (false);
            fconn.block (false);
        }
        else if (lastAspect)
            fixr->set_inconsistent (true);

        lastAspect = fixr->get_active ();
    }


  if (fixr->get_active() && !fixr->get_inconsistent()) {
 
//    int W = w->get_value ();
//    int H = h->get_value ();
    int W = nw;
    int H = nh;
    int X = nx;
    int Y = ny;
    if (W>=H)
      cropWidth2Resized (X, Y, W, H);
    else 
      cropHeight2Resized (X, Y, W, H);

    g_idle_add (refreshspins, new RefreshSpinHelper (this, true));
  }
}

void Crop::refreshSize () {

    if (!batchMode) {

        std::ostringstream ostrin;
        ostrin.precision (3);
    //    ostrin << h->get_value()/dpi->get_value() << " in x " << w->get_value()/dpi->get_value() << " in";;
        ostrin << nh/dpi->get_value() << " in x " << nw/dpi->get_value() << " in";;

        sizein->set_text (ostrin.str ());

        std::ostringstream ostrcm;
        ostrcm.precision (3);
    //    ostrcm << h->get_value()/dpi->get_value()*2.54 << " cm x " << w->get_value()/dpi->get_value()*2.54 << " cm";;
        ostrcm << nh/dpi->get_value()*2.54 << " cm x " << nw/dpi->get_value()*2.54 << " cm";;

        sizecm->set_text (ostrcm.str ());
    }
}

void Crop::setDimensions (int mw, int mh) {

  maxw = mw;
  maxh = mh;

  xconn.block (true);
  yconn.block (true);
  wconn.block (true);
  hconn.block (true);

  w->set_range (0, maxw);
  h->set_range (0, maxh);
  x->set_range (0, maxw);
  y->set_range (0, maxh);

  xconn.block (false);
  yconn.block (false);
  wconn.block (false);
  hconn.block (false);

  if (enabled->get_active()==false) {
    nx = 0;
    ny = 0;
    nw = mw;
    nh = mh;

    nsx = nx / lastScale;
    nsy = ny / lastScale;
    nsw = nw / lastScale;
    nsh = nh / lastScale;
    refreshSpins ();
  }  
  refreshSize ();
}

struct setdimparams {
    Crop* crop;
    int x;
    int y;
};

int setdim (void* data) {

    gdk_threads_enter ();
    setdimparams* params = (setdimparams*)data;
    params->crop->setDimensions (params->x, params->y);
    delete params;
    gdk_threads_leave ();
    return 0;
}

void Crop::sizeChanged (int x, int y, int ow, int oh) {

    setdimparams* params = new setdimparams;
    params->x = x;
    params->y = y;
    params->crop = this;
    g_idle_add (setdim, params);
}

bool Crop::refreshSpins (bool notify) {

  xconn.block (true);
  yconn.block (true);
  wconn.block (true);
  hconn.block (true);

  x->set_value (nx);
  y->set_value (ny);
  w->set_value (nw);
  h->set_value (nh);

  xDirty = true;
  yDirty = true;
  wDirty = true;
  hDirty = true;

  xconn.block (false);
  yconn.block (false);
  wconn.block (false);
  hconn.block (false);

  refreshSize ();
  if (notify)
    notifyListener ();

  return false;
}

void Crop::cropMoved (int &X, int &Y, int &W, int &H) {

//  W = w->get_value ();
//  H = h->get_value ();
  W = nw;
  H = nh;

  if (X+W>maxw)
    X = maxw-W;
  if (Y+H>maxh)
    Y = maxh-H;
  if (X<0)
    X = 0;
  if (Y<0)
    Y = 0;

  nx = X;
  ny = Y;
  nw = W;
  nh = H;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropWidth1Resized (int &X, int &Y, int &W, int &H) {

  if (W<0)
    W = 0;
  if (H<0)
    H = 0;

  if (X<0) {
    W += X;
    X = 0;
  }
  if (fixr->get_active()) {
    double r = getRatio();
    int W2max = (int)round(r*(maxh-Y));
    if (W>W2max) {
      X += W - W2max;
      W = W2max;  
    }
    H = (int)round(W / r);
  }

  nx = X;
  ny = Y;
  nw = W;
  nh = H;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropWidth2Resized (int &X, int &Y, int &W, int &H) {

//  X = x->get_value ();
//  Y = y->get_value ();
  X = nx;
  Y = ny;

  if (W<0)
    W = 0;
  if (H<0)
    H = 0;

  if (W>maxw-X)
    W = maxw-X;

  if (fixr->get_active()) {
    double r = getRatio();
    int W2max = (int)round(r*(maxh-Y));
    if (W>W2max)
      W = W2max;  
    H = (int)round(W / r);
  }

  nx = X;
  ny = Y;
  nw = W;
  nh = H;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropHeight1Resized (int &X, int &Y, int &W, int &H) {

  if (W<0)
    W = 0;
  if (H<0)
    H = 0;

  if (Y<0) {
    H += Y;
    Y = 0;
  }

  if (fixr->get_active()) {
    double r = getRatio();
    int H2max = (int)round((maxw-X) / r);
    if (H>H2max) {
      Y += H - H2max;
      H = H2max;  
    }
    W = (int)round(H * r);
  }

  nx = X;
  ny = Y;
  nw = W;
  nh = H;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropHeight2Resized (int &X, int &Y, int &W, int &H) {

//  X = x->get_value ();
//  Y = y->get_value ();
  X = nx;
  Y = ny;

  if (W<0)
    W = 0;
  if (H<0)
    H = 0;
  int H1max = maxh-Y;
  if (H>H1max)
    H = H1max;
  if (fixr->get_active()) {
    double r = getRatio ();
    int H2max = (int)round ((maxw-X) / r);
    if (H>H2max)
      H = H2max;  
    W = (int)round(H * r);
  }

  nx = X;
  ny = Y;
  nw = W;
  nh = H;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropInit (int &x, int &y, int &w, int &h) {

  nx = x;
  ny = y;
  nw = 1;
  nh = 1;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  w = 1; h = 1;
  
  econn.block (true);
  enabled->set_active (1);
  econn.block (false);
  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropResized (int &x, int &y, int& x2, int& y2) {

  if (x2<0)
    x2 = 0;
  if (y2<0)
    y2 = 0;
  if (x2>=maxw)
    x2 = maxw-1;
  if (y2>=maxh)
    y2 = maxh-1;

  int X, Y;
  int W;
  if (x<x2) {
    W = x2 - x + 1;
    X = x;
  }
  else {
    W = x - x2 + 1;
    X = x2;
  }
  int H;
  if (y<y2) {
    H = y2 - y + 1;
    Y = y;
  }
  else {
    H = y - y2 + 1;
    Y = y2;
  }

  if (W>maxw)
    W = maxw;
  if (H>maxh)
    H = maxh;

  if (fixr->get_active()) {
    double r = getRatio ();
    if (y<=y2) {
      int W2max = (int)round ((maxh-Y) * r);
      if (W>W2max)
        W = W2max;  
    }
    else {
      int W2max = (int)round (y * r);
      if (W>W2max)
        W = W2max;  
    }
    H = (int)round(W / r);
    if (x<x2)
      x2 = x + W - 1;
    else 
      x2 = x - W + 1;
    if (y<y2) 
      y2 = y + H - 1;
    else 
      y2 = y - H + 1;
  }

  if (x<x2) {
    W = x2 - x + 1;
    X = x;
  }
  else {
    W = x - x2 + 1;
    X = x2;
  }
  if (y<y2) {
    H = y2 - y + 1;
    Y = y;
  }
  else {
    H = y - y2 + 1;
    Y = y2;
  }

  nx = X;
  ny = Y;
  nw = W;
  nh = H;

  nsx = nx / lastScale;
  nsy = ny / lastScale;
  nsw = nw / lastScale;
  nsh = nh / lastScale;

  g_idle_add (refreshspins, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropManipReady () {

    g_idle_add (notifylistener, this);
}

double Crop::getRatio () {

  double r = -1.0;
  if (fixr->get_active()==false)
    return r;
  if (ratio->get_active_row_number()==0)
    r = 3.0/2.0;
  else if (ratio->get_active_row_number()==1)
    r = 4.0/3.0;
  else if (ratio->get_active_row_number()==2)
    r = 16.0/9.0;
  else if (ratio->get_active_row_number()==3)
    r = 16.0/10.0;
  else if (ratio->get_active_row_number()==4)
    r = 5.0/4.0;
  else if (ratio->get_active_row_number()==5)
    r = 2.0/1.0;
  else if (ratio->get_active_row_number()==6)
    r = 1.0/1.0;
  else if (ratio->get_active_row_number()==7)
    r = 1.414;
  if (orientation->get_active_row_number()==0)
    return r;
  else
    return 1.0 / r;
}

void Crop::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);

	ratio->append_text ("(Unchanged)");
	orientation->append_text ("(Unchanged)");
	guide->append_text ("(Unchanged)");
    removeIfThere (this, dpibox);
}
