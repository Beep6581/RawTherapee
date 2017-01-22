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
#include "crop.h"
#include "options.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

extern Options options;

class RefreshSpinHelper
{

public:
    Crop* crop;
    bool  notify;
    RefreshSpinHelper (Crop* _crop, bool _notify)
        : crop(_crop), notify(_notify) {}
};

Crop::Crop (): FoldableToolPanel(this, "crop", M("TP_CROP_LABEL"), false, true)
{

    clistener = nullptr;

    maxw = 3000;
    maxh = 2000;

    Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());

    hb1->pack_start (*Gtk::manage (new Gtk::Label (Glib::ustring(" ") + M("TP_CROP_X") + ": ")));
    x = Gtk::manage (new MySpinButton ());
    x->set_size_request (60, -1);
    hb1->pack_start (*x);

    hb1->pack_start (*Gtk::manage (new Gtk::Label (Glib::ustring("   ") + M("TP_CROP_Y") + ": ")));
    y = Gtk::manage (new MySpinButton ());
    y->set_size_request (60, -1);
    hb1->pack_start (*y);

    pack_start (*hb1, Gtk::PACK_SHRINK, 2);

    Gtk::HBox* hb2 = Gtk::manage (new Gtk::HBox ());

    hb2->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_W") + ": ")));
    w = Gtk::manage (new MySpinButton ());
    w->set_size_request (60, -1);
    hb2->pack_start (*w);

    hb2->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_H") + ": ")));
    h = Gtk::manage (new MySpinButton ());
    h->set_size_request (60, -1);
    hb2->pack_start (*h);

    pack_start (*hb2, Gtk::PACK_SHRINK, 4);

    selectCrop = Gtk::manage (new Gtk::Button (M("TP_CROP_SELECTCROP")));
    selectCrop->set_image (*Gtk::manage (new RTImage ("crop.png")));

    pack_start (*selectCrop, Gtk::PACK_SHRINK, 2);

    Gtk::HBox* hb3 = Gtk::manage (new Gtk::HBox ());

    fixr = Gtk::manage (new Gtk::CheckButton (M("TP_CROP_FIXRATIO")));
    fixr->set_active (1);

    hb3->pack_start (*fixr, Gtk::PACK_SHRINK, 4);

    ratio = Gtk::manage (new MyComboBoxText ());
    hb3->pack_start (*ratio, Gtk::PACK_EXPAND_WIDGET, 4);

    orientation = Gtk::manage (new MyComboBoxText ());
    hb3->pack_start (*orientation);

    pack_start (*hb3, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hb31 = Gtk::manage (new Gtk::HBox ());

    hb31->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_GUIDETYPE"))), Gtk::PACK_SHRINK, 4);
    guide = Gtk::manage (new MyComboBoxText ());
    hb31->pack_start (*guide);

    pack_start (*hb31, Gtk::PACK_SHRINK, 4);

    // ppibox START
    ppibox = Gtk::manage (new Gtk::VBox());
    ppibox->pack_start (*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_SHRINK, 2);

    Gtk::HBox* hb4 = Gtk::manage (new Gtk::HBox ());
    hb4->pack_start (*Gtk::manage (new Gtk::Label (M("TP_CROP_PPI"))));
    ppi = Gtk::manage (new MySpinButton ());
    ppi->set_size_request (60, -1);
    hb4->pack_start (*ppi);

    sizebox = Gtk::manage (new Gtk::VBox());

    sizecm = Gtk::manage (new Gtk::Label (M("GENERAL_NA") + " cm x " + M("GENERAL_NA") + " cm"));
    sizein = Gtk::manage (new Gtk::Label (M("GENERAL_NA") + " in x " + M("GENERAL_NA") + " in"));

    sizebox->pack_start (*sizecm, Gtk::PACK_SHRINK, 4);
    sizebox->pack_start (*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_SHRINK, 6);
    sizebox->pack_start (*sizein, Gtk::PACK_SHRINK, 4);
    sizebox->pack_start (*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_SHRINK, 6);
    sizebox->pack_start (*hb4, Gtk::PACK_SHRINK, 2);

    ppibox->pack_start (*sizebox, Gtk::PACK_SHRINK, 1);
    pack_start (*ppibox, Gtk::PACK_SHRINK, 0);

    ppi->set_value (300);
    // ppibox END

    /****************
    * Crop Ratio
    *****************/
    int NumberOfCropRatios = 26;    //!!! change this value when adding new crop ratios
    cropratio.resize (NumberOfCropRatios);

    cropratio[0].label  = "3:2";
    cropratio[0].value  = 3.0 / 2.0;
    cropratio[1].label  = "4:3";
    cropratio[1].value  = 4.0 / 3.0;
    cropratio[2].label  = "16:9";
    cropratio[2].value  = 16.0 / 9.0;
    cropratio[3].label  = "16:10";
    cropratio[3].value  = 16.0 / 10.0;
    cropratio[4].label  = "1:1";
    cropratio[4].value  = 1.0 / 1.0;
    cropratio[5].label  = "2:1";
    cropratio[5].value  = 2.0 / 1.0;
    cropratio[6].label  = "3:1";
    cropratio[6].value  = 3.0 / 1.0;
    cropratio[7].label  = "4:1";
    cropratio[7].value  = 4.0 / 1.0;
    cropratio[8].label  = "5:1";
    cropratio[8].value  = 5.0 / 1.0;
    cropratio[9].label  = "6:1";
    cropratio[9].value  = 6.0 / 1.0;
    cropratio[10].label = "7:1";
    cropratio[10].value = 7.0 / 1.0;
    cropratio[11].label = "4:5";
    cropratio[11].value = 4.0 / 5.0;
    cropratio[12].label = "5:7";
    cropratio[12].value = 5.0 / 7.0;
    cropratio[13].label = "6:7";
    cropratio[13].value = 6.0 / 7.0;
    cropratio[14].label = "6:17";
    cropratio[14].value = 6.0 / 17.0;
    cropratio[15].label = "24:65 - XPAN";
    cropratio[15].value = 24.0 / 65.0;
    cropratio[16].label = "1.414 - DIN EN ISO 216";
    cropratio[16].value = 1.414;
    cropratio[17].label = "3.5:5";
    cropratio[17].value = 3.5 / 5.0;
    cropratio[18].label = "8.5:11 - US Letter";
    cropratio[18].value = 8.5 / 11.0;
    cropratio[19].label = "9.5:12";
    cropratio[19].value = 9.5 / 12.0;
    cropratio[20].label = "10:12";
    cropratio[20].value = 10.0 / 12.0;
    cropratio[21].label = "11:14";
    cropratio[21].value = 11.0 / 14.0;
    cropratio[22].label = "11:17 - Tabloid";
    cropratio[22].value = 11.0 / 17.0;
    cropratio[23].label = "13:19";
    cropratio[23].value = 13.0 / 19.0;
    cropratio[24].label = "17:22";
    cropratio[24].value = 17.0 / 22.0;
    cropratio[25].label = "45:35 - ePassport";
    cropratio[25].value = 45.0 / 35.0;



    // populate the combobox
    for (int i = 0; i < NumberOfCropRatios; i++) {
        ratio->append (cropratio[i].label);
    }

    ratio->set_active (0);

    orientation->append (M("GENERAL_LANDSCAPE"));
    orientation->append (M("GENERAL_PORTRAIT"));
    orientation->append (M("GENERAL_ASIMAGE"));
    orientation->set_active (2);

    guide->append (M("TP_CROP_GTNONE"));
    guide->append (M("TP_CROP_GTFRAME"));
    guide->append (M("TP_CROP_GTRULETHIRDS"));
    guide->append (M("TP_CROP_GTDIAGONALS"));
    guide->append (M("TP_CROP_GTHARMMEANS"));
    guide->append (M("TP_CROP_GTGRID"));
    guide->append (M("TP_CROP_GTTRIANGLE1"));
    guide->append (M("TP_CROP_GTTRIANGLE2"));
    guide->append (M("TP_CROP_GTEPASSPORT"));
    guide->set_active (0);

    w->set_range (1, maxw);
    h->set_range (1, maxh);
    x->set_range (0, maxw - 1);
    y->set_range (0, maxh - 1);

    x->set_digits (0);
    x->set_increments (1, 100);
    x->set_value (0);

    y->set_digits (0);
    y->set_increments (1, 100);
    y->set_value (0);

    w->set_digits (0);
    w->set_increments (1, 100);
    w->set_value (200);

    h->set_digits (0);
    h->set_increments (1, 100);
    h->set_value (200);

    ppi->set_digits (0);
    ppi->set_increments (1, 100);
    ppi->set_range (50, 12000);
    ppi->set_value (300);

    xconn = x->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::positionChanged), true);
    yconn = y->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::positionChanged), true);
    wconn = w->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::widthChanged), true);
    hconn = h->signal_value_changed().connect ( sigc::mem_fun(*this, &Crop::heightChanged), true);
    fconn = fixr->signal_toggled().connect( sigc::mem_fun(*this, &Crop::ratioFixedChanged) );
    rconn = ratio->signal_changed().connect( sigc::mem_fun(*this, &Crop::ratioChanged) );
    oconn = orientation->signal_changed().connect( sigc::mem_fun(*this, &Crop::ratioChanged) );
    gconn = guide->signal_changed().connect( sigc::mem_fun(*this, &Crop::notifyListener) );
    selectCrop->signal_pressed().connect( sigc::mem_fun(*this, &Crop::selectPressed) );
    ppi->signal_value_changed().connect( sigc::mem_fun(*this, &Crop::refreshSize) );

    nx = ny = nw = nh = 0;
    lastRotationDeg = 0;
    show_all ();
}

void Crop::writeOptions ()
{

    options.cropPPI = (int)ppi->get_value ();
}

void Crop::readOptions ()
{

    disableListener ();

    ppi->set_value (options.cropPPI);

    enableListener ();
}

void Crop::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    xconn.block (true);
    yconn.block (true);
    wconn.block (true);
    hconn.block (true);
    rconn.block (true);
    fconn.block (true);
    oconn.block (true);
    gconn.block (true);

    setEnabled(pp->crop.enabled);

    // check if the new values are larger than the maximum
    double tmp, maxw, maxh;
    w->get_range (tmp, maxw);
    h->get_range (tmp, maxh);

    if (pp->crop.x + pp->crop.w > (int)maxw || pp->crop.y + pp->crop.h > (int)maxh) {
        setDimensions (pp->crop.x + pp->crop.w, pp->crop.y + pp->crop.h);
    }

    ratio->set_active_text (pp->crop.ratio);
    fixr->set_active (pp->crop.fixratio);

    const bool flip_orientation = pp->crop.fixratio && cropratio[ratio->get_active_row_number()].value < 1.0;

    if (pp->crop.orientation == "Landscape") {
        orientation->set_active (flip_orientation ? 1 : 0);
    } else if (pp->crop.orientation == "Portrait") {
        orientation->set_active (flip_orientation ? 0 : 1);
    } else {
        orientation->set_active (2);
    }

    if (pp->crop.guide == "None") {
        guide->set_active (0);
    } else if (pp->crop.guide == "Frame") {
        guide->set_active (1);
    } else if (pp->crop.guide == "Rule of thirds") {
        guide->set_active (2);
    } else if (pp->crop.guide == "Rule of diagonals") {
        guide->set_active (3);
    } else if (!strncmp(pp->crop.guide.data(), "Harmonic means", 14)) {
        guide->set_active (4);
    } else if (pp->crop.guide == "Grid") {
        guide->set_active (5);
    } else if (pp->crop.guide == "Golden Triangle 1") {
        guide->set_active (6);
    } else if (pp->crop.guide == "Golden Triangle 2") {
        guide->set_active (7);
    } else if (pp->crop.guide == "ePassport") {
        guide->set_active (8);
    }

    x->set_value(pp->crop.x);
    y->set_value(pp->crop.y);
    w->set_value(std::max(pp->crop.w, 1));
    h->set_value(std::max(pp->crop.h, 1));

    nx = pp->crop.x;
    ny = pp->crop.y;
    nw = pp->crop.w;
    nh = pp->crop.h;

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

        if (!pedited->crop.ratio) {
            ratio->set_active_text (M("GENERAL_UNCHANGED"));
        }

        if (!pedited->crop.orientation) {
            orientation->set_active_text (M("GENERAL_UNCHANGED"));
        }

        if (!pedited->crop.guide) {
            guide->set_active_text (M("GENERAL_UNCHANGED"));
        }

        set_inconsistent (multiImage && !pedited->crop.enabled);
        fixr->set_inconsistent (!pedited->crop.fixratio);
    }

    lastFixRatio = pp->crop.fixratio;

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

void Crop::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->crop.enabled = getEnabled ();
    pp->crop.x = nx;
    pp->crop.y = ny;
    pp->crop.w = nw;
    pp->crop.h = nh;
    pp->crop.fixratio = fixr->get_active ();
    pp->crop.ratio = ratio->get_active_text ();

    // for historical reasons we store orientation different if ratio is written as 2:3 instead of 3:2, but in GUI 'landscape' is always long side horizontal regardless of the ratio is written short or long side first.
    const bool flip_orientation = fixr->get_active() && cropratio[ratio->get_active_row_number()].value < 1.0;

    if (orientation->get_active_row_number() == 0) {
        pp->crop.orientation = flip_orientation ? "Portrait" : "Landscape";
    } else if (orientation->get_active_row_number() == 1) {
        pp->crop.orientation = flip_orientation ? "Landscape" : "Portrait";
    } else {
        pp->crop.orientation = "As Image";
    }

    if (guide->get_active_row_number() == 0) {
        pp->crop.guide = "None";
    } else if (guide->get_active_row_number() == 1) {
        pp->crop.guide = "Frame";
    } else if (guide->get_active_row_number() == 2) {
        pp->crop.guide = "Rule of thirds";
    } else if (guide->get_active_row_number() == 3) {
        pp->crop.guide = "Rule of diagonals";
    } else if (guide->get_active_row_number() == 4) {
        pp->crop.guide = "Harmonic means";
    } else if (guide->get_active_row_number() == 5) {
        pp->crop.guide = "Grid";
    } else if (guide->get_active_row_number() == 6) {
        pp->crop.guide = "Golden Triangle 1";
    } else if (guide->get_active_row_number() == 7) {
        pp->crop.guide = "Golden Triangle 2";
    } else if (guide->get_active_row_number() == 8) {
        pp->crop.guide = "ePassport";
    }

    if (pedited) {
        pedited->crop.enabled       = !get_inconsistent();
        pedited->crop.ratio         = ratio->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->crop.orientation   = orientation->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->crop.guide         = guide->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->crop.fixratio      = !fixr->get_inconsistent();
        pedited->crop.w             = wDirty;
        pedited->crop.h             = hDirty;
        pedited->crop.x             = xDirty;
        pedited->crop.y             = yDirty;
    }

}

void Crop::trim (ProcParams* pp, int ow, int oh)
{

    int xmin = pp->crop.x;
    int ymin = pp->crop.y;

    if (xmin > ow || ymin > oh) {
        // the crop is completely out of the image, so we disable the crop
        pp->crop.enabled = false;
        // and we set the values to the defaults
        pp->crop.x = 0;
        pp->crop.y = 0;
        pp->crop.w = ow;
        pp->crop.h = oh;
        // the ratio is now not guaranteed, so we set it off
        pp->crop.fixratio = false;
    } else {
        if ((xmin + pp->crop.w) > ow) {
            // crop overflow in the width dimension ; we trim it
            pp->crop.w = ow - xmin;
        }

        if ((ymin + pp->crop.h) > oh) {
            // crop overflow in the height dimension ; we trim it
            pp->crop.h = oh - ymin;
        }
    }
}

bool Crop::inImageArea (int x, int y)
{
    return x >= 0 && x < maxw && y >= 0 && y < maxh;
}

void Crop::selectPressed ()
{

    if (clistener) {
        clistener->cropSelectRequested ();
    }
}

void Crop::notifyListener ()
{

    if (listener && getEnabled ()) {
        if (nw == 1 && nh == 1) {
            setEnabled(false);
            nx = (int)x->get_value ();
            ny = (int)y->get_value ();
            nw = (int)w->get_value ();
            nh = (int)h->get_value ();
            listener->panelChanged (EvCrop, M("GENERAL_DISABLED"));
        } else {
            listener->panelChanged (EvCrop, Glib::ustring::compose ("%1=%2, %3=%4\n%5=%6, %7=%8", M("TP_CROP_X"), nx, M("TP_CROP_Y"), ny, M("TP_CROP_W"), nw, M("TP_CROP_H"), nh));
        }
    }
}

void Crop::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvCrop, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvCrop, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvCrop, M("GENERAL_DISABLED"));
        }
    }
}

int notifyListenerUI (void* data)
{
    (static_cast<Crop*>(data))->notifyListener ();
    return 0;
}

int refreshSpinsUI (void* data)
{
    RefreshSpinHelper* rsh = static_cast<RefreshSpinHelper*>(data);
    rsh->crop->refreshSpins (rsh->notify);
    delete rsh;
    return 0;
}

void Crop::hFlipCrop ()
{

    nx = maxw - nx - nw;
    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::vFlipCrop ()
{

    ny = maxh - ny - nh;
    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::rotateCrop (int deg, bool hflip, bool vflip)
{

    int rotation = (360 + deg - lastRotationDeg) % 360;

    if((hflip != vflip) && ((rotation % 180) == 90)) {
        rotation = (rotation + 180) % 360;
    }

    int tmp;

    switch (rotation) {
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

    lastRotationDeg = deg;
    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::positionChanged ()
{

    xDirty = true;
    yDirty = true;

    int X = (int)x->get_value ();
    int Y = (int)y->get_value ();
    int W = nw;
    int H = nh;
    cropMoved (X, Y, W, H);
    add_idle (notifyListenerUI, this);
}

void Crop::widthChanged ()
{

    wDirty = true;

    int X = nx;
    int Y = ny;
    int W = (int)w->get_value ();
    int H = nh;
    cropWidth2Resized (X, Y, W, H);
    add_idle (notifyListenerUI, this);
}

void Crop::heightChanged ()
{

    hDirty = true;

    int X = nx;
    int Y = ny;
    int W = nw;
    int H = (int)h->get_value ();
    cropHeight2Resized (X, Y, W, H);
    add_idle (notifyListenerUI, this);
}

// Fixed ratio toggle button
void Crop::ratioFixedChanged ()
{
    // Batch mode handling when enabling/disabling fixed crop
    if (batchMode && lastFixRatio != fixr->get_active ()) {
        if (fixr->get_inconsistent()) {
            fixr->set_inconsistent (false);
            fconn.block (true);
            fixr->set_active (false);
            fconn.block (false);
        } else if (lastFixRatio) {
            fixr->set_inconsistent (true);
        }
    }

    lastFixRatio = fixr->get_active ();
    adjustCropToRatio();
}

// change to orientation or ration
void Crop::ratioChanged ()
{
    if (!fixr->get_active ()) {
        fixr->set_active(true);    // will ajust ratio anyway
    } else {
        adjustCropToRatio();
    }
}

// Correct current crop if it doesn't fit
void Crop::adjustCropToRatio()
{
    if (fixr->get_active() && !fixr->get_inconsistent()) {

//        int W = w->get_value ();
//        int H = h->get_value ();
        int W = nw;
        int H = nh;
        int X = nx;
        int Y = ny;

        if (W >= H) {
            cropWidth2Resized (X, Y, W, H);
        } else {
            cropHeight2Resized (X, Y, W, H);
        }
    }

    // This will save the options
    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, true));
}

void Crop::refreshSize ()
{

    if (!batchMode) {

        std::ostringstream ostrin;
        ostrin.precision (3);
        //    ostrin << h->get_value()/ppi->get_value() << " in x " << w->get_value()/ppi->get_value() << " in";;
        ostrin << nh / ppi->get_value() << " in x " << nw / ppi->get_value() << " in";;

        sizein->set_text (ostrin.str ());

        std::ostringstream ostrcm;
        ostrcm.precision (3);
        //    ostrcm << h->get_value()/ppi->get_value()*2.54 << " cm x " << w->get_value()/ppi->get_value()*2.54 << " cm";;
        ostrcm << nh / ppi->get_value() * 2.54 << " cm x " << nw / ppi->get_value() * 2.54 << " cm";;

        sizecm->set_text (ostrcm.str ());
    }
}

/*
 * Set the maximum dimensions of the image. This method can be called with wrong values, then
 * called with the good ones !?
 */
void Crop::setDimensions (int mw, int mh)
{

    maxw = mw;
    maxh = mh;

    bool xconnWasBlocked = xconn.block (true);
    bool yconnWasBlocked = yconn.block (true);
    bool wconnWasBlocked = wconn.block (true);
    bool hconnWasBlocked = hconn.block (true);

    w->set_range (1, maxw);
    h->set_range (1, maxh);
    x->set_range (0, maxw - 1);
    y->set_range (0, maxh - 1);

    if (!xconnWasBlocked) {
        xconn.block (false);
    }

    if (!yconnWasBlocked) {
        yconn.block (false);
    }

    if (!wconnWasBlocked) {
        wconn.block (false);
    }

    if (!hconnWasBlocked) {
        hconn.block (false);
    }

    if (!getEnabled()) {
        nx = 0;
        ny = 0;
        nw = mw;
        nh = mh;

        refreshSpins ();
    }

    refreshSize ();
}

struct setdimparams {
    Crop* crop;
    int x;
    int y;
};

int sizeChangedUI (void* data)
{
    setdimparams* params = static_cast<setdimparams*>(data);
    params->crop->setDimensions (params->x, params->y);
    delete params;
    return 0;
}

void Crop::sizeChanged (int x, int y, int ow, int oh)
{

    setdimparams* params = new setdimparams;
    params->x = x;
    params->y = y;
    params->crop = this;
    add_idle (sizeChangedUI, params);
}

bool Crop::refreshSpins (bool notify)
{

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

    if (notify) {
        notifyListener ();
    }

    return false;
}

void Crop::cropMoved (int &X, int &Y, int &W, int &H)
{

//  W = w->get_value ();
//  H = h->get_value ();
    W = nw;
    H = nh;

    if (X + W > maxw) {
        X = maxw - W;
    }

    if (Y + H > maxh) {
        Y = maxh - H;
    }

    if (X < 0) {
        X = 0;
    }

    if (Y < 0) {
        Y = 0;
    }

    nx = X;
    ny = Y;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropWidth1Resized (int &X, int &Y, int &W, int &H)
{

    int oldXR = nx + nw;

    if (W < 0) {
        W = 0;
    }

    if (W > oldXR) {
        W = oldXR;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        H = (int)round(W / r);
        int Hmax = min(ny + nh, maxh - ny);

        if (H > Hmax) {
            H = Hmax;
            W = H * r;
        }

        ny = ny - (H - nh) / 2.0;

        if (ny < 0) {
            ny = 0;
        }

        if (ny + H > maxh) {
            ny = maxh - H;
        }
    }

    X = oldXR - W;
    Y = ny;
    nx = X;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropWidth2Resized (int &X, int &Y, int &W, int &H)
{

    if (W < 0) {
        W = 0;
    }

    if (W > maxw - nx) {
        W = maxw - nx;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        H = (int)round(W / r);
        int Hmax = min(ny + nh, maxh - ny);

        if (H > Hmax) {
            H = Hmax;
            W = H * r;
        }

        ny = ny - (H - nh) / 2.0;

        if (ny < 0) {
            ny = 0;
        }

        if (ny + H > maxh) {
            ny = maxh - H;
        }
    }

    X = nx;
    Y = ny;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropHeight1Resized (int &X, int &Y, int &W, int &H)
{

    int oldYB = ny + nh;

    if (H < 0) {
        H = 0;
    }

    if (H > oldYB) {
        H = oldYB;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        W = (int)round(H * r);
        int Wmax = min(nx + nw, maxw - nx);

        if (W > Wmax) {
            W = Wmax;
            H = W / r;
        }

        nx = nx - (W - nw) / 2.0;

        if (nx < 0) {
            nx = 0;
        }

        if (nx + W > maxw) {
            nx = maxw - W;
        }
    }

    X = nx;
    Y = oldYB - H;
    ny = Y;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropHeight2Resized (int &X, int &Y, int &W, int &H)
{

    if (H < 0) {
        H = 0;
    }

    if (H > maxh - ny) {
        H = maxh - ny;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        W = (int)round(H * r);
        int Wmax = min(nx + nw, maxw - nx);

        if (W > Wmax) {
            W = Wmax;
            H = W / r;
        }

        nx = nx - (W - nw) / 2.0; // nx must be floating point to avoid drifting

        if (nx < 0) {
            nx = 0;
        }

        if (nx + W > maxw) {
            nx = maxw - W;
        }
    }

    X = nx;
    Y = ny;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropTopLeftResized (int &X, int &Y, int &W, int &H)
{

    int oldXR = nx + nw; // right side
    int oldYB = ny + nh; // bottom side

    if (W < 0) {
        W = 0;
    }

    if (H < 0) {
        H = 0;
    }

    if (W > oldXR) {
        W = oldXR;
    }

    if (H > oldYB) {
        H = oldYB;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        W = (int)round(H * r);

        if (W > oldXR) {
            W = oldXR;
            H = (int)round(W / r);
        }
    }

    X = oldXR - W;
    Y = oldYB - H;
    nx = X;
    ny = Y;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropTopRightResized (int &X, int &Y, int &W, int &H)
{

    int oldYB = ny + nh;

    if (W < 0) {
        W = 0;
    }

    if (H < 0) {
        H = 0;
    }

    if (W > maxw - nx) {
        W = maxw - nx;
    }

    if (H > oldYB) {
        H = oldYB;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        W = (int)round(H * r);

        if (W > maxw - nx) {
            W = maxw - nx;
            H = (int)round(W / r);
        }
    }

    X = nx;
    Y = oldYB - H;
    ny = Y;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropBottomLeftResized (int &X, int &Y, int &W, int &H)
{

    int oldXR = nx + nw;

    if (W < 0) {
        W = 0;
    }

    if (H < 0) {
        H = 0;
    }

    if (W > oldXR) {
        W = oldXR;
    }

    if (H > maxh - ny) {
        H = maxh - ny;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        W = (int)round(H * r);

        if (W > oldXR) {
            W = oldXR;
            H = (int)round(W / r);
        }
    }

    X = oldXR - W;
    Y = ny;
    nx = X;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropBottomRightResized (int &X, int &Y, int &W, int &H)
{

    if (W < 0) {
        W = 0;
    }

    if (H < 0) {
        H = 0;
    }

    if (W > maxw - nx) {
        W = maxw - nx;
    }

    if (H > maxh - ny) {
        H = maxh - ny;
    }

    if (fixr->get_active()) {
        double r = getRatio();
        W = (int)round(H * r);

        if (W > maxw - nx) {
            W = maxw - nx;
            H = (int)round(W / r);
        }
    }

    X = nx;
    Y = ny;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropInit (int &x, int &y, int &w, int &h)
{

    nx = x;
    ny = y;
    nw = 1;
    nh = 1;

    w = 1;
    h = 1;

    setEnabled(true);
}

void Crop::cropResized (int &x, int &y, int& x2, int& y2)
{

    if (x2 < 0) {
        x2 = 0;
    }

    if (y2 < 0) {
        y2 = 0;
    }

    if (x2 >= maxw) {
        x2 = maxw - 1;
    }

    if (y2 >= maxh) {
        y2 = maxh - 1;
    }

    int X, Y;
    int W;

    if (x < x2) {
        W = x2 - x + 1;
        X = x;
    } else {
        W = x - x2 + 1;
        X = x2;
    }

    int H;

    if (y < y2) {
        H = y2 - y + 1;
        Y = y;
    } else {
        H = y - y2 + 1;
        Y = y2;
    }

    if (W > maxw) {
        W = maxw;
    }

    if (H > maxh) {
        H = maxh;
    }

    if (fixr->get_active()) {
        double r = getRatio ();

        if (y <= y2) {
            int W2max = (int)round ((maxh - Y) * r);

            if (W > W2max) {
                W = W2max;
            }
        } else {
            int W2max = (int)round (y * r);

            if (W > W2max) {
                W = W2max;
            }
        }

        H = (int)round(W / r);

        if (x < x2) {
            x2 = x + W - 1;
        } else {
            x2 = x - W + 1;
        }

        if (y < y2) {
            y2 = y + H - 1;
        } else {
            y2 = y - H + 1;
        }
    }

    if (x < x2) {
        W = x2 - x + 1;
        X = x;
    } else {
        W = x - x2 + 1;
        X = x2;
    }

    if (y < y2) {
        H = y2 - y + 1;
        Y = y;
    } else {
        H = y - y2 + 1;
        Y = y2;
    }

    nx = X;
    ny = Y;
    nw = W;
    nh = H;

    add_idle (refreshSpinsUI, new RefreshSpinHelper (this, false));
}

void Crop::cropManipReady ()
{

    add_idle (notifyListenerUI, this);
}

double Crop::getRatio ()
{

    double r = -1.0;

    if (!fixr->get_active()) {
        return r;
    }

    r = cropratio[ratio->get_active_row_number()].value;

    if (r < 1.0) {
        r = 1.0 / r;    // convert to long side first (eg 4:5 becomes 5:4)
    }

    if (orientation->get_active_row_number() == 0) {
        return r;
    } else if(orientation->get_active_row_number() == 1) {
        return 1.0 / r;
    } else {
        return maxh <= maxw ? r : 1.0 / r;
    }

}

void Crop::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    ratio->append (M("GENERAL_UNCHANGED"));
    orientation->append (M("GENERAL_UNCHANGED"));
    guide->append (M("GENERAL_UNCHANGED"));
    removeIfThere (this, ppibox);
}
