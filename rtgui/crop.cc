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
#include <vector>

#include "crop.h"

#include "aspectratios.h"
#include "options.h"
#include "rtimage.h"

#include "rtengine/procparams.h"
#include "rtengine/utils.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Crop::TOOL_NAME = "crop";

namespace
{

inline void get_custom_ratio(int w, int h, double &rw, double &rh)
{
    if (w < h) {
        double r = double(h) / double(w);
        int rr = r * 100 + 0.5;
        rw = 1.0;
        rh = rr / 100.0;
    } else {
        double r = double(w) / double(h);
        int rr = r * 100 + 0.5;
        rw = rr / 100.0;
        rh = 1.0;
    }
}

} // namespace

class Crop::CropRatios final
{
public:
    CropRatios() :
        ratios{
            {M("GENERAL_ASIMAGE"), 0.0},
            {M("GENERAL_CURRENT"), -1.0}
        }
    {
        fillAspectRatios(ratios);
    }

    std::vector<Glib::ustring> getLabels() const
    {
        std::vector<Glib::ustring> res;

        res.reserve(ratios.size());

        for (const auto& ratio : ratios) {
            res.push_back(ratio.label);
        }

        return res;
    }

    double getValue(std::size_t index) const
    {
        return
            index < ratios.size()
                ? ratios[index].value
                : ratios[0].value;
    }

    void updateCurrentRatio(double value)
    {
        ratios[1].value = value;
    }

private:
    std::vector<AspectRatio> ratios;
};

Crop::Crop():
    FoldableToolPanel(this, TOOL_NAME, M("TP_CROP_LABEL"), false, true),
    crop_ratios(new CropRatios),
    opt(0),
    wDirty(true),
    hDirty(true),
    xDirty(true),
    yDirty(true),
    lastFixRatio(true)
{

    clistener = nullptr;

    maxw = 3000;
    maxh = 2000;

    methodgrid = Gtk::manage(new Gtk::Grid());
    methodgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(methodgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Label* xlab = Gtk::manage (new Gtk::Label (M("TP_CROP_X") + ":"));
    setExpandAlignProperties(xlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    x = Gtk::manage (new MySpinButton ());
    setExpandAlignProperties(x, true, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    x->set_width_chars(6);

    Gtk::Label* ylab = Gtk::manage (new Gtk::Label (M("TP_CROP_Y") + ":"));
    setExpandAlignProperties(ylab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    y = Gtk::manage (new MySpinButton ());
    setExpandAlignProperties(y, true, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    y->set_width_chars(6);

    Gtk::Label* wlab = Gtk::manage (new Gtk::Label (M("TP_CROP_W") + ":"));
    setExpandAlignProperties(wlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    w = Gtk::manage (new MySpinButton ());
    setExpandAlignProperties(w, true, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    w->set_width_chars(6);

    Gtk::Label* hlab = Gtk::manage (new Gtk::Label (M("TP_CROP_H") + ":"));
    setExpandAlignProperties(hlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    h = Gtk::manage (new MySpinButton ());
    setExpandAlignProperties(h, true, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    h->set_width_chars(6);

    selectCrop = Gtk::manage (new Gtk::Button (M("TP_CROP_SELECTCROP")));
    setExpandAlignProperties(selectCrop, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    selectCrop->get_style_context()->add_class("independent");
    selectCrop->set_image (*Gtk::manage (new RTImage ("crop-small", Gtk::ICON_SIZE_BUTTON)));

    resetCrop = Gtk::manage (new Gtk::Button (M("TP_CROP_RESETCROP")));
    setExpandAlignProperties(resetCrop, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    resetCrop->get_style_context()->add_class("independent");
    resetCrop->set_image (*Gtk::manage (new RTImage ("undo-small", Gtk::ICON_SIZE_BUTTON)));

    methodgrid->attach (*xlab, 0, 0, 1, 1);
    methodgrid->attach (*x, 1, 0, 1, 1);
    methodgrid->attach (*ylab, 2, 0, 1, 1);
    methodgrid->attach (*y, 3, 0, 1, 1);
    methodgrid->attach (*wlab, 0, 1, 1, 1);
    methodgrid->attach (*w, 1, 1, 1, 1);
    methodgrid->attach (*hlab, 2, 1, 1, 1);
    methodgrid->attach (*h, 3, 1, 1, 1);
    methodgrid->attach (*selectCrop, 0, 2, 2, 1);
    methodgrid->attach (*resetCrop, 2, 2, 2, 1);
    pack_start (*methodgrid, Gtk::PACK_EXPAND_WIDGET, 0 );

    Gtk::Separator* methodseparator = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    methodseparator->get_style_context()->add_class("grid-row-separator");
    pack_start (*methodseparator, Gtk::PACK_SHRINK, 0);

    Gtk::Grid* settingsgrid = Gtk::manage(new Gtk::Grid());
    settingsgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(settingsgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    fixr = Gtk::manage (new Gtk::CheckButton (M("TP_CROP_FIXRATIO")));
    setExpandAlignProperties(fixr, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    fixr->set_active (1);

    Gtk::Grid* ratiogrid = Gtk::manage(new Gtk::Grid());
    ratiogrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(ratiogrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    ratio = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(ratio, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    orientation = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(orientation, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    customRatioLabel = Gtk::manage(new Gtk::Label(""));
    customRatioLabel->hide();
    setExpandAlignProperties(customRatioLabel, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);

    ratiogrid->set_column_homogeneous (true);
    ratiogrid->attach (*ratio, 0, 0, 1, 1);
    ratiogrid->attach (*customRatioLabel, 1, 0, 1, 1);
    ratiogrid->attach (*orientation, 1, 0, 1, 1);

    Gtk::Label* guidelab = Gtk::manage (new Gtk::Label (M("TP_CROP_GUIDETYPE")));
    setExpandAlignProperties(guidelab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    guide = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(guide, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    settingsgrid->attach (*fixr, 0, 0, 1, 1);
    settingsgrid->attach (*ratiogrid, 1, 0, 1, 1);
    settingsgrid->attach (*guidelab, 0, 1, 1, 1);
    settingsgrid->attach (*guide, 1, 1, 1, 1);
    pack_start (*settingsgrid, Gtk::PACK_SHRINK, 0 );


    // ppigrid START
    ppigrid = Gtk::manage(new Gtk::Grid());
    ppigrid->get_style_context()->add_class("grid-spacing");
    ppigrid->set_column_homogeneous (true);
    setExpandAlignProperties(ppigrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Separator* ppiseparator = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    ppiseparator->get_style_context()->add_class("grid-row-separator");

    Gtk::Grid* ppisubgrid = Gtk::manage(new Gtk::Grid());
    ppisubgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(ppisubgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Label* ppilab = Gtk::manage (new Gtk::Label (M("TP_CROP_PPI") + ":"));
    setExpandAlignProperties(ppilab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    ppi = Gtk::manage (new MySpinButton ());
    setExpandAlignProperties(ppi, true, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    ppi->set_width_chars(6);

    ppisubgrid->attach (*ppilab, 0, 0, 1, 1);
    ppisubgrid->attach (*ppi, 1, 0, 1, 1);

    sizecm = Gtk::manage (new Gtk::Label (M("GENERAL_NA") + " cm x " + M("GENERAL_NA") + " cm"));
    setExpandAlignProperties(sizecm, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);

    sizein = Gtk::manage (new Gtk::Label (M("GENERAL_NA") + " in x " + M("GENERAL_NA") + " in"));
    setExpandAlignProperties(sizein, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);

    ppigrid->attach (*ppiseparator, 0, 0, 2, 1);
    ppigrid->attach (*sizecm, 1, 1, 1, 1);
    ppigrid->attach (*sizein, 1, 2, 1, 1);
    ppigrid->attach (*ppisubgrid, 0, 1, 1, 2);
    pack_start (*ppigrid, Gtk::PACK_SHRINK, 0 );

    ppi->set_value (300);
    // ppigrid END

    // Populate the combobox
    for (const auto& label : crop_ratios->getLabels()) {
        ratio->append (label);
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
    guide->append (M("TP_CROP_GTCENTEREDSQUARE"));
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
    resetCrop->signal_pressed().connect( sigc::mem_fun(*this, &Crop::doresetCrop) );
    ppi->signal_value_changed().connect( sigc::mem_fun(*this, &Crop::refreshSize) );

    nx = ny = nw = nh = 0;
    lastRotationDeg = 0;

//GTK318
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
    methodgrid->set_row_spacing(4);
    methodgrid->set_column_spacing(4);
    settingsgrid->set_row_spacing(4);
    settingsgrid->set_column_spacing(4);
    ppigrid->set_row_spacing(4);
    ppigrid->set_column_spacing(4);
    ppisubgrid->set_row_spacing(4);
    ppisubgrid->set_column_spacing(4);
#endif
//GTK318

    show_all ();
}

Crop::~Crop()
{
    idle_register.destroy();
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

    const bool flip_orientation =
        pp->crop.fixratio
        && crop_ratios->getValue(ratio->get_active_row_number()) > 0
        && crop_ratios->getValue(ratio->get_active_row_number()) < 1.0;

    if (pp->crop.orientation == "Landscape") {
        orientation->set_active (flip_orientation ? 1 : 0);
    } else if (pp->crop.orientation == "Portrait") {
        orientation->set_active (flip_orientation ? 0 : 1);
    } else {
        orientation->set_active (2);
    }

    switch (pp->crop.guide) {
        case procparams::CropParams::Guide::NONE:
        case procparams::CropParams::Guide::FRAME:
        case procparams::CropParams::Guide::RULE_OF_THIRDS:
        case procparams::CropParams::Guide::RULE_OF_DIAGONALS:
        case procparams::CropParams::Guide::HARMONIC_MEANS:
        case procparams::CropParams::Guide::GRID:
        case procparams::CropParams::Guide::GOLDEN_TRIANGLE_1:
        case procparams::CropParams::Guide::GOLDEN_TRIANGLE_2:
        case procparams::CropParams::Guide::EPASSPORT:
        case procparams::CropParams::Guide::CENTERED_SQUARE: {
            guide->set_active(toUnderlying(pp->crop.guide));
        }
    }

    x->set_value(pp->crop.x);
    y->set_value(pp->crop.y);
    w->set_value(std::max(pp->crop.w, 1));
    h->set_value(std::max(pp->crop.h, 1));

    nx = pp->crop.x;
    ny = pp->crop.y;
    nw = pp->crop.w;
    nh = pp->crop.h;

    customRatioLabel->hide();
    orientation->show();
    if (pp->crop.ratio == "As Image") {
        ratio->set_active(0);
    } else if (pp->crop.ratio == "Current") {
        ratio->set_active(1);
        updateCurrentRatio();
        customRatioLabel->show();
        orientation->hide();
    } else {
        ratio->set_active_text (pp->crop.ratio);
    }
    fixr->set_active (pp->crop.fixratio);

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
    if (ratio->get_active_row_number() == 0) {
        pp->crop.ratio = "As Image";
    } else if (ratio->get_active_row_number() == 1) {
        pp->crop.ratio = "Current";
    } else {
        pp->crop.ratio = ratio->get_active_text ();
    }

    // for historical reasons we store orientation different if ratio is written as 2:3 instead of 3:2, but in GUI 'landscape' is always long side horizontal regardless of the ratio is written short or long side first.
    const bool flip_orientation =
        fixr->get_active()
        && crop_ratios->getValue(ratio->get_active_row_number()) > 0
        && crop_ratios->getValue(ratio->get_active_row_number()) < 1.0;

    if (orientation->get_active_row_number() == 0) {
        pp->crop.orientation = flip_orientation ? "Portrait" : "Landscape";
    } else if (orientation->get_active_row_number() == 1) {
        pp->crop.orientation = flip_orientation ? "Landscape" : "Portrait";
    } else {
        pp->crop.orientation = "As Image";
    }

    pp->crop.guide = procparams::CropParams::Guide(guide->get_active_row_number());

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

void Crop::doresetCrop ()
{
    xDirty = true;
    yDirty = true;
    wDirty = true;
    hDirty = true;

    // Reset ratio, ratio lock and orientation as well
    ratio->set_active(0);
    orientation->set_active(2);
    fixr->set_active(true);

    int X = 0;
    int Y = 0;
    int W = maxw;
    int H = maxh;
    cropResized (X, Y, W, H);
    idle_register.add(
        [this]() -> bool
        {
            notifyListener();
            return false;
        }
    );

    refreshSpins();
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

void Crop::hFlipCrop ()
{
    nx = maxw - nx - nw;
    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::vFlipCrop ()
{
    ny = maxh - ny - nh;
    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
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
    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
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
    idle_register.add(
        [this]() -> bool
        {
            notifyListener();
            return false;
        }
    );
}

void Crop::widthChanged ()
{

    wDirty = true;

    int X = nx;
    int Y = ny;
    int W = (int)w->get_value ();
    int H = nh;
    cropWidth2Resized (X, Y, W, H);
    idle_register.add(
        [this]() -> bool
        {
            notifyListener();
            return false;
        }
    );
}

void Crop::heightChanged ()
{

    hDirty = true;

    int X = nx;
    int Y = ny;
    int W = nw;
    int H = (int)h->get_value ();
    cropHeight2Resized (X, Y, W, H);
    idle_register.add(
        [this]() -> bool
        {
            notifyListener();
            return false;
        }
    );
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
    if (ratio->get_active_row_number() == 1) {
        orientation->hide();
        updateCurrentRatio();
        customRatioLabel->show();
    } else {
        orientation->show();
        customRatioLabel->hide();
    }

    if (!fixr->get_active ()) {
        fixr->set_active(true);    // will adjust ratio anyway
    } else {
        adjustCropToRatio();
    }
}

// Correct current crop if it doesn't fit
void Crop::adjustCropToRatio()
{
    if (fixr->get_active() && !fixr->get_inconsistent()) {
        int W1 = nw, W2 = nw;
        int H1 = nh, H2 = nh;
        int X1 = nx, X2 = nx;
        int Y1 = ny, Y2 = ny;

        float r = getRatio();

        H1 = round(W1 / r);
        Y1 = ny + (nh - H1)/2.0;
        if (Y1 < 0) {
            Y1 = 0;
        }
        if (H1 > maxh) {
            H1 = maxh;
            W1 = round(H1 * r);
            X1 = nx + (nw - W1)/2.0;
        }
        if (Y1+H1 > maxh) {
            Y1 = maxh - H1;
        }

        W2 = round(H2 * r);
        X2 = nx + (nw - W2)/2.0;
        if (X2 < 0) {
            X2 = 0;
        }
        if (W2 > maxw) {
            W2 = maxw;
            H2 = round(W2 / r);
            Y2 = ny + (nh - H2)/2.0;
        }
        if (X2+W2 > maxw) {
            X2 = maxw - W2;
        }

        if (W1 * H1 >= W2 * H2) {
            nx = X1;
            ny = Y1;
            nw = W1;
            nh = H1;
        } else {
            nx = X2;
            ny = Y2;
            nw = W2;
            nh = H2;
        }
    }

    // This will save the options
    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(true);
            return false;
        }
    );
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

void Crop::sizeChanged(int x, int y, int ow, int oh)
{
    idle_register.add(
        [this, x, y]() -> bool
        {
            setDimensions(x, y);
            return false;
        }
    );
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

    if (ratio->get_active_row_number() == 1 && !fixr->get_active()) {
        updateCurrentRatio();
    }

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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
//  Glib::signal_idle().connect (sigc::mem_fun(*this, &Crop::refreshSpins));
}

void Crop::cropWidth1Resized (int &X, int &Y, int &W, int &H, float custom_ratio)
{

    int oldXR = nx + nw;

    if (W < 0) {
        W = 0;
    }

    if (W > oldXR) {
        W = oldXR;
    }

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0.f ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropWidth2Resized (int &X, int &Y, int &W, int &H, float custom_ratio)
{

    if (W < 0) {
        W = 0;
    }

    if (W > maxw - nx) {
        W = maxw - nx;
    }

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropHeight1Resized (int &X, int &Y, int &W, int &H, float custom_ratio)
{

    int oldYB = ny + nh;

    if (H < 0) {
        H = 0;
    }

    if (H > oldYB) {
        H = oldYB;
    }

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropHeight2Resized (int &X, int &Y, int &W, int &H, float custom_ratio)
{

    if (H < 0) {
        H = 0;
    }

    if (H > maxh - ny) {
        H = maxh - ny;
    }

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropTopLeftResized (int &X, int &Y, int &W, int &H, float custom_ratio)
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

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropTopRightResized (int &X, int &Y, int &W, int &H, float custom_ratio)
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

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropBottomLeftResized (int &X, int &Y, int &W, int &H, float custom_ratio)
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

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropBottomRightResized (int &X, int &Y, int &W, int &H, float custom_ratio)
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

    if (fixr->get_active() || custom_ratio > 0) {
        double r = custom_ratio > 0 ? custom_ratio : static_cast<float>(getRatio());
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
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

    int W;

    if (x < x2) {
        W = x2 - x + 1;
    } else {
        W = x - x2 + 1;
    }


    int Y;
    if (y < y2) {
        Y = y;
    } else {
        Y = y2;
    }

    if (W > maxw) {
        W = maxw;
    }

    int H;
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

    int X;
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

    idle_register.add(
        [this]() -> bool
        {
            refreshSpins(false);
            return false;
        }
    );
}

void Crop::cropManipReady ()
{
    idle_register.add(
        [this]() -> bool
        {
            notifyListener();
            return false;
        }
    );
}

double Crop::getRatio () const
{
    double r = -1.0;

    if (!fixr->get_active()) {
        return r;
    }

    r = crop_ratios->getValue(ratio->get_active_row_number());
    if (!r) {
        r = maxh <= maxw ? float(maxh)/float(maxw) : float(maxw)/float(maxh);
    }

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
    removeIfThere (this, ppigrid);
    removeIfThere (methodgrid, selectCrop);
    removeIfThere (methodgrid, resetCrop);
}


void Crop::updateCurrentRatio()
{
    double rw, rh;
    get_custom_ratio(w->get_value(), h->get_value(), rw, rh);
    customRatioLabel->set_text(Glib::ustring::compose("%1:%2", rw, rh));
    crop_ratios->updateCurrentRatio(static_cast<double>(w->get_value()) / static_cast<double>(h->get_value()));
}
