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
#include "resize.h"
#include <iomanip>
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;

Resize::Resize () : FoldableToolPanel(this), maxw(100000), maxh(100000), lastEnabled(false) {

	cropw = 0;
	croph = 0;

    enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
    pack_start(*enabled);
    pack_start(*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_SHRINK, 2);

    Gtk::Table* combos = Gtk::manage (new Gtk::Table (2, 2));
    Gtk::Label *label = NULL;

    appliesTo = Gtk::manage (new MyComboBoxText ());
    appliesTo->append_text (M("TP_RESIZE_CROPPEDAREA"));
    appliesTo->append_text (M("TP_RESIZE_FULLIMAGE"));
    appliesTo->set_active (0);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_APPLIESTO")));
    label->set_alignment(0., 0.);
    combos->attach (*label, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*appliesTo, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

	method = Gtk::manage (new MyComboBoxText ());
	method->append_text (M("TP_RESIZE_NEAREST"));
	method->append_text (M("TP_RESIZE_BILINEAR"));
	method->append_text (M("TP_RESIZE_BICUBIC"));
	method->append_text (M("TP_RESIZE_BICUBICSF"));
	method->append_text (M("TP_RESIZE_BICUBICSH"));
	method->append_text (M("TP_RESIZE_LANCZOS"));
	method->set_active (0);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_METHOD")));
    label->set_alignment(0., 0.);
    combos->attach (*label, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*method, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

	spec = Gtk::manage (new MyComboBoxText ());
	spec->append_text (M("TP_RESIZE_SCALE"));
	spec->append_text (M("TP_RESIZE_WIDTH"));
	spec->append_text (M("TP_RESIZE_HEIGHT"));
	spec->append_text (M("TP_RESIZE_FITBOX"));
	spec->set_active (0);
    
    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_SPECIFY")));
    label->set_alignment(0., 0.);
	combos->attach (*label, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*spec, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

	pack_start (*combos, Gtk::PACK_SHRINK, 4);

    scale = new Adjuster (M("TP_RESIZE_SCALE"), 0.01, 4, 0.01, 1.);
    scale->setAdjusterListener (this);

    pack_start (*scale, Gtk::PACK_SHRINK, 4);

    sizeBox = Gtk::manage (new Gtk::VBox ());

    Gtk::HBox* sbox = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* wbox = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
    w = Gtk::manage (new MySpinButton ());
    h = Gtk::manage (new MySpinButton ());
    wbox->set_spacing(3);
    wbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_W"))), Gtk::PACK_SHRINK, 0);
    wbox->pack_start (*w);
    hbox->set_spacing(3);
    hbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_H"))), Gtk::PACK_SHRINK, 0);
    hbox->pack_start (*h);
    sbox->set_spacing(4);
    sbox->pack_start (*wbox);
    sbox->pack_start (*hbox);
    
    sizeBox->pack_start (*sbox, Gtk::PACK_SHRINK, 0);
    sizeBox->show_all ();
    sizeBox->reference ();
    
    w->set_digits (0);
    w->set_increments (1,100);
    w->set_value (800);
    w->set_range (32, 4*maxw);

    h->set_digits (0);
    h->set_increments (1,100);
    h->set_value (600);
    h->set_range (32, 4*maxh);

    wconn = w->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryWChanged), true);
    hconn = h->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryHChanged), true);
	aconn = appliesTo->signal_changed().connect ( sigc::mem_fun(*this, &Resize::appliesToChanged) );
	method->signal_changed().connect ( sigc::mem_fun(*this, &Resize::methodChanged) );
	sconn = spec->signal_changed().connect ( sigc::mem_fun(*this, &Resize::specChanged) );
	enaConn = enabled->signal_toggled().connect ( sigc::mem_fun(*this, &Resize::enabledToggled) );

    show_all();
}

Resize::~Resize () {

    delete scale;
    delete sizeBox;
}

void Resize::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    aconn.block (true);
    wconn.block (true);
    hconn.block (true);
    sconn.block (true);
    scale->block(true);

    scale->setValue (pp->resize.scale);
    w->set_value (pp->resize.width);
    h->set_value (pp->resize.height);
    enabled->set_active (pp->resize.enabled);
    spec->set_active (pp->resize.dataspec);
    updateGUI();

    appliesTo->set_active (0);
    if (pp->resize.appliesTo == "Cropped area")
    	appliesTo->set_active (0);
    else if (pp->resize.appliesTo == "Full image")
    	appliesTo->set_active (1);

    method->set_active (2);
    if (pp->resize.method == "Nearest")
        method->set_active (0);
    else if (pp->resize.method == "Bilinear")
        method->set_active (1);
    else if (pp->resize.method == "Bicubic")
        method->set_active (2);
    else if (pp->resize.method == "Bicubic (Softer)")
        method->set_active (3);
    else if (pp->resize.method == "Bicubic (Sharper)")
        method->set_active (4);
    else if (pp->resize.method == "Lanczos")
        method->set_active (5);
    else if (pp->resize.method == "Downscale (Better)" ||
             pp->resize.method == "Downscale (Faster)")
    {
        method->set_active (5);
    }

    wDirty = false;
    hDirty = false;

    if (pedited) {
        wDirty = pedited->resize.width;
        hDirty = pedited->resize.height;
        scale->setEditedState (pedited->resize.scale ? Edited : UnEdited);
        if (!pedited->resize.appliesTo)
            appliesTo->set_active (2);
        if (!pedited->resize.method)
            method->set_active (8);
        if (!pedited->resize.dataspec) 
            spec->set_active (4);
        enabled->set_inconsistent (!pedited->resize.enabled);
    }

    lastEnabled = pp->resize.enabled;

    scale->block(false);
    sconn.block (false);
    wconn.block (false);
    hconn.block (false);
    aconn.block (false);
    enableListener ();    
}

void Resize::write (ProcParams* pp, ParamsEdited* pedited) {
	int dataSpec = spec->get_active_row_number();

    pp->resize.scale  = scale->getValue();

    pp->resize.appliesTo = "Cropped area";
    if (appliesTo->get_active_row_number() == 0)
        pp->resize.appliesTo = "Cropped area";
    else if (appliesTo->get_active_row_number() == 1)
        pp->resize.appliesTo = "Full image";

    pp->resize.method = "Bicubic";
    if (method->get_active_row_number() == 0) 
        pp->resize.method = "Nearest";
    else if (method->get_active_row_number() == 1) 
        pp->resize.method = "Bilinear";
    else if (method->get_active_row_number() == 2) 
        pp->resize.method = "Bicubic";
    else if (method->get_active_row_number() == 3) 
        pp->resize.method = "Bicubic (Softer)";
    else if (method->get_active_row_number() == 4) 
        pp->resize.method = "Bicubic (Sharper)";
    else if (method->get_active_row_number() == 5) 
        pp->resize.method = "Lanczos";
        
    pp->resize.dataspec = dataSpec;
    pp->resize.width = w->get_value_as_int ();
    pp->resize.height = h->get_value_as_int ();
    pp->resize.enabled = enabled->get_active ();
    //printf("  L:%d   H:%d\n", pp->resize.width, pp->resize.height);

    if (pedited) {
        pedited->resize.enabled   = !enabled->get_inconsistent();
        pedited->resize.dataspec  = dataSpec != 4;
        pedited->resize.appliesTo = appliesTo->get_active_row_number() != 2;
        pedited->resize.method    = method->get_active_row_number() != 8;
        if (pedited->resize.dataspec) {
            pedited->resize.scale     = scale->getEditedState ();
            pedited->resize.width     = wDirty;
            pedited->resize.height    = hDirty;
        }
        else {
            pedited->resize.scale     = false;
            pedited->resize.width     = false;
            pedited->resize.height    = false;
        }
    }
}

void Resize::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    scale->setDefault (defParams->resize.scale);

    if (pedited) 
        scale->setDefaultEditedState (pedited->resize.scale ? Edited : UnEdited);
    else 
        scale->setDefaultEditedState (Irrelevant);
}

void Resize::adjusterChanged (Adjuster* a, double newval) {

    if (!batchMode) {
        wconn.block (true);
        hconn.block (true);
        h->set_value ((croph && appliesTo->get_active_row_number()==0 ? croph : maxh) * a->getValue ());
        w->set_value ((cropw && appliesTo->get_active_row_number()==0 ? cropw : maxw) * a->getValue ());
        wconn.block (false);
        hconn.block (false);
    }

    if (listener && (enabled->get_active () || batchMode)) 
        listener->panelChanged (EvResizeScale, Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(2), scale->getValue()));
}

int Resize::getComputedWidth() {

	if (cropw && appliesTo->get_active_row_number()==0)
		// we use the crop dimensions
		return (int)((double)(cropw) * (h->get_value()/(double)(croph)) + 0.5);
	else
		// we use the image dimensions
		return (int)((double)(maxw) * (h->get_value()/(double)(maxh)) + 0.5);
}

int Resize::getComputedHeight() {

	if (croph && appliesTo->get_active_row_number()==0)
		// we use the crop dimensions
		return (int)((double)(croph) * (w->get_value()/(double)(cropw)) + 0.5);
	else
		// we use the image dimensions
		return (int)((double)(maxh) * (w->get_value()/(double)(maxw)) + 0.5);
}

void Resize::appliesToChanged () {

	//printf("\nPASSAGE EN MODE \"%s\"\n\n", appliesTo->get_active_text().c_str());
	setDimensions();
    if (listener && (enabled->get_active () || batchMode)) {
    	//printf("Appel du listener\n");
        listener->panelChanged (EvResizeAppliesTo, appliesTo->get_active_text());
    }
}

void Resize::methodChanged () {

    if (listener && (enabled->get_active () || batchMode)) 
        listener->panelChanged (EvResizeMethod, method->get_active_text());
}

void Resize::update (bool isCropped, int cw, int ch, int ow, int oh) {

	// updating crop values now
    if (isCropped) {
		cropw = cw;
		croph = ch;
    }
    else {
		cropw = 0;
		croph = 0;
    }

    // updating the full image dimensions
    if (ow && oh) {
        maxw = ow;
        maxh = oh;
    }
    // updating the GUI synchronously
    setDimensions();
}

void Resize::sizeChanged (int mw, int mh, int ow, int oh) {

	// updating max values now
    maxw = ow;
    maxh = oh;

    // updating the GUI synchronously
    setDimensions();
}

void Resize::setDimensions () {

    int refw, refh;

    wconn.block (true);
    hconn.block (true);
    scale->block(true);

    if (appliesTo->get_active_row_number()==0 && cropw) {
    	// Applies to Cropped area
    	refw = cropw;
    	refh = croph;
    }
    else {
    	// Applies to Full image or crop is disabled
    	refw = maxw;
    	refh = maxh;
    }
    w->set_range (32, 4*refw);
    h->set_range (32, 4*refh);

	double tmpScale;
    switch (spec->get_active_row_number()) {
    case (0):	// Scale mode
        w->set_value((double)((int)( (double)(refw) * scale->getValue() + 0.5) ));
        h->set_value((double)((int)( (double)(refh) * scale->getValue() + 0.5) ));
        break;
    case (1):	// Width mode
		tmpScale = w->get_value() / (double)refw;
		scale->setValue (tmpScale);
		h->set_value((double)((int)( (double)(refh) * tmpScale + 0.5) ));
        break;
    case (2):	// Height mode
		tmpScale = h->get_value() / (double)refh;
		scale->setValue (tmpScale);
		w->set_value((double)((int)( (double)(refw) * tmpScale + 0.5) ));
		break;
    case (3): {	// Bounding box mode
    	double wSliderValue = w->get_value();
    	double hSliderValue = h->get_value();
    	if ( (wSliderValue/hSliderValue) < ((double)refw/(double)refh)) {
    		tmpScale = wSliderValue / (double)refw;
    	}
    	else {
    		tmpScale = hSliderValue / (double)refh;
    	}
		scale->setValue (tmpScale);
        break;
		}
    default:
    	break;
    }

    scale->block(false);
    wconn.block (false);
    hconn.block (false);
}

void Resize::fitBoxScale() {
	double tmpScale;
	double neww = w->get_value ();
	double newh = h->get_value ();

	if (cropw && appliesTo->get_active_row_number()==0) {
		// we use the crop dimensions
		if (((double)(cropw) / (double)(croph)) > (neww / newh)) {
			// the new scale is given by the image width
			tmpScale = neww / (double)(cropw);
		}
		else {
			// the new scale is given by the image height
			tmpScale = newh / (double)(croph);
		}
	}
	else {
		// we use the image dimensions
		if (((double)(maxw) / (double)(maxh)) > (neww / newh)) {
			// the new scale is given by the image width
			tmpScale = neww / (double)(maxw);
		}
		else {
			// the new scale is given by the image height
			tmpScale = newh / (double)(maxh);
		}
	}
	scale->setValue (tmpScale);
}

void Resize::entryWChanged () {

    wDirty = true;

    // updating width
    if (!batchMode) {
    	if (spec->get_active_row_number() == 3) {
    		// Fit box mode
    		fitBoxScale();
    	}
    	else {
    		// Other modes
			hconn.block (true);
			scale->block (true);

			h->set_value ((double)(getComputedHeight()));
			scale->setValue (w->get_value () / (cropw && appliesTo->get_active_row_number()==0 ? (double)cropw : (double)maxw));

			scale->block (false);
			hconn.block (false);
    	}
    }

    if (listener) {
    	if (spec->get_active_row_number() == 3)
    		notifyBBox();
    	else {
			if (enabled->get_active () || batchMode)
				listener->panelChanged (EvResizeWidth, Glib::ustring::format (w->get_value_as_int()));
    	}
    }
}

void Resize::entryHChanged () {

    hDirty = true;

    if (!batchMode && listener) {
    	if (spec->get_active_row_number() == 3) {
    		// Fit box mode
    		fitBoxScale();
    	}
    	else {
    		// Other modes
            wconn.block (true);
			scale->block (true);

            w->set_value ((double)(getComputedWidth()));
            scale->setValue (h->get_value () / (croph && appliesTo->get_active_row_number()==0 ? (double)croph : (double)maxh));

			scale->block (false);
            wconn.block (false);
    	}
    }

    if (listener) {
    	if (spec->get_active_row_number() == 3)
    		notifyBBox();
    	else {
    		if (enabled->get_active () || batchMode)
				listener->panelChanged (EvResizeHeight, Glib::ustring::format (h->get_value_as_int()));
    	}
    }
}

void Resize::specChanged () {

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        scale->sliderChanged ();
        break;
    case (1):
        // Width mode
        w->set_value((double)(getComputedWidth()));
        entryWChanged ();
        break;
    case (2):
        // Height mode
        h->set_value((double)(getComputedHeight()));
        entryHChanged ();
        break;
    case (3):
        // Bounding box mode
        notifyBBox();
        break;
    default:
        break;
    }
    updateGUI();
}

void Resize::updateGUI () {

    removeIfThere (this, scale, false);
    removeIfThere (this, sizeBox, false);

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        pack_start (*scale, Gtk::PACK_SHRINK, 4);
        break;
    case (1):
        // Width mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        w->set_sensitive (true);
        h->set_sensitive (false);
        break;
    case (2):
        // Height mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        w->set_sensitive (false);
        h->set_sensitive (true);
        break;
    case (3):
        // Bounding box mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        w->set_sensitive (true);
        h->set_sensitive (true);
        break;
    default:
        break;
    }
}

void Resize::notifyBBox() {
    if (listener && (enabled->get_active () || batchMode))
        listener->panelChanged (EvResizeBoundingBox, Glib::ustring::compose("(%1x%2)",(int)w->get_value(), (int)h->get_value() ));
}

void Resize::setBatchMode (bool batchMode) {

	method->append_text (M("GENERAL_UNCHANGED"));
	spec->append_text (M("GENERAL_UNCHANGED"));
    ToolPanel::setBatchMode (batchMode);
    scale->showEditedCB ();
}

void Resize::enabledToggled () {
    
    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaConn.block (true);
            enabled->set_active (false);
            enaConn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvResizeEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvResizeEnabled, M("GENERAL_DISABLED"));
    }
}

