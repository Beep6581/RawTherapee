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
#include <resize.h>
#include <iomanip>
#include <guiutils.h>

using namespace rtengine;
using namespace rtengine::procparams;

Resize::Resize () : maxw(100000), maxh(100000) {

    enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
    pack_start(*enabled);
    pack_start(*Gtk::manage (new  Gtk::HSeparator()), Gtk::PACK_SHRINK, 2);

    Gtk::Table* combos = Gtk::manage (new Gtk::Table (2, 2));

	method = Gtk::manage (new Gtk::ComboBoxText ());
	method->append_text (M("TP_RESIZE_NEAREST"));
	method->append_text (M("TP_RESIZE_BILINEAR"));
	method->append_text (M("TP_RESIZE_BICUBIC"));
	method->append_text (M("TP_RESIZE_BICUBICSF"));
	method->append_text (M("TP_RESIZE_BICUBICSH"));
	method->set_active (0);

    combos->attach (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_METHOD"))), 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*method, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

	spec = Gtk::manage (new Gtk::ComboBoxText ());
	spec->append_text ("Scale");
	spec->append_text ("Width");
    spec->append_text ("Height");
	method->set_active (0);
    
    combos->attach (*Gtk::manage (new Gtk::Label ("Specify:")), 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*spec, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

	pack_start (*combos, Gtk::PACK_SHRINK, 4);

    scale = new Adjuster (M("TP_RESIZE_SCALE"), 0.2, 4, 0.01, 1);
    scale->setAdjusterListener (this); 

    pack_start (*scale, Gtk::PACK_SHRINK, 4);

    sizeBox = Gtk::manage (new Gtk::VBox ());

    Gtk::HBox* sbox = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* wbox = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
    w = Gtk::manage (new Gtk::SpinButton ());
    h = Gtk::manage (new Gtk::SpinButton ());
    wbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_W"))), Gtk::PACK_SHRINK, 4);
    wbox->pack_start (*w);
    hbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_H"))), Gtk::PACK_SHRINK, 4);
    hbox->pack_start (*h);
    sbox->pack_start (*wbox);
    sbox->pack_start (*hbox);
    
    sizeBox->pack_start (*sbox, Gtk::PACK_SHRINK, 4);
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
	method->signal_changed().connect ( sigc::mem_fun(*this, &Resize::methodChanged) );
	spec->signal_changed().connect ( sigc::mem_fun(*this, &Resize::specChanged) );
	enaConn = enabled->signal_toggled().connect ( sigc::mem_fun(*this, &Resize::enabledToggled) );

    show_all();
}

Resize::~Resize () {

    delete scale;
    delete sizeBox;
}

void Resize::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    wconn.block (true);
    hconn.block (true);

    scale->setValue (pp->resize.scale);
    w->set_value (pp->resize.width);
    h->set_value (pp->resize.height);
    enabled->set_active (pp->resize.enabled);
    spec->set_active (pp->resize.dataspec);

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

    wDirty = false;
    hDirty = false;

    if (pedited) {
        wDirty = pedited->resize.width;
        hDirty = pedited->resize.height;
        scale->setEditedState (pedited->resize.scale ? Edited : UnEdited);
        if (!pedited->resize.method)
            method->set_active (5);
        if (!pedited->resize.dataspec) 
            spec->set_active (3);
        enabled->set_inconsistent (!pedited->resize.enabled);
    }

    lastEnabled = pp->resize.enabled;

    wconn.block (false);
    hconn.block (false);
    enableListener ();    
}

void Resize::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->resize.scale  = scale->getValue ();
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
        
    pp->resize.dataspec = spec->get_active_row_number();
    pp->resize.width = round (w->get_value ());
    pp->resize.height = round(h->get_value ());
    pp->resize.enabled = enabled->get_active ();

    if (pedited) {
        pedited->resize.enabled   = !enabled->get_inconsistent();
        pedited->resize.dataspec  = spec->get_active_row_number() != 3;
        pedited->resize.method    = method->get_active_row_number() != 5;
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
        h->set_value (maxh * a->getValue ());
        w->set_value (maxw * a->getValue ());
        wconn.block (false);
        hconn.block (false);
    }

    if (listener && (enabled->get_active () || batchMode)) 
        listener->panelChanged (EvResizeScale, Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(4), scale->getValue()));
}

void Resize::methodChanged () {

    if (listener && (enabled->get_active () || batchMode)) 
        listener->panelChanged (EvResizeMethod, method->get_active_text());
}

struct setrdimparams {
    Resize* resize;
    int mw;
    int mh;
    int ow;
    int oh;
};

int setrdim (void* data) {

    gdk_threads_enter ();
    setrdimparams* params = (setrdimparams*)data;
    params->resize->setDimensions (params->mw, params->mh, params->ow, params->oh);
    delete params;
    gdk_threads_leave ();
    return 0;
}

void Resize::sizeChanged (int mw, int mh, int ow, int oh) {

    setrdimparams* params = new setrdimparams;
    params->mw = mw;
    params->mh = mh;
    params->ow = ow;
    params->oh = oh;
    params->resize = this;
    g_idle_add (setrdim, params);
}

void Resize::setDimensions (int mw, int mh, int ow, int oh) {

    maxw = ow;
    maxh = oh;
	
    wconn.block (true);
    hconn.block (true);

    w->set_range (32, 4*maxw);
    h->set_range (32, 4*maxh);

    wconn.block (false);
    hconn.block (false);
}

void Resize::entryWChanged () {

    wDirty = true;

    if (!batchMode && listener) {
        hconn.block (true);
        h->set_value (w->get_value () * maxh / maxw);
        hconn.block (false);
        scale->setValue (w->get_value () / maxw);
    }

    if (listener && (enabled->get_active () || batchMode)) 
        listener->panelChanged (EvResizeWidth, Glib::ustring::format ((int)w->get_value()));

    
}

void Resize::entryHChanged () {

    hDirty = true;

    if (!batchMode && listener) {
        wconn.block (true);
        w->set_value (h->get_value () * maxw / maxh);
        wconn.block (false);
        scale->setValue (h->get_value () / maxh);
    }

    if (listener && (enabled->get_active () || batchMode)) 
        listener->panelChanged (EvResizeHeight, Glib::ustring::format ((int)h->get_value()));
}

void Resize::specChanged () {

    removeIfThere (this, scale, false);
    removeIfThere (this, sizeBox, false);

    if (spec->get_active_row_number() == 0) {
        pack_start (*scale, Gtk::PACK_SHRINK, 4);
        scale->sliderChanged ();
    }
    else if (spec->get_active_row_number() == 1) {
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        w->set_sensitive (true);
        h->set_sensitive (false);
        entryWChanged ();
    }
    else if (spec->get_active_row_number() == 2) {
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        h->set_sensitive (true);
        w->set_sensitive (false);
        entryHChanged ();
    }
}

void Resize::setBatchMode (bool batchMode) {

	method->append_text ("(Unchanged)");
	spec->append_text ("(Unchanged)");
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

