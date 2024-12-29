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
#include "lensgeom.h"

#include <iostream>

#include "eventmapper.h"
#include "guiutils.h"
#include "rtimage.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring LensGeometry::TOOL_NAME = "lensgeom";

LensGeometry::LensGeometry () : FoldableToolPanel(this, TOOL_NAME, M("TP_LENSGEOM_LABEL")), rlistener(nullptr), lastFill(false)
{
    auto m = ProcEventMapper::getInstance();
    EvTransScale  = m->newEvent(TRANSFORM, "HISTORY_MSG_TRANS_SCALE");
    EvTransMethod = m->newEvent(TRANSFORM, "HISTORY_MSG_TRANS_METHOD");

    Gtk::Box* hb1 = Gtk::manage (new Gtk::Box ());
    hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_RAW_DMETHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    method = Gtk::manage (new MyComboBoxText ());
    method->append(M("TP_LENSGEOM_LOG"));
    method->append(M("TP_LENSGEOM_LIN"));
    method->set_active(0);
    hb1->pack_end (*method, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb1, Gtk::PACK_SHRINK, 4);

    scale= Gtk::manage (new Adjuster (M("TP_LENSGEOM_SCALE"), 0.1, 10, 0.01, 1));
    scale->setAdjusterListener (this);
    scale->setLogScale(300, 0.1);
    pack_start (*scale);

    fill = Gtk::manage (new Gtk::CheckButton (M("TP_LENSGEOM_FILL")));
    pack_start (*fill);

    autoCrop = Gtk::manage (new Gtk::Button (M("TP_LENSGEOM_AUTOCROP")));
    autoCrop->set_image (*Gtk::manage (new RTImage ("crop-auto-small", Gtk::ICON_SIZE_BUTTON)));
    autoCrop->get_style_context()->add_class("independent");
    pack_start (*autoCrop, Gtk::PACK_SHRINK, 2);

    method->connect(method->signal_changed().connect(sigc::mem_fun(*this, &LensGeometry::methodChanged)));
    autoCrop->signal_pressed().connect(sigc::mem_fun(*this, &LensGeometry::autoCropPressed));
    fillConn = fill->signal_toggled().connect(sigc::mem_fun(*this, &LensGeometry::fillPressed));

    fill->set_active (true);
    scale->setEnabled(!fill->get_active());
    show_all ();
}

LensGeometry::~LensGeometry ()
{
    idle_register.destroy();
}

void LensGeometry::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
    method->block (true);
    method->set_active(pp->commonTrans.method == "log" ? 0 : 1);

    if (pedited) {
        if(!pedited->commonTrans.method) {
            method->set_active_text(M("GENERAL_UNCHANGED"));
        }

        fill->set_inconsistent (!pedited->commonTrans.autofill);
        scale->setEditedState (pedited->commonTrans.scale ? Edited : UnEdited);
    }

    fillConn.block (true);
    fill->set_active (pp->commonTrans.autofill);
    fillConn.block (false);
    autoCrop->set_sensitive (!pp->commonTrans.autofill);

    scale->setValue (pp->commonTrans.scale);

    lastFill = pp->commonTrans.autofill;

    method->block (false);
    scale->setEnabled(!fill->get_active());
    enableListener ();
}

void LensGeometry::write (ProcParams* pp, ParamsEdited* pedited)
{
    int currentRow = method->get_active_row_number();
    if( currentRow >= 0 && method->get_active_text() != M("GENERAL_UNCHANGED")) {
        pp->commonTrans.method = currentRow == 0 ? "log" : "lin";
    }
    pp->commonTrans.autofill = fill->get_active ();
    pp->commonTrans.scale = scale->getValue ();

    if (pedited) {
        pedited->commonTrans.method = method->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->commonTrans.autofill = !fill->get_inconsistent();
        pedited->commonTrans.scale = scale->getEditedState();
    }
}

void LensGeometry::autoCropPressed ()
{

    if (rlistener) {
        rlistener->autoCropRequested ();
    }
}

void LensGeometry::adjusterChanged(Adjuster *a, double newval)
{
    if (listener) {
        if (a == scale) {
            listener->panelChanged (EvTransScale,
                    Glib::ustring::format(scale->getValue()));
        }
        else {
            if (settings->verbose) {
                std::cout << "Unknown adjuster given in LensGeometry::adjusterChanged, file " << __FILE__ << " line " << __LINE__ << std::endl;
            }
        }
    }
}

void LensGeometry::fillPressed ()
{

    if (batchMode) {
        if (fill->get_inconsistent()) {
            fill->set_inconsistent (false);
            fillConn.block (true);
            fill->set_active (false);
            fillConn.block (false);
        } else if (lastFill) {
            fill->set_inconsistent (true);
        }

        lastFill = fill->get_active ();
    } else {
        autoCrop->set_sensitive (!fill->get_active());
    }

    if (listener) {
        if (fill->get_active ()) {
            listener->panelChanged (EvTransAutoFill, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvTransAutoFill, M("GENERAL_DISABLED"));
        }
    }
    scale->setEnabled(!fill->get_active());
}

void LensGeometry::methodChanged ()
{

    if (listener && method->get_active_row_number() >= 0) {
        listener->panelChanged(EvTransMethod, method->get_active_text());
    }
}

void LensGeometry::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    removeIfThere (this, autoCrop);
    scale->showEditedCB ();
}

