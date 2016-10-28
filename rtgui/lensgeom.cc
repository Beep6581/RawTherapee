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
#include "lensgeom.h"
#include "guiutils.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

LensGeometry::LensGeometry () : FoldableToolPanel(this, "lensgeom", M("TP_LENSGEOM_LABEL")), rlistener(nullptr)
{

    fill = Gtk::manage (new Gtk::CheckButton (M("TP_LENSGEOM_FILL")));
    pack_start (*fill);

    autoCrop = Gtk::manage (new Gtk::Button ());
    autoCrop->set_tooltip_text(M("TP_LENSGEOM_AUTOCROP"));
    autoCrop->set_image (*Gtk::manage (new RTImage ("crop-auto.png")));
    pack_start (*autoCrop, Gtk::PACK_SHRINK, 2);

    packBox = Gtk::manage (new ToolParamBlock ());
    pack_start (*packBox);

    autoCrop->signal_pressed().connect( sigc::mem_fun(*this, &LensGeometry::autoCropPressed) );
    fillConn = fill->signal_toggled().connect( sigc::mem_fun(*this, &LensGeometry::fillPressed) );

    fill->set_active (true);
    show_all ();
}

LensGeometry::~LensGeometry ()
{
    g_idle_remove_by_data(this);
}

void LensGeometry::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();

    if (pedited) {
        fill->set_inconsistent (!pedited->commonTrans.autofill);
    }

    fillConn.block (true);
    fill->set_active (pp->commonTrans.autofill);
    fillConn.block (false);
    autoCrop->set_sensitive (!pp->commonTrans.autofill);

    lastFill = pp->commonTrans.autofill;

    enableListener ();
}

void LensGeometry::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->commonTrans.autofill   = fill->get_active ();

    if (pedited) {
        pedited->commonTrans.autofill   = !fill->get_inconsistent();
    }
}

void LensGeometry::autoCropPressed ()
{

    if (rlistener) {
        rlistener->autoCropRequested ();
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
}

void LensGeometry::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    removeIfThere (this, autoCrop);
}

void LensGeometry::disableAutoFillIfActive ()
{
    g_idle_add(doDisableAutoFillIfActive, this);
}

int LensGeometry::doDisableAutoFillIfActive (void* data)
{
    GThreadLock lock; // Is this really needed?

    LensGeometry* const instance = static_cast<LensGeometry*>(data);

    if (!instance->batchMode) {
        if (instance->fill->get_active()) {
            instance->fillConn.block (true);
            instance->fill->set_active(false);

            if (instance->listener) {
                instance->listener->panelChanged (EvTransAutoFill, M("GENERAL_DISABLED"));
            }

            instance->fillConn.block (false);
        }
    }

    return 0;
}
