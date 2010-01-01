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
#include <colorshift.h>
#include <iomanip>

using namespace rtengine;
using namespace rtengine::procparams;

ColorShift::ColorShift () : ToolPanel(), aAdd(false), bAdd(false) {

  ashift   = Gtk::manage (new Adjuster (M("TP_COLORSHIFT_GREENMAGENTA"), -25, 25, 0.1, 0));
  pack_start (*ashift);

  bshift   = Gtk::manage (new Adjuster (M("TP_COLORSHIFT_BLUEYELLOW"), -25, 25, 0.1, 0));
  pack_start (*bshift);

  ashift->setAdjusterListener (this);
  bshift->setAdjusterListener (this);

  show_all ();
}

void ColorShift::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    if (pedited) {
        ashift->setEditedState (pedited->colorShift.a ? Edited : UnEdited);
        bshift->setEditedState (pedited->colorShift.b ? Edited : UnEdited);
    }

    ashift->setValue (pp->colorShift.a);
    bshift->setValue (pp->colorShift.b);
    
    enableListener ();
}

void ColorShift::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->colorShift.a = ashift->getValue ();
    pp->colorShift.b = bshift->getValue ();

    if (pedited) {
        pedited->colorShift.a = ashift->getEditedState ();
        pedited->colorShift.b = bshift->getEditedState ();
    }
}

void ColorShift::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    ashift->setDefault (defParams->colorShift.a);
    bshift->setDefault (defParams->colorShift.b);

    if (pedited) {
        ashift->setDefaultEditedState (pedited->colorShift.a ? Edited : UnEdited);
        bshift->setDefaultEditedState (pedited->colorShift.b ? Edited : UnEdited);
    }
    else {
        ashift->setDefaultEditedState (Irrelevant);
        bshift->setDefaultEditedState (Irrelevant);
    }
}

void ColorShift::adjusterChanged (Adjuster* a, double newval) {

    if (!listener)
        return;

    if (a==ashift) 
        listener->panelChanged (EvCShiftA, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
    else if (a==bshift) 
        listener->panelChanged (EvCShiftB, Glib::ustring::format (std::setw(2), std::fixed, std::setprecision(1), a->getValue()));
}

void ColorShift::setAdjusterBehavior (bool baadd, bool bbadd) {

    if (!aAdd && baadd || aAdd && !baadd) 
        ashift->setLimits (-25, 25, 0.1, 0);

    if (!bAdd && bbadd || bAdd && !bbadd) 
        bshift->setLimits (-25, 25, 0.1, 0);
    
    aAdd = baadd;
    bAdd = bbadd;
}

void ColorShift::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    ashift->showEditedCB ();
    bshift->showEditedCB ();
}
