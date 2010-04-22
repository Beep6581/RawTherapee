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
#include <whitebalance.h>
#include <iomanip>

#define MINTEMP 1200
#define MAXTEMP 12000
#define MINGREEN 0.02
#define MAXGREEN 5.0

extern Glib::ustring argv0;

using namespace rtengine;
using namespace rtengine::procparams;

WhiteBalance::WhiteBalance () : ToolPanel(), wbp(NULL), wblistener(NULL), tempAdd(false), greenAdd (false) {

  Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
  hbox->show ();
  Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_METHOD")));
  lab->show ();
  method = Gtk::manage (new Gtk::ComboBoxText ());
  method->show ();
  method->append_text (M("TP_WBALANCE_CAMERA"));
  method->append_text (M("TP_WBALANCE_AUTO"));
  method->append_text (M("TP_WBALANCE_CUSTOM"));
  method->set_active (0);
  hbox->pack_start (*lab, Gtk::PACK_SHRINK, 4);
  hbox->pack_start (*method);
  pack_start (*hbox, Gtk::PACK_SHRINK, 4);
  opt = 0;

  Gtk::HBox* spotbox = Gtk::manage (new Gtk::HBox ());
  spotbox->show ();

  spotbutton = Gtk::manage (new Gtk::Button (M("TP_WBALANCE_SPOTWB")));
  Gtk::Image* spotimg = Gtk::manage (new Gtk::Image (argv0+"/images/wbpicker16.png"));
  spotimg->show ();
  spotbutton->set_image (*spotimg);
  spotbutton->show ();

  spotbox->pack_start (*spotbutton);

  Gtk::Label* slab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_SIZE")));
  slab->show ();

  spotsize = Gtk::manage (new Gtk::ComboBoxText ());
  spotsize->show ();
  spotsize->append_text ("2");
  spotsize->append_text ("4");
  spotsize->append_text ("8");
  spotsize->append_text ("16");
  spotsize->append_text ("32");
  spotsize->set_active (2);

  spotbox->pack_end (*spotsize, Gtk::PACK_SHRINK, 4);
  spotbox->pack_end (*slab, Gtk::PACK_SHRINK, 4);

  pack_start (*spotbox, Gtk::PACK_SHRINK, 4);

  temp = Gtk::manage (new Adjuster (M("TP_WBALANCE_TEMPERATURE"), MINTEMP, MAXTEMP, 1, 4750));
  green = Gtk::manage (new Adjuster (M("TP_WBALANCE_GREEN"), MINGREEN, MAXGREEN, 0.001, 1.2));
  temp->show ();
  green->show ();

  pack_start (*temp);
  pack_start (*green);

  temp->setAdjusterListener (this);
  green->setAdjusterListener (this);

  spotbutton->signal_pressed().connect( sigc::mem_fun(*this, &WhiteBalance::spotPressed) );
  methconn = method->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::optChanged) );
  spotsize->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::spotSizeChanged) );
}

void WhiteBalance::adjusterChanged (Adjuster* a, double newval) {

    if (method->get_active_row_number()!=2) {
        disableListener ();    
        method->set_active (2);
        enableListener ();
    }

    if (listener) {
        if (a==temp) 
            listener->panelChanged (EvWBTemp, Glib::ustring::format ((int)a->getValue()));
        else if (a==green) 
            listener->panelChanged (EvWBGreen, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
    }
}

void WhiteBalance::optChanged () {

    if (opt!=method->get_active_row_number()) {
        opt = method->get_active_row_number();
        Glib::ustring meth = M("TP_WBALANCE_CUSTOM");
        if (opt==0 && wbp) {
            double ctemp; double cgreen;
            wbp->getCamWB (ctemp, cgreen);
            temp->setValue (tempAdd ? 0.0 : (int)ctemp);
            green->setValue (greenAdd ? 0.0 : cgreen);
            meth = M("TP_WBALANCE_CAMERA");
            if (batchMode) {
                temp->setEditedState (UnEdited);
                green->setEditedState (UnEdited);
            }
        }
        else if (opt==1 && wbp) {
            double ctemp; double cgreen;
            wbp->getAutoWB (ctemp, cgreen);
            temp->setValue (tempAdd ? 0.0 : (int)ctemp);
            green->setValue (greenAdd ? 0.0 : cgreen);
            meth = M("TP_WBALANCE_AUTO");
            if (batchMode) {
                temp->setEditedState (UnEdited);
                green->setEditedState (UnEdited);
            }
        }
        else if (opt==3) {
            meth = "(Unchanged)";
            temp->setEditedState (UnEdited);
            green->setEditedState (UnEdited);
        }
        if (listener) 
            listener->panelChanged (EvWBMethod, meth);
    }
}

void WhiteBalance::spotPressed () {

  if (wblistener)
    wblistener->spotWBRequested (getSize());
}

void WhiteBalance::spotSizeChanged () {

  if (wblistener)
    wblistener->spotWBRequested (getSize());
}

void WhiteBalance::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    methconn.block (true);

    if (pedited) {
        temp->setEditedState (UnEdited);
        green->setEditedState (UnEdited);
    }
    
    if (pedited && !pedited->wb.method) {
        method->set_active (3);
        opt = 3;
    }
    else {
        if (pp->wb.method == "Camera") {
            method->set_active (0);
            if (wbp) {
                double ctemp; double cgreen;
                wbp->getCamWB (ctemp, cgreen);
                temp->setValue (tempAdd ? 0.0 : (int)ctemp);
                green->setValue (greenAdd ? 0.0 : cgreen);
            }
            opt = 0;
        }
        else if (pp->wb.method == "Auto") {
            method->set_active (1);
            if (wbp) {
                double ctemp; double cgreen;
                wbp->getAutoWB (ctemp, cgreen);
                temp->setValue (tempAdd ? 0.0 : (int)ctemp);
                green->setValue (greenAdd ? 0.0 : cgreen);
            }
            opt = 1;
        }
        else if (pp->wb.method == "Custom") {
            method->set_active (2);
            temp->setValue (pp->wb.temperature);
            green->setValue (pp->wb.green);
            opt = 2;
            if (pedited) {
                temp->setEditedState (pedited->wb.temperature ? Edited : UnEdited);
                green->setEditedState (pedited->wb.green ? Edited : UnEdited);
            }
        }
    }
    methconn.block (false);
    enableListener ();
}

void WhiteBalance::write (ProcParams* pp, ParamsEdited* pedited) {
    
    if (pedited) {
        pedited->wb.temperature = temp->getEditedState ();
        pedited->wb.green = green->getEditedState ();
        pedited->wb.method = method->get_active_row_number()!=3;
    }

    if (method->get_active_row_number()==0)
        pp->wb.method = "Camera";
    else if (method->get_active_row_number()==1)
        pp->wb.method = "Auto";
    else if (method->get_active_row_number()>=2)
        pp->wb.method = "Custom";      
    
    pp->wb.temperature = (int)temp->getValue ();
    pp->wb.green = green->getValue ();
    
}

void WhiteBalance::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    if (wbp && defParams->wb.method == "Camera") {
        double ctemp; double cgreen;
        wbp->getCamWB (ctemp, cgreen);
        temp->setDefault (tempAdd ? 0 : (int)ctemp);
        green->setDefault (greenAdd ? 0 : cgreen);
    }
    else if (wbp && defParams->wb.method == "Auto") {
        double ctemp; double cgreen;
        wbp->getAutoWB (ctemp, cgreen);
        temp->setDefault (tempAdd ? 0 : (int)ctemp);
        green->setDefault (greenAdd ? 0 : cgreen);
    }
    else if (defParams->wb.method == "Custom") {
        temp->setDefault (defParams->wb.temperature);
        green->setDefault (defParams->wb.green);
    }
    if (pedited) {
        temp->setDefaultEditedState (pedited->wb.temperature ? Edited : UnEdited);
        green->setDefaultEditedState (pedited->wb.green ? Edited : UnEdited);
    }
    else {
        temp->setDefaultEditedState (Irrelevant);
        green->setDefaultEditedState (Irrelevant);
    }
}

void WhiteBalance::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    temp->showEditedCB ();
    green->showEditedCB ();
    method->append_text ("(Unchanged)");
}

int WhiteBalance::getSize () {

    return atoi(spotsize->get_active_text().c_str());
}

void WhiteBalance::setWB (int vtemp, double vgreen) {

    disableListener ();
    temp->setValue (vtemp);
    green->setValue (vgreen);
    method->set_active (2);
    temp->setEditedState (Edited);
    green->setEditedState (Edited);
    enableListener ();

    if (listener) 
        listener->panelChanged (EvWBTemp, Glib::ustring::compose("%1, %2", (int)temp->getValue(), Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), green->getValue())));
}
void WhiteBalance::setAdjusterBehavior (bool btempadd, bool bgreenadd) {

    if (!tempAdd && btempadd)
        temp->setLimits (-4000, +4000, 1, 0);
    else if (tempAdd && !btempadd)
        temp->setLimits (MINTEMP, MAXTEMP, 1, 4750);

    if (!greenAdd && bgreenadd)
        green->setLimits (-1.0, +1.0, 0.001, 0);
    else if (greenAdd && bgreenadd)
        green->setLimits (MINGREEN, MAXGREEN, 0.001, 1.2);
    
    tempAdd = btempadd;
    greenAdd = bgreenadd;
}

