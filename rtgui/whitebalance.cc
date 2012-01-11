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
#include "whitebalance.h"
#include <iomanip>
#include "rtimage.h"
#include "options.h"
#include "../rtengine/safegtk.h"

#define MINTEMP 2000   //1200
#define MAXTEMP 25000  //12000
#define MINGREEN 0.02
#define MAXGREEN 5.0

extern Glib::ustring argv0;

using namespace rtengine;
using namespace rtengine::procparams;

Glib::RefPtr<Gdk::Pixbuf> WhiteBalance::wbPixbufs[rtengine::procparams::WBT_CUSTOM+1];
/*
Glib::RefPtr<Gdk::Pixbuf> WhiteBalance::wbCameraPB, WhiteBalance::wbAutoPB, WhiteBalance::wbSunPB, WhiteBalance::wbTungstenPB,
                          WhiteBalance::wbCloudyPB, WhiteBalance::wbShadePB, WhiteBalance::wbFluorescentPB, WhiteBalance::wbLampPB,
                          WhiteBalance::wbFlashPB, WhiteBalance::wbLedPB, WhiteBalance::wbCustomPB;
*/

void WhiteBalance::init () {
    wbPixbufs[WBT_CAMERA]      = safe_create_from_file("wb-camera.png");
    wbPixbufs[WBT_AUTO]        = safe_create_from_file("wb-auto.png");
    wbPixbufs[WBT_DAYLIGHT]    = safe_create_from_file("wb-sun.png");
    wbPixbufs[WBT_CLOUDY]      = safe_create_from_file("wb-cloudy.png");
    wbPixbufs[WBT_SHADE]       = safe_create_from_file("wb-shade.png");
    wbPixbufs[WBT_TUNGSTEN]    = safe_create_from_file("wb-tungsten.png");
    wbPixbufs[WBT_FLUORESCENT] = safe_create_from_file("wb-fluorescent.png");
    wbPixbufs[WBT_LAMP]        = safe_create_from_file("wb-lamp.png");
    wbPixbufs[WBT_FLASH]       = safe_create_from_file("wb-flash.png");
    wbPixbufs[WBT_LED]         = safe_create_from_file("wb-led.png");
    wbPixbufs[WBT_CUSTOM]      = safe_create_from_file("wb-custom.png");
}

void WhiteBalance::cleanup () {
    for (unsigned int i=0; i<WBT_CUSTOM+1; i++) {
        wbPixbufs[i].reset();
    }
}

WhiteBalance::WhiteBalance () : Gtk::VBox(), FoldableToolPanel(this), wbp(NULL), wblistener(NULL) {

  Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
  hbox->show ();
  Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_METHOD")));
  lab->show ();

  // Create the Tree model
  refTreeModel = Gtk::TreeStore::create(methodColumns);
  // Create the Combobox
  method = Gtk::manage (new MyComboBox ());
  // Assign the model to the Combobox
  method->set_model(refTreeModel);

  custom_green = new double[WBParams::wbEntries.size()];
  enum WBTypes oldType = WBParams::wbEntries[0]->type;
  enum WBTypes currType;
  Gtk::TreeModel::Row row, childrow;
  for (unsigned int i=0; i<WBParams::wbEntries.size(); i++) {
      if (oldType != (currType = WBParams::wbEntries[i]->type)) {
          // New entry type
          if (currType == WBT_FLUORESCENT) {
              // Creating the Fluorescent subcategory header
              row = *(refTreeModel->append());
              row[methodColumns.colIcon] = wbPixbufs[currType];
              row[methodColumns.colLabel] = M("TP_WBALANCE_FLUO_HEADER");
              row[methodColumns.colId] = i+100;
          }
          if (currType == WBT_LAMP) {
              // Creating the Lamp subcategory header
              row = *(refTreeModel->append());
              row[methodColumns.colIcon] = wbPixbufs[currType];
              row[methodColumns.colLabel] = M("TP_WBALANCE_LAMP_HEADER");
              row[methodColumns.colId] = i+100;
          }
          if (currType == WBT_LED) {
              // Creating the LED subcategory header
              row = *(refTreeModel->append());
              row[methodColumns.colIcon] = wbPixbufs[currType];
              row[methodColumns.colLabel] = M("TP_WBALANCE_LED_HEADER");
              row[methodColumns.colId] = i+100;
          }
          if (currType == WBT_FLASH) {
              // Creating the Flash subcategory header
              row = *(refTreeModel->append());
              row[methodColumns.colIcon] = wbPixbufs[currType];
              row[methodColumns.colLabel] = M("TP_WBALANCE_FLASH_HEADER");
              row[methodColumns.colId] = i+100;
          }
      }
      if (currType == WBT_FLUORESCENT
       || currType == WBT_LAMP
       || currType == WBT_FLASH
       || currType == WBT_LED
      ) {
          childrow = *(refTreeModel->append(row.children()));
          childrow[methodColumns.colIcon] = wbPixbufs[currType];
          childrow[methodColumns.colLabel] = WBParams::wbEntries[i]->GUILabel;
          childrow[methodColumns.colId] = i;
      }
      else {
          row = *(refTreeModel->append());
          row[methodColumns.colIcon] = wbPixbufs[currType];
          row[methodColumns.colLabel] = WBParams::wbEntries[i]->GUILabel;
          row[methodColumns.colId] = i;
      }
      oldType = currType;

	  custom_green[i] = 1.0;
  }

  //Add the model columns to the Combo (which is a kind of view),
  //rendering them in the default way:
  method->pack_start(methodColumns.colIcon, false);
  method->pack_start(methodColumns.colLabel, true);

  method->set_active (0); // Camera
  method->show ();
  hbox->pack_start (*lab, Gtk::PACK_SHRINK, 4);
  hbox->pack_start (*method);
  pack_start (*hbox, Gtk::PACK_SHRINK, 4);
  opt = 0;

  Gtk::HBox* spotbox = Gtk::manage (new Gtk::HBox ());
  spotbox->show ();

  spotbutton = Gtk::manage (new Gtk::Button (M("TP_WBALANCE_SPOTWB")));
  Gtk::Image* spotimg = Gtk::manage (new RTImage ("gtk-color-picker-small.png"));
  spotimg->show ();
  spotbutton->set_image (*spotimg);
  spotbutton->show ();

  spotbox->pack_start (*spotbutton);

  Gtk::Label* slab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_SIZE")));
  slab->show ();

  spotsize = Gtk::manage (new MyComboBoxText ());
  spotsize->show ();
  spotsize->append_text ("2");  if (options.whiteBalanceSpotSize==2)  spotsize->set_active(0);
  spotsize->append_text ("4");  if (options.whiteBalanceSpotSize==4)  spotsize->set_active(1);
  spotsize->append_text ("8");  if (options.whiteBalanceSpotSize==8)  spotsize->set_active(2);
  spotsize->append_text ("16"); if (options.whiteBalanceSpotSize==16) spotsize->set_active(3);
  spotsize->append_text ("32"); if (options.whiteBalanceSpotSize==32) spotsize->set_active(4);

  spotbox->pack_end (*spotsize, Gtk::PACK_EXPAND_WIDGET, 4);
  spotbox->pack_end (*slab, Gtk::PACK_SHRINK, 4);

  pack_start (*spotbox, Gtk::PACK_SHRINK, 4);

  temp = Gtk::manage (new Adjuster (M("TP_WBALANCE_TEMPERATURE"), MINTEMP, MAXTEMP, 5, 4750));
  green = Gtk::manage (new Adjuster (M("TP_WBALANCE_GREEN"), MINGREEN, MAXGREEN, 0.001, 1.0));
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

WhiteBalance::~WhiteBalance () {
    delete custom_green;
}

void WhiteBalance::adjusterChanged (Adjuster* a, double newval) {

    int tVal = (int)temp->getValue();
    double gVal = green->getValue();

    Gtk::TreeModel::Row row = getActiveMethod();
    if (row == refTreeModel->children().end()) return;

    Glib::ustring colLabel = row[methodColumns.colLabel];
    WBEntry* ppMethod = findWBEntry (row[methodColumns.colLabel], WBLT_GUI);
    WBEntry* wbCustom = findWBEntry ("Custom", WBLT_PP);

    if (!ppMethod || ppMethod->ppLabel != wbCustom->ppLabel) {
        if (!ppMethod || a==temp || (ppMethod->type==WBT_CAMERA || ppMethod->type==WBT_AUTO) ) {
            methconn.block(true);
            opt = setActiveMethod(wbCustom->GUILabel);
            methconn.block(false);
        }
    }

    //cache custom WB setting to allow its recall
    if (a==temp)
        cache_customTemp (tVal);
    else
        cache_customGreen (gVal);

    if (listener) {
        if (a==temp)
            listener->panelChanged (EvWBTemp, Glib::ustring::format ((int)a->getValue()));
        else if (a==green)
            listener->panelChanged (EvWBGreen, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
    }
}

void WhiteBalance::optChanged () {

    Gtk::TreeModel::Row row = getActiveMethod();
    if (row == refTreeModel->children().end()) return;
    if (row[methodColumns.colId] >= 100) {
        // "Header" solutions are trapped ; the combo is then set to the previous value
        bool prevState = methconn.block(true);
        method->set_active(opt);
        methconn.block(prevState);
        return;
    }

    if (opt != row[methodColumns.colId]) {

        opt = row[methodColumns.colId];

        if (row[methodColumns.colLabel] == M("GENERAL_UNCHANGED")) {
            temp->setEditedState (UnEdited);
            green->setEditedState (UnEdited);
        }
        else {
            int methodId = findWBEntryId (row[methodColumns.colLabel], WBLT_GUI);
            WBEntry* currMethod = WBParams::wbEntries[methodId];

            switch (currMethod->type) {
            case WBT_CAMERA:
                if (wbp) {
                    double ctemp, cgreen;
                    wbp->getCamWB (ctemp, cgreen);
                    temp->setValue (temp->getAddMode() ? 0.0 : (int)ctemp);
                    green->setValue (green->getAddMode() ? 0.0 : cgreen);
                    if (batchMode) {
                        temp->setEditedState (UnEdited);
                        green->setEditedState (UnEdited);
                    }
                }
                break;
            case WBT_AUTO:
                if (wbp) {
                    double ctemp, cgreen;
                    wbp->getAutoWB (ctemp, cgreen);
                    if (ctemp != -1.0) {
                        temp->setValue (temp->getAddMode() ? 0.0 : (int)ctemp);
                        green->setValue (green->getAddMode() ? 0.0 : cgreen);
                    }
                    if (batchMode) {
                        temp->setEditedState (UnEdited);
                        green->setEditedState (UnEdited);
                    }
                }
                break;
            case WBT_CUSTOM:
                if (custom_temp>0){
                    temp->setValue (temp->getAddMode() ? 0.0 : custom_temp);
                }
                green->setValue (green->getAddMode() ? 0.0 : custom_green[methodId]);
                if (batchMode) {
                    temp->setEditedState (Edited);
                    green->setEditedState (Edited);
                }
                break;
            /* All other solution are the default cases
            case WBT_DAYLIGHT:
            case WBT_CLOUDY:
            case WBT_SHADE:
            case WBT_TUNGSTEN:
            case WBT_FLUORESCENT:
            case WBT_LAMP:
            case WBT_FLASH:
            case WBT_LED:*/
            default:
                // recall custom WB settings if it exists, set to 1.0 otherwise
                temp->setValue (temp->getAddMode() ? 0.0 : (double)(currMethod->temperature));
                green->setValue (green->getAddMode() ? 0.0 : custom_green[methodId]);
                if (batchMode) {
                    temp->setEditedState (Edited);
                    green->setEditedState (Edited);
                }
                break;
            }
        }

        if (listener)
            listener->panelChanged (EvWBMethod, row[methodColumns.colLabel]);
    }
}

void WhiteBalance::spotPressed () {

  if (wblistener)
    wblistener->spotWBRequested (getSize());
}

void WhiteBalance::spotSizeChanged () {
  options.whiteBalanceSpotSize=getSize();

  if (wblistener)
    wblistener->spotWBRequested (getSize());
}

void WhiteBalance::read (const ProcParams* pp, const ParamsEdited* pedited) {

    methconn.block (true);

    if (pedited) {
        // By default, temperature and green are said "UnEdited", but it may change later
        temp->setEditedState (UnEdited);
        green->setEditedState (UnEdited);
    }
    
    if (pedited && !pedited->wb.method) {
        opt = setActiveMethod(M("GENERAL_UNCHANGED"));
    }
    else {
        WBEntry* wbValues = findWBEntry(pp->wb.method, WBLT_PP);
        if (!wbValues)
            wbValues = findWBEntry("Camera", WBLT_PP);

        opt = setActiveMethod(wbValues->GUILabel);

        // temperature is reset to the associated temperature, or 0.0 if addMode is set.
        switch (wbValues->type) {
        case WBT_CUSTOM:
            temp->setValue (pp->wb.temperature);
            green->setValue (pp->wb.green);
            if (pedited) {
                // The user may have changed the temperature and green value
                temp->setEditedState (pedited->wb.temperature ? Edited : UnEdited);
                green->setEditedState (pedited->wb.green ? Edited : UnEdited);
            }
            //cache_customWB (pp->wb.temperature, pp->wb.green);
            break;
        case WBT_CAMERA:
            if (wbp) {
                double ctemp; double cgreen;
                wbp->getCamWB (ctemp, cgreen);

                // Set the camera's temperature value, or 0.0 if in ADD mode
                temp->setValue (temp->getAddMode() ? 0.0 : ctemp);
                // Set the camera's green value, or 0.0 if in ADD mode
                green->setValue (green->getAddMode() ? 0.0 : cgreen);

                //cache_customWB ((int)ctemp, cgreen); // this will be used to set initial Custom WB setting
            }
            break;
        case WBT_AUTO:
            if (wbp) {
                double ctemp; double cgreen;
                wbp->getAutoWB (ctemp, cgreen);

                if (ctemp != -1.0) {
                    // Set the automatics temperature value, or 0.0 if in ADD mode
                    temp->setValue (temp->getAddMode() ? 0.0 : ctemp);
                    // Set the automatics green value, or 0.0 if in ADD mode
                    green->setValue (green->getAddMode() ? 0.0 : cgreen);
                }

                //cache_customWB ((int)ctemp, cgreen); // this will be used to set initial Custom WB setting
            }
            break;
        /*
        All those types are the "default" case:
        case WBT_DAYLIGHT:
        case WBT_CLOUDY:
        case WBT_SHADE:
        case WBT_TUNGSTEN:
        case WBT_FLUORESCENT:
        case WBT_LAMP:
        case WBT_FLASH:
        case WBT_LED:
        */
        default:
            // Set the associated temperature, or 0.0 if in ADD mode
            temp->setValue(temp->getAddMode() ? 0.0 : (double)wbValues->temperature);
            // Set the stored temperature, or 0.0 if in ADD mode
            green->setValue(green->getAddMode() ? 0.0 : pp->wb.green);

            // The user may have changed the green value even for predefined WB values
            if (pedited) {
                green->setEditedState (pedited->wb.green ? Edited : UnEdited);
            }
            //cache_customGreen (pp->wb.green);
            break;
        }
    }
    methconn.block (false);
}

void WhiteBalance::write (ProcParams* pp, ParamsEdited* pedited) {
    
    Gtk::TreeModel::Row row = getActiveMethod();

    if (pedited) {
        pedited->wb.temperature = temp->getEditedState ();
        pedited->wb.green = green->getEditedState ();
        pedited->wb.method = row[methodColumns.colLabel]!=M("GENERAL_UNCHANGED");
    }

    WBEntry* ppMethod = findWBEntry (row[methodColumns.colLabel], WBLT_GUI);

    if (ppMethod)
        pp->wb.method = ppMethod->ppLabel;
    pp->wb.temperature = temp->getIntValue ();
    pp->wb.green = green->getValue ();
    
}

void WhiteBalance::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    if (wbp && defParams->wb.method == "Camera") {
        double ctemp; double cgreen;
        wbp->getCamWB (ctemp, cgreen);
        temp->setDefault (temp->getAddMode() ? 0 : (int)ctemp);
        green->setDefault (green->getAddMode() ? 0 : cgreen);
    }
    else if (wbp && defParams->wb.method == "Auto") {
        // this setDefaults method is called too early ; the wbp has been set,
        // but wbp is not ready to provide!
        double ctemp; double cgreen;
        wbp->getAutoWB (ctemp, cgreen);
        if (ctemp != -1.0) {
            temp->setDefault (temp->getAddMode() ? 0 : (int)ctemp);
            green->setDefault (green->getAddMode() ? 0 : cgreen);
        }
        else {
            // 6504 & 1.0 = same values as in ProcParams::setDefaults
            temp->setDefault (temp->getAddMode() ? 0 : 6504);
            green->setDefault (green->getAddMode() ? 0 : 1.0);
        }
    }
    else {
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
    Gtk::TreeModel::Row row = *(refTreeModel->append());
    row[methodColumns.colId] = WBParams::wbEntries.size();
    row[methodColumns.colLabel] = M("GENERAL_UNCHANGED");

}

int WhiteBalance::getSize () {

    return atoi(spotsize->get_active_text().c_str());
}

void WhiteBalance::setWB (int vtemp, double vgreen) {

    methconn.block(true);
    WBEntry *wbValues = findWBEntry("Custom", WBLT_PP);
    temp->setValue (vtemp);
    green->setValue (vgreen);
    opt = setActiveMethod(wbValues->GUILabel);
    cache_customWB (vtemp,vgreen); // sequence in which this call is made is important; must be before "method->set_active (2);"
    temp->setEditedState (Edited);
    green->setEditedState (Edited);
    methconn.block(false);

    if (listener) 
        listener->panelChanged (EvWBTemp, Glib::ustring::compose("%1, %2", (int)temp->getValue(), Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), green->getValue())));
}

void WhiteBalance::setAdjusterBehavior (bool tempadd, bool greenadd) {

	temp->setAddMode(tempadd);
	green->setAddMode(greenadd);
}

void WhiteBalance::trimValues (rtengine::procparams::ProcParams* pp) {

	temp->trimValue(pp->wb.temperature);
	green->trimValue(pp->wb.green);
}

inline void WhiteBalance::cache_customTemp(int temp) {
    custom_temp = temp;
}

void WhiteBalance::cache_customGreen(double green) {
    Gtk::TreeModel::Row row = getActiveMethod();
    if (row == refTreeModel->children().end()) return;

    custom_green[row[methodColumns.colId]] = green;
    //printf("WhiteBalance::cache_customWB(%d, %f): the \"green\" value of \"%s\" has been set to: %f\n", temp, green, row[methodColumns.colLabel], custom_green[row[methodColumns.colId]]);
}

void WhiteBalance::cache_customWB(int temp, double green) {
    cache_customTemp (temp);
    cache_customGreen (green);
}

int WhiteBalance::findWBEntryId (Glib::ustring label, enum WB_LabelType lblType) {
    for (unsigned int i=0; i<WBParams::wbEntries.size(); i++) {
        if (label == (lblType == WBLT_GUI ? WBParams::wbEntries[i]->GUILabel : WBParams::wbEntries[i]->ppLabel))
            return i;
    }
    return -1;
}

WBEntry* WhiteBalance::findWBEntry (Glib::ustring label, enum WB_LabelType lblType) {
    for (unsigned int i=0; i<WBParams::wbEntries.size(); i++) {
        if (label == (lblType == WBLT_GUI ? WBParams::wbEntries[i]->GUILabel : WBParams::wbEntries[i]->ppLabel))
            return WBParams::wbEntries[i];
    }
    return NULL;
}

int WhiteBalance::_setActiveMethod(Glib::ustring &label, Gtk::TreeModel::Children &children) {
    int found = -1;
    for(Gtk::TreeModel::Children::iterator iter = children.begin(); iter != children.end() && found==-1; ++iter) {
      Gtk::TreeModel::Row row = *iter;
      if (row[methodColumns.colLabel] == label) {
          method->set_active(iter);
          found = method->get_active_row_number();
      }
      if (found !=-1)
          return found;

      Gtk::TreeModel::Children childs = row.children();
      if (childs.size()) {
          found = _setActiveMethod(label, childs);
          if (found !=-1)
              return found;
      }
    }
    // Entry not found
    return -1;
}

int WhiteBalance::setActiveMethod(Glib::ustring label) {
    Gtk::TreeModel::Children children = refTreeModel->children();
    return _setActiveMethod(label, children);
}

inline Gtk::TreeRow WhiteBalance::getActiveMethod () {
    return *(method->get_active());
}
