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
#include "whitebalance.h"

#include <iomanip>

#include "rtimage.h"
#include "options.h"

#define MINTEMP 1500   //1200
#define MAXTEMP 60000  //12000
#define CENTERTEMP 4750
#define MINGREEN 0.02
#define MAXGREEN 10.0
#define MINEQUAL 0.8
#define MAXEQUAL 1.5

using namespace rtengine;
using namespace rtengine::procparams;

Glib::RefPtr<Gdk::Pixbuf> WhiteBalance::wbPixbufs[toUnderlying(WBEntry::Type::CUSTOM) + 1];

void WhiteBalance::init ()
{
    wbPixbufs[toUnderlying(WBEntry::Type::CAMERA)]      = RTImage::createPixbufFromFile ("wb-camera-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::AUTO)]        = RTImage::createPixbufFromFile ("wb-auto-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::DAYLIGHT)]    = RTImage::createPixbufFromFile ("wb-sun-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::CLOUDY)]      = RTImage::createPixbufFromFile ("wb-cloudy-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::SHADE)]       = RTImage::createPixbufFromFile ("wb-shade-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::WATER)]       = RTImage::createPixbufFromFile ("wb-water-small.png");
  //wbPixbufs[toUnderlying(WBEntry::Type::WATER2)]      = RTImage::createPixbufFromFile ("wb-water-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::TUNGSTEN)]    = RTImage::createPixbufFromFile ("wb-tungsten-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::FLUORESCENT)] = RTImage::createPixbufFromFile ("wb-fluorescent-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::LAMP)]        = RTImage::createPixbufFromFile ("wb-lamp-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::FLASH)]       = RTImage::createPixbufFromFile ("wb-flash-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::LED)]         = RTImage::createPixbufFromFile ("wb-led-small.png");
    wbPixbufs[toUnderlying(WBEntry::Type::CUSTOM)]      = RTImage::createPixbufFromFile ("wb-custom-small.png");
}

void WhiteBalance::cleanup ()
{
    for (int i = 0; i < toUnderlying(WBEntry::Type::CUSTOM) + 1; i++) {
        wbPixbufs[i].reset();
    }
}

static double wbSlider2Temp(double sval)
{

    // slider range: 0 - 10000
    double temp;

    if (sval <= 5000) {
        // linear below center-temp
        temp = MINTEMP + (sval / 5000.0) * (CENTERTEMP - MINTEMP);
    } else {
        const double slope = (double)(CENTERTEMP - MINTEMP) / (MAXTEMP - CENTERTEMP);
        double x = (sval - 5000) / 5000; // x 0..1
        double y = x * slope + (1.0 - slope) * pow(x, 4.0);
        //double y = pow(x, 4.0);
        temp = CENTERTEMP + y * (MAXTEMP - CENTERTEMP);
    }

    if (temp < MINTEMP) {
        temp = MINTEMP;
    }

    if (temp > MAXTEMP) {
        temp = MAXTEMP;
    }

    return temp;
}

static double wbTemp2Slider(double temp)
{

    double sval;

    if (temp <= CENTERTEMP) {
        sval = ((temp - MINTEMP) / (CENTERTEMP - MINTEMP)) * 5000.0;
    } else {
        const double slope = (double)(CENTERTEMP - MINTEMP) / (MAXTEMP - CENTERTEMP);
        const double y = (temp - CENTERTEMP) / (MAXTEMP - CENTERTEMP);
        double x = pow(y, 0.25); // rough guess of x, will be a little lower
        double k = 0.1;
        bool add = true;

        // the y=f(x) function is a mess to invert, therefore we have this trial-refinement loop instead.
        // from tests, worst case is about 20 iterations, ie no problem
        for (;;) {
            double y1 = x * slope + (1.0 - slope) * pow(x, 4.0);

            if (5000 * fabs(y1 - y) < 0.1) {
                break;
            }

            if (y1 < y) {
                if (!add) {
                    k /= 2;
                }

                x += k;
                add = true;
            } else {
                if (add) {
                    k /= 2;
                }

                x -= k;
                add = false;
            }
        }

        sval = 5000.0 + x * 5000.0;
    }

    if (sval < 0) {
        sval = 0;
    }

    if (sval > 10000) {
        sval = 10000;
    }

    return sval;
}

WhiteBalance::WhiteBalance () : FoldableToolPanel(this, "whitebalance", M("TP_WBALANCE_LABEL"), true, true), wbp(nullptr), wblistener(nullptr)
{

    Gtk::Grid* methodgrid = Gtk::manage(new Gtk::Grid());
    methodgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(methodgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Label* lab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_METHOD") + ":"));
    setExpandAlignProperties(lab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    // Create the Tree model
    refTreeModel = Gtk::TreeStore::create(methodColumns);
    // Create the Combobox
    method = Gtk::manage (new MyComboBox ());
    setExpandAlignProperties(method, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    // Assign the model to the Combobox
    method->set_model(refTreeModel);

    WBEntry::Type oldType = WBParams::getWbEntries()[0].type;
    WBEntry::Type currType;
    Gtk::TreeModel::Row row, childrow;

    for (unsigned int i = 0; i < WBParams::getWbEntries().size(); i++) {
        if (oldType != (currType = WBParams::getWbEntries()[i].type)) {
            // New entry type
            if (currType == WBEntry::Type::FLUORESCENT) {
                // Creating the Fluorescent subcategory header
                row = *(refTreeModel->append());
                row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
                row[methodColumns.colLabel] = M("TP_WBALANCE_FLUO_HEADER");
                row[methodColumns.colId] = i + 100;
            }

            if (currType == WBEntry::Type::AUTO) {
                // Creating the auto category
                row = *(refTreeModel->append());
                row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
                row[methodColumns.colLabel] = M("TP_WBALANCE_AUTO_HEADER");
                row[methodColumns.colId] = i + 100;
            }

            if (currType == WBEntry::Type::WATER) {
                // Creating the under water subcategory header
                row = *(refTreeModel->append());
                row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
                row[methodColumns.colLabel] = M("TP_WBALANCE_WATER_HEADER");
                row[methodColumns.colId] = i + 100;
            }

            if (currType == WBEntry::Type::LAMP) {
                // Creating the Lamp subcategory header
                row = *(refTreeModel->append());
                row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
                row[methodColumns.colLabel] = M("TP_WBALANCE_LAMP_HEADER");
                row[methodColumns.colId] = i + 100;
            }

            if (currType == WBEntry::Type::LED) {
                // Creating the LED subcategory header
                row = *(refTreeModel->append());
                row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
                row[methodColumns.colLabel] = M("TP_WBALANCE_LED_HEADER");
                row[methodColumns.colId] = i + 100;
            }

            if (currType == WBEntry::Type::FLASH) {
                // Creating the Flash subcategory header
                row = *(refTreeModel->append());
                row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
                row[methodColumns.colLabel] = M("TP_WBALANCE_FLASH_HEADER");
                row[methodColumns.colId] = i + 100;
            }
        }

        if (currType == WBEntry::Type::FLUORESCENT
                || currType == WBEntry::Type::LAMP
                || currType == WBEntry::Type::WATER
                || currType == WBEntry::Type::FLASH
                || currType == WBEntry::Type::LED
                || currType == WBEntry::Type::AUTO
           ) {
            childrow = *(refTreeModel->append(row.children()));
            childrow[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
            childrow[methodColumns.colLabel] = WBParams::getWbEntries()[i].GUILabel;
            childrow[methodColumns.colId] = i;
        } else {
            row = *(refTreeModel->append());
            row[methodColumns.colIcon] = wbPixbufs[toUnderlying(currType)];
            row[methodColumns.colLabel] = WBParams::getWbEntries()[i].GUILabel;
            row[methodColumns.colId] = i;
        }

        oldType = currType;

        custom_green = 1.0;
        custom_equal = 1.0;
    }

    //Add the model columns to the Combo (which is a kind of view),
    //rendering them in the default way:
    method->pack_start(methodColumns.colIcon, false);
    method->pack_start(methodColumns.colLabel, true);

    std::vector<Gtk::CellRenderer*> cells = method->get_cells();
    Gtk::CellRendererText* cellRenderer = dynamic_cast<Gtk::CellRendererText*>(cells.at(1));
    cellRenderer->property_ellipsize() = Pango::ELLIPSIZE_MIDDLE;

    method->set_active (0); // Camera
    methodgrid->attach (*lab, 0, 0, 1, 1);
    methodgrid->attach (*method, 1, 0, 1, 1);
    pack_start (*methodgrid, Gtk::PACK_SHRINK, 0 );
    opt = 0;

    Gtk::Grid* spotgrid = Gtk::manage(new Gtk::Grid());
    spotgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(spotgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    spotbutton = Gtk::manage (new Gtk::Button (M("TP_WBALANCE_PICKER")));
    setExpandAlignProperties(spotbutton, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    spotbutton->get_style_context()->add_class("independent");
    spotbutton->set_tooltip_text(M("TP_WBALANCE_SPOTWB"));
    spotbutton->set_image (*Gtk::manage (new RTImage ("color-picker-small.png")));

    Gtk::Label* slab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_SIZE")));
    setExpandAlignProperties(slab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    Gtk::Grid* wbsizehelper = Gtk::manage(new Gtk::Grid());
    wbsizehelper->set_name("WB-Size-Helper");
    setExpandAlignProperties(wbsizehelper, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    spotsize = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(spotsize, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    spotsize->append ("2");

    if (options.whiteBalanceSpotSize == 2) {
        spotsize->set_active(0);
    }

    spotsize->append ("4");

    if (options.whiteBalanceSpotSize == 4) {
        spotsize->set_active(1);
    }

    spotsize->append ("8");

    if (options.whiteBalanceSpotSize == 8) {
        spotsize->set_active(2);
    }

    spotsize->append ("16");

    if (options.whiteBalanceSpotSize == 16) {
        spotsize->set_active(3);
    }

    spotsize->append ("32");

    if (options.whiteBalanceSpotSize == 32) {
        spotsize->set_active(4);
    }

    wbsizehelper->attach (*spotsize, 0, 0, 1, 1);

    spotgrid->attach (*spotbutton, 0, 0, 1, 1);
    spotgrid->attach (*slab, 1, 0, 1, 1);
    spotgrid->attach (*wbsizehelper, 2, 0, 1, 1);
    pack_start (*spotgrid, Gtk::PACK_SHRINK, 0 );

    Gtk::Separator *separator = Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    separator->get_style_context()->add_class("grid-row-separator");
    pack_start (*separator, Gtk::PACK_SHRINK, 0);

    Gtk::Image* itempL =  Gtk::manage (new RTImage ("circle-blue-small.png"));
    Gtk::Image* itempR =  Gtk::manage (new RTImage ("circle-yellow-small.png"));
    Gtk::Image* igreenL = Gtk::manage (new RTImage ("circle-magenta-small.png"));
    Gtk::Image* igreenR = Gtk::manage (new RTImage ("circle-green-small.png"));
    Gtk::Image* iblueredL = Gtk::manage (new RTImage ("circle-blue-small.png"));
    Gtk::Image* iblueredR = Gtk::manage (new RTImage ("circle-red-small.png"));
    Gtk::Image* itempbiasL =  Gtk::manage (new RTImage ("circle-blue-small.png"));
    Gtk::Image* itempbiasR =  Gtk::manage (new RTImage ("circle-yellow-small.png"));

    StudLabel = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    StudLabel->set_tooltip_text(M("TP_WBALANCE_STUDLABEL_TOOLTIP"));

    temp = Gtk::manage (new Adjuster (M("TP_WBALANCE_TEMPERATURE"), MINTEMP, MAXTEMP, 5, CENTERTEMP, itempL, itempR, &wbSlider2Temp, &wbTemp2Slider));
    green = Gtk::manage (new Adjuster (M("TP_WBALANCE_GREEN"), MINGREEN, MAXGREEN, 0.001, 1.0, igreenL, igreenR));
    equal = Gtk::manage (new Adjuster (M("TP_WBALANCE_EQBLUERED"), MINEQUAL, MAXEQUAL, 0.001, 1.0, iblueredL, iblueredR));
    tempBias = Gtk::manage (new Adjuster(M("TP_WBALANCE_TEMPBIAS"), -0.5, 0.5, 0.01, 0.0, itempbiasL, itempbiasR));
    cache_customTemp (0);
    cache_customGreen (0);
    cache_customEqual (0);
    equal->set_tooltip_markup (M("TP_WBALANCE_EQBLUERED_TOOLTIP"));
    tempBias->set_tooltip_markup (M("TP_WBALANCE_TEMPBIAS_TOOLTIP"));
    temp->show ();
    green->show ();
    equal->show ();
    tempBias->show ();

    /*  Gtk::Box* boxgreen = Gtk::manage (new Gtk::Box ());
    boxgreen->show ();

    boxgreen->pack_start(*igreenL);
    boxgreen->pack_start(*green);
    boxgreen->pack_start(*igreenR);*/
    pack_start(*StudLabel);

    pack_start (*temp);
    //pack_start (*boxgreen);
    pack_start (*green);
    pack_start (*equal);
    pack_start (*tempBias);

    temp->setAdjusterListener (this);
    green->setAdjusterListener (this);
    equal->setAdjusterListener (this);
    tempBias->setAdjusterListener (this);

    spotbutton->signal_pressed().connect( sigc::mem_fun(*this, &WhiteBalance::spotPressed) );
    methconn = method->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::optChanged) );
    spotsize->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::spotSizeChanged) );
}

WhiteBalance::~WhiteBalance()
{
    idle_register.destroy();
}

void WhiteBalance::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvWBEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvWBEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvWBEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void WhiteBalance::adjusterChanged(Adjuster* a, double newval)
{
    int tVal = (int)temp->getValue();
    double gVal = green->getValue();
    double eVal = equal->getValue();
    Gtk::TreeModel::Row row = getActiveMethod();

    if (row == refTreeModel->children().end()) {
        return;
    }

    const std::pair<bool, const WBEntry&> ppMethod = findWBEntry (row[methodColumns.colLabel], WBLT_GUI);
    const std::pair<bool, const WBEntry&> wbCustom = findWBEntry ("Custom", WBLT_PP);

    if (
        !ppMethod.first
        || (
            ppMethod.second.ppLabel != wbCustom.second.ppLabel
            && !(
                (
                    a == equal
                    || a == tempBias
                )
                && ppMethod.second.type == WBEntry::Type::AUTO
            )
        )
    ) {
        methconn.block(true);
        opt = setActiveMethod(wbCustom.second.GUILabel);
        tempBias->set_sensitive(false);

        cache_customWB (tVal, gVal);
        if (a != equal) {
            cache_customEqual(eVal);
        }
        methconn.block(false);
    }

    //cache custom WB setting to allow its recall
    if (a == temp) {
        cache_customTemp (tVal);
    } else if (a == green) {
        cache_customGreen (gVal);
    } else if (a == equal) {
        cache_customEqual (eVal);
    }

        // Recomputing AutoWB if it's the current method will happen in improccoordinator.cc

    if (listener && getEnabled()) {
        if (a == temp) {
            listener->panelChanged (EvWBTemp, Glib::ustring::format ((int)a->getValue()));
        } else if (a == green) {
            listener->panelChanged (EvWBGreen, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
        } else if (a == equal) {
            listener->panelChanged (EvWBequal, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
        } else if (a == tempBias) {
            listener->panelChanged (EvWBtempBias, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(2), a->getValue()));
        }
    }
}

void WhiteBalance::optChanged ()
{
    Gtk::TreeModel::Row row = getActiveMethod();

    if (row == refTreeModel->children().end()) {
        return;
    }

    if (row[methodColumns.colId] >= 100) {
        // "Header" solutions are trapped ; the combo is then set to the previous value
        bool prevState = methconn.block(true);
        method->set_active(opt);
        methconn.block(prevState);
        return;
    }
    StudLabel->hide();

    if (opt != row[methodColumns.colId]) {

        opt = row[methodColumns.colId];

        if (row[methodColumns.colLabel] == M("GENERAL_UNCHANGED")) {
            temp->setEditedState (UnEdited);
            green->setEditedState (UnEdited);
            equal->setEditedState (UnEdited);
            tempBias->setEditedState (UnEdited);
        } else {
            unsigned int methodId = findWBEntryId (row[methodColumns.colLabel], WBLT_GUI);
            const WBEntry& currMethod = WBParams::getWbEntries()[methodId];

            tempBias->set_sensitive(currMethod.type == WBEntry::Type::AUTO);
            bool autit = (currMethod.ppLabel == "autitcgreen");
            if (autit) {
                StudLabel->show();
            } else {
                StudLabel->hide();
            }

            switch (currMethod.type) {
            case WBEntry::Type::CAMERA:
                if (wbp) {
                    double ctemp, cgreen;
                    wbp->getCamWB (ctemp, cgreen);
                    temp->setValue (temp->getAddMode() ? 0.0 : (int)ctemp);
                    green->setValue (green->getAddMode() ? 0.0 : cgreen);
                    equal->setValue (equal->getAddMode() ? 0.0 : 1.0);

                    if (batchMode) {
                        temp->setEditedState (UnEdited);
                        green->setEditedState (UnEdited);
                        equal->setEditedState (UnEdited);
                    }
                }

                break;

            case WBEntry::Type::AUTO:
                if (wbp) {
                    if (batchMode) {
                        temp->setEditedState (UnEdited);
                        green->setEditedState (UnEdited);
                        // equal remain as is
                    }

                    // Recomputing AutoWB will happen in improccoordinator.cc
                }

                break;

            case WBEntry::Type::CUSTOM:
                if (custom_temp > 0) {
                    temp->setValue (temp->getAddMode() ? 0.0 : custom_temp);
                    green->setValue (green->getAddMode() ? 0.0 : custom_green);
                    equal->setValue (equal->getAddMode() ? 0.0 : custom_equal);
                } else {
                    cache_customTemp (temp->getValue());
                    cache_customGreen (green->getValue());
                    cache_customEqual (equal->getValue());
                }

                if (batchMode) {
                    temp->setEditedState (Edited);
                    green->setEditedState (Edited);
                    equal->setEditedState (Edited);
                }

                break;

            /* All other solution are the default cases
            case WBEntry::Type::DAYLIGHT:
            case WBEntry::Type::CLOUDY:
            case WBEntry::Type::SHADE:
            case WBEntry::Type::TUNGSTEN:
            case WBEntry::Type::FLUORESCENT:
            case WBEntry::Type::LAMP:
            case WBEntry::Type::FLASH:
            case WBEntry::Type::LED:*/
            default:
                // recall custom WB settings if it exists, set to 1.0 otherwise
                temp->setValue  ( temp->getAddMode() ? 0.0 : (double)(currMethod.temperature));
                green->setValue (green->getAddMode() ? 0.0 : (double)(currMethod.green));
                equal->setValue (equal->getAddMode() ? 0.0 : (double)(currMethod.equal));

                if (batchMode) {
                    temp->setEditedState (Edited);
                    green->setEditedState (Edited);
                    equal->setEditedState (Edited);
                }

                break;
            }
        }

        if (listener && getEnabled()) {
            listener->panelChanged (EvWBMethod, row[methodColumns.colLabel]);
        }
    }
}

void WhiteBalance::spotPressed ()
{
    StudLabel->hide();

    if (wblistener) {
        wblistener->spotWBRequested (getSize());
    }
}

void WhiteBalance::spotSizeChanged ()
{
    options.whiteBalanceSpotSize = getSize();

    if (wblistener) {
        wblistener->spotWBRequested (getSize());
    }
}

void WhiteBalance::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    methconn.block (true);
    equal->setValue (pp->wb.equal);
    tempBias->setValue (pp->wb.tempBias);
    tempBias->set_sensitive(true);

    if (pedited) {
        // By default, temperature and green are said "UnEdited", but it may change later
        temp->setEditedState (UnEdited);
        green->setEditedState (UnEdited);
        equal->setEditedState (pedited->wb.equal ? Edited : UnEdited);
        tempBias->setEditedState (pedited->wb.tempBias ? Edited : UnEdited);
    }

    if (pedited && !pedited->wb.method) {
        opt = setActiveMethod(M("GENERAL_UNCHANGED"));
    } else {
        const WBEntry& wbValues =
            [this, pp]() -> const WBEntry&
            {
                const std::pair<bool, const WBEntry&> res = findWBEntry(pp->wb.method, WBLT_PP);
                return
                    !res.first
                        ? findWBEntry("Camera", WBLT_PP).second
                        : res.second;
            }();

        opt = setActiveMethod(wbValues.GUILabel);

        // temperature is reset to the associated temperature, or 0.0 if addMode is set.
        switch (wbValues.type) {
        case WBEntry::Type::CUSTOM:
            temp->setValue (temp->getAddMode() ? 0.0 : pp->wb.temperature);
            green->setValue (green->getAddMode() ? 0.0 : pp->wb.green);
            equal->setValue (equal->getAddMode() ? 0.0 : pp->wb.equal);
            tempBias->setValue (tempBias->getAddMode() ? 0.0 : pp->wb.tempBias);
            cache_customTemp (pp->wb.temperature);
            cache_customGreen (pp->wb.green);
            cache_customEqual (pp->wb.equal);

            if (pedited) {
                // The user may have changed the temperature and green value
                temp->setEditedState (pedited->wb.temperature ? Edited : UnEdited);
                green->setEditedState (pedited->wb.green ? Edited : UnEdited);
            }

            break;

        case WBEntry::Type::CAMERA:
            if (wbp) {
                double ctemp = -1.0;
                double cgreen = -1.0;
                wbp->getCamWB (ctemp, cgreen);

                if (ctemp != -1.0) {
                    // Set the camera's temperature value, or 0.0 if in ADD mode
                    temp->setValue (temp->getAddMode() ? 0.0 : ctemp);
                    // Set the camera's green value, or 0.0 if in ADD mode
                    green->setValue (green->getAddMode() ? 0.0 : cgreen);
                    equal->setValue (equal->getAddMode() ? 0.0 : 1.);
                } else {
                    temp->setValue (temp->getAddMode() ? 0.0 : pp->wb.temperature);
                    green->setValue (green->getAddMode() ? 0.0 : pp->wb.green);
                    equal->setValue (equal->getAddMode() ? 0.0 : pp->wb.equal);
                }
                tempBias->setValue (equal->getAddMode() ? 0.0 : pp->wb.tempBias);
            }

            break;

        case WBEntry::Type::AUTO:
            // the equalizer's value is restored for the AutoWB
            equal->setValue (equal->getAddMode() ? 0.0 : pp->wb.equal);
            tempBias->setValue (tempBias->getAddMode() ? 0.0 : pp->wb.tempBias);

            // set default values first if in ADD mode, otherwise keep the current ones
            if (temp->getAddMode() ) {
                temp->setValue (0.0);
            }

            if (green->getAddMode()) {
                green->setValue (0.0);
            }

            // Recomputing AutoWB will happen in improccoordinator.cc

            break;

        /*
        All those types are the "default" case:
        case WBEntry::Type::DAYLIGHT:
        case WBEntry::Type::CLOUDY:
        case WBEntry::Type::SHADE:
        case WBEntry::Type::TUNGSTEN:
        case WBEntry::Type::FLUORESCENT:
        case WBEntry::Type::LAMP:
        case WBEntry::Type::FLASH:
        case WBEntry::Type::LED:
        */
        default:
            // Set the associated temperature, or 0.0 if in ADD mode
            temp->setValue(temp->getAddMode() ? 0.0 : (double)wbValues.temperature);
            // Set the stored temperature, or 0.0 if in ADD mode
            green->setValue(green->getAddMode() ? 0.0 : pp->wb.green);
            equal->setValue(equal->getAddMode() ? 0.0 : pp->wb.equal);
            tempBias->setValue(equal->getAddMode() ? 0.0 : pp->wb.tempBias);

            // The user may have changed the green value even for predefined WB values
            if (pedited) {
                green->setEditedState (pedited->wb.green ? Edited : UnEdited);
                equal->setEditedState (pedited->wb.equal ? Edited : UnEdited);
                tempBias->setEditedState (pedited->wb.tempBias ? Edited : UnEdited);
            }

            //cache_customGreen (pp->wb.green);
            break;
        }

        tempBias->set_sensitive(wbValues.type == WBEntry::Type::AUTO);
        bool autit = (wbValues.ppLabel == "autitcgreen");
        if (autit) {
            StudLabel->show();
        } else {
            StudLabel->hide();
        }
        
    }

    setEnabled(pp->wb.enabled);
    if (pedited) {
        set_inconsistent(multiImage && !pedited->wb.enabled);
    }


    methconn.block (false);
    enableListener ();
}

void WhiteBalance::write (ProcParams* pp, ParamsEdited* pedited)
{

    Gtk::TreeModel::Row row = getActiveMethod();

    if (pedited) {
        pedited->wb.temperature = temp->getEditedState ();
        pedited->wb.green = green->getEditedState ();
        pedited->wb.equal = equal->getEditedState ();
        pedited->wb.tempBias = tempBias->getEditedState ();
        pedited->wb.method = row[methodColumns.colLabel] != M("GENERAL_UNCHANGED");
        pedited->wb.enabled = !get_inconsistent();
    }

    pp->wb.enabled = getEnabled();

    const std::pair<bool, const WBEntry&> ppMethod = findWBEntry (row[methodColumns.colLabel], WBLT_GUI);

    if (ppMethod.first) {
        pp->wb.method = ppMethod.second.ppLabel;
    }

    pp->wb.temperature = temp->getIntValue ();
    pp->wb.green = green->getValue ();
    pp->wb.equal = equal->getValue ();
    pp->wb.tempBias = tempBias->getValue ();
}

void WhiteBalance::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    equal->setDefault (defParams->wb.equal);
    tempBias->setDefault (defParams->wb.tempBias);

    if (wbp && defParams->wb.method == "Camera") {
        double ctemp;
        double cgreen;
        wbp->getCamWB (ctemp, cgreen);

        // FIXME: Seems to be always -1.0, called too early? Broken!
        if (ctemp != -1.0) {
            temp->setDefault (temp->getAddMode() ? 0 : (int)ctemp);
            green->setDefault (green->getAddMode() ? 0 : cgreen);
        }
    } else {
        temp->setDefault (defParams->wb.temperature);
        green->setDefault (defParams->wb.green);
    }
    // Recomputing AutoWB if it's the current method will happen in improccoordinator.cc

    if (pedited) {
        temp->setDefaultEditedState (pedited->wb.temperature ? Edited : UnEdited);
        green->setDefaultEditedState (pedited->wb.green ? Edited : UnEdited);
        equal->setDefaultEditedState (pedited->wb.equal ? Edited : UnEdited);
        tempBias->setDefaultEditedState (pedited->wb.tempBias ? Edited : UnEdited);
    } else {
        temp->setDefaultEditedState (Irrelevant);
        green->setDefaultEditedState (Irrelevant);
        equal->setDefaultEditedState (Irrelevant);
        tempBias->setDefaultEditedState (Irrelevant);
    }
}

void WhiteBalance::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    temp->showEditedCB ();
    green->showEditedCB ();
    equal->showEditedCB ();
    tempBias->showEditedCB ();
    Gtk::TreeModel::Row row = *(refTreeModel->append());
    row[methodColumns.colId] = WBParams::getWbEntries().size();
    row[methodColumns.colLabel] = M("GENERAL_UNCHANGED");

}

int WhiteBalance::getSize ()
{

    return atoi(spotsize->get_active_text().c_str());
}

void WhiteBalance::setWB (int vtemp, double vgreen)
{

    methconn.block(true);
    const std::pair<bool, const WBEntry&> wbValues = findWBEntry("Custom", WBLT_PP);
    setEnabled(true);
    temp->setValue (vtemp);
    green->setValue (vgreen);
    opt = setActiveMethod(wbValues.second.GUILabel);
    cache_customWB (vtemp, vgreen); // sequence in which this call is made is important; must be before "method->set_active (2);"
    cache_customEqual(equal->getValue());
    temp->setEditedState (Edited);
    green->setEditedState (Edited);
    methconn.block(false);

    if (listener) {
        listener->panelChanged (EvWBTemp, Glib::ustring::compose("%1, %2", (int)temp->getValue(), Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), green->getValue())));
    }

}

void WhiteBalance::setAdjusterBehavior (bool tempadd, bool greenadd, bool equaladd, bool tempbiasadd)
{

    temp->setAddMode(tempadd);
    green->setAddMode(greenadd);
    equal->setAddMode(equaladd);
    tempBias->setAddMode(tempbiasadd);
}

void WhiteBalance::trimValues (rtengine::procparams::ProcParams* pp)
{

    temp->trimValue(pp->wb.temperature);
    green->trimValue(pp->wb.green);
    equal->trimValue(pp->wb.equal);
    tempBias->trimValue(pp->wb.tempBias);
}

inline void WhiteBalance::cache_customTemp(int temp)
{
    custom_temp = temp;
}

void WhiteBalance::cache_customGreen(double green)
{
    custom_green = green;
}
void WhiteBalance::cache_customEqual(double equal)
{
    custom_equal = equal;
}

void WhiteBalance::cache_customWB(int temp, double green)
{
    cache_customTemp (temp);
    cache_customGreen (green);
}

unsigned int WhiteBalance::findWBEntryId (const Glib::ustring& label, enum WB_LabelType lblType)
{
    for (unsigned int i = 0; i < WBParams::getWbEntries().size(); i++) {
        if (label == (lblType == WBLT_GUI ? WBParams::getWbEntries()[i].GUILabel : WBParams::getWbEntries()[i].ppLabel)) {
            return i;
        }
    }

    return 0; // default to camera wb
}

std::pair<bool, const WBEntry&> WhiteBalance::findWBEntry(const Glib::ustring& label, enum WB_LabelType lblType)
{
    for (unsigned int i = 0; i < WBParams::getWbEntries().size(); ++i) {
        if (label == (lblType == WBLT_GUI ? WBParams::getWbEntries()[i].GUILabel : WBParams::getWbEntries()[i].ppLabel)) {
            return {true, WBParams::getWbEntries()[i]};
        }
    }

    return {false, WBParams::getWbEntries()[0]};
}

int WhiteBalance::_setActiveMethod(Glib::ustring &label, Gtk::TreeModel::Children &children)
{
    int found = -1;

    for(Gtk::TreeModel::Children::iterator iter = children.begin(); iter != children.end() && found == -1; ++iter) {
        Gtk::TreeModel::Row row = *iter;

        if (row[methodColumns.colLabel] == label) {
            method->set_active(iter);
            found = row[methodColumns.colId];
        }

        if (found != -1) {
            return found;
        }

        Gtk::TreeModel::Children childs = row.children();

        if (childs.size()) {
            found = _setActiveMethod(label, childs);

            if (found != -1) {
                return found;
            }
        }
    }

    // Entry not found
    return -1;
}

int WhiteBalance::setActiveMethod(Glib::ustring label)
{
    Gtk::TreeModel::Children children = refTreeModel->children();
    return _setActiveMethod(label, children);
}

inline Gtk::TreeRow WhiteBalance::getActiveMethod ()
{
    return *(method->get_active());
}

void WhiteBalance::WBChanged(double temperature, double greenVal, float studgood)
{
    idle_register.add(
        [this, temperature, greenVal, studgood]() -> bool
        {
            disableListener();
            setEnabled(true);
            temp->setValue(temperature);
            green->setValue(greenVal);
            StudLabel->set_text(
                Glib::ustring::compose(M("TP_WBALANCE_STUDLABEL"),
                                   Glib::ustring::format(std::fixed, std::setprecision(4), studgood))
            );            
            temp->setDefault(temperature);
            green->setDefault(greenVal);
            enableListener();

            return false;
        }
    );
}
