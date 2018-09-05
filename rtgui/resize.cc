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
#include "guiutils.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

Resize::Resize () : FoldableToolPanel(this, "resize", M("TP_RESIZE_LABEL"), false, true), maxw(100000), maxh(100000)
{
    auto m = ProcEventMapper::getInstance();
    EvResizeAllowUpscaling = m->newEvent(RESIZE, "HISTORY_MSG_RESIZE_ALLOWUPSCALING");
    EvResizePPI = m->newEvent(RESIZE, "HISTORY_MSG_RESIZE_PIXELSPERINCH");

    cropw = 0;
    croph = 0;

    Gtk::Table* paramTable = Gtk::manage (new Gtk::Table (2, 2));

    appliesTo = Gtk::manage (new MyComboBoxText ());
    appliesTo->append (M("TP_RESIZE_CROPPEDAREA"));
    appliesTo->append (M("TP_RESIZE_FULLIMAGE"));
    appliesTo->set_active (0);

    Gtk::Label *appliesToLbl = Gtk::manage (new Gtk::Label (M("TP_RESIZE_APPLIESTO")));
    appliesToLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    paramTable->attach (*appliesToLbl, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    paramTable->attach (*appliesTo, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    // See Resize::methodChanged() when adding a new method.
    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_RESIZE_LANCZOS"));
    method->append (M("TP_RESIZE_NEAREST"));
    method->set_active (0);

    Gtk::Label *methodLbl = Gtk::manage (new Gtk::Label (M("TP_RESIZE_METHOD")));
    methodLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    paramTable->attach (*methodLbl, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    paramTable->attach (*method, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    spec = Gtk::manage (new MyComboBoxText ());
    spec->append (M("TP_RESIZE_SCALE"));
    spec->append (M("TP_RESIZE_WIDTH"));
    spec->append (M("TP_RESIZE_HEIGHT"));
    spec->append (M("TP_RESIZE_FITBOX"));
    spec->set_active (0);

    Gtk::Label *specifyLbl = Gtk::manage (new Gtk::Label (M("TP_RESIZE_SPECIFY")));
    specifyLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    paramTable->attach (*specifyLbl, 0, 1, 2, 3, Gtk::FILL, Gtk::SHRINK, 2, 2);
    paramTable->attach (*spec, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    allowUpscaling = Gtk::manage(new Gtk::CheckButton(M("TP_RESIZE_ALLOW_UPSCALING"))); // TODO Wrong capitalization used, should be sentence-case.
    allowUpscaling->signal_toggled().connect(sigc::mem_fun(*this, &Resize::allowUpscalingChanged));
    paramTable->attach (*allowUpscaling, 1, 2, 3, 4, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    pack_start (*paramTable, Gtk::PACK_SHRINK, 0);

    scale = new Adjuster (M("TP_RESIZE_SCALE"), 0.01, MAX_SCALE, 0.01, 1.);
    scale->setAdjusterListener (this);

    pack_start (*scale, Gtk::PACK_SHRINK, 4);

    sizeVB = Gtk::manage (new Gtk::VBox ());

    uom = Gtk::manage(new MyComboBoxText());
    uom->append(M("GENERAL_UOM_CENTIMETERS"));
    uom->append(M("GENERAL_UOM_INCHES"));
    uom->set_active(0);

    ppiSB = Gtk::manage(new MySpinButton());
    wPx = Gtk::manage(new MySpinButton());
    hPx = Gtk::manage(new MySpinButton());
    wPhys = Gtk::manage(new MySpinButton());
    hPhys = Gtk::manage(new MySpinButton());

    ppiSB->set_hexpand(true);
    wPx->set_hexpand(true);
    hPx->set_hexpand(true);
    wPhys->set_hexpand(true);
    hPhys->set_hexpand(true);

    Gtk::Label* wPxLbl = Gtk::manage(new Gtk::Label(M("GENERAL_WIDTH"))); // TODO delete all TP_RESIZE_W keys
    Gtk::Label* hPxLbl = Gtk::manage(new Gtk::Label(M("GENERAL_HEIGHT"))); // TODO delete all TP_RESIZE_H keys
    Gtk::Label* ppiLbl = Gtk::manage(new Gtk::Label(M("GENERAL_PIXELSPERINCH")));
    Gtk::Label* uomLbl = Gtk::manage(new Gtk::Label(M("GENERAL_UOM")));
    Gtk::Label* wUomLbl = Gtk::manage(new Gtk::Label(M("GENERAL_WIDTH")));
    Gtk::Label* hUomLbl = Gtk::manage(new Gtk::Label(M("GENERAL_HEIGHT")));

    ppiLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    uomLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    wUomLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    hUomLbl->set_alignment(Gtk::ALIGN_END, Gtk::ALIGN_CENTER);

    Gtk::Grid* dimensionGrid = Gtk::manage(new Gtk::Grid());
    dimensionGrid->get_style_context()->add_class("grid-spacing");

    dimensionGrid->attach(*wPxLbl, 0, 0, 1, 1);
    dimensionGrid->attach(*wPx, 1, 0, 1, 1);
    dimensionGrid->attach(*hPxLbl, 2, 0, 1, 1);
    dimensionGrid->attach(*hPx, 3, 0, 1, 1);

    MyExpander* physicalME = Gtk::manage(new MyExpander(false, M("Physical Dimensions")));

    Gtk::Grid* physicalGrid = Gtk::manage(new Gtk::Grid());
    physicalGrid->get_style_context()->add_class("grid-spacing");

    physicalGrid->attach(*ppiLbl, 0, 0, 2, 1);
    physicalGrid->attach(*ppiSB, 2, 0, 2, 1);
    physicalGrid->attach(*uomLbl, 0, 1, 2, 1);
    physicalGrid->attach(*uom, 2, 1, 2, 1);
    physicalGrid->attach(*wUomLbl, 0, 2, 1, 1);
    physicalGrid->attach(*wPhys, 1, 2, 1, 1);
    physicalGrid->attach(*hUomLbl, 2, 2, 1, 1);
    physicalGrid->attach(*hPhys, 3, 2, 1, 1);

    physicalME->add(*physicalGrid, false);
    physicalME->setLevel(2);

    sizeVB->pack_start(*dimensionGrid, Gtk::PACK_SHRINK, 0);
    sizeVB->pack_start(*physicalME, Gtk::PACK_SHRINK, 0);

    sizeVB->show_all();
    sizeVB->reference();
    physicalME->set_expanded(false);

    wPx->set_digits(0);
    wPx->set_increments(1, 100);
    wPx->set_range(32, MAX_SCALE * maxw);
    wPx->set_value(800);

    hPx->set_digits(0);
    hPx->set_increments(1, 100);
    hPx->set_range(32, MAX_SCALE * maxh);
    hPx->set_value(600);

    ppiSB->set_digits(0);
    ppiSB->set_increments(1, 100);
    ppiSB->set_range(1, 1200);
    ppiSB->set_value(300);

    wPhys->set_digits(1);
    wPhys->set_increments(0.1, 1);
    wPhys->set_range(0.01, 100);
    wPhys->set_value(10);

    hPhys->set_digits(1);
    hPhys->set_increments(0.1, 1);
    hPhys->set_range(0.01, 100);
    hPhys->set_value(15);

    wPxConn = wPx->signal_value_changed().connect(sigc::mem_fun(*this, &Resize::entryWPxChanged), true);
    hPxConn = hPx->signal_value_changed().connect(sigc::mem_fun(*this, &Resize::entryHPxChanged), true);
    wPhysConn = wPhys->signal_value_changed().connect(sigc::mem_fun(*this, &Resize::entryWPhysChanged), true);
    hPhysConn = hPhys->signal_value_changed().connect(sigc::mem_fun(*this, &Resize::entryHPhysChanged), true);
    ppiConn = ppiSB->signal_value_changed().connect(sigc::mem_fun(*this, &Resize::ppiChanged));
    uomConn = uom->signal_changed().connect(sigc::mem_fun(*this, &Resize::uomChanged));
    appliesToConn = appliesTo->signal_changed().connect(sigc::mem_fun(*this, &Resize::appliesToChanged));
    method->signal_changed().connect(sigc::mem_fun(*this, &Resize::methodChanged));
    scaleConn = spec->signal_changed().connect(sigc::mem_fun(*this, &Resize::specChanged));

    packBox = Gtk::manage (new ToolParamBlock ());
    pack_end (*packBox);
    packBox->hide();
    packBox->set_tooltip_markup (M("TP_PRSHARPENING_TOOLTIP"));

    show_all();
}

Resize::~Resize ()
{
    idle_register.destroy();
    delete scale;
    delete sizeVB;
}

void Resize::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
    appliesToConn.block (true);
    ppiConn.block(true);
    wPxConn.block (true);
    hPxConn.block (true);
    /*
    wPhysConn.block(true);
    hPhysConn.block(true);
    */
    scaleConn.block (true);
    scale->block(true);

    ppiSB->set_value(pp->resize.ppi);
    scale->setValue (pp->resize.scale);
    wPx->set_value (pp->resize.width);
    hPx->set_value (pp->resize.height);
    setEnabled (pp->resize.enabled);
    spec->set_active (pp->resize.dataspec);
    allowUpscaling->set_active(pp->resize.allowUpscaling);

    updatePhysDimensions();
    /*
    double wInches = wPx->get_value() / ppiSB->get_value_as_int();
    double hInches = hPx->get_value() / ppiSB->get_value_as_int();
    if (uom->get_active_row_number() == 0) {
        wPhys->set_value(wInches * 2.54);
        hPhys->set_value(hInches * 2.54);
    } else {
        wPhys->set_value(wInches);
        hPhys->set_value(hInches);
    }
    */

    updateGUI();

    appliesTo->set_active (0);

    if (pp->resize.appliesTo == "Cropped area") {
        appliesTo->set_active (0);
    } else if (pp->resize.appliesTo == "Full image") {
        appliesTo->set_active (1);
    }

    if (pp->resize.method == "Lanczos") {
        method->set_active (0);
    } else if (pp->resize.method == "Nearest") {
        method->set_active (1);
    } else {
        method->set_active (0);
    }

    wDirty = false;
    hDirty = false;

    if (pedited) {
        wDirty = pedited->resize.width;
        // TODO Do I need to add ppiSB here?
        hDirty = pedited->resize.height;
        scale->setEditedState (pedited->resize.scale ? Edited : UnEdited);

        if (!pedited->resize.appliesTo) {
            appliesTo->set_active (2);
        }

        if (!pedited->resize.method) {
            method->set_active (3);
        }

        if (!pedited->resize.dataspec) {
            spec->set_active (4);
        }

        allowUpscaling->set_inconsistent(!pedited->resize.allowUpscaling);
        set_inconsistent (multiImage && !pedited->resize.enabled);
    }

    scale->block(false);
    scaleConn.block (false);
    ppiConn.block(false);
    wPxConn.block (false);
    hPxConn.block (false);
    /*
    wPhysConn.block(false);
    hPhysConn.block(false);
    */
    appliesToConn.block (false);
    enableListener ();
}

void Resize::write (ProcParams* pp, ParamsEdited* pedited)
{
    int dataSpec = spec->get_active_row_number();

    pp->resize.scale  = scale->getValue();

    pp->resize.appliesTo = "Cropped area";

    if (appliesTo->get_active_row_number() == 0) {
        pp->resize.appliesTo = "Cropped area";
    } else if (appliesTo->get_active_row_number() == 1) {
        pp->resize.appliesTo = "Full image";
    }

    pp->resize.method = "Lanczos";

    if (method->get_active_row_number() == 0) {
        pp->resize.method = "Lanczos";
    } else if (method->get_active_row_number() == 1) {
        pp->resize.method = "Nearest";
    }

    pp->resize.ppi = ppiSB->get_value_as_int();
    pp->resize.dataspec = dataSpec;
    pp->resize.width = wPx->get_value_as_int ();
    pp->resize.height = hPx->get_value_as_int ();
    pp->resize.enabled = getEnabled ();

    pp->resize.allowUpscaling = allowUpscaling->get_active();

    if (pedited) {
        pedited->resize.enabled   = !get_inconsistent();
        pedited->resize.dataspec  = dataSpec != MAX_SCALE;
        pedited->resize.appliesTo = appliesTo->get_active_row_number() != 2;
        pedited->resize.method    = method->get_active_row_number() != 3;

        // TODO Do I need to add ppiSB here?

        if (pedited->resize.dataspec) {
            pedited->resize.scale     = scale->getEditedState ();
            pedited->resize.width     = wDirty;
            pedited->resize.height    = hDirty;
        } else {
            pedited->resize.scale     = false;
            pedited->resize.width     = false;
            pedited->resize.height    = false;
        }
        pedited->resize.allowUpscaling = !allowUpscaling->get_inconsistent();
    }
}

void Resize::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    scale->setDefault (defParams->resize.scale);

    if (pedited) {
        scale->setDefaultEditedState (pedited->resize.scale ? Edited : UnEdited);
    } else {
        scale->setDefaultEditedState (Irrelevant);
    }
}

void Resize::adjusterChanged(Adjuster* a, double newval)
{
    if (!batchMode) {
        wPxConn.block (true);
        hPxConn.block (true);
        wPx->set_value ((cropw && appliesTo->get_active_row_number() == 0 ? cropw : maxw) * a->getValue ());
        hPx->set_value ((croph && appliesTo->get_active_row_number() == 0 ? croph : maxh) * a->getValue ());
        wPxConn.block (false);
        hPxConn.block (false);
    }

    if (listener && (getEnabled () || batchMode)) {
        listener->panelChanged (EvResizeScale, Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(2), scale->getValue()));
    }
}

void Resize::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

int Resize::getComputedWidth()
{

    if (cropw && appliesTo->get_active_row_number() == 0)
        // use the crop dimensions
    {
        return (int)((double)(cropw) * (hPx->get_value() / (double)(croph)) + 0.5);
    } else
        // use the image dimensions
    {
        return (int)((double)(maxw) * (hPx->get_value() / (double)(maxh)) + 0.5);
    }
}

int Resize::getComputedHeight()
{

    if (croph && appliesTo->get_active_row_number() == 0)
        // use the crop dimensions
    {
        return (int)((double)(croph) * (wPx->get_value() / (double)(cropw)) + 0.5);
    } else
        // use the image dimensions
    {
        return (int)((double)(maxh) * (wPx->get_value() / (double)(maxw)) + 0.5);
    }
}

void Resize::appliesToChanged ()
{

    //printf("\nPASSAGE EN MODE \"%s\"\n\n", appliesTo->get_active_text().c_str());
    setDimensions();

    if (listener && (getEnabled () || batchMode)) {
        //printf("Appel du listener\n");
        listener->panelChanged (EvResizeAppliesTo, appliesTo->get_active_text());
    }
}

void Resize::methodChanged ()
{

    if (listener && (getEnabled () || batchMode)) {
        listener->panelChanged (EvResizeMethod, method->get_active_text());
    }

    // Post-resize Sharpening assumes the image is in Lab space, and currently Lanczos is the only method which uses that space, and Lanczos is on row 0.
    if (method->get_active_row_number() == 0) {
        packBox->set_sensitive(true);
    } else {
        packBox->set_sensitive(false);
    }
}

void Resize::update (bool isCropped, int cw, int ch, int ow, int oh)
{

    // updating crop values now
    if (isCropped) {
        cropw = cw;
        croph = ch;
    } else {
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

void Resize::sizeChanged(int mw, int mh, int ow, int oh)
{
    // updating max values now
    maxw = ow;
    maxh = oh;

    // updating the GUI synchronously
    setDimensions();
}

void Resize::setDimensions ()
{
    const auto func = [](gpointer data) -> gboolean {
        Resize* const self = static_cast<Resize*>(data);

        self->wPxConn.block(true);
        self->hPxConn.block(true);
        self->scale->block(true);

        int refw, refh;

        if (self->appliesTo->get_active_row_number() == 0 && self->cropw) {
            // Applies to Cropped area
            refw = self->cropw;
            refh = self->croph;
        } else {
            // Applies to Full image or crop is disabled
            refw = self->maxw;
            refh = self->maxh;
        }

        self->wPx->set_range(32, MAX_SCALE * refw);
        self->hPx->set_range(32, MAX_SCALE * refh);

        switch (self->spec->get_active_row_number()) {
            case 0: {
                // Scale mode
                self->wPx->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * self->scale->getValue() + 0.5)));
                self->hPx->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * self->scale->getValue() + 0.5)));
                break;
            }

            case 1: {
                // Width mode
                const double tmp_scale = self->wPx->get_value() / static_cast<double>(refw);
                self->scale->setValue(tmp_scale);
                self->hPx->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * tmp_scale + 0.5)));
                break;
            }

            case 2: {
                // Height mode
                const double tmp_scale = self->hPx->get_value() / static_cast<double>(refh);
                self->scale->setValue(tmp_scale);
                self->wPx->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * tmp_scale + 0.5)));
                break;
            }

            case 3: {
                // Bounding box mode
                const double tmp_scale =
                    self->wPx->get_value() / self->hPx->get_value() < static_cast<double>(refw) / static_cast<double>(refh)
                        ? self->wPx->get_value() / static_cast<double>(refw)
                        : self->hPx->get_value() / static_cast<double>(refh);

                self->scale->setValue(tmp_scale);
                break;
            }

            default: {
                break;
            }
        }

        self->scale->block(false);
        self->wPxConn.block(false);
        self->hPxConn.block(false);

        return FALSE;
    };

    idle_register.add(func, this);
}

void Resize::fitBoxScale()
{
    double tmpScale;
    double neww = wPx->get_value ();
    double newh = hPx->get_value ();

    if (cropw && appliesTo->get_active_row_number() == 0) {
        // use the crop dimensions
        if (((double)(cropw) / (double)(croph)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(cropw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(croph);
        }
    } else {
        // use the image dimensions
        if (((double)(maxw) / (double)(maxh)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(maxw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(maxh);
        }
    }

    scale->setValue (tmpScale);
}

void Resize::entryWPxChanged ()
{

    wDirty = true;

    if (!batchMode) {
        if (spec->get_active_row_number() == 3) {
            // Fit box mode
            fitBoxScale();
        } else {
            // Other modes
            hPxConn.block (true);
            scale->block (true);

            hPx->set_value ((double)(getComputedHeight()));
            scale->setValue (wPx->get_value () / (cropw && appliesTo->get_active_row_number() == 0 ? (double)cropw : (double)maxw));

            scale->block (false);
            hPxConn.block (false);
        }

        updatePhysDimensions();

    }

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            notifyBBox();
        } else {
            if (getEnabled () || batchMode) {
                listener->panelChanged (EvResizeWidth, Glib::ustring::format (wPx->get_value_as_int()));
            }
        }
    }
}

void Resize::entryHPxChanged ()
{

    hDirty = true;

    if (!batchMode && listener) {
        if (spec->get_active_row_number() == 3) {
            // Fit box mode
            fitBoxScale();
        } else {
            // Other modes
            wPxConn.block (true);
            scale->block (true);

            wPx->set_value ((double)(getComputedWidth()));
            scale->setValue (hPx->get_value () / (croph && appliesTo->get_active_row_number() == 0 ? (double)croph : (double)maxh));

            scale->block (false);
            wPxConn.block (false);
        }

        updatePhysDimensions();

    }

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            notifyBBox();
        } else {
            if (getEnabled () || batchMode) {
                listener->panelChanged (EvResizeHeight, Glib::ustring::format (hPx->get_value_as_int()));
            }
        }
    }
}

void Resize::entryWPhysChanged()
{
    wDirty = true;

    if (!batchMode) {
        if (spec->get_active_row_number() != 0) {
            // Anything other than scale mode
            if (uom->get_active_row_number() == 0) {
                wPx->set_value((wPhys->get_value() / 2.54) * ppiSB->get_value_as_int());
            } else {
                wPx->set_value(wPhys->get_value() * ppiSB->get_value_as_int());
            }
            fitBoxScale();
        }
    }

    if (listener) {
        if (spec->get_active_row_number() != 0) {
            notifyBBox();
        } else {
            if (getEnabled () || batchMode) {
                listener->panelChanged (EvResizeWidth, Glib::ustring::format (wPx->get_value_as_int()));
            }
        }
    }
}

void Resize::entryHPhysChanged()
{
    hDirty = true;

    if (!batchMode) {
        if (spec->get_active_row_number() != 0) {
            // Anything other than scale mode
            if (uom->get_active_row_number() == 0) {
                hPx->set_value((hPhys->get_value() / 2.54) * ppiSB->get_value_as_int());
            } else {
                hPx->set_value(hPhys->get_value() * ppiSB->get_value_as_int());
            }
            fitBoxScale();
        }
    }

    if (listener) {
        if (spec->get_active_row_number() != 0) {
            notifyBBox();
        } else {
            if (getEnabled () || batchMode) {
                listener->panelChanged (EvResizeHeight, Glib::ustring::format (hPx->get_value_as_int()));
            }
        }
    }
}

void Resize::ppiChanged()
{
    updatePxDimensions();
}

void Resize::uomChanged()
{
    updatePxDimensions();
}

void Resize::updatePxDimensions()
{
    wPxConn.block(true);
    hPxConn.block(true);

    if (uom->get_active_row_number() == 0) {
        wPx->set_value((wPhys->get_value() / 2.54) * ppiSB->get_value_as_int());
        hPx->set_value((hPhys->get_value() / 2.54) * ppiSB->get_value_as_int());
    } else {
        wPx->set_value(wPhys->get_value() * ppiSB->get_value_as_int());
        hPx->set_value(hPhys->get_value() * ppiSB->get_value_as_int());
    }

    wPxConn.block(false);
    hPxConn.block(false);
}

void Resize::updatePhysDimensions()
{
    wPhysConn.block(true);
    hPhysConn.block(true);

    double wInches = wPx->get_value() / ppiSB->get_value_as_int();
    double hInches = hPx->get_value() / ppiSB->get_value_as_int();
    if (uom->get_active_row_number() == 0) {
        wPhys->set_value(wInches * 2.54);
        hPhys->set_value(hInches * 2.54);
    } else {
        wPhys->set_value(wInches);
        hPhys->set_value(hInches);
    }

    wPhysConn.block(false);
    hPhysConn.block(false);
}

void Resize::specChanged ()
{

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        scale->sliderChanged ();
        break;

    case (1):
        // Width mode
        wPx->set_value((double)(getComputedWidth()));
        entryWPxChanged ();
        break;

    case (2):
        // Height mode
        hPx->set_value((double)(getComputedHeight()));
        entryHPxChanged ();
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

void Resize::updateGUI ()
{

    removeIfThere (this, scale, false);
    removeIfThere (this, sizeVB, false);

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        pack_start (*scale, Gtk::PACK_SHRINK, 4);
        break;

    case (1):
        // Width mode
        pack_start (*sizeVB, Gtk::PACK_SHRINK, 4);
        wPx->set_sensitive (true);
        hPx->set_sensitive (false);
        wPhys->set_sensitive (true);
        hPhys->set_sensitive (false);
        break;

    case (2):
        // Height mode
        pack_start (*sizeVB, Gtk::PACK_SHRINK, 4);
        wPx->set_sensitive (false);
        hPx->set_sensitive (true);
        wPhys->set_sensitive (false);
        hPhys->set_sensitive (true);
        break;

    case (3):
        // Bounding box mode
        pack_start (*sizeVB, Gtk::PACK_SHRINK, 4);
        wPx->set_sensitive (true);
        hPx->set_sensitive (true);
        wPhys->set_sensitive (true);
        hPhys->set_sensitive (true);
        break;

    default:
        break;
    }
}

void Resize::notifyBBox()
{
    if (listener && (getEnabled () || batchMode)) {
        listener->panelChanged (EvResizeBoundingBox, Glib::ustring::compose("(%1x%2)", (int)wPx->get_value(), (int)hPx->get_value() ));
    }
}

void Resize::setBatchMode (bool batchMode)
{

    method->append (M("GENERAL_UNCHANGED"));
    spec->append (M("GENERAL_UNCHANGED"));
    ToolPanel::setBatchMode (batchMode);
    scale->showEditedCB ();
}

void Resize::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Resize::allowUpscalingChanged()
{

    if (listener) {
        if (allowUpscaling->get_inconsistent()) {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_UNCHANGED"));
        } else if (allowUpscaling->get_active()) {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_DISABLED"));
        }
    }
}


void Resize::setAdjusterBehavior (bool scaleadd)
{

    scale->setAddMode(scaleadd);
}

void Resize::trimValues (rtengine::procparams::ProcParams* pp)
{

    scale->trimValue(pp->resize.scale);
}
