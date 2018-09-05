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

    cropw = 0;
    croph = 0;

    Gtk::Table* paramTable = Gtk::manage (new Gtk::Table (2, 2));

    appliesTo = Gtk::manage (new MyComboBoxText ());
    appliesTo->append (M("TP_RESIZE_CROPPEDAREA"));
    appliesTo->append (M("TP_RESIZE_FULLIMAGE"));
    appliesTo->set_active (0);

    Gtk::Label *label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_APPLIESTO")));
    label->set_alignment(0., 0.);
    paramTable->attach (*label, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    paramTable->attach (*appliesTo, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    // See Resize::methodChanged() when adding a new method.
    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_RESIZE_LANCZOS"));
    method->append (M("TP_RESIZE_NEAREST"));
    method->set_active (0);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_METHOD")));
    label->set_alignment(0., 0.);
    paramTable->attach (*label, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    paramTable->attach (*method, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    spec = Gtk::manage (new MyComboBoxText ());
    spec->append (M("TP_RESIZE_SCALE"));
    spec->append (M("TP_RESIZE_WIDTH"));
    spec->append (M("TP_RESIZE_HEIGHT"));
    spec->append (M("TP_RESIZE_FITBOX"));
    spec->set_active (0);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_SPECIFY")));
    label->set_alignment(0., 0.);
    paramTable->attach (*label, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    paramTable->attach (*spec, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    pack_start (*paramTable, Gtk::PACK_SHRINK, 4);

    scale = new Adjuster (M("TP_RESIZE_SCALE"), 0.01, MAX_SCALE, 0.01, 1.);
    scale->setAdjusterListener (this);

    pack_start (*scale, Gtk::PACK_SHRINK, 4);

    sizeVB = Gtk::manage (new Gtk::VBox ());

    Gtk::HBox* sHB = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* wPxHB = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* hPxHB = Gtk::manage (new Gtk::HBox ());
    wPx = Gtk::manage (new MySpinButton ());
    hPx = Gtk::manage (new MySpinButton ());
    wPxHB->set_spacing(3);
    wPxHB->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_W"))), Gtk::PACK_SHRINK, 0);
    wPxHB->pack_start (*wPx);
    hPxHB->set_spacing(3);
    hPxHB->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_H"))), Gtk::PACK_SHRINK, 0);
    hPxHB->pack_start (*hPx);
    sHB->set_spacing(4);
    sHB->pack_start (*wPxHB);
    sHB->pack_start (*hPxHB);

    sizeVB->pack_start (*sHB, Gtk::PACK_SHRINK, 0);

    allowUpscaling = Gtk::manage(new Gtk::CheckButton(M("TP_RESIZE_ALLOW_UPSCALING")));
    sizeVB->pack_start(*allowUpscaling);
    allowUpscaling->signal_toggled().connect(sigc::mem_fun(*this, &Resize::allowUpscalingChanged));

    sizeVB->show_all ();
    sizeVB->reference ();

    wPx->set_digits (0);
    wPx->set_increments (1, 100);
    wPx->set_value (800);
    wPx->set_range (32, MAX_SCALE * maxw);

    hPx->set_digits (0);
    hPx->set_increments (1, 100);
    hPx->set_value (600);
    hPx->set_range (32, MAX_SCALE * maxh);

    wPxConn = wPx->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryWChanged), true);
    hPxConn = hPx->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryHChanged), true);
    appliesToConn = appliesTo->signal_changed().connect ( sigc::mem_fun(*this, &Resize::appliesToChanged) );
    method->signal_changed().connect ( sigc::mem_fun(*this, &Resize::methodChanged) );
    scaleConn = spec->signal_changed().connect ( sigc::mem_fun(*this, &Resize::specChanged) );

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
    wPxConn.block (true);
    hPxConn.block (true);
    scaleConn.block (true);
    scale->block(true);

    scale->setValue (pp->resize.scale);
    wPx->set_value (pp->resize.width);
    hPx->set_value (pp->resize.height);
    setEnabled (pp->resize.enabled);
    spec->set_active (pp->resize.dataspec);
    allowUpscaling->set_active(pp->resize.allowUpscaling);
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
    wPxConn.block (false);
    hPxConn.block (false);
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

    pp->resize.dataspec = dataSpec;
    pp->resize.width = wPx->get_value_as_int ();
    pp->resize.height = hPx->get_value_as_int ();
    pp->resize.enabled = getEnabled ();
    //printf("  L:%d   hPx:%d\n", pp->resize.width, pp->resize.height);

    pp->resize.allowUpscaling = allowUpscaling->get_active();

    if (pedited) {
        pedited->resize.enabled   = !get_inconsistent();
        pedited->resize.dataspec  = dataSpec != MAX_SCALE;
        pedited->resize.appliesTo = appliesTo->get_active_row_number() != 2;
        pedited->resize.method    = method->get_active_row_number() != 3;

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
        // we use the crop dimensions
    {
        return (int)((double)(cropw) * (hPx->get_value() / (double)(croph)) + 0.5);
    } else
        // we use the image dimensions
    {
        return (int)((double)(maxw) * (hPx->get_value() / (double)(maxh)) + 0.5);
    }
}

int Resize::getComputedHeight()
{

    if (croph && appliesTo->get_active_row_number() == 0)
        // we use the crop dimensions
    {
        return (int)((double)(croph) * (wPx->get_value() / (double)(cropw)) + 0.5);
    } else
        // we use the image dimensions
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
        // we use the crop dimensions
        if (((double)(cropw) / (double)(croph)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(cropw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(croph);
        }
    } else {
        // we use the image dimensions
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

void Resize::entryWChanged ()
{

    wDirty = true;

    // updating width
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

void Resize::entryHChanged ()
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
        entryWChanged ();
        break;

    case (2):
        // Height mode
        hPx->set_value((double)(getComputedHeight()));
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
        break;

    case (2):
        // Height mode
        pack_start (*sizeVB, Gtk::PACK_SHRINK, 4);
        wPx->set_sensitive (false);
        hPx->set_sensitive (true);
        break;

    case (3):
        // Bounding box mode
        pack_start (*sizeVB, Gtk::PACK_SHRINK, 4);
        wPx->set_sensitive (true);
        hPx->set_sensitive (true);
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
