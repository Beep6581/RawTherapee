/*
 *  Adapted from ART.
 *
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
#include "eventmapper.h"
#include "toneequalizer.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring ToneEqualizer::TOOL_NAME = "toneequalizer";

ToneEqualizer::ToneEqualizer(): FoldableToolPanel(this, TOOL_NAME, M("TP_TONE_EQUALIZER_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(AUTOEXP, "HISTORY_MSG_TONE_EQUALIZER_ENABLED");
    EvBands = m->newEvent(AUTOEXP, "HISTORY_MSG_TONE_EQUALIZER_BANDS");
    EvRegularization = m->newEvent(AUTOEXP, "HISTORY_MSG_TONE_EQUALIZER_REGULARIZATION");
    EvColormap = m->newEvent(AUTOEXP, "HISTORY_MSG_TONE_EQUALIZER_SHOW_COLOR_MAP");
    EvPivot = m->newEvent(AUTOEXP, "HISTORY_MSG_TONE_EQUALIZER_PIVOT");

    std::array<const char *, 5> images = {
        "purple",
        "blue",
        "gray",
        "yellow",
        "red"
    };
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i] = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_BAND_" + std::to_string(i)), -100, 100, 1, 0, Gtk::manage(new RTImage(Glib::ustring("circle-") + images[i] + "-small"))));
        bands[i]->setAdjusterListener(this);
        pack_start(*bands[i]);
        bands[i]->showIcons(false);
    }

    pivot = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_PIVOT"), -12, 12, 0.05, 0));
    pivot->setLogScale(64, 0, true);
    pivot->setAdjusterListener(this);
    pack_start(*pivot);

    pack_start(*Gtk::manage(new Gtk::HSeparator()));
    regularization = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_DETAIL"), -5, 5, 1, 0));
    regularization->setAdjusterListener(this);
    pack_start(*regularization);

    show_colormap = Gtk::manage(new CheckBox(M("TP_TONE_EQUALIZER_SHOW_COLOR_MAP"), multiImage));
    pack_start(*show_colormap);
    show_colormap->setCheckBoxListener(this);

    show_all_children ();
}


void ToneEqualizer::read(const ProcParams *pp, const ParamsEdited* pedited)
{
    disableListener();

    if (pedited) {
        set_inconsistent(multiImage && !pedited->toneEqualizer.enabled);
        for (size_t i = 0; i < bands.size(); ++i) {
            bands[i]->setEditedState(pedited->toneEqualizer.bands[i] ? Edited : UnEdited);
        }
        regularization->setEditedState(pedited->toneEqualizer.regularization ? Edited : UnEdited);
        pivot->setEditedState(pedited->toneEqualizer.pivot ? Edited : UnEdited);
        show_colormap->setEdited(pedited->toneEqualizer.show_colormap ? Edited : UnEdited);
    }

    setEnabled(pp->toneEqualizer.enabled);

    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->setValue(pp->toneEqualizer.bands[i]);
        bands[i]->showIcons(pp->toneEqualizer.show_colormap);
    }
    regularization->setValue(pp->toneEqualizer.regularization);

    pivot->setValue(pp->toneEqualizer.pivot);
    show_colormap->setValue(pp->toneEqualizer.show_colormap);

    enableListener();
}


void ToneEqualizer::write(ProcParams *pp, ParamsEdited* pedited)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        pp->toneEqualizer.bands[i] = bands[i]->getValue();
    }
    pp->toneEqualizer.enabled = getEnabled();
    pp->toneEqualizer.regularization = regularization->getValue();
    pp->toneEqualizer.show_colormap = show_colormap->getLastActive();
    pp->toneEqualizer.pivot = pivot->getValue();

    if (pedited) {
        auto &edited = pedited->toneEqualizer;
        edited.enabled = !get_inconsistent();
        for (size_t i = 0; i < bands.size(); ++i) {
            edited.bands[i] = bands[i]->getEditedState();
        }
        edited.regularization = regularization->getEditedState();
        edited.pivot = pivot->getEditedState();
        edited.show_colormap = show_colormap->getEdited();
    }
}


void ToneEqualizer::setDefaults(const ProcParams *defParams, const ParamsEdited* pedited)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->setDefault(defParams->toneEqualizer.bands[i]);
    }
    regularization->setDefault(defParams->toneEqualizer.regularization);

    pivot->setDefault(defParams->toneEqualizer.pivot);
    inital_params = defParams->toneEqualizer;

    if (pedited) {
        auto &edited = pedited->toneEqualizer;
        for (size_t i = 0; i < bands.size(); ++i) {
            bands[i]->setDefaultEditedState(edited.bands[i] ? Edited : UnEdited);
        }
        regularization->setDefaultEditedState(edited.regularization ? Edited : UnEdited);
        pivot->setDefaultEditedState(edited.pivot ? Edited : UnEdited);
    } else {
        for (auto band : bands) {
            band->setDefaultEditedState(Irrelevant);
        }
        regularization->setDefaultEditedState(Irrelevant);
        pivot->setDefaultEditedState(Irrelevant);
    }
}


void ToneEqualizer::adjusterChanged(Adjuster *a, double newval)
{
    if (listener && getEnabled()) {
        if (a == regularization) {
            listener->panelChanged(EvRegularization, Glib::ustring::format(a->getValue()));
        } else if (a == pivot) {
            listener->panelChanged(EvPivot, Glib::ustring::format(a->getValue()));
        } else {
            Glib::ustring s;
            for (size_t i = 0; i < bands.size(); ++i) {
                s += Glib::ustring::format((int)bands[i]->getValue()) + " ";
            }
            listener->panelChanged(EvBands, s);
        }
    }
}


void ToneEqualizer::adjusterAutoToggled(Adjuster *a)
{
}


void ToneEqualizer::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void ToneEqualizer::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    if (batchMode) {
        for (auto band : bands) {
            band->showEditedCB();
        }
        regularization->showEditedCB();
        pivot->showEditedCB();
    }
}


void ToneEqualizer::setAdjusterBehavior(bool bands_add, bool regularization_add, bool pivot_add)
{
    for (auto band : bands) {
        band->setAddMode(bands_add);
    }
    regularization->setAddMode(regularization_add);
    pivot->setAddMode(pivot_add);
}


void ToneEqualizer::checkBoxToggled(CheckBox *c, CheckValue newval)
{
    if (c == show_colormap) {
        colormapToggled();
    }
}


void ToneEqualizer::colormapToggled()
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->showIcons(show_colormap->getLastActive());
    }
    if (listener && getEnabled()) {
        listener->panelChanged(EvColormap, show_colormap->getLastActive() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void ToneEqualizer::trimValues(rtengine::procparams::ProcParams *pp)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->trimValue(pp->toneEqualizer.bands[i]);
    }
    regularization->trimValue(pp->toneEqualizer.regularization);
    pivot->trimValue(pp->toneEqualizer.pivot);
}

