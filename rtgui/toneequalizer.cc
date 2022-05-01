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


ToneEqualizer::ToneEqualizer(): FoldableToolPanel(this, "toneequalizer", M("TP_TONE_EQUALIZER_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_ENABLED");
    EvBands = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_BANDS");
    EvRegularization = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_REGULARIZATION");
    EvColormap = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_SHOW_COLOR_MAP");
    EvPivot = m->newEvent(RGBCURVE, "HISTORY_MSG_TONE_EQUALIZER_PIVOT");

    std::array<const char *, 5> images = {
        "purple",
        "blue",
        "gray",
        "yellow",
        "red"
    };
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i] = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_BAND_" + std::to_string(i)), -100, 100, 1, 0, Gtk::manage(new RTImage(Glib::ustring("circle-") + images[i] + "-small.png"))));
        bands[i]->setAdjusterListener(this);
        pack_start(*bands[i]);
        bands[i]->showIcons(false);
    }

    pivot = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_PIVOT"), -12, 12, 0.05, 0));
    pivot->setLogScale(64, 0, true);
    pivot->setAdjusterListener(this);
    pack_start(*pivot);

    pack_start(*Gtk::manage(new Gtk::HSeparator()));
    regularization = Gtk::manage(new Adjuster(M("TP_TONE_EQUALIZER_DETAIL"), -5, 5, 1, 5));
    regularization->setAdjusterListener(this);
    pack_start(*regularization);

    show_colormap = Gtk::manage(new Gtk::CheckButton(M("TP_TONE_EQUALIZER_SHOW_COLOR_MAP")));
    pack_start(*show_colormap);
    show_colormap->signal_toggled().connect(sigc::mem_fun(this, &ToneEqualizer::colormapToggled));

    show_all_children ();
}


void ToneEqualizer::read(const ProcParams *pp, const ParamsEdited* pedited)
{
    disableListener();

    setEnabled(pp->toneEqualizer.enabled);

    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->setValue(pp->toneEqualizer.bands[i]);
        bands[i]->showIcons(pp->toneEqualizer.show_colormap);
    }
    regularization->setValue(pp->toneEqualizer.regularization);

    pivot->setValue(pp->toneEqualizer.pivot);
    show_colormap->set_active(pp->toneEqualizer.show_colormap);

    enableListener();
}


void ToneEqualizer::write(ProcParams *pp, ParamsEdited* pedited)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        pp->toneEqualizer.bands[i] = bands[i]->getValue();
    }
    pp->toneEqualizer.enabled = getEnabled();
    pp->toneEqualizer.regularization = regularization->getValue();
    pp->toneEqualizer.show_colormap = show_colormap->get_active();
    pp->toneEqualizer.pivot = pivot->getValue();
}


void ToneEqualizer::setDefaults(const ProcParams *defParams, const ParamsEdited* pedited)
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->setDefault(defParams->toneEqualizer.bands[i]);
    }
    regularization->setDefault(defParams->toneEqualizer.regularization);

    pivot->setDefault(defParams->toneEqualizer.pivot);
    inital_params = defParams->toneEqualizer;
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


void ToneEqualizer::colormapToggled()
{
    for (size_t i = 0; i < bands.size(); ++i) {
        bands[i]->showIcons(show_colormap->get_active());
    }
    if (listener && getEnabled()) {
        listener->panelChanged(EvColormap, show_colormap->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
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

