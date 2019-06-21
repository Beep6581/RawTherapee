/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Romei <aldrop8@gmail.com>
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
#include <iomanip>

#include "filmnegative.h"

#include "editwidgets.h"
#include "eventmapper.h"
#include "options.h"
#include "rtimage.h"

#include "../rtengine/procparams.h"

namespace
{

Adjuster* createExponentAdjuster(AdjusterListener* listener, const Glib::ustring& label, double defaultVal)
{
    Adjuster* const adj = Gtk::manage(new Adjuster(label, 0.3, 6, 0.001, defaultVal)); // exponent
    adj->setAdjusterListener(listener);

    if (adj->delay < options.adjusterMaxDelay) {
        adj->delay = options.adjusterMaxDelay;
    }

    adj->show();
    return adj;
}

}

FilmNegative::FilmNegative() :
    FoldableToolPanel(this, "filmnegative", M("TP_FILMNEGATIVE_LABEL"), false, true),
    EditSubscriber(ET_OBJECTS),
    evFilmNegativeExponents(ProcEventMapper::getInstance()->newEvent(FIRST, "HISTORY_MSG_FILMNEGATIVE_EXPONENTS")),
    evFilmNegativeEnabled(ProcEventMapper::getInstance()->newEvent(FIRST, "HISTORY_MSG_FILMNEGATIVE_ENABLED")),
    fnp(nullptr),
    redExp(createExponentAdjuster(this, M("TP_FILMNEGATIVE_RED"), 2.72)),
    greenExp(createExponentAdjuster(this, M("TP_FILMNEGATIVE_GREEN"), 2.0)),
    blueExp(createExponentAdjuster(this, M("TP_FILMNEGATIVE_BLUE"), 1.72)),
    spotgrid(Gtk::manage(new Gtk::Grid())),
    spotbutton(Gtk::manage(new Gtk::ToggleButton(M("TP_FILMNEGATIVE_PICK")))),
    redRatio(redExp->getValue() / greenExp->getValue()),
    blueRatio(blueExp->getValue() / greenExp->getValue())
{
    spotgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(spotgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    setExpandAlignProperties(spotbutton, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    spotbutton->get_style_context()->add_class("independent");
    spotbutton->set_tooltip_text(M("TP_FILMNEGATIVE_GUESS_TOOLTIP"));
    spotbutton->set_image (*Gtk::manage (new RTImage ("color-picker-small.png")));

    // TODO make spot size configurable ?

    // Gtk::Label* slab = Gtk::manage (new Gtk::Label (M("TP_WBALANCE_SIZE")));
    // setExpandAlignProperties(slab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    // Gtk::Grid* wbsizehelper = Gtk::manage(new Gtk::Grid());
    // wbsizehelper->set_name("WB-Size-Helper");
    // setExpandAlignProperties(wbsizehelper, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    // spotsize = Gtk::manage (new MyComboBoxText ());
    // setExpandAlignProperties(spotsize, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    // spotsize->append ("2");
    // spotsize->set_active(0);
    // spotsize->append ("4");

    spotgrid->attach (*spotbutton, 0, 1, 1, 1);
    // spotgrid->attach (*slab, 1, 0, 1, 1);
    // spotgrid->attach (*wbsizehelper, 2, 0, 1, 1);

    pack_start(*redExp, Gtk::PACK_SHRINK, 0);
    pack_start(*greenExp, Gtk::PACK_SHRINK, 0);
    pack_start(*blueExp, Gtk::PACK_SHRINK, 0);
    pack_start(*spotgrid, Gtk::PACK_SHRINK, 0);

    spotbutton->signal_toggled().connect(sigc::mem_fun(*this, &FilmNegative::editToggled));
    // spotsize->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::spotSizeChanged) );

    // Editing geometry; create the spot rectangle
    Rectangle* const spotRect = new Rectangle();
    spotRect->filled = false;
    
    visibleGeometry.push_back(spotRect);

    // Stick a dummy rectangle over the whole image in mouseOverGeometry.
    // This is to make sure the getCursor() call is fired everywhere.
    Rectangle* const imgRect = new Rectangle();
    imgRect->filled = true;

    mouseOverGeometry.push_back(imgRect);
}

FilmNegative::~FilmNegative()
{
    for (auto geometry : visibleGeometry) {
        delete geometry;
    }

    for (auto geometry : mouseOverGeometry) {
        delete geometry;
    }
}

void FilmNegative::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();

    if (pedited) {
        redExp->setEditedState(pedited->filmNegative.redExp ? Edited : UnEdited);
        greenExp->setEditedState(pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueExp->setEditedState(pedited->filmNegative.blueExp ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->filmNegative.enabled);
    }

    setEnabled(pp->filmNegative.enabled);
    redExp->setValue(pp->filmNegative.redExp);
    greenExp->setValue(pp->filmNegative.greenExp);
    blueExp->setValue(pp->filmNegative.blueExp);

    enableListener();
}

void FilmNegative::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->filmNegative.redExp = redExp->getValue();
    pp->filmNegative.greenExp = greenExp->getValue();
    pp->filmNegative.blueExp = blueExp->getValue();
    pp->filmNegative.enabled = getEnabled();

    if (pedited) {
        pedited->filmNegative.redExp = redExp->getEditedState();
        pedited->filmNegative.greenExp = greenExp->getEditedState();
        pedited->filmNegative.blueExp = blueExp->getEditedState();
        pedited->filmNegative.enabled = !get_inconsistent();
    }
}

void FilmNegative::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    redExp->setValue(defParams->filmNegative.redExp);
    greenExp->setValue(defParams->filmNegative.greenExp);
    blueExp->setValue(defParams->filmNegative.blueExp);

    if (pedited) {
        redExp->setDefaultEditedState(pedited->filmNegative.redExp ? Edited : UnEdited);
        greenExp->setDefaultEditedState(pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueExp->setDefaultEditedState(pedited->filmNegative.blueExp ? Edited : UnEdited);
    } else {
        redExp->setDefaultEditedState(Irrelevant);
        greenExp->setDefaultEditedState(Irrelevant);
        blueExp->setDefaultEditedState(Irrelevant);
    }
}

void FilmNegative::setBatchMode(bool batchMode)
{
    if (batchMode) {
        spotConn.disconnect();
        removeIfThere(this, spotgrid, false);
        ToolPanel::setBatchMode(batchMode);
        redExp->showEditedCB();
        greenExp->showEditedCB();
        blueExp->showEditedCB();
    }
}

void FilmNegative::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        if (a == redExp || a == greenExp || a == blueExp) {
            disableListener();
            if (a == greenExp) {
                redExp->setValue(a->getValue() * redRatio);
                blueExp->setValue(a->getValue() * blueRatio);
            }
            else if (a == redExp) {
                redRatio = newval / greenExp->getValue();
            }
            else if (a == blueExp) {
                blueRatio = newval / greenExp->getValue();
            }
            enableListener();

            if (getEnabled()) {
                listener->panelChanged(
                    evFilmNegativeExponents,
                    Glib::ustring::compose(
                        "R=%1 ; G=%2 ; B=%3",
                        redExp->getTextValue(),
                        greenExp->getTextValue(),
                        blueExp->getTextValue()
                    )
                );
            }
        }
    }
}

void FilmNegative::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void FilmNegative::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(evFilmNegativeEnabled, M("GENERAL_UNCHANGED"));
        }
        else if (getEnabled()) {
            listener->panelChanged(evFilmNegativeEnabled, M("GENERAL_ENABLED"));
        }
        else {
            listener->panelChanged(evFilmNegativeEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void FilmNegative::setFilmNegProvider(FilmNegProvider* provider)
{
    fnp = provider;
}

void FilmNegative::setEditProvider(EditDataProvider* provider)
{
    EditSubscriber::setEditProvider(provider);
}

CursorShape FilmNegative::getCursor(int objectID) const
{
   return CSSpotWB;
}

bool FilmNegative::mouseOver(int modifierKey)
{
    EditDataProvider* const provider = getEditProvider();
    Rectangle* const spotRect = static_cast<Rectangle*>(visibleGeometry.at(0));
    spotRect->setXYWH(provider->posImage.x - 16, provider->posImage.y - 16, 32, 32);

    return true;
}

bool FilmNegative::button1Pressed(int modifierKey)
{
    EditDataProvider* const provider = getEditProvider();

    EditSubscriber::action = EditSubscriber::Action::NONE;

    if (listener) {
        refSpotCoords.push_back(provider->posImage);

        if (refSpotCoords.size() == 2) {
            // User has selected 2 reference gray spots. Calculating new exponents
            // from channel values and updating parameters.

            std::array<float, 3> newExps;
            if (fnp->getFilmNegativeExponents(refSpotCoords[0], refSpotCoords[1], newExps)) {
                disableListener();
                redExp->setValue(newExps[0]);
                greenExp->setValue(newExps[1]);
                blueExp->setValue(newExps[2]);
                redRatio = redExp->getValue() / greenExp->getValue();
                blueRatio = blueExp->getValue() / greenExp->getValue();
                enableListener();

                if (listener && getEnabled()) {
                    listener->panelChanged(
                        evFilmNegativeExponents,
                        Glib::ustring::compose(
                            "R=%1 ; G=%2 ; B=%3",
                            redExp->getTextValue(),
                            greenExp->getTextValue(),
                            blueExp->getTextValue()
                        )
                    );
                }
            }

            switchOffEditMode();
        }
    }

    return true;
}

bool FilmNegative::button1Released()
{
    EditSubscriber::action = EditSubscriber::Action::NONE;
    return true;
}

void FilmNegative::switchOffEditMode()
{
    refSpotCoords.clear();
    unsubscribe();
    spotbutton->set_active(false);
}

void FilmNegative::editToggled()
{
    if (spotbutton->get_active()) {
        subscribe();

        int w, h;
        getEditProvider()->getImageSize(w, h);

        // Stick a dummy rectangle over the whole image in mouseOverGeometry.
        // This is to make sure the getCursor() call is fired everywhere.
        Rectangle* const imgRect = static_cast<Rectangle*>(mouseOverGeometry.at(0));
        imgRect->setXYWH(0, 0, w, h);
    } else {
        refSpotCoords.clear();
        unsubscribe();
    }
}
