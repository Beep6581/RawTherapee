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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
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

Adjuster* createExponentAdjuster(AdjusterListener* listener, const Glib::ustring& label, double minV, double maxV, double defaultVal)
{
    Adjuster* const adj = Gtk::manage(new Adjuster(label, minV, maxV, 0.001, defaultVal));
    adj->setAdjusterListener(listener);
    adj->setLogScale(6, 1, true);

    if (adj->delay < options.adjusterMaxDelay) {
        adj->delay = options.adjusterMaxDelay;
    }

    adj->show();
    return adj;
}

Glib::ustring fmt(const float v)
{
    if (v <= 0.f) {
        return "-";
    } else {
        return Glib::ustring::format(std::fixed, std::setprecision(1), v);
    }    
}

}

FilmNegative::FilmNegative() :
    FoldableToolPanel(this, "filmnegative", M("TP_FILMNEGATIVE_LABEL"), false, true),
    EditSubscriber(ET_OBJECTS),
    evFilmNegativeExponents(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_VALUES")),
    evFilmNegativeEnabled(ProcEventMapper::getInstance()->newEvent(ALL, "HISTORY_MSG_FILMNEGATIVE_ENABLED")),
    evFilmBaseValues(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_FILMBASE")),
    evFilmNegativeBalance(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_BALANCE")),
    evOldFilmNegativeExponents(ProcEventMapper::getInstance()->newEvent(ALL, "HISTORY_MSG_FILMNEGATIVE_VALUES")),
    filmBaseGreenValue(0.f),
    fnp(nullptr),
    greenExp(createExponentAdjuster(this, M("TP_FILMNEGATIVE_GREEN"), 0.3, 4, 1.5)),  // master exponent (green channel)
    redRatio(createExponentAdjuster(this, M("TP_FILMNEGATIVE_RED"), 0.3, 3, (2.04 / 1.5))), // ratio of red exponent to master exponent
    blueRatio(createExponentAdjuster(this, M("TP_FILMNEGATIVE_BLUE"), 0.3, 3, (1.29 / 1.5))), // ratio of blue exponent to master exponent
    spotgrid(Gtk::manage(new Gtk::Grid())),
    spotbutton(Gtk::manage(new Gtk::ToggleButton(M("TP_FILMNEGATIVE_PICK")))),
    filmBaseLabel(Gtk::manage(new Gtk::Label(M("TP_FILMNEGATIVE_BASE_LABEL")))), //, Gtk::ALIGN_CENTER))),
    filmBaseSpotButton(Gtk::manage(new Gtk::ToggleButton(M("TP_FILMNEGATIVE_FILMBASE_PICK")))),
    redBalance(createExponentAdjuster(this, M("TP_FILMNEGATIVE_REDBALANCE"), 0.05, 20, 1)),  // red balance
    blueBalance(createExponentAdjuster(this, M("TP_FILMNEGATIVE_BLUEBALANCE"), 0.05, 20, 1)) // blue balance
{
    spotgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(spotgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    setExpandAlignProperties(spotbutton, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    spotbutton->get_style_context()->add_class("independent");
    spotbutton->set_tooltip_text(M("TP_FILMNEGATIVE_GUESS_TOOLTIP"));
    spotbutton->set_image(*Gtk::manage(new RTImage("color-picker-small.png")));

    filmBaseSpotButton->set_tooltip_text(M("TP_FILMNEGATIVE_FILMBASE_TOOLTIP"));
//    setExpandAlignProperties(filmBaseValuesLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    setExpandAlignProperties(filmBaseLabel, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    filmBaseLabel->set_justify(Gtk::Justification::JUSTIFY_CENTER);

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

    spotgrid->attach(*spotbutton, 0, 1, 1, 1);
    // spotgrid->attach (*slab, 1, 0, 1, 1);
    // spotgrid->attach (*wbsizehelper, 2, 0, 1, 1);

    pack_start(*greenExp, Gtk::PACK_SHRINK, 0);
    pack_start(*redRatio, Gtk::PACK_SHRINK, 0);
    pack_start(*blueRatio, Gtk::PACK_SHRINK, 0);
    pack_start(*spotgrid, Gtk::PACK_SHRINK, 0);

//    pack_start(*oldMethod, Gtk::PACK_SHRINK, 0);

    Gtk::HSeparator* const sep = Gtk::manage(new  Gtk::HSeparator());
    sep->get_style_context()->add_class("grid-row-separator");
    pack_start(*sep, Gtk::PACK_SHRINK, 0);

//    Gtk::Grid* const fbGrid = Gtk::manage(new Gtk::Grid());
//    fbGrid->attach(*filmBaseLabel, 0, 0, 1, 1);
//    fbGrid->attach(*filmBaseValuesLabel, 1, 0, 1, 1);
//    pack_start(*fbGrid, Gtk::PACK_SHRINK, 0);
    pack_start(*filmBaseLabel, Gtk::PACK_SHRINK, 0);

    pack_start(*filmBaseSpotButton, Gtk::PACK_SHRINK, 0);

    pack_start(*redBalance, Gtk::PACK_SHRINK, 0);
    pack_start(*blueBalance, Gtk::PACK_SHRINK, 0);

    spotbutton->signal_toggled().connect(sigc::mem_fun(*this, &FilmNegative::editToggled));
    // spotsize->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::spotSizeChanged) );

    filmBaseSpotButton->signal_toggled().connect(sigc::mem_fun(*this, &FilmNegative::baseSpotToggled));

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
    idle_register.destroy();

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
        redRatio->setEditedState(pedited->filmNegative.redRatio ? Edited : UnEdited);
        greenExp->setEditedState(pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueRatio->setEditedState(pedited->filmNegative.blueRatio ? Edited : UnEdited);
        redBalance->setEditedState(pedited->filmNegative.redBalance ? Edited : UnEdited);
        blueBalance->setEditedState(pedited->filmNegative.blueBalance ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->filmNegative.enabled);
    }

    setEnabled(pp->filmNegative.enabled);
    redRatio->setValue(pp->filmNegative.redRatio);
    greenExp->setValue(pp->filmNegative.greenExp);
    blueRatio->setValue(pp->filmNegative.blueRatio);

    filmBaseGreenValue = pp->filmNegative.greenBase;

    if (filmBaseGreenValue == -1) {

        filmBaseLabel->set_markup(batchMode
            ? M("TP_FILMNEGATIVE_UPGRADE_LABEL_BATCH")
            : M("TP_FILMNEGATIVE_UPGRADE_LABEL"));

        filmBaseSpotButton->set_label(M("TP_FILMNEGATIVE_UPGRADE_BUTTON"));
        filmBaseSpotButton->set_tooltip_markup(M("TP_FILMNEGATIVE_UPGRADE_TOOLTIP"));

        filmBaseSpotButton->set_sensitive(true);

        redBalance->set_sensitive(false);
        blueBalance->set_sensitive(false);
    } else {
        // If base values are not set in params, estimated values will be passed in later
        // (after processing) via FilmNegListener
        filmBaseLabel->set_markup(
            Glib::ustring::compose(M("TP_FILMNEGATIVE_BASE_LABEL"), fmt(filmBaseGreenValue)));            
        filmBaseSpotButton->set_label(M("TP_FILMNEGATIVE_FILMBASE_PICK"));
        filmBaseSpotButton->set_sensitive(true);
        filmBaseSpotButton->set_tooltip_text(M("TP_FILMNEGATIVE_FILMBASE_TOOLTIP"));

        redBalance->set_sensitive(true);
        blueBalance->set_sensitive(true);
        redBalance->setValue(pp->filmNegative.redBalance);
        blueBalance->setValue(pp->filmNegative.blueBalance);
    }

    enableListener();
}

void FilmNegative::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->filmNegative.redRatio = redRatio->getValue();
    pp->filmNegative.greenExp = greenExp->getValue();
    pp->filmNegative.blueRatio = blueRatio->getValue();

    pp->filmNegative.enabled = getEnabled();

    if (pedited) {
        pedited->filmNegative.redRatio = redRatio->getEditedState();
        pedited->filmNegative.greenExp = greenExp->getEditedState();
        pedited->filmNegative.blueRatio = blueRatio->getEditedState();
        pedited->filmNegative.greenBase = filmBaseGreenValue != pp->filmNegative.greenBase;
        pedited->filmNegative.redBalance = redBalance->getEditedState();
        pedited->filmNegative.blueBalance = blueBalance->getEditedState();
        pedited->filmNegative.enabled = !get_inconsistent();
    }

    pp->filmNegative.greenBase = filmBaseGreenValue;
    
    pp->filmNegative.redBalance = redBalance->getValue();
    pp->filmNegative.blueBalance = blueBalance->getValue();

}

void FilmNegative::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    redRatio->setValue(defParams->filmNegative.redRatio);
    greenExp->setValue(defParams->filmNegative.greenExp);
    blueRatio->setValue(defParams->filmNegative.blueRatio);
    redBalance->setValue(defParams->filmNegative.redBalance);
    blueBalance->setValue(defParams->filmNegative.blueBalance);

    if (pedited) {
        redRatio->setDefaultEditedState(pedited->filmNegative.redRatio ? Edited : UnEdited);
        greenExp->setDefaultEditedState(pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueRatio->setDefaultEditedState(pedited->filmNegative.blueRatio ? Edited : UnEdited);
        redBalance->setDefaultEditedState(pedited->filmNegative.redBalance ? Edited : UnEdited);
        blueBalance->setDefaultEditedState(pedited->filmNegative.blueBalance ? Edited : UnEdited);
    } else {
        redRatio->setDefaultEditedState(Irrelevant);
        greenExp->setDefaultEditedState(Irrelevant);
        blueRatio->setDefaultEditedState(Irrelevant);
        redBalance->setDefaultEditedState(Irrelevant);
        blueBalance->setDefaultEditedState(Irrelevant);
    }
}

void FilmNegative::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    if (batchMode) {
        removeIfThere(this, spotgrid, false);
        removeIfThere(this, filmBaseSpotButton, false);
        redRatio->showEditedCB();
        greenExp->showEditedCB();
        blueRatio->showEditedCB();
        redBalance->showEditedCB();
        blueBalance->showEditedCB();
    }
}

void FilmNegative::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == redRatio || a == greenExp || a == blueRatio) {
            // Exponents were changed; if we are in legacy mode,
            // trigger ALL processing, otherwise ALLNORAW
            listener->panelChanged(
                (filmBaseGreenValue == -1) ? evOldFilmNegativeExponents : evFilmNegativeExponents,
                Glib::ustring::compose(
                    "Ref=%1\nR=%2\nB=%3",
                    greenExp->getValue(),
                    redRatio->getValue(),
                    blueRatio->getValue()
                )
            );
        } else if (a == redBalance || a == blueBalance) {
            listener->panelChanged(
                evFilmNegativeBalance,
                Glib::ustring::compose(
                    "B=%1 G=%2",
                    redBalance->getValue(),
                    blueBalance->getValue()
                )
            );
        }
    }
}

void FilmNegative::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(evFilmNegativeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(evFilmNegativeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(evFilmNegativeEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void FilmNegative::filmBaseValuesChanged(float greenBase, float redBal, float blueBal)
{

    idle_register.add(
        [this, greenBase, redBal, blueBal]() -> bool
        {
            filmBaseGreenValue = greenBase;

            disableListener();

            filmBaseLabel->set_markup(
                Glib::ustring::compose(M("TP_FILMNEGATIVE_BASE_LABEL"), fmt(filmBaseGreenValue)));            
            filmBaseSpotButton->set_label(M("TP_FILMNEGATIVE_FILMBASE_PICK"));
            filmBaseSpotButton->set_sensitive(true);
            filmBaseSpotButton->set_tooltip_text(M("TP_FILMNEGATIVE_FILMBASE_TOOLTIP"));

            if (filmBaseGreenValue == -1) {
                redBalance->set_sensitive(false);
                blueBalance->set_sensitive(false);
            } else {
                redBalance->set_sensitive(true);
                blueBalance->set_sensitive(true);
                filmBaseSpotButton->set_sensitive(true);
                redBalance->setValue(redBal);
                blueBalance->setValue(blueBal);
            }
            enableListener();
            return false;
        }
    );

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
        if (spotbutton->get_active()) {

            refSpotCoords.push_back(provider->posImage);

            float rBal, bBal;

            if (refSpotCoords.size() == 2) {
                // User has selected 2 reference gray spots. Calculating new exponents
                // from channel values and updating parameters.
                std::array<float, 3> newExps = { };

                if (fnp->getFilmNegativeExponents(refSpotCoords[0], refSpotCoords[1], newExps, rBal, bBal)) {
                    disableListener();
                    // Leaving green exponent unchanged, setting red and blue exponents based on
                    // the ratios between newly calculated exponents.
                    redRatio->setValue(newExps[0] / newExps[1]);
                    blueRatio->setValue(newExps[2] / newExps[1]);

                    if (filmBaseGreenValue != -1) {
                        redBalance->setValue(rBal);
                        blueBalance->setValue(bBal);
                    }

                    enableListener();

                }

                switchOffEditMode();

                if (listener && getEnabled()) {
                    if (filmBaseGreenValue == -1) {
                        // Legacy mode: trigger entire raw processing
                        listener->panelChanged(
                            evOldFilmNegativeExponents,
                            Glib::ustring::compose(
                                "Ref=%1\nR=%2\nB=%3",
                                greenExp->getValue(),
                                redRatio->getValue(),
                                blueRatio->getValue()
                            )
                        );
                    } else {
                        listener->panelChanged(
                            evFilmNegativeExponents,
                            Glib::ustring::compose(
                                "Ref=%1\nR=%2\nB=%3\nBalance: R=%4 B=%5",
                                greenExp->getValue(),
                                redRatio->getValue(),
                                blueRatio->getValue(),
                                redBalance->getValue(),
                                blueBalance->getValue()
                            )
                        );
                    }

                }

            } else if (filmBaseGreenValue != -1) {

                fnp->getFilmNegativeBalance(provider->posImage, 32, rBal, bBal);

                disableListener();
                redBalance->setValue(rBal);
                blueBalance->setValue(bBal);
                enableListener();

                // const Glib::ustring vs = formatBaseValues(filmBaseValues);

                // filmBaseValuesLabel->set_text(vs);

                // if (listener && getEnabled()) {
                //     listener->panelChanged(evFilmBaseValues, vs);
                // }

                if (listener && getEnabled()) {
                    listener->panelChanged(
                        evFilmNegativeBalance,
                        Glib::ustring::compose(
                            "R=%1 B=%2",
                            redBalance->getValue(),
                            blueBalance->getValue()
                        )
                    );
                }

            }


        } else if (filmBaseSpotButton->get_active()) {

            filmBaseGreenValue = fnp->getFilmBaseGreen(provider->posImage, 32);

            filmBaseLabel->set_text(
                Glib::ustring::compose(M("TP_FILMNEGATIVE_BASE_LABEL"), fmt(filmBaseGreenValue)));

            if (listener && getEnabled()) {
                listener->panelChanged(evFilmBaseValues, fmt(filmBaseGreenValue));
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
    filmBaseSpotButton->set_active(false);
}

void FilmNegative::editToggled()
{
    if (spotbutton->get_active()) {

        filmBaseSpotButton->set_active(false);
        refSpotCoords.clear();

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


void FilmNegative::baseSpotToggled()
{
    if (filmBaseSpotButton->get_active()) {

        spotbutton->set_active(false);
        refSpotCoords.clear();

        if (filmBaseGreenValue == -1) {
            filmBaseSpotButton->set_active(false);
            filmBaseSpotButton->set_sensitive(false);

            filmBaseGreenValue = 0.f;
            listener->panelChanged(evFilmNegativeEnabled, M("HISTORY_MSG_FILMNEGATIVE_UPGRADE"));
        } else {

            subscribe();

            int w, h;
            getEditProvider()->getImageSize(w, h);

            // Stick a dummy rectangle over the whole image in mouseOverGeometry.
            // This is to make sure the getCursor() call is fired everywhere.
            Rectangle* const imgRect = static_cast<Rectangle*>(mouseOverGeometry.at(0));
            imgRect->setXYWH(0, 0, w, h);
        }
    } else {
        refSpotCoords.clear();
        unsubscribe();
    }
}
