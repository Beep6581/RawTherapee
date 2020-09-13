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

    adj->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    adj->show();
    return adj;
}

Adjuster* createLevelAdjuster(AdjusterListener* listener, const Glib::ustring& label)
{
    Adjuster* const adj = Gtk::manage(new Adjuster(label, 1.0, 65535.0, 1.0, 65535.0 / 512.0));
    adj->setAdjusterListener(listener);
    adj->setLogScale(6, 4000.0, true);

    adj->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    adj->show();
    return adj;
}

Adjuster* createBalanceAdjuster(AdjusterListener* listener, const Glib::ustring& label, double minV, double maxV, double defaultVal)
{
    Adjuster* const adj = Gtk::manage(new Adjuster(label, minV, maxV, 0.01, defaultVal));
    adj->setAdjusterListener(listener);
    adj->setLogScale(6, 1, true);

    adj->setDelay(std::max(options.adjusterMinDelay, options.adjusterMaxDelay));

    adj->show();
    return adj;
}


Glib::ustring fmt(const RGB& rgb)
{
    if (rgb.r <= 0.f && rgb.g <= 0.f && rgb.b <= 0.f) {
        return "- - -";
    } else {
        return Glib::ustring::format(std::fixed, std::setprecision(1), rgb.r) + " " +
               Glib::ustring::format(std::fixed, std::setprecision(1), rgb.g) + " " +
               Glib::ustring::format(std::fixed, std::setprecision(1), rgb.b);
    }    
}


RGB getFilmNegativeExponents(const RGB &ref1, const RGB &ref2) // , const RGB &clearValsOut, const RGB &denseValsOut)
{
    using rtengine::settings;

    RGB clearVals = ref1;
    RGB denseVals = ref2;

    // Detect which one is the dense spot, based on green channel
    if (clearVals.g < denseVals.g) {
        std::swap(clearVals, denseVals);
        //std::swap(clearValsOut, denseValsOut);
    }

    if (settings->verbose) {
        printf("Clear input values: R=%g G=%g B=%g\n", static_cast<double>(clearVals.r), static_cast<double>(clearVals.g), static_cast<double>(clearVals.b));
        printf("Dense input values: R=%g G=%g B=%g\n", static_cast<double>(denseVals.r), static_cast<double>(denseVals.g), static_cast<double>(denseVals.b));

        // printf("Clear output values: R=%g G=%g B=%g\n", static_cast<double>(clearValsOut.r), static_cast<double>(clearValsOut.g), static_cast<double>(clearValsOut.b));
        // printf("Dense output values: R=%g G=%g B=%g\n", static_cast<double>(denseValsOut.r), static_cast<double>(denseValsOut.g), static_cast<double>(denseValsOut.b));
    }

    const float denseGreenRatio = clearVals.g / denseVals.g;

    // Calculate logarithms in arbitrary base
    const auto logBase =
        [](float base, float num) -> float
        {
            return std::log(num) / std::log(base);
        };

    // const auto ratio =
    //     [](float a, float b) -> float
    //     {
    //         return a > b ? a / b : b / a;
    //     };

    RGB newExps;
    newExps.r = logBase(clearVals.r / denseVals.r, denseGreenRatio );
    newExps.g = 1.f; // logBase(ratio(clearVals.g, denseVals.g), ratio(denseValsOut.g, clearValsOut.g) );
    newExps.b = logBase(clearVals.b / denseVals.b, denseGreenRatio );



    if (settings->verbose) {
        printf("New exponents:  R=%g G=%g B=%g\n", static_cast<double>(newExps.r), static_cast<double>(newExps.g), static_cast<double>(newExps.b));
    }

    // // Re-adjust color balance based on dense spot values and new exponents
    // calcBalance(rtengine::max(static_cast<float>(params->filmNegative.refInput.g), 1.f),
    //     -newExps[0], -newExps[1], -newExps[2],
    //     denseVals[0], denseVals[1], denseVals[2],
    //     rBal, bBal);

    return newExps;

}

}

FilmNegative::FilmNegative() :
    FoldableToolPanel(this, "filmnegative", M("TP_FILMNEGATIVE_LABEL"), false, true),
    EditSubscriber(ET_OBJECTS),
    evFilmNegativeExponents(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_VALUES")),
    evFilmNegativeEnabled(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_ENABLED")),
    evFilmNegativeRefSpot(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_REF_SPOT")),
    evFilmNegativeBalance(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_BALANCE")),
    evFilmNegativeColorSpace(ProcEventMapper::getInstance()->newEvent(ALLNORAW, "HISTORY_MSG_FILMNEGATIVE_COLORSPACE")),
    refInputValues({0.f, 0.f, 0.f}),
    paramsUpgraded(false),
    fnp(nullptr),
    colorSpace(Gtk::manage(new MyComboBoxText())),
    greenExp(createExponentAdjuster(this, M("TP_FILMNEGATIVE_GREEN"), 0.3, 4, 1.5)),  // master exponent (green channel)
    redRatio(createExponentAdjuster(this, M("TP_FILMNEGATIVE_RED"), 0.3, 3, (2.04 / 1.5))), // ratio of red exponent to master exponent
    blueRatio(createExponentAdjuster(this, M("TP_FILMNEGATIVE_BLUE"), 0.3, 3, (1.29 / 1.5))), // ratio of blue exponent to master exponent
    spotgrid(Gtk::manage(new Gtk::Grid())),
    spotbutton(Gtk::manage(new Gtk::ToggleButton(M("TP_FILMNEGATIVE_PICK")))),
    refInputLabel(Gtk::manage(new Gtk::Label(M("TP_FILMNEGATIVE_REF_LABEL")))), //, Gtk::ALIGN_CENTER))),
    refSpotButton(Gtk::manage(new Gtk::ToggleButton(M("TP_FILMNEGATIVE_REF_PICK")))),
    outputLevel(createLevelAdjuster(this, M("TP_FILMNEGATIVE_OUT_LEVEL"))),  // ref level
    greenBalance(createBalanceAdjuster(this, M("TP_FILMNEGATIVE_GREENBALANCE"), 0.1, 10, 1.0)),  // green balance
    blueBalance(createBalanceAdjuster(this, M("TP_FILMNEGATIVE_BLUEBALANCE"), 0.1, 10, 1.0))  // blue balance
{
    spotgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(spotgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    setExpandAlignProperties(spotbutton, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    spotbutton->get_style_context()->add_class("independent");
    spotbutton->set_tooltip_text(M("TP_FILMNEGATIVE_GUESS_TOOLTIP"));
    spotbutton->set_image(*Gtk::manage(new RTImage("color-picker-small.png")));

    refSpotButton->set_tooltip_text(M("TP_FILMNEGATIVE_REF_TOOLTIP"));

    setExpandAlignProperties(refInputLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
//    refInputLabel->set_justify(Gtk::Justification::JUSTIFY_CENTER);
//    refInputLabel->set_line_wrap(true);

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

    colorSpace->append(M("TP_FILMNEGATIVE_COLORSPACE_CAMERA"));
    colorSpace->append(M("TP_FILMNEGATIVE_COLORSPACE_WORKING"));
    colorSpace->append(M("TP_FILMNEGATIVE_COLORSPACE_BUILTIN"));
    setExpandAlignProperties(colorSpace, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Grid* csGrid = Gtk::manage(new Gtk::Grid());
    Gtk::Label* csLabel = Gtk::manage (new Gtk::Label (M("TP_FILMNEGATIVE_COLORSPACE")));
    csGrid->attach(*csLabel, 0, 0, 1, 1);
    csGrid->attach(*colorSpace, 1, 0, 1, 1);

    pack_start(*csGrid);

    colorSpace->set_active((int)rtengine::procparams::FilmNegativeParams::ColorSpace::BUILTIN);
    colorSpace->signal_changed().connect(sigc::mem_fun(*this, &FilmNegative::colorSpaceChanged));
    colorSpace->show();

    pack_start(*greenExp, Gtk::PACK_SHRINK, 0);
    pack_start(*redRatio, Gtk::PACK_SHRINK, 0);
    pack_start(*blueRatio, Gtk::PACK_SHRINK, 0);
    pack_start(*spotgrid, Gtk::PACK_SHRINK, 0);

//    pack_start(*oldMethod, Gtk::PACK_SHRINK, 0);

    Gtk::HSeparator* const sep = Gtk::manage(new  Gtk::HSeparator());
    sep->get_style_context()->add_class("grid-row-separator");
    pack_start(*sep, Gtk::PACK_SHRINK, 0);

//    Gtk::Grid* const fbGrid = Gtk::manage(new Gtk::Grid());
//    fbGrid->attach(*refInputLabel, 0, 0, 1, 1);
//    fbGrid->attach(*filmBaseValuesLabel, 1, 0, 1, 1);
//    pack_start(*fbGrid, Gtk::PACK_SHRINK, 0);
    pack_start(*refInputLabel, Gtk::PACK_SHRINK, 0);

    pack_start(*outputLevel, Gtk::PACK_SHRINK, 0);
    pack_start(*greenBalance, Gtk::PACK_SHRINK, 0);
    pack_start(*blueBalance, Gtk::PACK_SHRINK, 0);

    pack_start(*refSpotButton, Gtk::PACK_SHRINK, 0);

    spotbutton->signal_toggled().connect(sigc::mem_fun(*this, &FilmNegative::editToggled));
    // spotsize->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::spotSizeChanged) );

    refSpotButton->signal_toggled().connect(sigc::mem_fun(*this, &FilmNegative::baseSpotToggled));

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
        outputLevel->setEditedState(pedited->filmNegative.refOutput  ? Edited : UnEdited);
        blueBalance->setEditedState(pedited->filmNegative.refOutput ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->filmNegative.enabled);
    }

    setEnabled(pp->filmNegative.enabled);

    colorSpace->set_active((int)pp->filmNegative.colorSpace);
    redRatio->setValue(pp->filmNegative.redRatio);
    greenExp->setValue(pp->filmNegative.greenExp);
    blueRatio->setValue(pp->filmNegative.blueRatio);

    refInputValues = pp->filmNegative.refInput;

    // If base values are not set in params, estimated values will be passed in later
    // (after processing) via FilmNegListener
    refInputLabel->set_markup(
        Glib::ustring::compose(M("TP_FILMNEGATIVE_REF_LABEL"), fmt(refInputValues)));

    outputLevel->setValue(pp->filmNegative.refOutput.r);
    greenBalance->setValue(pp->filmNegative.refOutput.g / pp->filmNegative.refOutput.r);
    blueBalance->setValue(pp->filmNegative.refOutput.b / pp->filmNegative.refOutput.r);

    enableListener();
}

void FilmNegative::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    if (colorSpace->get_active_row_number() != 3) {  // UNCHANGED entry, see setBatchMode
        pp->filmNegative.colorSpace = rtengine::procparams::FilmNegativeParams::ColorSpace(colorSpace->get_active_row_number());
    }

    pp->filmNegative.redRatio = redRatio->getValue();
    pp->filmNegative.greenExp = greenExp->getValue();
    pp->filmNegative.blueRatio = blueRatio->getValue();

    pp->filmNegative.enabled = getEnabled();

    if (pedited) {
        pedited->filmNegative.colorSpace = colorSpace->get_active_row_number() != 3; // UNCHANGED entry, see setBatchMode
        pedited->filmNegative.redRatio = redRatio->getEditedState();
        pedited->filmNegative.greenExp = greenExp->getEditedState();
        pedited->filmNegative.blueRatio = blueRatio->getEditedState();
        pedited->filmNegative.refInput = refInputValues != pp->filmNegative.refInput;
        pedited->filmNegative.refOutput = outputLevel->getEditedState() ||  greenBalance->getEditedState() || blueBalance->getEditedState();
        pedited->filmNegative.enabled = !get_inconsistent();
    }

    pp->filmNegative.refInput = refInputValues;
    
    pp->filmNegative.refOutput.r = outputLevel->getValue();
    pp->filmNegative.refOutput.g = greenBalance->getValue() * outputLevel->getValue();
    pp->filmNegative.refOutput.b = blueBalance->getValue() * outputLevel->getValue();

    if (paramsUpgraded) {
        pp->filmNegative.backCompat = BackCompat::CURRENT;
    }

}

void FilmNegative::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    redRatio->setValue(defParams->filmNegative.redRatio);
    greenExp->setValue(defParams->filmNegative.greenExp);
    blueRatio->setValue(defParams->filmNegative.blueRatio);

    outputLevel->setValue(defParams->filmNegative.refOutput.r);
    greenBalance->setValue(defParams->filmNegative.refOutput.g / defParams->filmNegative.refOutput.r);
    blueBalance->setValue(defParams->filmNegative.refOutput.b / defParams->filmNegative.refOutput.r);

    if (pedited) {
        redRatio->setDefaultEditedState(pedited->filmNegative.redRatio ? Edited : UnEdited);
        greenExp->setDefaultEditedState(pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueRatio->setDefaultEditedState(pedited->filmNegative.blueRatio ? Edited : UnEdited);

        outputLevel->setDefaultEditedState(pedited->filmNegative.refOutput ? Edited : UnEdited);
        greenBalance->setDefaultEditedState(pedited->filmNegative.refOutput ? Edited : UnEdited);
        blueBalance->setDefaultEditedState(pedited->filmNegative.refOutput ? Edited : UnEdited);
    } else {
        redRatio->setDefaultEditedState(Irrelevant);
        greenExp->setDefaultEditedState(Irrelevant);
        blueRatio->setDefaultEditedState(Irrelevant);
        outputLevel->setDefaultEditedState(Irrelevant);
        greenBalance->setDefaultEditedState(Irrelevant);
        blueBalance->setDefaultEditedState(Irrelevant);
    }
}

void FilmNegative::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    if (batchMode) {
        removeIfThere(this, spotgrid, false);
        removeIfThere(this, refSpotButton, false);
        colorSpace->append(M("GENERAL_UNCHANGED"));
        colorSpace->set_active_text(M("GENERAL_UNCHANGED"));
        redRatio->showEditedCB();
        greenExp->showEditedCB();
        blueRatio->showEditedCB();
        outputLevel->showEditedCB();
        greenBalance->showEditedCB();
        blueBalance->showEditedCB();
    }
}

void FilmNegative::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == redRatio || a == greenExp || a == blueRatio) {
            listener->panelChanged(
                evFilmNegativeExponents,
                Glib::ustring::compose(
                    "Ref=%1\nR=%2\nB=%3",
                    greenExp->getValue(),
                    redRatio->getValue(),
                    blueRatio->getValue()
                )
            );
        } else if (a == outputLevel || a == greenBalance || a == blueBalance) {

            listener->panelChanged(
                evFilmNegativeBalance,
                Glib::ustring::compose(
                    "Lev=%1 G=%2 B=%3",
                    outputLevel->getValue(),
                    greenBalance->getValue(),
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

void FilmNegative::colorSpaceChanged()
{
    if (listener) {
        listener->panelChanged(evFilmNegativeColorSpace, colorSpace->get_active_text());
    }
}

void FilmNegative::filmBaseValuesChanged(const RGB &refInput, const RGB &refOutput)
{

    idle_register.add(
        [this, refInput, refOutput]() -> bool
        {
            refInputValues = refInput;
            paramsUpgraded = true;

            disableListener();

            refInputLabel->set_markup(
                Glib::ustring::compose(M("TP_FILMNEGATIVE_REF_LABEL"), fmt(refInputValues)));            
            
            outputLevel->setValue(refOutput.r);
            greenBalance->setValue(refOutput.g / refOutput.r);
            blueBalance->setValue(refOutput.b / refOutput.r);

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

            if (refSpotCoords.size() == 2) {
                // User has selected 2 reference gray spots. Calculating new exponents
                // from channel values and updating parameters.

                RGB ref1, ref2, dummy;
                if (fnp->getFilmNegativeSpot(refSpotCoords[0], 32, ref1, dummy) &&
                    fnp->getFilmNegativeSpot(refSpotCoords[1], 32, ref2, dummy)) {

                    disableListener();

                    RGB newExps = getFilmNegativeExponents(ref1, ref2);

                    // Leaving green exponent unchanged, setting red and blue exponents based on
                    // the ratios between newly calculated exponents.
                    redRatio->setValue(newExps.r / newExps.g);
                    blueRatio->setValue(newExps.b / newExps.g);

                    enableListener();

                    if (getEnabled()) {
                        listener->panelChanged(
                            evFilmNegativeExponents,
                            Glib::ustring::compose(
                                "Ref=%1\nR=%2\nB=%3",
                                greenExp->getValue(),
                                redRatio->getValue(),
                                blueRatio->getValue()
                            )
                        );
                    }

                }

                switchOffEditMode();

            }


        } else if (refSpotButton->get_active()) {

            RGB filmBaseOut;
            fnp->getFilmNegativeSpot(provider->posImage, 32, refInputValues, filmBaseOut);
            
            disableListener();

            outputLevel->setValue(filmBaseOut.r); // rtengine::MAXVALF / 512.f);
            greenBalance->setValue(1.0); // filmBaseOut.g / filmBaseOut.r);
            blueBalance->setValue(1.0); // filmBaseOut.g / filmBaseOut.r);

            refInputLabel->set_text(
                Glib::ustring::compose(M("TP_FILMNEGATIVE_REF_LABEL"), fmt(refInputValues)));

            enableListener();

            listener->panelChanged(
                evFilmNegativeRefSpot,
                Glib::ustring::compose(
                    "%1, %2, %3",
                    round(refInputValues.r), round(refInputValues.g), round(refInputValues.b)
                )
            );

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
    refSpotButton->set_active(false);
}

void FilmNegative::editToggled()
{
    if (spotbutton->get_active()) {

        refSpotButton->set_active(false);
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
    if (refSpotButton->get_active()) {

        spotbutton->set_active(false);
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
