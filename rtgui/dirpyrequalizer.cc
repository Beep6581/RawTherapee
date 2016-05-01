/*
 *  This file is part of RawTherapee.
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
 *
 *  © 2010 Emil Martinec <ejmartin@uchicago.edu>
 */

#include "dirpyrequalizer.h"

using namespace rtengine;
using namespace rtengine::procparams;

DirPyrEqualizer::DirPyrEqualizer () : FoldableToolPanel(this, "dirpyrequalizer", M("TP_DIRPYREQUALIZER_LABEL"), true, true)
{

    std::vector<GradientMilestone> milestones;

    float r, g, b;
    Color::hsv2rgb01(0.7500, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.    , r, g, b) ); // hsv: 0.75   rad: -0.9
    Color::hsv2rgb01(0.8560, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.1470, r, g, b) ); // hsv: 0.856  rad: -0.4
    Color::hsv2rgb01(0.9200, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.2353, r, g, b) ); // hsv: 0.92   rad: -0.1
    Color::hsv2rgb01(0.9300, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.2647, r, g, b) ); // hsv: 0.93   rad:  0
    Color::hsv2rgb01(0.9600, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.3380, r, g, b) ); // hsv: 0.96   rad:  0.25
    Color::hsv2rgb01(1.0000, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.4412, r, g, b) ); // hsv: 1.     rad:  0.6
    Color::hsv2rgb01(0.0675, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.6176, r, g, b) ); // hsv: 0.0675 rad:  1.2
    Color::hsv2rgb01(0.0900, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.6764, r, g, b) ); // hsv: 0.09   rad:  1.4
    Color::hsv2rgb01(0.1700, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.7647, r, g, b) ); // hsv: 0.17   rad:  1.7
    Color::hsv2rgb01(0.2650, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(0.8824, r, g, b) ); // hsv: 0.265  rad:  2.1
    Color::hsv2rgb01(0.3240, 0.5, 0.5, r, g, b);
    milestones.push_back( GradientMilestone(1.    , r, g, b) ); // hsv: 0.324  rad:  2.5

    Gtk::VBox * cbVBox = Gtk::manage ( new Gtk::VBox());
    cbVBox->set_border_width(4);
    cbVBox->set_spacing(2);

    cdbox = Gtk::manage (new Gtk::HBox ());
    labmcd = Gtk::manage (new Gtk::Label (M("TP_CBDL_METHOD") + ":"));
    cdbox->pack_start (*labmcd, Gtk::PACK_SHRINK, 1);

    cbdlMethod = Gtk::manage (new MyComboBoxText ());
    cbdlMethod->append_text (M("TP_CBDL_BEF"));
    cbdlMethod->append_text (M("TP_CBDL_AFT"));
    cbdlMethod->set_active(0);
    cbdlMethodConn = cbdlMethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrEqualizer::cbdlMethodChanged) );
    cbdlMethod->set_tooltip_markup (M("TP_CBDL_METHOD_TOOLTIP"));
    cdbox->pack_start(*cbdlMethod);
    cbVBox->pack_start(*cdbox);
    pack_start(*cbVBox);

    setEnabledTooltipMarkup(M("TP_SHARPENING_TOOLTIP"));

    Gtk::HBox * buttonBox1 = Gtk::manage (new Gtk::HBox(true, 10));
    pack_start(*buttonBox1);

    Gtk::Button * lumacontrastMinusButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")));
    buttonBox1->pack_start(*lumacontrastMinusButton);
    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumacontrastMinusPressed));

    Gtk::Button * lumaneutralButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")));
    buttonBox1->pack_start(*lumaneutralButton);
    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumaneutralPressed));

    Gtk::Button * lumacontrastPlusButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")));
    buttonBox1->pack_start(*lumacontrastPlusButton);
    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumacontrastPlusPressed));

    buttonBox1->show_all_children();

    Gtk::HSeparator *separator2 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator2, Gtk::PACK_SHRINK, 2);

    for(int i = 0; i < 6; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);

        if     (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if(i == 5) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage ( new Adjuster (ss, 0, 4, 0.01, 1.0) );
        multiplier[i]->setAdjusterListener(this);
        pack_start(*multiplier[i]);
    }

    Gtk::HSeparator *separator3 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator3, Gtk::PACK_SHRINK, 2);

    threshold = Gtk::manage ( new Adjuster (M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 1, 0.01, 0.2) );
    threshold->setAdjusterListener(this);
    pack_start(*threshold);

    Gtk::HSeparator *separator4 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator4, Gtk::PACK_SHRINK, 2);
    /*
        algoHBox = Gtk::manage (new Gtk::HBox ());
        algoHBox->set_border_width (0);
        algoHBox->set_spacing (2);
        algoHBox->set_tooltip_markup (M("TP_DIRPYREQUALIZER_ALGO_TOOLTIP"));
    */
//    alLabel = Gtk::manage (new Gtk::Label (M("TP_DIRPYREQUALIZER_ALGO")+":"));
//  algoHBox->pack_start (*alLabel, Gtk::PACK_SHRINK);
    /*
        algo = Gtk::manage (new MyComboBoxText ());
        algo->append_text (M("TP_DIRPYREQUALIZER_ALGO_FI"));
        algo->append_text (M("TP_DIRPYREQUALIZER_ALGO_LA"));
        algo->set_active (1);
    //  algoHBox->pack_start (*algo);
    //  pack_start(*algoHBox);
        algoconn = algo->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrEqualizer::algoChanged) );
    */
    hueskin = Gtk::manage (new ThresholdAdjuster (M("TP_DIRPYREQUALIZER_HUESKIN"), -40., 210., -5., 25., 170., 120., 0, false));        //default (b_l 0, t_l 30, b_r 170, t_r 120);
    hueskin->set_tooltip_markup (M("TP_DIRPYREQUALIZER_HUESKIN_TOOLTIP"));

    hueskin->setBgGradient(milestones);
    pack_start(*hueskin);

    skinprotect = Gtk::manage ( new Adjuster (M("TP_DIRPYREQUALIZER_SKIN"), -100, 100, 1, 0.) );
    skinprotect->setAdjusterListener(this);
    pack_start(*skinprotect);
    skinprotect->set_tooltip_markup (M("TP_DIRPYREQUALIZER_SKIN_TOOLTIP"));

    gamutlab = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYREQUALIZER_ARTIF")));
    gamutlab->set_active (true);
    pack_start(*gamutlab);
    gamutlabConn = gamutlab->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrEqualizer::gamutlabToggled) );
    gamutlab->set_tooltip_markup (M("TP_DIRPYREQUALIZER_TOOLTIP"));

    hueskin->setAdjusterListener (this);

    show_all_children ();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}

DirPyrEqualizer::~DirPyrEqualizer ()
{

}

void DirPyrEqualizer::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
    cbdlMethodConn.block(true);

    if (pedited) {

        set_inconsistent (multiImage && !pedited->dirpyrequalizer.enabled);
        gamutlab->set_inconsistent (!pedited->dirpyrequalizer.gamutlab);

        if (!pedited->dirpyrequalizer.cbdlMethod) {
            cbdlMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        for(int i = 0; i < 6; i++) {
            multiplier[i]->setEditedState (pedited->dirpyrequalizer.mult[i] ? Edited : UnEdited);
        }

        threshold->setEditedState (pedited->dirpyrequalizer.threshold ? Edited : UnEdited);
        skinprotect->setEditedState (pedited->dirpyrequalizer.skinprotect ? Edited : UnEdited);
        hueskin->setEditedState     (pedited->dirpyrequalizer.hueskin ? Edited : UnEdited);
    }

    setEnabled(pp->dirpyrequalizer.enabled);

    /*
        algoconn.block(true);
        if (pedited && !pedited->dirpyrequalizer.algo)
            algo->set_active (2);
        else if (pp->dirpyrequalizer.algo=="FI")
            algo->set_active (0);
        else if (pp->dirpyrequalizer.algo=="LA")
            algo->set_active (1);
        algoconn.block(false);
        algoChanged();
    */
    gamutlabConn.block (true);
    gamutlab->set_active (pp->dirpyrequalizer.gamutlab);
    gamutlabConn.block (false);
    lastgamutlab = pp->dirpyrequalizer.gamutlab;

    for (int i = 0; i < 6; i++) {
        multiplier[i]->setValue(pp->dirpyrequalizer.mult[i]);
    }

    threshold->setValue(pp->dirpyrequalizer.threshold);
    skinprotect->setValue(pp->dirpyrequalizer.skinprotect);
    hueskin->setValue<int>(pp->dirpyrequalizer.hueskin);

    if (pp->dirpyrequalizer.cbdlMethod == "bef") {
        cbdlMethod->set_active (0);
    } else if (pp->dirpyrequalizer.cbdlMethod == "aft") {
        cbdlMethod->set_active (1);
    }

    cbdlMethodChanged ();
    cbdlMethodConn.block(false);

    enableListener ();
}

void DirPyrEqualizer::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->dirpyrequalizer.enabled = getEnabled();
    pp->dirpyrequalizer.gamutlab = gamutlab->get_active ();
    pp->dirpyrequalizer.hueskin        = hueskin->getValue<int> ();

    for (int i = 0; i < 6; i++) {
        pp->dirpyrequalizer.mult[i] = multiplier[i]->getValue();
    }

    pp->dirpyrequalizer.threshold = threshold->getValue();
    pp->dirpyrequalizer.skinprotect = skinprotect->getValue();

    if (pedited) {

        pedited->dirpyrequalizer.enabled =  !get_inconsistent();
        pedited->dirpyrequalizer.hueskin        = hueskin->getEditedState ();
        pedited->dirpyrequalizer.cbdlMethod    = cbdlMethod->get_active_text() != M("GENERAL_UNCHANGED");

        for(int i = 0; i < 6; i++) {
            pedited->dirpyrequalizer.mult[i] = multiplier[i]->getEditedState();
        }

        pedited->dirpyrequalizer.threshold = threshold->getEditedState();
        pedited->dirpyrequalizer.skinprotect = skinprotect->getEditedState();
//       pedited->dirpyrequalizer.algo          = algo->get_active_text()!=M("GENERAL_UNCHANGED");
    }


    if (cbdlMethod->get_active_row_number() == 0) {
        pp->dirpyrequalizer.cbdlMethod = "bef";
    } else if (cbdlMethod->get_active_row_number() == 1) {
        pp->dirpyrequalizer.cbdlMethod = "aft";
    }

    /*    if (algo->get_active_row_number()==0)
            pp->dirpyrequalizer.algo = "FI";
        else if (algo->get_active_row_number()==1)
            pp->dirpyrequalizer.algo = "LA";
            */
}
/*
void DirPyrEqualizer::algoChanged () {
    if (listener && (multiImage||enabled->get_active()) ) {
        listener->panelChanged (EvDirPyrEqualizeralg, algo->get_active_text ());}
}
*/
void DirPyrEqualizer::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    for (int i = 0; i < 6; i++) {
        multiplier[i]->setDefault(defParams->dirpyrequalizer.mult[i]);
    }

    threshold->setDefault(defParams->dirpyrequalizer.threshold);
    hueskin->setDefault<int> (defParams->dirpyrequalizer.hueskin);

    if (pedited) {
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefaultEditedState(pedited->dirpyrequalizer.mult[i] ? Edited : UnEdited);
        }

        threshold->setDefaultEditedState(pedited->dirpyrequalizer.threshold ? Edited : UnEdited);
        skinprotect->setDefaultEditedState(pedited->dirpyrequalizer.skinprotect ? Edited : UnEdited);
        hueskin->setDefaultEditedState  (pedited->dirpyrequalizer.hueskin ? Edited : UnEdited);
    } else {
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefaultEditedState(Irrelevant);
        }

        threshold->setDefaultEditedState(Irrelevant);
        skinprotect->setDefaultEditedState(Irrelevant);
        hueskin->setDefaultEditedState (Irrelevant);
    }
}

void DirPyrEqualizer::adjusterChanged (ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight)
{
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDirPyrEqualizerHueskin, hueskin->getHistoryString());
    }
}


void DirPyrEqualizer::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);

    for (int i = 0; i < 6; i++) {
        multiplier[i]->showEditedCB();
    }

    threshold->showEditedCB();
    skinprotect->showEditedCB();
    hueskin->showEditedCB ();
//   algo->append_text (M("GENERAL_UNCHANGED"));
}

void DirPyrEqualizer::cbdlMethodChanged()
{

    if (listener) {
        listener->panelChanged (EvcbdlMethod, cbdlMethod->get_active_text ());
    }
}



void DirPyrEqualizer::adjusterChanged (Adjuster* a, double newval)
{

    if (listener && getEnabled()) {
        if (a == threshold) {
            listener->panelChanged (EvDirPyrEqualizerThreshold,
                                    Glib::ustring::compose("%1",
                                            Glib::ustring::format(std::fixed, std::setprecision(2), threshold->getValue()))
                                   );
        } else if (a == skinprotect) {
            listener->panelChanged (EvDirPyrEqualizerSkin,
                                    Glib::ustring::compose("%1",
                                            Glib::ustring::format(std::fixed, std::setprecision(2), skinprotect->getValue()))
                                   );
        } else {
            listener->panelChanged (EvDirPyrEqualizer,
                                    Glib::ustring::compose("%1, %2, %3, %4, %5, %6",
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[0]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[1]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[2]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[3]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[4]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[5]->getValue()))
                                   );
        }
    }
}

void DirPyrEqualizer::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void DirPyrEqualizer::gamutlabToggled ()
{

    if (batchMode) {
        if (gamutlab->get_inconsistent()) {
            gamutlab->set_inconsistent (false);
            gamutlabConn.block (true);
            gamutlab->set_active (false);
            gamutlabConn.block (false);
        } else if (lastgamutlab) {
            gamutlab->set_inconsistent (true);
        }

        lastgamutlab = gamutlab->get_active ();
    }

    if (listener) {
        if (gamutlab->get_active ()) {
            listener->panelChanged (EvDirPyrEqlgamutlab, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDirPyrEqlgamutlab, M("GENERAL_DISABLED"));
        }
    }
}

void DirPyrEqualizer::lumaneutralPressed ()
{

    for (int i = 0; i < 6; i++) {
        multiplier[i]->setValue(1.0);
        adjusterChanged(multiplier[i], 1.0);
    }
}


void DirPyrEqualizer::lumacontrastPlusPressed ()
{

    for (int i = 0; i < 6; i++) {
        float inc = 0.05 * (6 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}


void DirPyrEqualizer::lumacontrastMinusPressed ()
{

    for (int i = 0; i < 6; i++) {
        float inc = -0.05 * (6 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}

void DirPyrEqualizer::setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool skinadd)
{

    for (int i = 0; i < 6; i++) {
        multiplier[i]->setAddMode(multiplieradd);
    }

    threshold->setAddMode(thresholdadd);
    skinprotect->setAddMode(skinadd);
}

void DirPyrEqualizer::trimValues (rtengine::procparams::ProcParams* pp)
{

    for (int i = 0; i < 6; i++) {
        multiplier[i]->trimValue(pp->dirpyrequalizer.mult[i]);
    }

    threshold->trimValue(pp->dirpyrequalizer.threshold);
    skinprotect->trimValue(pp->dirpyrequalizer.skinprotect);
}
