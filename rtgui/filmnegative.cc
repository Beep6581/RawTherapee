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
#include "filmnegative.h"

#include <iomanip>

#include "rtimage.h"
#include "options.h"
#include "editwidgets.h"
#include "eventmapper.h"


using namespace rtengine;
using namespace rtengine::procparams;


FilmNegative::FilmNegative () : FoldableToolPanel(this, "filmnegative", M("TP_FILMNEGATIVE_LABEL"), false, true), EditSubscriber(ET_OBJECTS)
{

    auto mkExponentAdjuster = [this](Glib::ustring label, double defaultVal) {
        Adjuster *adj = Gtk::manage(new Adjuster (label, 0.3, 6, 0.001, defaultVal)); //exponent
        adj->setAdjusterListener (this);
    
        if (adj->delay < options.adjusterMaxDelay) {
            adj->delay = options.adjusterMaxDelay;
        }
    
        adj->show();
        return adj;
    };

    redExp   = mkExponentAdjuster(M("TP_FILMNEGATIVE_RED"), 2.72);
    greenExp = mkExponentAdjuster(M("TP_FILMNEGATIVE_GREEN"), 2.0);
    blueExp  = mkExponentAdjuster(M("TP_FILMNEGATIVE_BLUE"), 1.72);

    redExp->set_sensitive(false);
    blueExp->set_sensitive(false);

    redRatio = redExp->getValue() / greenExp->getValue();
    blueRatio = blueExp->getValue() / greenExp->getValue();

    auto m = ProcEventMapper::getInstance();
    EvFilmNegativeEnabled = m->newEvent(FIRST, "HISTORY_MSG_FILMNEGATIVE_ENABLED");
    EvFilmNegativeExponents = m->newEvent(FIRST, "HISTORY_MSG_FILMNEGATIVE_EXPONENTS");

    lockChannels = Gtk::manage (new Gtk::CheckButton (M("TP_FILMNEGATIVE_LOCKCHANNELS")));
    lockChannels->set_tooltip_text(M("TP_FILMNEGATIVE_LOCKCHANNELS_TOOLTIP"));
    lockChannels->set_active (true);


    spotgrid = Gtk::manage(new Gtk::Grid());
    spotgrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(spotgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    spotbutton = Gtk::manage (new Gtk::ToggleButton (M("TP_FILMNEGATIVE_PICK")));
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
//    spotgrid->attach (*slab, 1, 0, 1, 1);
    // spotgrid->attach (*wbsizehelper, 2, 0, 1, 1);

    pack_start (*lockChannels, Gtk::PACK_SHRINK, 0);
    pack_start (*redExp, Gtk::PACK_SHRINK, 0);
    pack_start (*greenExp, Gtk::PACK_SHRINK, 0);
    pack_start (*blueExp, Gtk::PACK_SHRINK, 0);
    pack_start (*spotgrid, Gtk::PACK_SHRINK, 0 );

    lockChannels->signal_toggled().connect( sigc::mem_fun(*this, &FilmNegative::lockChannelsToggled) );
    spotbutton->signal_toggled().connect( sigc::mem_fun(*this, &FilmNegative::editToggled) );
//    spotsize->signal_changed().connect( sigc::mem_fun(*this, &WhiteBalance::spotSizeChanged) );


    // Editing geometry; create the spot rectangle
    Rectangle *spotRect = new Rectangle();
    spotRect->filled = false;
    
    EditSubscriber::visibleGeometry.push_back( spotRect );

    // Stick a dummy rectangle over the whole image in mouseOverGeometry.
    // This is to make sure the getCursor call is fired everywhere.
    Rectangle *imgRect = new Rectangle();
    imgRect->filled = true;

    EditSubscriber::mouseOverGeometry.push_back( imgRect );

}

FilmNegative::~FilmNegative()
{
//    idle_register.destroy();

    for (std::vector<Geometry*>::const_iterator i = visibleGeometry.begin(); i != visibleGeometry.end(); ++i) {
        delete *i;
    }

    for (std::vector<Geometry*>::const_iterator i = mouseOverGeometry.begin(); i != mouseOverGeometry.end(); ++i) {
        delete *i;
    }

}



void FilmNegative::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvFilmNegativeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvFilmNegativeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvFilmNegativeEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void FilmNegative::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        if(a == redExp || a == greenExp || a == blueExp) {
            disableListener();
            if(a == greenExp) {
                redExp->setValue(a->getValue() * redRatio);
                blueExp->setValue(a->getValue() * blueRatio);
            } else if(a == redExp) {
                redRatio = newval / greenExp->getValue();
            } else if(a == blueExp) {
                blueRatio = newval / greenExp->getValue();
            }
            enableListener();

            if(getEnabled()) {
                listener->panelChanged (EvFilmNegativeExponents, Glib::ustring::compose (
                    "R=%1 ; G=%2 ; B=%3", redExp->getTextValue(), greenExp->getTextValue(), blueExp->getTextValue()));
            }
        }
    }
}

void FilmNegative::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void FilmNegative::setEditProvider (EditDataProvider* provider)
{
    EditSubscriber::setEditProvider(provider);
}

void FilmNegative::editToggled ()
{
    if (spotbutton->get_active()) {
        subscribe();

        int w, h;
        getEditProvider()->getImageSize(w, h);

        // Stick a dummy rectangle over the whole image in mouseOverGeometry.
        // This is to make sure the getCursor call is fired everywhere.
        const auto imgRect = static_cast<Rectangle*>(mouseOverGeometry.at(0));
        imgRect->setXYWH(0, 0, w, h);

    } else {
        this->refSpotCoords.clear();        
        unsubscribe();
    }
}

void FilmNegative::lockChannelsToggled ()
{
    bool unlocked = !lockChannels->get_active();
    redExp->set_sensitive(unlocked);
    blueExp->set_sensitive(unlocked);
}

void FilmNegative::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if(pedited) {
        redExp->setEditedState(pedited->filmNegative.redExp ? Edited : UnEdited);
        greenExp->setEditedState(pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueExp->setEditedState(pedited->filmNegative.blueExp ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->filmNegative.enabled);
    }

    setEnabled(pp->filmNegative.enabled);
    redExp->setValue(pp->filmNegative.redExp);
    greenExp->setValue(pp->filmNegative.greenExp);
    blueExp->setValue(pp->filmNegative.blueExp);

    enableListener ();
}

void FilmNegative::write (ProcParams* pp, ParamsEdited* pedited)
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

void FilmNegative::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
    redExp->setValue(defParams->filmNegative.redExp);
    greenExp->setValue(defParams->filmNegative.greenExp);
    blueExp->setValue(defParams->filmNegative.blueExp);

    if (pedited) {
        redExp->setDefaultEditedState (pedited->filmNegative.redExp ? Edited : UnEdited);
        greenExp->setDefaultEditedState (pedited->filmNegative.greenExp ? Edited : UnEdited);
        blueExp->setDefaultEditedState (pedited->filmNegative.blueExp ? Edited : UnEdited);
    } else {
        redExp->setDefaultEditedState (Irrelevant);
        greenExp->setDefaultEditedState (Irrelevant);
        blueExp->setDefaultEditedState (Irrelevant);
    }
}

void FilmNegative::setBatchMode (bool batchMode)
{
    if(batchMode) {
        spotConn.disconnect();
        lockChannelsConn.disconnect();
        removeIfThere(this, spotgrid, false);
        removeIfThere(this, lockChannels, false);
        redExp->set_sensitive(true);
        blueExp->set_sensitive(true);
        ToolPanel::setBatchMode (batchMode);
        redExp->showEditedCB ();
        greenExp->showEditedCB ();
        blueExp->showEditedCB ();
    }
}

bool FilmNegative::mouseOver(int modifierKey)
{
    EditDataProvider *provider = getEditProvider();
    const auto spotRect = static_cast<Rectangle*>(visibleGeometry.at(0));
    spotRect->setXYWH(provider->posImage.x - 16, provider->posImage.y - 16, 32, 32);

    return true;
}

bool FilmNegative::button1Pressed(int modifierKey)
{
    EditDataProvider *provider = getEditProvider();

    if(provider) { // debug. remove me
        printf("x=%d y=%d pv1=%f pv2=%f pv3=%f\n", provider->posImage.x, provider->posImage.y, provider->getPipetteVal1(), provider->getPipetteVal2(), provider->getPipetteVal3());
    }
    
    EditSubscriber::action = EditSubscriber::Action::NONE;

    if (listener) {

        refSpotCoords.push_back(provider->posImage);

        if(refSpotCoords.size() == 2) {
            
            // User has selected 2 reference gray spots. Calculating new exponents
            // from channel values and updating parameters.

            float newExps[3];
            if(fnp->getFilmNegativeExponents(refSpotCoords[0], refSpotCoords[1], newExps)) {
                disableListener();
                redExp->setValue(newExps[0]);
                greenExp->setValue(newExps[1]);
                blueExp->setValue(newExps[2]);
                redRatio = redExp->getValue() / greenExp->getValue();
                blueRatio = blueExp->getValue() / greenExp->getValue();
                enableListener();

                if (listener && getEnabled()) {
                    listener->panelChanged (EvFilmNegativeExponents, Glib::ustring::compose (
                        "R=%1 ; G=%2 ; B=%3", redExp->getTextValue(), greenExp->getTextValue(), blueExp->getTextValue()));
                }
            }

            switchOffEditMode();
        }
    }

    return true;
}

bool FilmNegative::button1Released ()
{
    EditDataProvider *provider = getEditProvider();

    if(provider) { // debug. remove me
        printf("x=%d y=%d pv1=%f pv2=%f pv3=%f\n", provider->posImage.x, provider->posImage.y, provider->getPipetteVal1(), provider->getPipetteVal2(), provider->getPipetteVal3());
    }

    EditSubscriber::action = EditSubscriber::Action::NONE;
    return true;
}

// TODO remove me ; couldn't make Action::PICKING work
bool FilmNegative::pick1 (bool picked) {
    EditDataProvider *provider = getEditProvider();
    if(provider) { // debug. remove me
        printf("Picked pick=%d x=%d y=%d pv1=%f pv2=%f pv3=%f\n", picked, provider->posImage.x, provider->posImage.y, provider->getPipetteVal1(), provider->getPipetteVal2(), provider->getPipetteVal3());
    }
    return true;
}

CursorShape FilmNegative::getCursor(int objectID) const
{
   return CSSpotWB;
}

void FilmNegative::switchOffEditMode ()
{
    refSpotCoords.clear();
    unsubscribe();
    spotbutton->set_active(false);
}