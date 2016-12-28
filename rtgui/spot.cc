/*
 *  This file is part of RawTherapee.
 */
#include "spot.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;

#define STATIC_VISIBLE_OBJ_NBR 6
#define STATIC_MO_OBJ_NBR 6

Spot::Spot() : FoldableToolPanel (this, "spot", M ("TP_SPOT_LABEL"), true, true), EditSubscriber (ET_OBJECTS), lastObject (-1), activeSpot (-1),
    sourceIcon ("spot-normal.png", "spot-active.png", "spot-active.png", "spot-prelight.png", "", Geometry::DP_CENTERCENTER), editedCheckBox(NULL)
{
    countLabel = Gtk::manage (new Gtk::Label (Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), 0)));

    edit = Gtk::manage (new Gtk::ToggleButton());
    edit->add (*Gtk::manage (new RTImage ("editmodehand.png")));
    editConn = edit->signal_toggled().connect ( sigc::mem_fun (*this, &Spot::editToggled) );

    reset = Gtk::manage (new Gtk::Button ());
    reset->add (*Gtk::manage (new RTImage ("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png")));
    reset->set_relief (Gtk::RELIEF_NONE);
    reset->set_border_width (0);
    reset->signal_clicked().connect ( sigc::mem_fun (*this, &Spot::resetPressed) );

    labelBox = Gtk::manage (new Gtk::HBox());
    labelBox->set_spacing (2);
    labelBox->pack_start (*countLabel, false, false, 0);
    labelBox->pack_end (*edit, false, false, 0);
    labelBox->pack_end (*reset, false, false, 0);
    pack_start (*labelBox);

    sourceIcon.datum = Geometry::IMAGE;
    sourceIcon.setActive (false);
    sourceIcon.state = Geometry::ACTIVE;
    sourceCircle.datum = Geometry::IMAGE;
    sourceCircle.setActive (false);
    sourceCircle.radiusInImageSpace = true;
    sourceMODisc.datum = Geometry::IMAGE;
    sourceMODisc.setActive (false);
    sourceMODisc.radiusInImageSpace = true;
    sourceMODisc.filled = true;
    sourceMODisc.innerLineWidth = 0.;
    targetCircle.datum = Geometry::IMAGE;
    targetCircle.setActive (false);
    targetCircle.radiusInImageSpace = true;
    targetMODisc.datum = Geometry::IMAGE;
    targetMODisc.setActive (false);
    targetMODisc.radiusInImageSpace = true;
    targetMODisc.filled = true;
    targetMODisc.innerLineWidth = 0.;
    sourceFeatherCircle.datum = Geometry::IMAGE;
    sourceFeatherCircle.setActive (false);
    sourceFeatherCircle.radiusInImageSpace = true;
    targetFeatherCircle.datum = Geometry::IMAGE;
    targetFeatherCircle.setActive (false);
    targetFeatherCircle.radiusInImageSpace = true;
    link.datum = Geometry::IMAGE;
    link.setActive (false);

    show_all();
}

Spot::~Spot()
{
    // delete all dynamically allocated geometry
    if (EditSubscriber::visibleGeometry.size()) {
        for (size_t i = 0; i < EditSubscriber::visibleGeometry.size() - STATIC_VISIBLE_OBJ_NBR; ++i) { // static visible geometry at the end if the list
            delete EditSubscriber::visibleGeometry.at (i);
        }
    }
    // We do not delete the mouseOverGeometry, because the referenced objects are either
    // shared with visibleGeometry or instantiated by the class's ctor
}

void Spot::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    size_t oldSize = spots.size();
    spots = pp->spot.entries;

    if (pedited) {
        set_inconsistent (multiImage && !pedited->spot.enabled);
    }

    setEnabled (pp->spot.enabled);
    lastEnabled = pp->spot.enabled;

    if (spots.size() != oldSize) {
        createGeometry();
    }

    updateGeometry();

    enableListener ();
}

void Spot::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->spot.enabled = getEnabled();
    pp->spot.entries = spots;

    if (pedited) {
        pedited->spot.enabled = !get_inconsistent();
        pedited->spot.entries = !editedCheckBox->get_active();
    }
}

void Spot::resetPressed()
{
    if (batchMode) {
        // no need to handle the Geometry in batch mode, since point editing is disabled
        spots.clear();
        editedConn.block (true);
        editedCheckBox->set_active (true);
        editedConn.block (false);

        if (listener) {
            listener->panelChanged (EvSpotEntry, Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), spots.size()));
        }
    } else {
        if (!spots.empty()) {
            spots.clear();
            activeSpot = -1;
            lastObject = -1;
            createGeometry();
            updateGeometry();

            if (listener) {
                listener->panelChanged (EvSpotEntry, Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), 0));
            }
        }
    }
}

void Spot::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);

    if (batchMode) {
        removeIfThere (labelBox, edit, false);

        if (!editedCheckBox) {
            removeIfThere (labelBox, countLabel, false);
            countLabel = NULL;
            editedCheckBox = Gtk::manage (new Gtk::CheckButton (Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), 0)));
            labelBox->pack_start (*editedCheckBox, Gtk::PACK_SHRINK, 2);
            labelBox->reorder_child (*editedCheckBox, 0);
            editedConn = editedCheckBox->signal_toggled().connect ( sigc::mem_fun (*this, &Spot::editedToggled) );
            editedCheckBox->show();
        }
    }
}

void Spot::editedToggled ()
{
    if (listener) {
        listener->panelChanged (EvSpotEntry, !editedCheckBox->get_active() ? M ("GENERAL_UNCHANGED") : Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), spots.size()));
    }
}


void Spot::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvSpotEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvSpotEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvSpotEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Spot::setEditProvider (EditDataProvider* provider)
{
    EditSubscriber::setEditProvider (provider);
}

void Spot::editToggled ()
{
    if (edit->get_active()) {
        subscribe();
    } else {
        unsubscribe();
    }
}

Geometry* Spot::getVisibleGeometryFromMO (int MOID)
{
    if (MOID == -1) {
        return NULL;
    }

    if (MOID == 0) {
        return getActiveSpotIcon();
    }

    if (MOID == 1) { // sourceMODisc
        return &sourceIcon;
    }

    return EditSubscriber::mouseOverGeometry.at (MOID);
}

void Spot::createGeometry ()
{
    int nbrEntry = spots.size();

    if (!batchMode) {
        countLabel->set_text(Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), nbrEntry));
    }

    //printf("CreateGeometry(%d)\n", nbrEntry);
    // delete all dynamically allocated geometry
    if (EditSubscriber::visibleGeometry.size() > STATIC_VISIBLE_OBJ_NBR)
        for (size_t i = 0; i < EditSubscriber::visibleGeometry.size() - STATIC_VISIBLE_OBJ_NBR; ++i) { // static visible geometry at the end if the list
            delete EditSubscriber::visibleGeometry.at (i);
        }

    size_t i = 0, j = 0;
    EditSubscriber::mouseOverGeometry.resize (STATIC_MO_OBJ_NBR + nbrEntry);
    EditSubscriber::visibleGeometry.resize (nbrEntry + STATIC_VISIBLE_OBJ_NBR);

    EditSubscriber::mouseOverGeometry.at (i++) = &targetMODisc;        // STATIC_MO_OBJ_NBR + 0
    EditSubscriber::mouseOverGeometry.at (i++) = &sourceMODisc;        // STATIC_MO_OBJ_NBR + 1
    EditSubscriber::mouseOverGeometry.at (i++) = &targetCircle;        // STATIC_MO_OBJ_NBR + 2
    EditSubscriber::mouseOverGeometry.at (i++) = &sourceCircle;        // STATIC_MO_OBJ_NBR + 3
    EditSubscriber::mouseOverGeometry.at (i++) = &targetFeatherCircle; // STATIC_MO_OBJ_NBR + 4
    EditSubscriber::mouseOverGeometry.at (i++) = &sourceFeatherCircle; // STATIC_MO_OBJ_NBR + 5

    // recreate all spots geometry
    Cairo::RefPtr<Cairo::ImageSurface> normalImg   = sourceIcon.getNormalImg();
    Cairo::RefPtr<Cairo::ImageSurface> prelightImg = sourceIcon.getPrelightImg();
    Cairo::RefPtr<Cairo::ImageSurface> activeImg   = sourceIcon.getActiveImg();

    for (; j < EditSubscriber::visibleGeometry.size() - STATIC_VISIBLE_OBJ_NBR; ++i, ++j) {
        EditSubscriber::mouseOverGeometry.at (i) = EditSubscriber::visibleGeometry.at (j) = new OPIcon (normalImg, activeImg, prelightImg, Cairo::RefPtr<Cairo::ImageSurface>(NULL), Cairo::RefPtr<Cairo::ImageSurface>(NULL), Geometry::DP_CENTERCENTER);
        EditSubscriber::visibleGeometry.at (j)->setActive (true);
        EditSubscriber::visibleGeometry.at (j)->datum = Geometry::IMAGE;
        EditSubscriber::visibleGeometry.at (j)->state = Geometry::NORMAL;
        //printf("mouseOverGeometry.at(%d) = %p\n", (unsigned int)i, (void*)EditSubscriber::mouseOverGeometry.at(i));
    }

    EditSubscriber::visibleGeometry.at (j++) = &sourceIcon;          // STATIC_VISIBLE_OBJ_NBR + 0
    EditSubscriber::visibleGeometry.at (j++) = &sourceFeatherCircle; // STATIC_VISIBLE_OBJ_NBR + 1
    EditSubscriber::visibleGeometry.at (j++) = &link;                // STATIC_VISIBLE_OBJ_NBR + 2
    EditSubscriber::visibleGeometry.at (j++) = &sourceCircle;        // STATIC_VISIBLE_OBJ_NBR + 3
    EditSubscriber::visibleGeometry.at (j++) = &targetFeatherCircle; // STATIC_VISIBLE_OBJ_NBR + 4
    EditSubscriber::visibleGeometry.at (j++) = &targetCircle;        // STATIC_VISIBLE_OBJ_NBR + 5
}

void Spot::updateGeometry()
{
    EditDataProvider* dataProvider = getEditProvider();

    if (dataProvider) {
        int imW, imH;
        dataProvider->getImageSize (imW, imH);

        if (activeSpot > -1) {
            // Target point circle
            targetCircle.center = spots.at (activeSpot).targetPos;
            targetCircle.radius = spots.at (activeSpot).radius;
            targetCircle.setActive (true);

            // Target point Mouse Over disc
            targetMODisc.center = targetCircle.center;
            targetMODisc.radius = targetCircle.radius;
            targetMODisc.setActive (true);

            // Source point Icon
            sourceIcon.position = spots.at (activeSpot).sourcePos;
            sourceIcon.setActive (true);

            // Source point circle
            sourceCircle.center = spots.at (activeSpot).sourcePos;
            sourceCircle.radius = spots.at (activeSpot).radius;
            sourceCircle.setActive (true);

            // Source point Mouse Over disc
            sourceMODisc.center = sourceCircle.center;
            sourceMODisc.radius = sourceCircle.radius;
            sourceMODisc.setActive (true);

            // Target point feather circle
            targetFeatherCircle.center = spots.at (activeSpot).targetPos;
            targetFeatherCircle.radius = float (spots.at (activeSpot).radius) * (1.f + spots.at (activeSpot).feather);
            targetFeatherCircle.radiusInImageSpace = true;
            targetFeatherCircle.setActive (true);

            // Source point feather circle
            sourceFeatherCircle.center = spots.at (activeSpot).sourcePos;
            sourceFeatherCircle.radius = targetFeatherCircle.radius;
            sourceFeatherCircle.setActive (true);

            // Link line
            PolarCoord p;
            p = targetCircle.center - sourceCircle.center;

            if (p.radius > sourceCircle.radius + targetCircle.radius) {
                PolarCoord p2 (sourceCircle.radius, p.angle);
                Coord p3;
                p3 = p2;
                link.begin = sourceCircle.center + p3;
                p2.set (targetCircle.radius, p.angle + 180);
                p3 = p2;
                link.end = targetCircle.center + p3;
                link.setActive (true);
            } else {
                link.setActive (false);
            }
        } else {
            targetCircle.setActive (false);
            targetMODisc.setActive (false);
            sourceIcon.setActive (false);
            sourceCircle.setActive (false);
            sourceMODisc.setActive (false);
            targetFeatherCircle.setActive (false);
            sourceFeatherCircle.setActive (false);
            link.setActive (false);
        }

        for (size_t i = 0; i < spots.size(); ++i) {
            // Target point icon
            OPIcon* geom = static_cast<OPIcon*> (EditSubscriber::visibleGeometry.at (i));
            geom->position = spots.at (i).targetPos;
            geom->setActive (true);

            if (int (i) == activeSpot) {
                geom->setHoverable (false);
            }
        }
    }
}

OPIcon *Spot::getActiveSpotIcon()
{
    if (activeSpot > -1) {
        return static_cast<OPIcon*> (EditSubscriber::visibleGeometry.at (activeSpot));
    }

    return NULL;
}

void Spot::addNewEntry()
{
    EditDataProvider* editProvider = getEditProvider();
    // we create a new entry
    SpotEntry se;
    se.targetPos = editProvider->posImage;
    se.sourcePos = se.targetPos;
    spots.push_back (se); // this make a copy of se ...
    activeSpot = spots.size() - 1;
    lastObject = 1;

    //printf("ActiveSpot = %d\n", activeSpot);

    createGeometry();
    updateGeometry();
    EditSubscriber::visibleGeometry.at (activeSpot)->state = Geometry::ACTIVE;
    sourceIcon.state = Geometry::DRAGGED;
    // TODO: find a way to disable the active spot's Mouse Over geometry but still displaying its location...

    if (listener) {
        listener->panelChanged (EvSpotEntry, M ("TP_SPOT_ENTRYCHANGED"));
    }
}

void Spot::deleteSelectedEntry()
{
    // delete the activeSpot
    if (activeSpot > -1) {
        std::vector<rtengine::SpotEntry>::iterator i = spots.begin();
        for (int j = 0; j < activeSpot; ++j) {
            ++i;
        }
        spots.erase (i);
    }

    lastObject = -1;
    activeSpot = -1;

    createGeometry();
    updateGeometry();

    if (listener) {
        listener->panelChanged (EvSpotEntry, M ("TP_SPOT_ENTRYCHANGED"));
    }
}

// TODO
CursorShape Spot::getCursor (const int objectID)
{
    return CSOpenHand;
}

bool Spot::mouseOver (const int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider && editProvider->object != lastObject) {
        if (lastObject > -1) {
            if (EditSubscriber::mouseOverGeometry.at (lastObject) == &targetMODisc) {
                getVisibleGeometryFromMO (lastObject)->state = Geometry::ACTIVE;
            } else {
                getVisibleGeometryFromMO (lastObject)->state = Geometry::NORMAL;
            }

            sourceIcon.state = Geometry::ACTIVE;
        }

        if (editProvider->object > -1) {
            getVisibleGeometryFromMO (editProvider->object)->state = Geometry::PRELIGHT;

            if (editProvider->object >= STATIC_MO_OBJ_NBR) {
                // a Spot is being edited
                int oldActiveSpot = activeSpot;
                activeSpot = editProvider->object - STATIC_MO_OBJ_NBR;

                if (activeSpot != oldActiveSpot) {
                    if (oldActiveSpot > -1) {
                        EditSubscriber::visibleGeometry.at (oldActiveSpot)->state = Geometry::NORMAL;
                        EditSubscriber::mouseOverGeometry.at (oldActiveSpot + STATIC_MO_OBJ_NBR)->state = Geometry::NORMAL;
                    }

                    EditSubscriber::visibleGeometry.at (activeSpot)->state = Geometry::PRELIGHT;
                    EditSubscriber::mouseOverGeometry.at (activeSpot + STATIC_MO_OBJ_NBR)->state = Geometry::PRELIGHT;
                    //printf("ActiveSpot = %d (was %d before)\n", activeSpot, oldActiveSpot);
                }
            }
        }

        lastObject = editProvider->object;

        if (lastObject > -1 && EditSubscriber::mouseOverGeometry.at (lastObject) == getActiveSpotIcon()) {
            lastObject = 0;  // targetMODisc
        }

        updateGeometry();
        return true;
    }

    return false;
}

// Create a new Target and Source point or start the drag of a Target point under the cursor
bool Spot::button1Pressed (const int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider) {
        if (lastObject == -1 && (modifierKey & GDK_CONTROL_MASK)) {
            addNewEntry();
            EditSubscriber::action = ES_ACTION_DRAGGING;
            return true;
        } else if (lastObject > -1) {
            getVisibleGeometryFromMO (lastObject)->state = Geometry::DRAGGED;
            EditSubscriber::action = ES_ACTION_DRAGGING;
            return true;
        }
    }

    return false;
}

// End the drag of a Target point
bool Spot::button1Released()
{
    Geometry *loGeom = getVisibleGeometryFromMO (lastObject);

    if (!loGeom) {
        EditSubscriber::action = ES_ACTION_NONE;
        return false;
    }

    loGeom->state = Geometry::PRELIGHT;
    EditSubscriber::action = ES_ACTION_NONE;
    updateGeometry();
    return true;
}

// Delete a point
bool Spot::button2Pressed (const int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (!editProvider || lastObject == -1 || activeSpot==-1) {
        return false;
    }

    if (!(modifierKey & (GDK_SHIFT_MASK|GDK_SHIFT_MASK))) {
        EditSubscriber::action = ES_ACTION_PICKING;
    }

    return false;
}

// Create a new Target and Source point or start the drag of a Target point under the cursor
bool Spot::button3Pressed (const int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (!editProvider || lastObject == -1 || activeSpot==-1) {
        return false;
    }

    if ((modifierKey & GDK_CONTROL_MASK) && (EditSubscriber::mouseOverGeometry.at (lastObject) == &targetMODisc || lastObject >= STATIC_MO_OBJ_NBR)) {
        lastObject = 1;  // sourceMODisc
        sourceIcon.state = Geometry::DRAGGED;
        EditSubscriber::action = ES_ACTION_DRAGGING;
        return true;
    } else if (!(modifierKey & (GDK_SHIFT_MASK|GDK_SHIFT_MASK))) {
        EditSubscriber::action = ES_ACTION_PICKING;
    }

    return false;
}

bool Spot::button3Released()
{
    Geometry *loGeom = getVisibleGeometryFromMO (lastObject);

    if (!loGeom) {
        EditSubscriber::action = ES_ACTION_NONE;
        return false;
    }

    lastObject = -1;
    sourceIcon.state = Geometry::ACTIVE;
    updateGeometry();
    EditSubscriber::action = ES_ACTION_NONE;
    return true;

    return false;
}

bool Spot::drag1 (const int modifierKey)
{
    EditDataProvider *editProvider = getEditProvider();
    int imW, imH;
    editProvider->getImageSize (imW, imH);
    bool modified = false;

    //printf("Drag1 / LastObject=%d\n", lastObject);

    Geometry *loGeom = EditSubscriber::mouseOverGeometry.at (lastObject);

    if (loGeom == &sourceMODisc) {
        //printf("sourceMODisc / deltaPrevImage = %d / %d\n", editProvider->deltaPrevImage.x, editProvider->deltaPrevImage.y);
        rtengine::Coord currPos = spots.at (activeSpot).sourcePos;
        spots.at (activeSpot).sourcePos += editProvider->deltaPrevImage;
        spots.at (activeSpot).sourcePos.clip (imW, imH);

        if (spots.at (activeSpot).sourcePos != currPos) {
            modified = true;
        }

        EditSubscriber::mouseOverGeometry.at (activeSpot + STATIC_MO_OBJ_NBR)->state = Geometry::DRAGGED;
    } else if (loGeom == &targetMODisc || lastObject >= STATIC_MO_OBJ_NBR) {
        //printf("targetMODisc / deltaPrevImage = %d / %d\n", editProvider->deltaPrevImage.x, editProvider->deltaPrevImage.y);
        rtengine::Coord currPos = spots.at (activeSpot).targetPos;
        spots.at (activeSpot).targetPos += editProvider->deltaPrevImage;
        spots.at (activeSpot).targetPos.clip (imW, imH);

        if (spots.at (activeSpot).targetPos != currPos) {
            modified = true;
        }
    } else if (loGeom == &sourceCircle) {
        //printf("sourceCircle / deltaPrevImage = %d / %d\n", editProvider->deltaImage.x, editProvider->deltaImage.y);
        int lastRadius = spots.at (activeSpot).radius;
        rtengine::Coord currPos = editProvider->posImage + editProvider->deltaImage;
        rtengine::PolarCoord currPolar (currPos - spots.at (activeSpot).sourcePos);
        spots.at (activeSpot).radius = LIM<int> (int (currPolar.radius), SpotParams::minRadius, SpotParams::maxRadius);

        if (spots.at (activeSpot).radius != lastRadius) {
            modified = true;
        }
    } else if (loGeom == &targetCircle) {
        //printf("targetCircle / deltaPrevImage = %d / %d\n", editProvider->deltaImage.x, editProvider->deltaImage.y);
        int lastRadius = spots.at (activeSpot).radius;
        rtengine::Coord currPos = editProvider->posImage + editProvider->deltaImage;
        rtengine::PolarCoord currPolar (currPos - spots.at (activeSpot).targetPos);
        spots.at (activeSpot).radius = LIM<int> (int (currPolar.radius), SpotParams::minRadius, SpotParams::maxRadius);

        if (spots.at (activeSpot).radius != lastRadius) {
            modified = true;
        }
    } else if (loGeom == &sourceFeatherCircle) {
        //printf("sourceFeatherCircle / deltaPrevImage = %d / %d\n", editProvider->deltaImage.x, editProvider->deltaImage.y);
        float currFeather = spots.at (activeSpot).feather;
        rtengine::Coord currPos = editProvider->posImage + editProvider->deltaImage;
        rtengine::PolarCoord currPolar (currPos -spots.at (activeSpot).sourcePos);
        spots.at (activeSpot).feather = LIM01<float> ((currPolar.radius - double (spots.at (activeSpot).radius)) / double (spots.at (activeSpot).radius));

        if (spots.at (activeSpot).feather != currFeather) {
            modified = true;
        }
    } else if (loGeom == &targetFeatherCircle) {
        //printf("targetFeatherCircle / deltaPrevImage = %d / %d\n", editProvider->deltaImage.x, editProvider->deltaImage.y);
        float currFeather = spots.at (activeSpot).feather;
        rtengine::Coord currPos = editProvider->posImage + editProvider->deltaImage;
        rtengine::PolarCoord currPolar (currPos - spots.at (activeSpot).targetPos);
        spots.at (activeSpot).feather = LIM01<float> ((currPolar.radius - double (spots.at (activeSpot).radius)) / double (spots.at (activeSpot).radius));

        if (spots.at (activeSpot).feather != currFeather) {
            modified = true;
        }
    }

    if (listener && modified) {
        updateGeometry();
        listener->panelChanged (EvSpotEntry, M ("TP_SPOT_ENTRYCHANGED"));
    }

    return modified;
}

bool Spot::drag3 (const int modifierKey)
{
    EditDataProvider *editProvider = getEditProvider();
    int imW, imH;
    editProvider->getImageSize (imW, imH);
    bool modified = false;

    Geometry *loGeom = EditSubscriber::mouseOverGeometry.at (lastObject);

    if (loGeom == &sourceMODisc) {
        rtengine::Coord currPos = spots.at (activeSpot).sourcePos;
        spots.at (activeSpot).sourcePos += editProvider->deltaPrevImage;
        spots.at (activeSpot).sourcePos.clip (imW, imH);

        if (spots.at (activeSpot).sourcePos != currPos) {
            modified = true;
        }
    }

    if (listener) {
        updateGeometry();
        listener->panelChanged (EvSpotEntry, M ("TP_SPOT_ENTRYCHANGED"));
    }

    return modified;
}

bool Spot::pick2(const bool picked)
{
    return pick3(picked);
}

bool Spot::pick3 (const bool picked)
{
    EditDataProvider* editProvider = getEditProvider();

    if (!picked) {
        if (editProvider->object != lastObject) {
            return false;
        }
    }

    // Object is picked, we delete it
    deleteSelectedEntry();
    EditSubscriber::action = ES_ACTION_NONE;
    updateGeometry();
    return true;
}


void Spot::switchOffEditMode ()
{
    if (edit->get_active()) {
        // switching off the toggle button
        bool wasBlocked = editConn.block (true);
        edit->set_active (false);

        if (!wasBlocked) {
            editConn.block (false);
        }
    }

    EditSubscriber::switchOffEditMode();  // disconnect
}

