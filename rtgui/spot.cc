/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

#include "editcallbacks.h"
#include "spot.h"
#include "rtimage.h"
#include <iomanip>
#include "rtengine/rt_math.h"
#include "guiutils.h"
#include "eventmapper.h"
#include "rtengine/refreshmap.h"

using namespace rtengine;
using namespace rtengine::procparams;

namespace
{

enum GeometryIndex {
    MO_TARGET_DISK,
    MO_SOURCE_DISC,
    MO_TARGET_CIRCLE,
    MO_SOURCE_CIRCLE,
    MO_TARGET_FEATHER_CIRCLE,
    MO_SOURCE_FEATHER_CIRCLE,
    MO_OBJECT_COUNT,

    VISIBLE_SOURCE_ICON = 0,
    VISIBLE_SOURCE_FEATHER_CIRCLE,
    VISIBLE_LINK,
    VISIBLE_SOURCE_CIRCLE,
    VISIBLE_TARGET_FEATHER_CIRCLE,
    VISIBLE_TARGET_CIRCLE,
    VISIBLE_OBJECT_COUNT
};

}

const Glib::ustring Spot::TOOL_NAME = "spot";

Spot::Spot() :
    FoldableToolPanel(this, TOOL_NAME, M ("TP_SPOT_LABEL"), true, true),
    EditSubscriber(ET_OBJECTS),
    draggedSide(DraggedSide::NONE),
    lastObject(-1),
    activeSpot(-1),
    sourceIcon("spot-normal", "spot-active", "spot-prelight", "", "", Geometry::DP_CENTERCENTER),
    editedCheckBox(nullptr)
{
    countLabel = Gtk::manage (new Gtk::Label (Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), 0)));

    edit = Gtk::manage (new Gtk::ToggleButton());
    edit->add (*Gtk::manage (new RTImage ("edit-point")));
    editConn = edit->signal_toggled().connect ( sigc::mem_fun (*this, &Spot::editToggled) );
    edit->set_tooltip_text(M("TP_SPOT_HINT"));

    reset = Gtk::manage (new Gtk::Button ());
    reset->add (*Gtk::manage (new RTImage ("undo-small")));
    reset->set_relief (Gtk::RELIEF_NONE);
    reset->set_border_width (0);
    reset->signal_clicked().connect ( sigc::mem_fun (*this, &Spot::resetPressed) );

    spotSize = Gtk::manage(new Adjuster(M("TP_SPOT_DEFAULT_SIZE"), SpotParams::minRadius, SpotParams::maxRadius, 1, 25));

    labelBox = Gtk::manage (new Gtk::Box());
    labelBox->set_spacing (2);
    labelBox->pack_start (*countLabel, false, false, 0);
    labelBox->pack_end (*edit, false, false, 0);
    labelBox->pack_end (*reset, false, false, 0);
    labelBox->pack_end (*spotSize, false, false, 0);
    pack_start (*labelBox);

    sourceIcon.datum = Geometry::IMAGE;
    sourceIcon.setActive (false);
    sourceIcon.state = Geometry::ACTIVE;
    sourceCircle.datum = Geometry::IMAGE;
    sourceCircle.setActive (false);
    sourceCircle.radiusInImageSpace = true;
    sourceCircle.setDashed(true);
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
    sourceFeatherCircle.setDashed(true);
    sourceFeatherCircle.innerLineWidth = 0.7;
    targetFeatherCircle.datum = Geometry::IMAGE;
    targetFeatherCircle.setActive (false);
    targetFeatherCircle.radiusInImageSpace = true;
    targetFeatherCircle.innerLineWidth = 0.7;
    link.datum = Geometry::IMAGE;
    link.setActive (false);

    auto m = ProcEventMapper::getInstance();
    EvSpotEnabled = m->newEvent(ALLNORAW, "HISTORY_MSG_SPOT");
    EvSpotEnabledOPA = m->newEvent(SPOTADJUST, "HISTORY_MSG_SPOT");
    EvSpotEntry = m->newEvent(SPOTADJUST, "HISTORY_MSG_SPOT_ENTRY");
    EvSpotEntryOPA = m->newEvent(SPOTADJUST, "HISTORY_MSG_SPOT_ENTRY");

    show_all();
}

Spot::~Spot()
{
    // delete all dynamically allocated geometry
    if (EditSubscriber::visibleGeometry.size()) {
        for (size_t i = 0; i < EditSubscriber::visibleGeometry.size() - VISIBLE_OBJECT_COUNT; ++i) { // static visible geometry at the end of the list
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
    activeSpot = -1;
    lastObject = -1;

    if (batchMode) {
        editedCheckBox->set_label(Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), spots.size()));
    }
    else {
        if (spots.size() != oldSize) {
            createGeometry();
        }

        updateGeometry();
    }

    enableListener ();
}

void Spot::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->spot.enabled = getEnabled();
    pp->spot.entries = spots;

    if (pedited) {
        pedited->spot.enabled = !get_inconsistent();
        pedited->spot.entries = editedCheckBox->get_active();
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

        editedCheckBox->set_label(Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), spots.size()));

        if (listener) {
            listener->panelChanged (EvSpotEntry, Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), 0));
        }
    } else {
        if (!spots.empty()) {
            EditSubscriber::action = EditSubscriber::Action::NONE;
            spots.clear();
            activeSpot = -1;
            lastObject = -1;
            createGeometry();
            updateGeometry();

            if (listener) {
                listener->panelChanged (edit->get_active() ? EvSpotEntryOPA : EvSpotEntry, Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), 0));
            }
        }
    }
}

/**
 * Release anything that's currently being dragged.
 */
void Spot::releaseEdit()
{
    Geometry *loGeom = getVisibleGeometryFromMO (lastObject);

    EditSubscriber::action = EditSubscriber::Action::NONE;

    if (!loGeom) {
        return;
    }

    loGeom->state = Geometry::NORMAL;
    sourceIcon.state = Geometry::NORMAL;
    draggedSide = DraggedSide::NONE;
    updateGeometry();
}

void Spot::setBatchMode (bool batchMode)
{
    ToolPanel::setBatchMode (batchMode);

    if (batchMode) {
        removeIfThere (labelBox, edit, false);

        if (!editedCheckBox) {
            removeIfThere (labelBox, countLabel, false);
            countLabel = nullptr;
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
            listener->panelChanged (edit->get_active() ? EvSpotEnabledOPA : EvSpotEnabled, M ("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (edit->get_active() ? EvSpotEnabledOPA : EvSpotEnabled, M ("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (edit->get_active() ? EvSpotEnabledOPA : EvSpotEnabled, M ("GENERAL_DISABLED"));
        }
    }
}

void Spot::setEditProvider (EditDataProvider* provider)
{
    EditSubscriber::setEditProvider (provider);
}

void Spot::editToggled ()
{
    if (listener) {
        if (edit->get_active()) {
            listener->setTweakOperator(this);
            listener->refreshPreview(EvSpotEnabledOPA); // reprocess the preview w/o creating History entry
            subscribe();
        } else {
            releaseEdit();
            unsubscribe();
            listener->unsetTweakOperator(this);
            listener->refreshPreview(EvSpotEnabled); // reprocess the preview w/o creating History entry
        }
    }
}

Geometry* Spot::getVisibleGeometryFromMO (int MOID)
{
    if (MOID == -1) {
        return nullptr;
    }

    if (MOID == MO_TARGET_DISK) {
        return getActiveSpotIcon();
    }

    if (MOID == MO_SOURCE_DISC) {
        return &sourceIcon;
    }

    if (MOID > MO_OBJECT_COUNT) {
        return EditSubscriber::visibleGeometry.at(MOID - MO_OBJECT_COUNT);
    }

    return EditSubscriber::mouseOverGeometry.at (MOID);
}

void Spot::createGeometry ()
{
    int nbrEntry = spots.size();

    if (!batchMode) {
        countLabel->set_text (Glib::ustring::compose (M ("TP_SPOT_COUNTLABEL"), nbrEntry));
    }

    //printf("CreateGeometry(%d)\n", nbrEntry);
    // delete all dynamically allocated geometry
    if (EditSubscriber::visibleGeometry.size() > VISIBLE_OBJECT_COUNT)
        for (size_t i = 0; i < EditSubscriber::visibleGeometry.size() - VISIBLE_OBJECT_COUNT; ++i) { // static visible geometry at the end of the list
            delete EditSubscriber::visibleGeometry.at (i);
        }

    // mouse over geometry starts with the static geometry, then the spot's icon geometry
    EditSubscriber::mouseOverGeometry.resize (MO_OBJECT_COUNT + nbrEntry);
    // visible geometry starts with the spot's icon geometry, then the static geometry
    EditSubscriber::visibleGeometry.resize (nbrEntry + VISIBLE_OBJECT_COUNT);

    size_t i = 0, j = 0;
    assert(i == MO_TARGET_DISK);
    EditSubscriber::mouseOverGeometry.at (i++) = &targetMODisc;        // MO_OBJECT_COUNT + 0
    assert(i == MO_SOURCE_DISC);
    EditSubscriber::mouseOverGeometry.at (i++) = &sourceMODisc;        // MO_OBJECT_COUNT + 1
    assert(i == MO_TARGET_CIRCLE);
    EditSubscriber::mouseOverGeometry.at (i++) = &targetCircle;        // MO_OBJECT_COUNT + 2
    assert(i == MO_SOURCE_CIRCLE);
    EditSubscriber::mouseOverGeometry.at (i++) = &sourceCircle;        // MO_OBJECT_COUNT + 3
    assert(i == MO_TARGET_FEATHER_CIRCLE);
    EditSubscriber::mouseOverGeometry.at (i++) = &targetFeatherCircle; // MO_OBJECT_COUNT + 4
    assert(i == MO_SOURCE_FEATHER_CIRCLE);
    EditSubscriber::mouseOverGeometry.at (i++) = &sourceFeatherCircle; // MO_OBJECT_COUNT + 5

    // recreate all spots geometry
    std::shared_ptr<RTSurface> normalImg   = sourceIcon.getNormalImg();
    std::shared_ptr<RTSurface> prelightImg = sourceIcon.getPrelightImg();
    std::shared_ptr<RTSurface> activeImg   = sourceIcon.getActiveImg();

    for (; j < EditSubscriber::visibleGeometry.size() - VISIBLE_OBJECT_COUNT; ++i, ++j) {
        EditSubscriber::mouseOverGeometry.at (i) = EditSubscriber::visibleGeometry.at (j) = new OPIcon (normalImg, activeImg, prelightImg, nullptr, nullptr, Geometry::DP_CENTERCENTER);
        EditSubscriber::visibleGeometry.at (j)->setActive (true);
        EditSubscriber::visibleGeometry.at (j)->datum = Geometry::IMAGE;
        EditSubscriber::visibleGeometry.at (j)->state = Geometry::NORMAL;
        //printf("mouseOverGeometry.at(%d) = %p\n", (unsigned int)i, (void*)EditSubscriber::mouseOverGeometry.at(i));
    }

    int visibleOffset = j;
    assert(j - visibleOffset == VISIBLE_SOURCE_ICON);
    EditSubscriber::visibleGeometry.at (j++) = &sourceIcon;          // VISIBLE_OBJECT_COUNT + 0
    assert(j - visibleOffset == VISIBLE_SOURCE_FEATHER_CIRCLE);
    EditSubscriber::visibleGeometry.at (j++) = &sourceFeatherCircle; // VISIBLE_OBJECT_COUNT + 1
    assert(j - visibleOffset == VISIBLE_LINK);
    EditSubscriber::visibleGeometry.at (j++) = &link;                // VISIBLE_OBJECT_COUNT + 2
    assert(j - visibleOffset == VISIBLE_SOURCE_CIRCLE);
    EditSubscriber::visibleGeometry.at (j++) = &sourceCircle;        // VISIBLE_OBJECT_COUNT + 3
    assert(j - visibleOffset == VISIBLE_TARGET_FEATHER_CIRCLE);
    EditSubscriber::visibleGeometry.at (j++) = &targetFeatherCircle; // VISIBLE_OBJECT_COUNT + 4
    assert(j - visibleOffset == VISIBLE_TARGET_CIRCLE);
    EditSubscriber::visibleGeometry.at (j++) = &targetCircle;        // VISIBLE_OBJECT_COUNT + 5
    static_cast<void>(visibleOffset);
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
                link.setActive (draggedSide == DraggedSide::NONE);
            } else {
                link.setActive (false);
            }

            sourceCircle.setVisible(draggedSide != DraggedSide::SOURCE);
            targetCircle.setVisible(draggedSide != DraggedSide::TARGET);
        } else {
            targetCircle.state = Geometry::NORMAL;
            sourceCircle.state = Geometry::NORMAL;
            targetFeatherCircle.state = Geometry::NORMAL;
            sourceFeatherCircle.state = Geometry::NORMAL;

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

    return nullptr;
}

void Spot::addNewEntry()
{
    EditDataProvider* editProvider = getEditProvider();
    // we create a new entry
    SpotEntry se;
    se.radius = spotSize->getIntValue();
    se.targetPos = editProvider->posImage;
    se.sourcePos = se.targetPos;
    spots.push_back (se); // this make a copy of se ...
    activeSpot = spots.size() - 1;
    lastObject = MO_SOURCE_DISC;

    //printf("ActiveSpot = %d\n", activeSpot);

    createGeometry();
    updateGeometry();
    EditSubscriber::visibleGeometry.at (activeSpot)->state = Geometry::ACTIVE;
    sourceIcon.state = Geometry::DRAGGED;
    // TODO: find a way to disable the active spot's Mouse Over geometry but still displaying its location...

    if (listener) {
        listener->panelChanged (EvSpotEntryOPA, M ("TP_SPOT_ENTRYCHANGED"));
    }
}

void Spot::deleteSelectedEntry()
{
    // delete the activeSpot
    if (activeSpot > -1) {
        std::vector<rtengine::procparams::SpotEntry>::iterator i = spots.begin();

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

CursorShape Spot::getCursor (int objectID, int xPos, int yPos) const
{
    const EditDataProvider* editProvider = getEditProvider();
    if (editProvider && activeSpot > -1) {
        if (draggedSide != DraggedSide::NONE) {
            return CSEmpty;
        }

        if (objectID == MO_TARGET_DISK || objectID == MO_SOURCE_DISC) {
            return CSMove2D;
        }
        if (objectID >= MO_TARGET_CIRCLE && objectID <= MO_SOURCE_FEATHER_CIRCLE) {
            Coord delta(Coord(xPos, yPos) - ((objectID == MO_SOURCE_CIRCLE || objectID == MO_SOURCE_FEATHER_CIRCLE) ? spots.at(activeSpot).sourcePos : spots.at(activeSpot).targetPos));
            PolarCoord polarPos(delta);
            if (polarPos.angle < 0.) {
                polarPos.angle += 180.;
            }
            if (polarPos.angle < 22.5 || polarPos.angle >= 157.5) {
                return CSMove1DH;
            }
            if (polarPos.angle < 67.5) {
                return CSResizeBottomRight;
            }
            if (polarPos.angle < 112.5) {
                return CSMove1DV;
            }
            return CSResizeBottomLeft;
        }
    }
    return CSCrosshair;
}

bool Spot::mouseOver (int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider && editProvider->getObject() != lastObject) {
        if (lastObject > -1) {
            if (EditSubscriber::mouseOverGeometry.at (lastObject) == &targetMODisc) {
                getVisibleGeometryFromMO (lastObject)->state = Geometry::ACTIVE;
            } else {
                getVisibleGeometryFromMO (lastObject)->state = Geometry::NORMAL;
            }

            sourceIcon.state = Geometry::ACTIVE;
        }

        if (editProvider->getObject() > -1) {
            getVisibleGeometryFromMO (editProvider->getObject())->state = Geometry::PRELIGHT;

            if (editProvider->getObject() >= MO_OBJECT_COUNT) {
                // a Spot is being edited
                int oldActiveSpot = activeSpot;
                activeSpot = editProvider->getObject() - MO_OBJECT_COUNT;

                if (activeSpot != oldActiveSpot) {
                    if (oldActiveSpot > -1) {
                        EditSubscriber::visibleGeometry.at (oldActiveSpot)->state = Geometry::NORMAL;
                        EditSubscriber::mouseOverGeometry.at (oldActiveSpot + MO_OBJECT_COUNT)->state = Geometry::NORMAL;
                    }

                    EditSubscriber::visibleGeometry.at (activeSpot)->state = Geometry::PRELIGHT;
                    EditSubscriber::mouseOverGeometry.at (activeSpot + MO_OBJECT_COUNT)->state = Geometry::PRELIGHT;
                    //printf("ActiveSpot = %d (was %d before)\n", activeSpot, oldActiveSpot);
                }
            }
        }

        lastObject = editProvider->getObject();

        if (lastObject > -1 && EditSubscriber::mouseOverGeometry.at (lastObject) == getActiveSpotIcon()) {
            lastObject = MO_TARGET_DISK;
        }

        updateGeometry();
        return true;
    }

    return false;
}

// Create a new Target and Source point or start the drag of a Target point under the cursor
bool Spot::button1Pressed (int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (editProvider) {
        if (lastObject == -1 && (modifierKey & GDK_CONTROL_MASK)) {
            int imW, imH;
            const auto startPos = editProvider->posImage;
            editProvider->getImageSize(imW, imH);
            if (startPos.x < 0 || startPos.y < 0 || startPos.x > imW || startPos.y > imH) {
                return false; // Outside of image area.
            }
            draggedSide = DraggedSide::SOURCE;
            addNewEntry();
            EditSubscriber::action = EditSubscriber::Action::DRAGGING;
            return true;
        } else if (lastObject > -1) {
            draggedSide = lastObject == MO_TARGET_DISK ? DraggedSide::TARGET : lastObject == MO_SOURCE_DISC ? DraggedSide::SOURCE : DraggedSide::NONE;
            getVisibleGeometryFromMO (lastObject)->state = Geometry::DRAGGED;
            EditSubscriber::action = EditSubscriber::Action::DRAGGING;
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
        EditSubscriber::action = EditSubscriber::Action::NONE;
        return false;
    }

    loGeom->state = Geometry::PRELIGHT;
    EditSubscriber::action = EditSubscriber::Action::NONE;
    draggedSide = DraggedSide::NONE;
    updateGeometry();
    return true;
}

// Delete a point
bool Spot::button2Pressed (int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (!editProvider || lastObject == -1 || activeSpot == -1) {
        return false;
    }

    if (! (modifierKey & (GDK_SHIFT_MASK | GDK_SHIFT_MASK))) {
        EditSubscriber::action = EditSubscriber::Action::PICKING;
    }

    return false;
}

// Create a new Target and Source point or start the drag of a Target point under the cursor
bool Spot::button3Pressed (int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();

    if (!editProvider || lastObject == -1 || activeSpot == -1) {
        return false;
    }

    if ((modifierKey & GDK_CONTROL_MASK) && (EditSubscriber::mouseOverGeometry.at (lastObject) == &targetMODisc || lastObject >= MO_OBJECT_COUNT)) {
        lastObject = MO_SOURCE_DISC;
        sourceIcon.state = Geometry::DRAGGED;
        EditSubscriber::action = EditSubscriber::Action::DRAGGING;
        draggedSide = DraggedSide::SOURCE;
        return true;
    } else if (! (modifierKey & (GDK_SHIFT_MASK | GDK_SHIFT_MASK))) {
        EditSubscriber::action = EditSubscriber::Action::PICKING;
        return true;
    }

    return false;
}

bool Spot::button3Released()
{
    Geometry *loGeom = getVisibleGeometryFromMO (lastObject);

    if (!loGeom) {
        EditSubscriber::action = EditSubscriber::Action::NONE;
        return false;
    }

    lastObject = -1;
    sourceIcon.state = Geometry::ACTIVE;
    draggedSide = DraggedSide::NONE;
    updateGeometry();
    EditSubscriber::action = EditSubscriber::Action::NONE;
    return true;
}

bool Spot::drag1 (int modifierKey)
{
    if (EditSubscriber::action != EditSubscriber::Action::DRAGGING) {
        return false;
    }

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

        EditSubscriber::mouseOverGeometry.at (activeSpot + MO_OBJECT_COUNT)->state = Geometry::DRAGGED;
    } else if (loGeom == &targetMODisc || lastObject >= MO_OBJECT_COUNT) {
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
        rtengine::PolarCoord currPolar (currPos - spots.at (activeSpot).sourcePos);
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

bool Spot::drag3 (int modifierKey)
{
    if (EditSubscriber::action != EditSubscriber::Action::DRAGGING) {
        return false;
    }

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

bool Spot::pick2 (bool picked)
{
    return pick3 (picked);
}

bool Spot::pick3 (bool picked)
{
    EditDataProvider* editProvider = getEditProvider();

    if (!picked) {
        if (editProvider->getObject() != lastObject) {
            return false;
        }
    }

    // Object is picked, we delete it
    deleteSelectedEntry();
    EditSubscriber::action = EditSubscriber::Action::NONE;
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
    listener->unsetTweakOperator(this);
    listener->refreshPreview(EvSpotEnabled); // reprocess the preview w/o creating History entry
}


void Spot::tweakParams(procparams::ProcParams& pparams)
{
    //params->raw.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
    //params->raw.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);

    // -> disabling all transform
    //pparams.coarse = CoarseTransformParams();
    pparams.lensProf = LensProfParams();
    pparams.cacorrection = CACorrParams();
    pparams.distortion = DistortionParams();
    pparams.rotate = RotateParams();
    pparams.perspective = PerspectiveParams();
    pparams.vignetting = VignettingParams();

    // -> disabling standard crop
    pparams.crop.enabled = false;

    // -> disabling time consuming and unnecessary tool
    pparams.sh.enabled = false;
    pparams.blackwhite.enabled = false;
    pparams.dehaze.enabled = false;
    pparams.wavelet.enabled = false;
    pparams.filmSimulation.enabled = false;
    pparams.sharpenEdge.enabled = false;
    pparams.sharpenMicro.enabled = false;
    pparams.sharpening.enabled = false;
    pparams.softlight.enabled = false;
    pparams.gradient.enabled = false;
    pparams.pcvignette.enabled = false;
    pparams.colorappearance.enabled = false;
    pparams.locallab.enabled = false;
   // pparams.toneCurve.hrenabled = false;  // not sure for this one, it could be useful for ExpComp w/o performance penalty
    pparams.toneCurve.histmatching = false;
}
