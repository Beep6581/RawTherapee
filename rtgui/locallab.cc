/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>frame
 *
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
 *  2017 Jacques Desmis <jdesmis@gmail.com>
 *  2019 Pierre Cabrera <pierre.cab@gmail.com>
 */
#include "locallab.h"

#include "options.h"
#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace procparams;

extern Options options;

Locallab::Locallab():
    FoldableToolPanel(this, "locallab", M("TP_LOCALLAB_LABEL"), false, true),

    // Spot control panel widget
    expsettings(Gtk::manage(new ControlSpotPanel()))
{
    const bool showtooltip = options.showtooltip;

    // Create panel widget to receive Locallab GUI elements
    ToolVBox* const panel = Gtk::manage(new ToolVBox());

    // Add spot control panel to panel widget
    expsettings->setLevel(2); // TODO Move this to controlspotpanel.cc
    panel->pack_start(*expsettings->getExpander(), false, false);

    // Create Locallab tools
    expcolor = Gtk::manage(new LocallabColor());
    expexpose = Gtk::manage(new LocallabExposure());
    expshadhigh = Gtk::manage(new LocallabShadow());
    expvibrance = Gtk::manage(new LocallabVibrance());
    expsoft = Gtk::manage(new LocallabSoft());
    expblur = Gtk::manage(new LocallabBlur());
    exptonemap = Gtk::manage(new LocallabTone());
    expreti = Gtk::manage(new LocallabRetinex());
    expsharp = Gtk::manage(new LocallabSharp());
    expcontrast = Gtk::manage(new LocallabContrast());
    expcbdl = Gtk::manage(new LocallabCBDL());
    expdenoi = Gtk::manage(new LocallabDenoise());

    // Add Locallab tools to panel widget
    addTool(panel, expcolor);
    addTool(panel, expexpose);
    addTool(panel, expshadhigh);
    addTool(panel, expvibrance);
    addTool(panel, expsoft);
    addTool(panel, expblur);
    addTool(panel, exptonemap);
    addTool(panel, expreti);
    addTool(panel, expsharp);
    addTool(panel, expcontrast);
    addTool(panel, expcbdl);
    addTool(panel, expdenoi);

    // Add panel widget to Locallab GUI
    pack_start(*panel);

    // Show all widgets
    show_all();

    // By default, if no photo is loaded, all Locallab tools are removed and it's not possible to add them
    // (to be necessary called after "show_all" function)
    setParamEditable(false);
}

void Locallab::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // printf("Locallab read\n");

    // Disable all listeners
    disableListener();

    // Update Locallab activation state
    setEnabled(pp->locallab.enabled);

    // Transmit Locallab activation state to Locallab tools
    for (auto tool : locallabTools) {
        tool->isLocallabActivated(exp->getEnabled());
    }

    // TODO Manage it with read function in controlspotpanel.cc
    // Delete all existent spots
    std::vector<int>* const list = expsettings->getSpotIdList();

    for (size_t i = 0; i < list->size(); i++) {
        expsettings->deleteControlSpot(list->at(i));
    }

    // TODO Manage it with read function in controlspotpanel.cc
    // Add existent spots based on pp
    ControlSpotPanel::SpotRow* const r = new ControlSpotPanel::SpotRow();

    for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
        r->id = pp->locallab.spots.at(i).id;
        r->name = pp->locallab.spots.at(i).name;
        r->isvisible = pp->locallab.spots.at(i).isvisible;

        if (pp->locallab.spots.at(i).shape == "ELI") {
            r->shape = 0;
        } else {
            r->shape = 1;
        }

        if (pp->locallab.spots.at(i).spotMethod == "norm") {
            r->spotMethod = 0;
        } else {
            r->spotMethod = 1;
        }

        r->sensiexclu = pp->locallab.spots.at(i).sensiexclu;
        r->structexclu = pp->locallab.spots.at(i).structexclu;
        r->struc = pp->locallab.spots.at(i).struc;

        if (pp->locallab.spots.at(i).shapeMethod == "IND") {
            r->shapeMethod = 0;
        } else if (pp->locallab.spots.at(i).shapeMethod == "SYM") {
            r->shapeMethod = 1;
        } else if (pp->locallab.spots.at(i).shapeMethod == "INDSL") {
            r->shapeMethod = 2;
        } else {
            r->shapeMethod = 3;
        }

        r->locX = pp->locallab.spots.at(i).locX;
        r->locXL = pp->locallab.spots.at(i).locXL;
        r->locY = pp->locallab.spots.at(i).locY;
        r->locYT = pp->locallab.spots.at(i).locYT;
        r->centerX = pp->locallab.spots.at(i).centerX;
        r->centerY = pp->locallab.spots.at(i).centerY;
        r->circrad = pp->locallab.spots.at(i).circrad;

        if (pp->locallab.spots.at(i).qualityMethod == "enh") {
            r->qualityMethod = 0;
        } else {
            r->qualityMethod = 1;
        }

        r->transit = pp->locallab.spots.at(i).transit;
        r->thresh = pp->locallab.spots.at(i).thresh;
        r->iter = pp->locallab.spots.at(i).iter;
        r->balan = pp->locallab.spots.at(i).balan;
        r->transitweak = pp->locallab.spots.at(i).transitweak;
        r->transitgrad = pp->locallab.spots.at(i).transitgrad;
        r->avoid = pp->locallab.spots.at(i).avoid;

        expsettings->addControlSpot(r);
    }

    // Select active spot
    if (pp->locallab.nbspot > 0) {
        expsettings->setSelectedSpot(pp->locallab.spots.at(pp->locallab.selspot).id);
    }

    // Update each Locallab tools GUI
    for (auto tool : locallabTools) {
        tool->read(pp, pedited);
    }

    // Specific case: if there is no spot, GUI isn't anymore editable (i.e. Locallab tool cannot be managed)
    if (pp->locallab.nbspot > 0) {
        setParamEditable(true);
    } else {
        setParamEditable(false);
    }

    // Enable all listeners
    enableListener();

    // Open/re-open all Locallab tools expanders
    openAllTools();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void Locallab::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    // Update Locallab activation state
    pp->locallab.enabled = getEnabled();

    // Transmit Locallab activation state to Locallab tools (in case of updated)
    for (auto tool : locallabTools) {
        tool->isLocallabActivated(exp->getEnabled());
    }

    const int spotPanelEvent = expsettings->getEventType();
    int spotId;
    ControlSpotPanel::SpotRow* r;
    LocallabParams::LocallabSpot* newSpot;

    int imW, imH; // Size of image
    int prW, prH; // Size of preview area
    int prX, prY; // Coord of preview area center
    EditDataProvider* const provider = expsettings->getEditProvider();

    switch (spotPanelEvent) {
        case (ControlSpotPanel::SpotCreation): // Spot creation event
            // Spot creation (default initialization)
            newSpot = new LocallabParams::LocallabSpot();
            spotId = expsettings->getNewId();
            r = new ControlSpotPanel::SpotRow();
            r->id = newSpot->id = spotId;
            r->name = newSpot->name = M("TP_LOCALLAB_SPOTNAME") + std::to_string(spotId);
            r->isvisible = newSpot->isvisible;

            if (newSpot->shape == "ELI") {
                r->shape = 0;
            } else {
                r->shape = 1;
            }

            if (newSpot->spotMethod == "norm") {
                r->spotMethod = 0;
            } else {
                r->spotMethod = 1;
            }

            r->sensiexclu = newSpot->sensiexclu;
            r->structexclu = newSpot->structexclu;
            r->struc = newSpot->struc;

            if (newSpot->shapeMethod == "IND") {
                r->shapeMethod = 0;
            } else if (newSpot->shapeMethod == "SYM") {
                r->shapeMethod = 1;
            } else if (newSpot->shapeMethod == "INDSL") {
                r->shapeMethod = 2;
            } else {
                r->shapeMethod = 3;
            }

            // Calculate spot size and center position according to preview area
            if (provider && !batchMode) {
                provider->getImageSize(imW, imH);
                provider->getPreviewCenterPos(prX, prY);
                provider->getPreviewSize(prW, prH);

                if (imW && imH) { // Image loaded
                    // Spot center position computation
                    newSpot->centerX = rtengine::LIM(int(int((double)prX - (double)imW / 2.) * 2000. / (double)imW), -1000, 1000);
                    newSpot->centerY = rtengine::LIM(int(int((double)prY - (double)imH / 2.) * 2000. / (double)imH), -1000, 1000);
                    // Ellipse/rectangle size computation
                    newSpot->locX = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->locX);
                    newSpot->locXL = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->locXL);
                    newSpot->locY = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->locY);
                    newSpot->locYT = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->locYT);
                }
            }

            r->locX = newSpot->locX;
            r->locXL = newSpot->locXL;
            r->locY = newSpot->locY;
            r->locYT = newSpot->locYT;
            r->centerX = newSpot->centerX;
            r->centerY = newSpot->centerY;

            r->circrad = newSpot->circrad;

            if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 0;
            } else {
                r->qualityMethod = 1;
            }

            r->transit = newSpot->transit;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            r->balan = newSpot->balan;
            r->transitweak = newSpot->transitweak;
            r->transitgrad = newSpot->transitgrad;
            r->avoid = newSpot->avoid;
            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.nbspot++;
            pp->locallab.selspot = pp->locallab.nbspot - 1;
            pp->locallab.spots.push_back(*newSpot);

            // New created spot selection
            expsettings->setSelectedSpot(spotId);

            // Update Locallab tools GUI with new created spot
            disableListener();

            for (auto tool : locallabTools) {
                tool->read(pp, pedited);
            }

            enableListener();

            if (pp->locallab.nbspot == 1) {
                setParamEditable(true);
            }

            // Update default values according to selected spot
            setDefaults(pp, pedited);

            // Open/re-open all Locallab tools expanders
            openAllTools();

            // Note: No need to manage pedited as batch mode is deactivated for Locallab

            break;

        case (ControlSpotPanel::SpotDeletion): // Spot deletion event
            // Get deleted spot index in ProcParams and update it
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
                if (pp->locallab.spots.at(i).id == spotId) {
                    // ProcParams update
                    pp->locallab.nbspot--;
                    pp->locallab.spots.erase(pp->locallab.spots.begin() + i);
                    expsettings->deleteControlSpot(spotId);

                    // Select the first remaining spot before deleted one
                    if (pp->locallab.nbspot > 0) {
                        for (int j = i - 1; j >= 0; j--) { // procparams spots uses zero-based index whereas spot ids use one-based index
                            if (expsettings->setSelectedSpot(j + 1)) { // True if an existing spot has been selected on controlspotpanel
                                pp->locallab.selspot = j;

                                break;
                            }
                        }
                    } else {
                        // Reset selspot
                        pp->locallab.selspot = 0;
                    }

                    // Update Locallab tools GUI with selected spot
                    disableListener();

                    for (auto tool : locallabTools) {
                        tool->read(pp, pedited);
                    }

                    enableListener();

                    if (pp->locallab.nbspot == 0) {
                        setParamEditable(false);
                    }

                    // Update default values according to selected spot
                    setDefaults(pp, pedited);

                    // Open/re-open all Locallab tools expanders
                    openAllTools();

                    // Note: No need to manage pedited as batch mode is deactivated for Locallab

                    break;
                }
            }

            break;

        case (ControlSpotPanel::SpotSelection):  // Spot selection event
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
                if (pp->locallab.spots.at(i).id == spotId) {
                    pp->locallab.selspot = i;
                    break;
                }
            }

            // Update control spots and Locallab tools GUI with selected spot
            expsettings->setSelectedSpot(spotId);
            disableListener();

            for (auto tool : locallabTools) {
                tool->read(pp, pedited);
            }

            enableListener();

            // Update locallab tools mask background
            if (pp->locallab.selspot < (int)maskBackRef.size()) {
                const double huer = maskBackRef.at(pp->locallab.selspot).huer;
                const double lumar = maskBackRef.at(pp->locallab.selspot).lumar;
                const double chromar = maskBackRef.at(pp->locallab.selspot).chromar;

                for (auto tool : locallabTools) {
                    tool->refChanged(huer, lumar, chromar);
                }
            }

            // Update default values according to selected spot
            setDefaults(pp, pedited);

            // Open/re-open all Locallab tools expanders
            openAllTools();

            // Note: No need to manage pedited as batch mode is deactivated for Locallab

            break;

        case (ControlSpotPanel::SpotDuplication): // Spot duplication event
            newSpot = nullptr;
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot && i < (int)pp->locallab.spots.size(); i++) {
                if (pp->locallab.spots.at(i).id == spotId) {
                    newSpot = new LocallabParams::LocallabSpot(pp->locallab.spots.at(i));
                    break;
                }
            }

            if (!newSpot) {
                break;
            }

            // Spot creation (initialization at currently selected spot)
            spotId = expsettings->getNewId();
            r = new ControlSpotPanel::SpotRow();
            r->id = newSpot->id = spotId;
            r->name = newSpot->name = newSpot->name + " - " + M("TP_LOCALLAB_DUPLSPOTNAME");
            r->isvisible = newSpot->isvisible;

            if (newSpot->shape == "ELI") {
                r->shape = 0;
            } else {
                r->shape = 1;
            }

            if (newSpot->spotMethod == "norm") {
                r->spotMethod = 0;
            } else {
                r->spotMethod = 1;
            }

            r->sensiexclu = newSpot->sensiexclu;
            r->structexclu = newSpot->structexclu;
            r->struc = newSpot->struc;

            if (newSpot->shapeMethod == "IND") {
                r->shapeMethod = 0;
            } else if (newSpot->shapeMethod == "SYM") {
                r->shapeMethod = 1;
            } else if (newSpot->shapeMethod == "INDSL") {
                r->shapeMethod = 2;
            } else {
                r->shapeMethod = 3;
            }

            // Calculate spot size and center position according to preview area
            if (provider && !batchMode) {
                provider->getImageSize(imW, imH);
                provider->getPreviewCenterPos(prX, prY);
                provider->getPreviewSize(prW, prH);

                if (imW && imH) { // Image loaded
                    // Spot center position computation
                    newSpot->centerX = rtengine::LIM(int(int((double)prX - (double)imW / 2.) * 2000. / (double)imW), -1000, 1000);
                    newSpot->centerY = rtengine::LIM(int(int((double)prY - (double)imH / 2.) * 2000. / (double)imH), -1000, 1000);
                    // Ellipse/rectangle size computation
                    newSpot->locX = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->locX);
                    newSpot->locXL = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->locXL);
                    newSpot->locY = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->locY);
                    newSpot->locYT = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->locYT);
                }
            }

            r->locX = newSpot->locX;
            r->locXL = newSpot->locXL;
            r->locY = newSpot->locY;
            r->locYT = newSpot->locYT;
            r->centerX = newSpot->centerX;
            r->centerY = newSpot->centerY;

            r->circrad = newSpot->circrad;

            if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 0;
            } else {
                r->qualityMethod = 1;
            }

            r->transit = newSpot->transit;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            r->balan = newSpot->balan;
            r->transitweak = newSpot->transitweak;
            r->transitgrad = newSpot->transitgrad;
            r->avoid = newSpot->avoid;
            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.nbspot++;
            pp->locallab.selspot = pp->locallab.nbspot - 1;
            pp->locallab.spots.push_back(*newSpot);

            // New created spot selection
            expsettings->setSelectedSpot(spotId);

            // Update Locallab tools GUI with new created spot
            disableListener();

            for (auto tool : locallabTools) {
                tool->read(pp, pedited);
            }

            enableListener();

            // Update default values according to selected spot
            setDefaults(pp, pedited);

            // Open/re-open all Locallab tools expanders
            openAllTools();

            // Note: No need to manage pedited as batch mode is deactivated for Locallab

            break;

        case (ControlSpotPanel::SpotAllVisibilityChanged): // Event when updating visibility of all spots
            r = expsettings->getSpot(expsettings->getSelectedSpot());

            // ProcParams update
            for (size_t i = 0; i < pp->locallab.spots.size(); i++) {
                pp->locallab.spots.at(i).isvisible = r->isvisible;
            }

            // Note: No need to manage pedited as batch mode is deactivated for Locallab

            break;

        default: // Spot or locallab GUI updated
            if (pp->locallab.nbspot > 0) {
                r = expsettings->getSpot(expsettings->getSelectedSpot());

                // ProcParams update
                if (pp->locallab.selspot < (int)pp->locallab.spots.size()) {
                    // Control spot settings
                    pp->locallab.spots.at(pp->locallab.selspot).name = r->name;
                    pp->locallab.spots.at(pp->locallab.selspot).isvisible = r->isvisible;

                    if (r->shape == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).shape = "ELI";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).shape = "RECT";
                    }

                    if (r->spotMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "norm";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "exc";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).sensiexclu = r->sensiexclu;
                    pp->locallab.spots.at(pp->locallab.selspot).structexclu = r->structexclu;
                    pp->locallab.spots.at(pp->locallab.selspot).struc = r->struc;

                    if (r->shapeMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "IND";
                    } else if (r->shapeMethod == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "SYM";
                    } else if (r->shapeMethod == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "INDSL";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "SYMSL";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).locX = r->locX;
                    pp->locallab.spots.at(pp->locallab.selspot).locXL = r->locXL;
                    pp->locallab.spots.at(pp->locallab.selspot).locY = r->locY;
                    pp->locallab.spots.at(pp->locallab.selspot).locYT = r->locYT;
                    pp->locallab.spots.at(pp->locallab.selspot).centerX = r->centerX;
                    pp->locallab.spots.at(pp->locallab.selspot).centerY = r->centerY;
                    pp->locallab.spots.at(pp->locallab.selspot).circrad = r->circrad;

                    if (r->qualityMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enh";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enhden";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).transit = r->transit;
                    pp->locallab.spots.at(pp->locallab.selspot).thresh = r->thresh;
                    pp->locallab.spots.at(pp->locallab.selspot).iter = r->iter;
                    pp->locallab.spots.at(pp->locallab.selspot).balan = r->balan;
                    pp->locallab.spots.at(pp->locallab.selspot).transitweak = r->transitweak;
                    pp->locallab.spots.at(pp->locallab.selspot).transitgrad = r->transitgrad;
                    pp->locallab.spots.at(pp->locallab.selspot).avoid = r->avoid;
                }

                for (auto tool : locallabTools) {
                    tool->write(pp, pedited);
                }

                // Note: No need to manage pedited as batch mode is deactivated for Locallab
            }
    }
}

/*
 * Note:
 * By default, this function is called when a new image/profile is loaded (after read function). In this case,
 * if there is at least one spot, default values are set to selected spot ones.
 * To keep having default values according to selected spot, this function shall also be called in the following
 * situations (after having called write function for controlspotpanel):
 * - After spot creation
 * - After spot deletion
 * - After spot selection
 * - After spot duplication
 */
void Locallab::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    // Set default values in spot panel control
    expsettings->setDefaults(defParams, pedited);

    // Set defaut values in Locallab tools
    for (auto tool : locallabTools) {
        tool->setDefaults(defParams, pedited);
    }
}

void Locallab::setListener(ToolPanelListener* tpl)
{
    this->listener = tpl;

    // Set listener for spot control panel
    expsettings->setListener(tpl);

    // Set listener for locallab tools
    for (auto tool : locallabTools) {
        tool->setListener(tpl);
    }
}

void Locallab::refChanged(const std::vector<locallabRef> &ref, int selspot)
{
    // Saving transmitted mask background data
    maskBackRef = ref;

    // Update locallab tools mask background
    if (selspot < (int)maskBackRef.size()) {
        const double huer = maskBackRef.at(selspot).huer;
        const double lumar = maskBackRef.at(selspot).lumar;
        const double chromar = maskBackRef.at(selspot).chromar;

        for (auto tool : locallabTools) {
            tool->refChanged(huer, lumar, chromar);
        }
    }
}

void Locallab::resetMaskVisibility()
{
    // Indicate to spot control panel that no more mask preview is active
    expsettings->setMaskPrevActive(false);

    // Reset mask preview for all Locallab tools
    for (auto tool : locallabTools) {
        tool->resetMaskView();
    }
}

Locallab::llMaskVisibility Locallab::getMaskVisibility() const
{
    // Get mask preview from Locallab tools
    int colorMask, expMask, SHMask, softMask, tmMask, retiMask, cbMask;

    for (auto tool : locallabTools) {
        tool->getMaskView(colorMask, expMask, SHMask, softMask, tmMask, retiMask, cbMask);
    }

    // Indicate to spot control panel if one mask preview is active
    const bool isMaskActive = (colorMask == 0) || (expMask == 0) || (SHMask == 0) ||
                              (softMask == 0) || (tmMask == 0) || (retiMask == 0) || (cbMask == 0);
    expsettings->setMaskPrevActive(isMaskActive);

    return {colorMask, expMask, SHMask, softMask, tmMask, retiMask, cbMask};
}

void Locallab::setEditProvider(EditDataProvider * provider)
{
    expsettings->setEditProvider(provider);
}

void Locallab::subscribe()
{
    expsettings->subscribe();
}

void Locallab::unsubscribe()
{
    expsettings->unsubscribe();
}

void Locallab::enabledChanged()
{
    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::autoOpenCurve()
{
    // TODO Actually autoOpenCurve only considers linearity state of selected spot curve
}

void Locallab::foldAllButOne(LocallabTool* except)
{
    for (auto tool : locallabTools) {
        if (tool != except) {
            // All other tool expanders are fold
            tool->setExpanded(false);
        } else {
            // If fold, selected tool expander is unfold
            if (!tool->getExpanded()) {
                tool->setExpanded(true);
            }
        }
    }
}

void Locallab::addTool(Gtk::Box* where, LocallabTool* tool)
{
    where->pack_start(*tool->getExpander(), false, false);
    locallabTools.push_back(tool);
    tool->setLocallabToolListener(this);
}

void Locallab::openAllTools()
{
    for (auto tool : locallabTools) {
        tool->setExpanded(true);
    }
}

void Locallab::setParamEditable(bool cond)
{
    // Update params editable state for controlspotpanel
    expsettings->setParamEditable(cond); // TODO Move this code to controlspotpanel.cc when there is zero spot

    // Enable/disable possibility to add Locallab tool
    // TODO To implement

    // Remove all Locallab tool only if cond is false
    // TODO To be managed in locallabtools first
    /*
    if (!cond) {
        for (auto tool : locallabTools) {
            tool->addLocallabTool(false);
        }
    }
    */
}

void Locallab::resetOtherMaskView(LocallabTool* current)
{
    // Reset mask view GUI for all other Locallab tools except current
    for (auto tool : locallabTools) {
        if (tool != current) {
            tool->resetMaskView();
        }
    }
}
