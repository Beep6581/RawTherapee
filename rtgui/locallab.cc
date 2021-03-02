/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as publishfed by
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

/* ==== LocallabToolList ==== */
LocallabToolList::LocallabToolList():
    // Tool list GUI elements
    list(Gtk::manage(new MyComboBox())),
    listTreeModel(Gtk::ListStore::create(toolRow)),

    // Tool list listener
    listListener(nullptr)
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    list->set_model(listTreeModel);
    list->pack_start(toolRow.name);
    listConn = list->signal_changed().connect(sigc::mem_fun(*this, &LocallabToolList::toolRowSelected));
    list->set_tooltip_text(M("TP_LOCALLAB_LIST_TOOLTIP"));
    // Append title row to list
    // Important: Title row shall always be the first one
    const auto titleRow = *(listTreeModel->append());
    titleRow[toolRow.id] = 0;
    titleRow[toolRow.name] = M("TP_LOCALLAB_LIST_NAME");
    listConn.block(true);
    list->set_active(titleRow);
    listConn.block(false);

    // Add ComboBox to LocallabToolList widget
    add(*list);
}

void LocallabToolList::addToolRow(const Glib::ustring &toolname, const int id)
{
    // Disable event management
    listConn.block(true);

    // Add tool name according to id
    Gtk::TreeIter insertAfter;

    for (auto &r : listTreeModel->children()) {
        if (r[toolRow.id] < id) {
            insertAfter = *r; // Tool name shall be added at least after this row
        } else {
            break; // Tool name shall be added before this row
        }
    }

    // Note: There is always at list one row (i.e. title one)

    const auto newRow = *(listTreeModel->insert_after(insertAfter));
    newRow[toolRow.id] = id;
    newRow[toolRow.name] = toolname;

    // Select title row (i.e. always first row)
    list->set_active(0);

    // Enable event management
    listConn.block(false);
}

void LocallabToolList::removeToolRow(const Glib::ustring &toolname)
{
    // Disable event management
    listConn.block(true);

    // Remove tool name row
    for (auto &r : listTreeModel->children()) {
        if (r[toolRow.name] == toolname) {
            listTreeModel->erase(*r);
            break;
        }
    }

    // Select title row (i.e. always first row)
    list->set_active(0);

    // Enable event management
    listConn.block(false);
}

void LocallabToolList::removeAllTool()
{
    // Disable event management
    listConn.block(true);

    // Remove all tools
    listTreeModel->clear();

    // Add title row again
    const auto titleRow = *(listTreeModel->append());
    titleRow[toolRow.id] = 0;
    titleRow[toolRow.name] = M("TP_LOCALLAB_LIST_NAME");

    // Select title row (i.e. always first row)
    list->set_active(0);

    // Enable event management
    listConn.block(false);
}

void LocallabToolList::toolRowSelected()
{
    // Get selected tool name
    const auto selRow = *(list->get_active());
    const Glib::ustring toolname = selRow[toolRow.name];

    // Remove selected tool name for ComboBox
    removeToolRow(toolname);

    // Warm tool list listener
    if (listListener) {
        listListener->locallabToolToAdd(toolname);
    }
}

/* ==== Locallab ==== */
Locallab::Locallab():
    FoldableToolPanel(this, "locallab", M("TP_LOCALLAB_LABEL"), false, true),

    // Spot control panel widget
    expsettings(Gtk::manage(new ControlSpotPanel())),

    // Tool list widget
    toollist(Gtk::manage(new LocallabToolList())),

    // Create Locallab tools
    expcolor(Gtk::manage(new LocallabColor())),
    expexpose(Gtk::manage(new LocallabExposure())),
    expshadhigh(Gtk::manage(new LocallabShadow())),
    expvibrance(Gtk::manage(new LocallabVibrance())),
    expsoft(Gtk::manage(new LocallabSoft())),
    expblur(Gtk::manage(new LocallabBlur())),
    exptonemap(Gtk::manage(new LocallabTone())),
    expreti(Gtk::manage(new LocallabRetinex())),
    expsharp(Gtk::manage(new LocallabSharp())),
    expcontrast(Gtk::manage(new LocallabContrast())),
    expcbdl(Gtk::manage(new LocallabCBDL())),
    explog(Gtk::manage(new LocallabLog())),
    expmask(Gtk::manage(new LocallabMask())),

    // Other widgets
    resetshowButton(Gtk::manage(new Gtk::Button(M("TP_LOCALLAB_RESETSHOW"))))
{
    set_orientation(Gtk::ORIENTATION_VERTICAL);
    
    // Create panel widget to receive Locallab GUI elements
    ToolVBox* const panel = Gtk::manage(new ToolVBox());
    panel->set_spacing(2);

    // Add spot control panel to panel widget
    expsettings->setControlPanelListener(this);
    expsettings->setLevel(2);
    panel->pack_start(*expsettings->getExpander(), false, false);

    // Add separator
    Gtk::Separator* const separator = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
    panel->pack_start(*separator, false, false);

    // Add tool list widget
    toollist->setLocallabToolListListener(this);
    panel->pack_start(*toollist, false, false);

    // Add Locallab tools to panel widget
    ToolVBox* const toolpanel = Gtk::manage(new ToolVBox());
    toolpanel->set_name("LocallabToolPanel");
    addTool(toolpanel, expcolor);
    addTool(toolpanel, expshadhigh);
    addTool(toolpanel, expvibrance);
    addTool(toolpanel, explog);
    addTool(toolpanel, expexpose);
    addTool(toolpanel, expmask);
    addTool(toolpanel, expsoft);
    addTool(toolpanel, expblur);
    addTool(toolpanel, exptonemap);
    addTool(toolpanel, expreti);
    addTool(toolpanel, expsharp);
    addTool(toolpanel, expcontrast);
    addTool(toolpanel, expcbdl);
    panel->pack_start(*toolpanel, false, false);

    // Add separator
 //   Gtk::Separator* const separator2 = Gtk::manage(new Gtk::Separator(Gtk::ORIENTATION_HORIZONTAL));
 //   panel->pack_start(*separator2, false, false);

    // Add mask reset button to panel widget
    resetshowButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::resetshowPressed));
   // panel->pack_start(*resetshowButton);

    // Add panel widget to Locallab GUI
    pack_start(*panel);

    // Show all widgets
    show_all();

    // Update Locallab tools advice tooltips visibility based on saved option
    for (auto tool : locallabTools) {
        tool->updateAdviceTooltips(options.showtooltip);
    }

    // By default, if no photo is loaded, all Locallab tools are removed and it's not possible to add them
    // (to be necessary called after "show_all" function)
    setParamEditable(false);
}

void Locallab::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
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
    const int spotNb = expsettings->getSpotNumber();

    for (int i = spotNb - 1; i >= 0; i--) {
        expsettings->deleteControlSpot(i);
    }

    // TODO Manage it with read function in controlspotpanel.cc
    // Add existent spots based on pp
    ControlSpotPanel::SpotRow* const r = new ControlSpotPanel::SpotRow();

    for (int i = 0; i < (int)pp->locallab.spots.size(); i++) {
        r->name = pp->locallab.spots.at(i).name;
        r->isvisible = pp->locallab.spots.at(i).isvisible;

        if (pp->locallab.spots.at(i).shape == "ELI") {
            r->shape = 0;
        } else {
            r->shape = 1;
        }

        if (pp->locallab.spots.at(i).prevMethod == "hide") {
            r->prevMethod = 0;
        } else {
            r->prevMethod = 1;
        }

        if (pp->locallab.spots.at(i).spotMethod == "norm") {
            r->spotMethod = 0;
        } else if(pp->locallab.spots.at(i).spotMethod == "exc"){
            r->spotMethod = 1;
        } else if (pp->locallab.spots.at(i).spotMethod == "full"){
            r->spotMethod = 2;
        }

        r->sensiexclu = pp->locallab.spots.at(i).sensiexclu;
        r->structexclu = pp->locallab.spots.at(i).structexclu;

        if (pp->locallab.spots.at(i).shapeMethod == "IND") {
            r->shapeMethod = 0;
        } else if (pp->locallab.spots.at(i).shapeMethod == "SYM") {
            r->shapeMethod = 1;
        } else if (pp->locallab.spots.at(i).shapeMethod == "INDSL") {
            r->shapeMethod = 2;
        } else {
            r->shapeMethod = 3;
        }

        r->locX = pp->locallab.spots.at(i).loc.at(0);
        r->locXL = pp->locallab.spots.at(i).loc.at(1);
        r->locY = pp->locallab.spots.at(i).loc.at(2);
        r->locYT = pp->locallab.spots.at(i).loc.at(3);
        r->centerX = pp->locallab.spots.at(i).centerX;
        r->centerY = pp->locallab.spots.at(i).centerY;
        r->circrad = pp->locallab.spots.at(i).circrad;

        if (pp->locallab.spots.at(i).qualityMethod == "enh") {
            r->qualityMethod = 0;
        } else {
            r->qualityMethod = 1;
        }

        r->transit = pp->locallab.spots.at(i).transit;
        r->transitweak = pp->locallab.spots.at(i).transitweak;
        r->transitgrad = pp->locallab.spots.at(i).transitgrad;
        r->feather = pp->locallab.spots.at(i).feather;
        r->struc = pp->locallab.spots.at(i).struc;
        r->thresh = pp->locallab.spots.at(i).thresh;
        r->iter = pp->locallab.spots.at(i).iter;
        r->balan = pp->locallab.spots.at(i).balan;
        r->balanh = pp->locallab.spots.at(i).balanh;
        r->colorde = pp->locallab.spots.at(i).colorde;
        r->colorscope = pp->locallab.spots.at(i).colorscope;
        r->avoidrad = pp->locallab.spots.at(i).avoidrad;
        r->hishow = pp->locallab.spots.at(i).hishow;
        r->activ = pp->locallab.spots.at(i).activ;
        r->avoid = pp->locallab.spots.at(i).avoid;
        r->avoidmun = pp->locallab.spots.at(i).avoidmun;
        r->blwh = pp->locallab.spots.at(i).blwh;
        r->recurs = pp->locallab.spots.at(i).recurs;
        r->laplac = true; //pp->locallab.spots.at(i).laplac;
        r->deltae = pp->locallab.spots.at(i).deltae;
        r->scopemask = pp->locallab.spots.at(i).scopemask;
        r->shortc = pp->locallab.spots.at(i).shortc;
        r->lumask = pp->locallab.spots.at(i).lumask;
        r->savrest = pp->locallab.spots.at(i).savrest;

        if (pp->locallab.spots.at(i).complexMethod == "sim") {
            r->complexMethod = 0;
        } else  if (pp->locallab.spots.at(i).complexMethod == "mod") {
            r->complexMethod = 1;
        } else  if (pp->locallab.spots.at(i).complexMethod == "all") {
            r->complexMethod = 2;
        }

        if (pp->locallab.spots.at(i).wavMethod == "D2") {
            r->wavMethod = 0;
        } else  if (pp->locallab.spots.at(i).wavMethod == "D4") {
            r->wavMethod = 1;
        } else  if (pp->locallab.spots.at(i).wavMethod == "D6") {
            r->wavMethod = 2;
        } else  if (pp->locallab.spots.at(i).wavMethod == "D10") {
            r->wavMethod = 3;
        } else  if (pp->locallab.spots.at(i).wavMethod == "D14") {
            r->wavMethod = 4;
        }

        expsettings->addControlSpot(r);
    }

    // Select active spot
    if (pp->locallab.spots.size() > 0) {
        expsettings->setSelectedSpot(pp->locallab.selspot);
    }

    // Update each Locallab tools GUI
    for (auto tool : locallabTools) {
        tool->read(pp, pedited);
    }

    // Update tool list widget
    int toolNb = 0;
    toollist->removeAllTool(); // Reset Locallab list firstly

    for (auto tool : locallabTools) {
        toolNb++;

        if (!tool->isLocallabToolAdded()) {
            toollist->addToolRow(tool->getToolName(), toolNb);
        }
    }

    // Specific case: if there is no spot, GUI isn't anymore editable (i.e. Locallab tool cannot be managed)
    if (pp->locallab.spots.size() > 0) {
        setParamEditable(true);
    } else {
        setParamEditable(false);
    }

    // Enable all listeners
    enableListener();

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
    int spotIndex;
    ControlSpotPanel::SpotRow* r;
    rtengine::procparams::LocallabParams::LocallabSpot* newSpot;

    int imW, imH; // Size of image
    int prW, prH; // Size of preview area
    int prX, prY; // Coord of preview area center
    EditDataProvider* const provider = expsettings->getEditProvider();

    int toolNb;

    switch (spotPanelEvent) {
        case (ControlSpotPanel::SpotCreation): // Spot creation event
            // Spot creation (default initialization)
            newSpot = new LocallabParams::LocallabSpot();
            r = new ControlSpotPanel::SpotRow();
            r->name = newSpot->name = M("TP_LOCALLAB_SPOTNAME");
            r->isvisible = newSpot->isvisible;

            if (newSpot->shape == "ELI") {
                r->shape = 0;
            } else {
                r->shape = 1;
            }

            if (newSpot->prevMethod == "hide") {
                r->prevMethod = 0;
            } else {
                r->prevMethod = 1;
            }


            if (newSpot->spotMethod == "norm") {
                r->spotMethod = 0;
            } else if(newSpot->spotMethod == "exc") {
                r->spotMethod = 1;
            } else if(newSpot->spotMethod == "full") {
                r->spotMethod = 2;
            }

            r->sensiexclu = newSpot->sensiexclu;
            r->structexclu = newSpot->structexclu;

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
                    newSpot->loc.at(0) = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->loc.at(0));
                    newSpot->loc.at(1) = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->loc.at(1));
                    newSpot->loc.at(2) = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->loc.at(2));
                    newSpot->loc.at(3) = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->loc.at(3));
                }
            }

            r->locX = newSpot->loc.at(0);
            r->locXL = newSpot->loc.at(1);
            r->locY = newSpot->loc.at(2);
            r->locYT = newSpot->loc.at(3);
            r->centerX = newSpot->centerX;
            r->centerY = newSpot->centerY;

            r->circrad = newSpot->circrad;

            if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 0;
            } else {
                r->qualityMethod = 1;
            }

            r->transit = newSpot->transit;
            r->transitweak = newSpot->transitweak;
            r->transitgrad = newSpot->transitgrad;
            r->feather = newSpot->feather;
            r->struc = newSpot->struc;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            r->balan = newSpot->balan;
            r->balanh = newSpot->balanh;
            r->colorde = newSpot->colorde;
            r->colorscope = newSpot->colorscope;
            r->avoidrad = newSpot->avoidrad;
            r->hishow = newSpot->hishow;
            r->activ = newSpot->activ;
            r->avoid = newSpot->avoid;
            r->avoidmun = newSpot->avoidmun;
            r->blwh = newSpot->blwh;
            r->recurs = newSpot->recurs;
            r->laplac = newSpot->laplac;
            r->deltae = newSpot->deltae;
            r->scopemask = newSpot->scopemask;
            r->shortc = newSpot->shortc;
            r->lumask = newSpot->lumask;
            r->savrest = newSpot->savrest;

            if (newSpot->complexMethod == "sim") {
                r->complexMethod = 0;
            } else  if (newSpot->complexMethod == "mod") {
                r->complexMethod = 1;
            } else  if (newSpot->complexMethod == "all") {
                r->complexMethod = 2;
            }

            if (newSpot->wavMethod == "D2") {
                r->wavMethod = 0;
            } else if (newSpot->wavMethod == "D4") {
                r->wavMethod = 1;
            } else if (newSpot->wavMethod == "D6") {
                r->wavMethod = 2;
            } else if (newSpot->wavMethod == "D10") {
                r->wavMethod = 3;
            } else if (newSpot->wavMethod == "D14") {
                r->wavMethod = 4;
            }

            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.spots.push_back(*newSpot);
            pp->locallab.selspot = pp->locallab.spots.size() - 1;

            // New created spot selection
            expsettings->setSelectedSpot(pp->locallab.selspot);

            // Update Locallab tools GUI with new created spot
            disableListener();

            for (auto tool : locallabTools) {
                tool->read(pp, pedited);
            }

            enableListener();

            // Update tool list widget
            toolNb = 0;
            toollist->removeAllTool(); // Reset Locallab list firstly

            for (auto tool : locallabTools) {
                toolNb++;

                if (!tool->isLocallabToolAdded()) {
                    toollist->addToolRow(tool->getToolName(), toolNb);
                }
            }

            if (pp->locallab.spots.size() == 1) {
                setParamEditable(true);
            }

            // Update default values according to selected spot
            setDefaults(pp, pedited);

            // Note: No need to manage pedited as batch mode is deactivated for Locallab

            break;

        case (ControlSpotPanel::SpotDeletion): // Spot deletion event
            // Get deleted spot index in ProcParams and update it
            spotIndex = expsettings->getSelectedSpot();

            for (int i = 0; i < (int)pp->locallab.spots.size(); i++) {
                if (i == spotIndex) {
                    // ProcParams update
                    pp->locallab.spots.erase(pp->locallab.spots.begin() + i);
                    expsettings->deleteControlSpot(spotIndex);

                    // Select the first remaining spot before deleted one
                    if (pp->locallab.spots.size() > 0) {
                        for (int j = i - 1; j >= 0; j--) {
                            if (expsettings->setSelectedSpot(j)) { // True if an existing spot has been selected on controlspotpanel
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

                    // Update tool list widget
                    toolNb = 0;
                    toollist->removeAllTool(); // Reset Locallab list firstly

                    for (auto tool : locallabTools) {
                        toolNb++;

                        if (!tool->isLocallabToolAdded()) {
                            toollist->addToolRow(tool->getToolName(), toolNb);
                        }
                    }

                    if (pp->locallab.spots.size() == 0) {
                        setParamEditable(false);
                    }

                    // Update default values according to selected spot
                    setDefaults(pp, pedited);

                    // Note: No need to manage pedited as batch mode is deactivated for Locallab

                    break;
                }
            }

            break;

        case (ControlSpotPanel::SpotSelection):  // Spot selection event
            pp->locallab.selspot = expsettings->getSelectedSpot();

            // Update control spots and Locallab tools GUI with selected spot
            expsettings->setSelectedSpot(pp->locallab.selspot);
            disableListener();

            for (auto tool : locallabTools) {
                tool->read(pp, pedited);
            }

            enableListener();

            // Update tool list widget
            toolNb = 0;
            toollist->removeAllTool(); // Reset Locallab list firstly

            for (auto tool : locallabTools) {
                toolNb++;

                if (!tool->isLocallabToolAdded()) {
                    toollist->addToolRow(tool->getToolName(), toolNb);
                }
            }

            // Update locallab tools mask background
            if (pp->locallab.selspot < (int)maskBackRef.size()) {
                const double huer = maskBackRef.at(pp->locallab.selspot).huer;
                const double lumar = maskBackRef.at(pp->locallab.selspot).lumar;
                const double chromar = maskBackRef.at(pp->locallab.selspot).chromar;

                for (auto tool : locallabTools) {
                    tool->refChanged(huer, lumar, chromar);
                }
            }

            // Update Locallab Retinex tool min/max
            if (pp->locallab.selspot < (int)retiMinMax.size()) {
                const double cdma = retiMinMax.at(pp->locallab.selspot).cdma;
                const double cdmin = retiMinMax.at(pp->locallab.selspot).cdmin;
                const double mini = retiMinMax.at(pp->locallab.selspot).mini;
                const double maxi = retiMinMax.at(pp->locallab.selspot).maxi;
                const double Tmean = retiMinMax.at(pp->locallab.selspot).Tmean;
                const double Tsigma = retiMinMax.at(pp->locallab.selspot).Tsigma;
                const double Tmin = retiMinMax.at(pp->locallab.selspot).Tmin;
                const double Tmax = retiMinMax.at(pp->locallab.selspot).Tmax;

                expreti->updateMinMax(cdma, cdmin, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
            }

            // Update default values according to selected spot
            setDefaults(pp, pedited);

            // Note: No need to manage pedited as batch mode is deactivated for Locallab

            break;

        case (ControlSpotPanel::SpotDuplication): // Spot duplication event
            newSpot = nullptr;
            spotIndex = expsettings->getSelectedSpot();

            for (int i = 0; i < (int)pp->locallab.spots.size(); i++) {
                if (i == spotIndex) {
                    newSpot = new LocallabParams::LocallabSpot(pp->locallab.spots.at(i));
                    break;
                }
            }

            if (!newSpot) {
                break;
            }

            // Spot creation (initialization at currently selected spot)
            r = new ControlSpotPanel::SpotRow();
            r->name = newSpot->name = newSpot->name + " - " + M("TP_LOCALLAB_DUPLSPOTNAME");
            r->isvisible = newSpot->isvisible;

            if (newSpot->shape == "ELI") {
                r->shape = 0;
            } else {
                r->shape = 1;
            }

            if (newSpot->prevMethod == "hide") {
                r->prevMethod = 0;
            } else {
                r->prevMethod = 1;
            }

            if (newSpot->spotMethod == "norm") {
                r->spotMethod = 0;
            } else if (newSpot->spotMethod == "exc") {
                r->spotMethod = 1;
            } else if (newSpot->spotMethod == "full") {
                r->spotMethod = 2;
            }

            r->sensiexclu = newSpot->sensiexclu;
            r->structexclu = newSpot->structexclu;

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
                    newSpot->loc.at(0) = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->loc.at(0));
                    newSpot->loc.at(1) = rtengine::LIM(int(((double)prW / 2. - 5.) * 2000. / (double)imW), 2, newSpot->loc.at(1));
                    newSpot->loc.at(2) = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->loc.at(2));
                    newSpot->loc.at(3) = rtengine::LIM(int(((double)prH / 2. - 5.) * 2000. / (double)imH), 2, newSpot->loc.at(3));
                }
            }

            r->locX = newSpot->loc.at(0);
            r->locXL = newSpot->loc.at(1);
            r->locY = newSpot->loc.at(2);
            r->locYT = newSpot->loc.at(3);
            r->centerX = newSpot->centerX;
            r->centerY = newSpot->centerY;

            r->circrad = newSpot->circrad;

            if (newSpot->qualityMethod == "enh") {
                r->qualityMethod = 0;
            } else {
                r->qualityMethod = 1;
            }

            r->transit = newSpot->transit;
            r->transitweak = newSpot->transitweak;
            r->transitgrad = newSpot->transitgrad;
            r->feather = newSpot->feather;
            r->struc = newSpot->struc;
            r->thresh = newSpot->thresh;
            r->iter = newSpot->iter;
            r->balan = newSpot->balan;
            r->balanh = newSpot->balanh;
            r->colorde = newSpot->colorde;
            r->colorscope = newSpot->colorscope;
            r->avoidrad = newSpot->avoidrad;
            r->activ = newSpot->activ;
            r->avoid = newSpot->avoid;
            r->avoidmun = newSpot->avoidmun;
            r->blwh = newSpot->blwh;
            r->recurs = newSpot->recurs;
            r->laplac = newSpot->laplac;
            r->deltae = newSpot->deltae;
            r->scopemask = newSpot->scopemask;
            r->shortc = newSpot->shortc;
            r->lumask = newSpot->lumask;
            r->savrest = newSpot->savrest;

            if (newSpot->complexMethod == "sim") {
                r->complexMethod = 0;
            } else  if (newSpot->complexMethod == "mod") {
                r->complexMethod = 1;
            } else  if (newSpot->complexMethod == "all") {
                r->complexMethod = 2;
            }

            if (newSpot->wavMethod == "D2") {
                r->wavMethod = 0;
            } else if (newSpot->wavMethod == "D4") {
                r->wavMethod = 1;
            } else if (newSpot->wavMethod == "D6") {
                r->wavMethod = 2;
            } else if (newSpot->wavMethod == "D10") {
                r->wavMethod = 3;
            } else if (newSpot->wavMethod == "D14") {
                r->wavMethod = 4;
            }

            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.spots.push_back(*newSpot);
            pp->locallab.selspot = pp->locallab.spots.size() - 1;


            // New created spot selection
            expsettings->setSelectedSpot(pp->locallab.selspot);

            // Update Locallab tools GUI with new created spot
            disableListener();

            for (auto tool : locallabTools) {
                tool->read(pp, pedited);
            }

            enableListener();

            // Update tool list widget
            toolNb = 0;
            toollist->removeAllTool(); // Reset Locallab list firstly

            for (auto tool : locallabTools) {
                toolNb++;

                if (!tool->isLocallabToolAdded()) {
                    toollist->addToolRow(tool->getToolName(), toolNb);
                }
            }

            // Update default values according to selected spot
            setDefaults(pp, pedited);

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
            if (pp->locallab.spots.size() > 0) {
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

                    if (r->prevMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).prevMethod = "hide";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).prevMethod = "show";
                    }


                    if (r->spotMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "norm";
                    } else if (r->spotMethod == 1){
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "exc";
                    } else if (r->spotMethod == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).spotMethod = "full";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).sensiexclu = r->sensiexclu;
                    pp->locallab.spots.at(pp->locallab.selspot).structexclu = r->structexclu;

                    if (r->shapeMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "IND";
                    } else if (r->shapeMethod == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "SYM";
                    } else if (r->shapeMethod == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "INDSL";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).shapeMethod = "SYMSL";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).loc.at(0) = r->locX;
                    pp->locallab.spots.at(pp->locallab.selspot).loc.at(1) = r->locXL;
                    pp->locallab.spots.at(pp->locallab.selspot).loc.at(2) = r->locY;
                    pp->locallab.spots.at(pp->locallab.selspot).loc.at(3) = r->locYT;
                    pp->locallab.spots.at(pp->locallab.selspot).centerX = r->centerX;
                    pp->locallab.spots.at(pp->locallab.selspot).centerY = r->centerY;
                    pp->locallab.spots.at(pp->locallab.selspot).circrad = r->circrad;

                    if (r->qualityMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enh";
                    } else {
                        pp->locallab.spots.at(pp->locallab.selspot).qualityMethod = "enhden";
                    }

                    pp->locallab.spots.at(pp->locallab.selspot).transit = r->transit;
                    pp->locallab.spots.at(pp->locallab.selspot).transitweak = r->transitweak;
                    pp->locallab.spots.at(pp->locallab.selspot).transitgrad = r->transitgrad;
                    pp->locallab.spots.at(pp->locallab.selspot).feather = r->feather;
                    pp->locallab.spots.at(pp->locallab.selspot).struc = r->struc;
                    pp->locallab.spots.at(pp->locallab.selspot).thresh = r->thresh;
                    pp->locallab.spots.at(pp->locallab.selspot).iter = r->iter;
                    pp->locallab.spots.at(pp->locallab.selspot).balan = r->balan;
                    pp->locallab.spots.at(pp->locallab.selspot).balanh = r->balanh;
                    pp->locallab.spots.at(pp->locallab.selspot).colorde = r->colorde;
                    pp->locallab.spots.at(pp->locallab.selspot).colorscope = r->colorscope;
                    pp->locallab.spots.at(pp->locallab.selspot).avoidrad = r->avoidrad;
                    pp->locallab.spots.at(pp->locallab.selspot).hishow = r->hishow;
                    pp->locallab.spots.at(pp->locallab.selspot).activ = r->activ;
                    pp->locallab.spots.at(pp->locallab.selspot).avoid = r->avoid;
                    pp->locallab.spots.at(pp->locallab.selspot).avoidmun = r->avoidmun;
                    pp->locallab.spots.at(pp->locallab.selspot).blwh = r->blwh;
                    pp->locallab.spots.at(pp->locallab.selspot).recurs = r->recurs;
                    pp->locallab.spots.at(pp->locallab.selspot).laplac = r->laplac;
                    pp->locallab.spots.at(pp->locallab.selspot).deltae = r->deltae;
                    pp->locallab.spots.at(pp->locallab.selspot).scopemask = r->scopemask;
                    pp->locallab.spots.at(pp->locallab.selspot).shortc = r->shortc;
                    pp->locallab.spots.at(pp->locallab.selspot).lumask = r->lumask;
                    pp->locallab.spots.at(pp->locallab.selspot).savrest = r->savrest;

                    if (r->complexMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).complexMethod = "sim";
                    } else if (r->complexMethod == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).complexMethod = "mod";
                    } else if (r->complexMethod == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).complexMethod = "all";
                    }

                    if (r->wavMethod == 0) {
                        pp->locallab.spots.at(pp->locallab.selspot).wavMethod = "D2";
                    } else if (r->wavMethod == 1) {
                        pp->locallab.spots.at(pp->locallab.selspot).wavMethod = "D4";
                    } else if (r->wavMethod == 2) {
                        pp->locallab.spots.at(pp->locallab.selspot).wavMethod = "D6";
                    } else if (r->wavMethod == 3) {
                        pp->locallab.spots.at(pp->locallab.selspot).wavMethod = "D10";
                    } else if (r->wavMethod == 4) {
                        pp->locallab.spots.at(pp->locallab.selspot).wavMethod = "D14";
                    }
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

    // Set default values in Locallab tools
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

void Locallab::minmaxChanged(const std::vector<locallabRetiMinMax> &minmax, int selspot)
{
    // Saving transmitted min/max data
    retiMinMax = minmax;

    // Update Locallab Retinex tool min/max
    if (selspot < (int)retiMinMax.size()) {
        const double cdma = retiMinMax.at(selspot).cdma;
        const double cdmin = retiMinMax.at(selspot).cdmin;
        const double mini = retiMinMax.at(selspot).mini;
        const double maxi = retiMinMax.at(selspot).maxi;
        const double Tmean = retiMinMax.at(selspot).Tmean;
        const double Tsigma = retiMinMax.at(selspot).Tsigma;
        const double Tmin = retiMinMax.at(selspot).Tmin;
        const double Tmax = retiMinMax.at(selspot).Tmax;

        expreti->updateMinMax(cdma, cdmin, mini, maxi, Tmean, Tsigma, Tmin, Tmax);
    }
}

void Locallab::logencodChanged(const float blackev, const float whiteev, const float sourceg, const float sourceab, const float targetg)
{
    // Update Locallab Log Encoding accordingly
    explog->updateAutocompute(blackev, whiteev, sourceg, sourceab, targetg);
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

    // Reset deltaE preview
    expsettings->resetDeltaEPreview();

    // Reset mask preview for all Locallab tools
    for (auto tool : locallabTools) {
        tool->resetMaskView();
    }
}

Locallab::llMaskVisibility Locallab::getMaskVisibility() const
{
    // Get deltaE preview state
    const bool prevDeltaE = expsettings->isDeltaEPrevActive();

    // Get mask preview from Locallab tools
    int colorMask, colorMaskinv, expMask, expMaskinv, shMask, shMaskinv, vibMask, softMask, blMask, tmMask, retiMask, sharMask, lcMask, cbMask, logMask, maskMask;

    for (auto tool : locallabTools) {
        tool->getMaskView(colorMask, colorMaskinv, expMask, expMaskinv, shMask, shMaskinv, vibMask, softMask, blMask, tmMask, retiMask, sharMask, lcMask, cbMask, logMask, maskMask);
    }

    // Indicate to spot control panel if one mask preview is active
    const bool isMaskActive = (colorMask == 0) || (colorMaskinv == 0) || (expMask == 0) || (expMaskinv == 0) ||
                              (shMask == 0) || (shMaskinv == 0) || (vibMask == 0) || (softMask == 0) ||
                              (blMask == 0) || (tmMask == 0) || (retiMask == 0) || (sharMask == 0) ||
                              (lcMask == 0) || (cbMask == 0) || (logMask == 0) || (maskMask == 0);
    expsettings->setMaskPrevActive(isMaskActive);

    return {prevDeltaE, colorMask, colorMaskinv, expMask, expMaskinv, shMask, shMaskinv, vibMask, softMask, blMask, tmMask, retiMask, sharMask, lcMask, cbMask, logMask, maskMask};
}

void Locallab::resetshowPressed()
{
    // Raise event to reset mask
    if (listener) {
        listener->panelChanged(Evlocallabshowreset, "");
    }
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

void Locallab::openAllTools()
{
    // Set default visibility for settings panel sub-expanders
    expsettings->setDefaultExpanderVisibility();

    for (auto tool : locallabTools) {
        tool->setExpanded(true);

        // Set default visibility for tool sub-expanders
        tool->setDefaultExpanderVisibility();
    }
}

void Locallab::updateShowtooltipVisibility(bool showtooltip)
{
    for (auto tool : locallabTools) {
        tool->updateAdviceTooltips(showtooltip);
    }
}

void Locallab::addTool(Gtk::Box* where, LocallabTool* tool)
{
    tool->getExpander()->setLevel(3);
    where->pack_start(*tool->getExpander(), false, false);
    locallabTools.push_back(tool);
    tool->setLocallabToolListener(this);
}

void Locallab::setParamEditable(bool cond)
{
    // Update params editable state for controlspotpanel
    expsettings->setParamEditable(cond); // TODO Move this code to controlspotpanel.cc when there is zero spot

    // Enable/disable possibility to add Locallab tool
    toollist->set_sensitive(cond);

    // Remove all Locallab tool (without raising event) only if cond is false
    if (!cond) {
        for (auto tool : locallabTools) {
            tool->removeLocallabTool(false);
        }
    }
}

void Locallab::resetToolMaskView()
{
    // Reset mask view GUI for all other Locallab tools
    for (auto tool : locallabTools) {
        tool->resetMaskView();
    }
}

void Locallab::resetOtherMaskView(LocallabTool* current)
{
    // Reset deltaE preview
    expsettings->resetDeltaEPreview();

    // Reset mask view GUI for all other Locallab tools except current
    for (auto tool : locallabTools) {
        if (tool != current) {
            tool->resetMaskView();
        }
    }
}

void Locallab::toolRemoved(LocallabTool* current)
{
    // Update tool list widget according to removed tool
    int toolNb = 0;

    for (auto tool : locallabTools) {
        toolNb++;

        if (tool == current) {
            toollist->addToolRow(tool->getToolName(), toolNb);
        }
    }
}

void Locallab::locallabToolToAdd(const Glib::ustring &toolname)
{
    for (auto tool : locallabTools) {
        if (tool->getToolName() == toolname) {
            // Set expanders visibility default state when adding tool
            tool->setExpanded(true);
            tool->setDefaultExpanderVisibility();

            // Add tool
            tool->addLocallabTool(true);
        }
    }
}
