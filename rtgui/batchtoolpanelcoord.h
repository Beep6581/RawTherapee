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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "fileselectionchangelistener.h"
#include "paramsedited.h"
#include "thumbnaillistener.h"
#include "toolpanelcoord.h"

#include "rtengine/procevents.h"
#include "rtengine/procparams.h"

class FilePanel;
class Thumbnail;
class BatchToolPanelCoordinator final :
    public ToolPanelCoordinator,
    public FileSelectionChangeListener,
    public BatchPParamsChangeListener,
    public ThumbnailListener
{
protected:
    rtengine::procparams::ProcParams pparams;
    ParamsEdited pparamsEdited;
    std::vector<Thumbnail*> selected;
    std::vector<Glib::ustring> selFileNames;
    std::vector<rtengine::procparams::ProcParams> initialPP;
    bool somethingChanged;
    bool blockedUpdate;
    FilePanel* parent;

    void closeSession (bool save = true);
    void initSession ();

public:

    explicit BatchToolPanelCoordinator (FilePanel* parent);

    // FileSelectionChangeListener interface
    void selectionChanged (const std::vector<Thumbnail*>& selected) override;

    // toolpanellistener interface
    void panelChanged(const rtengine::ProcEvent& event, const Glib::ustring& descr) override;
    void setTweakOperator (rtengine::TweakOperator *tOperator) override;
    void unsetTweakOperator (rtengine::TweakOperator *tOperator) override;

    // profilechangelistener interface
    void profileChange(
        const rtengine::procparams::PartialProfile* nparams,
        const rtengine::ProcEvent& event,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr,
        bool fromLastSave = false
    ) override;

    // wbprovider interface
    void getAutoWB (double& temp, double& green, double equal, rtengine::StandardObserver observer, double tempBias) override;
    void getCamWB (double& temp, double& green, rtengine::StandardObserver observer) override;

    // thumbnaillistener interface
    void procParamsChanged (Thumbnail* thm, int whoChangedIt, bool upgradeHint) override;

    // batchpparamschangelistener interface
    void beginBatchPParamsChange(int numberOfEntries) override;
    void endBatchPParamsChange() override;

    // imageareatoollistener interface
    void spotWBselected (int x, int y, Thumbnail* thm = nullptr) override;
    void cropSelectionReady () override;
    void rotateSelectionReady (double rotate_deg, Thumbnail* thm = nullptr) override;
    CropGUIListener* startCropEditing (Thumbnail* thm = nullptr) override;

    void optionsChanged ();
};
