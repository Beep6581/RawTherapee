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
#ifndef __BATCHTOOLPANELCCORD__
#define __BATCHTOOLPANELCCORD__

#include "thumbnail.h"
#include "toolpanelcoord.h"
#include "fileselectionchangelistener.h"
#include "../rtengine/rtengine.h"
#include "paramsedited.h"
#include "thumbnaillistener.h"

class FilePanel;
class BatchToolPanelCoordinator :
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
    void selectionChanged (const std::vector<Thumbnail*>& selected);

    // toolpanellistener interface
    void panelChanged   (rtengine::ProcEvent event, const Glib::ustring& descr);

    // profilechangelistener interface
    void profileChange  (const rtengine::procparams::PartialProfile* nparams, rtengine::ProcEvent event, const Glib::ustring& descr, const ParamsEdited* paramsEdited = nullptr);

    // wbprovider interface
    void getAutoWB (double& temp, double& green, double equal, double tempBias);
    void getCamWB (double& temp, double& green);

    // thumbnaillistener interface
    void procParamsChanged (Thumbnail* thm, int whoChangedIt);

    // batchpparamschangelistener interface
    void beginBatchPParamsChange(int numberOfEntries);
    void endBatchPParamsChange();

    // imageareatoollistener interface
    void spotWBselected (int x, int y, Thumbnail* thm = nullptr);
    void cropSelectionReady ();
    void rotateSelectionReady (double rotate_deg, Thumbnail* thm = nullptr);
    CropGUIListener* startCropEditing (Thumbnail* thm = nullptr);

    void optionsChanged ();
};

#endif
