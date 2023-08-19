/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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

#include <gtkmm.h>
#include "toolpanel.h"
#include "exifpanel.h"
#include "iptcpanel.h"

class MetaDataPanel: public Gtk::Box, public ToolPanel {
private:
    rtengine::ProcEvent EvMetaDataMode;
    MyComboBoxText *metadataMode;
    Gtk::Notebook *tagsNotebook;
    ExifPanel *exifpanel;
    IPTCPanel *iptcpanel;

    void metaDataModeChanged();

public:
    MetaDataPanel();
    ~MetaDataPanel() override;
    
    void setBatchMode(bool batchMode) override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;

    void setImageData(const rtengine::FramesMetaData* id);
    void setListener(ToolPanelListener *tpl) override;

    void setProgressListener(rtengine::ProgressListener *pl);
};

