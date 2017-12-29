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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "metadatapanel.h"
#include "eventmapper.h"
#include "../rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;


MetaDataPanel::MetaDataPanel()
{
    EvMetaDataMode = ProcEventMapper::getInstance()->newEvent(M_VOID, "HISTORY_MSG_METADATA_MODE");

    Gtk::HBox *box = Gtk::manage(new Gtk::HBox());
    box->pack_start(*Gtk::manage(new Gtk::Label(M("TP_METADATA_MODE") + ": ")), Gtk::PACK_SHRINK, 4);
    metadataMode = Gtk::manage(new MyComboBoxText());
    metadataMode->append(M("TP_METADATA_TUNNEL"));
    metadataMode->append(M("TP_METADATA_EDIT"));
    metadataMode->append(M("TP_METADATA_STRIP"));
    metadataMode->set_active(0);
    box->pack_end(*metadataMode, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start(*box, Gtk::PACK_SHRINK, 4);

    metadataMode->signal_changed().connect(sigc::mem_fun(*this, &MetaDataPanel::metaDataModeChanged));

    tagsNotebook = Gtk::manage(new Gtk::Notebook());
    exifpanel = Gtk::manage(new ExifPanel());
    iptcpanel = Gtk::manage(new IPTCPanel());
    tagsNotebook->set_name("MetaPanelNotebook");
    tagsNotebook->append_page(*exifpanel, M("MAIN_TAB_EXIF"));
    tagsNotebook->append_page(*iptcpanel, M("MAIN_TAB_IPTC"));

    pack_end(*tagsNotebook);
}


void MetaDataPanel::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    metadataMode->append(M("GENERAL_UNCHANGED"));
    tagsNotebook->remove_page(-1);
    tagsNotebook->remove_page(-1);
}


void MetaDataPanel::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();
    metadataMode->set_active(int(pp->metadata.mode));
    if (pedited) {
        if (!pedited->metadata.mode) {
            metadataMode->set_active(3);
        }
    }

    exifpanel->read(pp, pedited);
    iptcpanel->read(pp, pedited);
    
    enableListener();
}


void MetaDataPanel::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->metadata.mode = static_cast<MetaDataParams::Mode>(max(metadataMode->get_active_row_number(), 2));
    if (pedited) {
        pedited->metadata.mode = metadataMode->get_active_row_number() != 3;
    }

    exifpanel->write(pp, pedited);
    iptcpanel->write(pp, pedited);
}


void MetaDataPanel::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    exifpanel->setDefaults(defParams, pedited);
    iptcpanel->setDefaults(defParams, pedited);
}


void MetaDataPanel::setImageData(const rtengine::FramesMetaData* id)
{
    exifpanel->setImageData(id);
    iptcpanel->setImageData(id);
}


void MetaDataPanel::setListener(ToolPanelListener *tpl)
{
    ToolPanel::setListener(tpl);
    exifpanel->setListener(tpl);
    iptcpanel->setListener(tpl);
}


void MetaDataPanel::metaDataModeChanged()
{
    if (listener) {
        listener->panelChanged(EvMetaDataMode, M("HISTORY_CHANGED"));
    }
}
