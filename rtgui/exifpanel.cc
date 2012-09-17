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
#include "exifpanel.h"
#include "../rtengine/safegtk.h"
#include "guiutils.h"
#include "rtimage.h"
#include <cstdio>

using namespace rtengine;
using namespace rtengine::procparams;

ExifPanel::ExifPanel () : idata(NULL) {

   exifTree = Gtk::manage(new Gtk::TreeView());
   scrolledWindow = Gtk::manage(new Gtk::ScrolledWindow());

   exifTree->set_headers_visible(false);
   exifTree->set_rules_hint(false);
   exifTree->set_reorderable(false);
   exifTree->set_enable_search(true);
   exifTree->get_selection()->set_mode (Gtk::SELECTION_MULTIPLE);
   scrolledWindow->set_border_width(2);
   scrolledWindow->set_shadow_type(Gtk::SHADOW_NONE);
   scrolledWindow->set_policy(Gtk::POLICY_ALWAYS, Gtk::POLICY_ALWAYS);
   scrolledWindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);
   //scrolledWindow->add(*exifTree);

   Gtk::VBox* vboxMain = Gtk::manage (new Gtk::VBox ());
   Gtk::HSeparator* hseptPE = Gtk::manage (new Gtk::HSeparator ());
   Gtk::Image* peImg = Gtk::manage (new RTImage("PanelEnding.png"));
   // add panel ending and exifTree to vboxMain
   vboxMain->pack_end(*peImg, Gtk::PACK_SHRINK, 0);
   vboxMain->pack_end(*hseptPE, Gtk::PACK_SHRINK, 4);
   vboxMain->pack_start(*exifTree, Gtk::PACK_EXPAND_WIDGET,0);
   scrolledWindow->add(*vboxMain);
   
   exifTreeModel = Gtk::TreeStore::create(exifColumns);
   exifTree->set_model (exifTreeModel);

   Gtk::TreeView::Column *viewcol = Gtk::manage(new Gtk::TreeView::Column ("Field Name"));
   Gtk::CellRendererPixbuf* render_pb = Gtk::manage(new Gtk::CellRendererPixbuf ());
   Gtk::CellRendererText *render_txt = Gtk::manage(new Gtk::CellRendererText());
   viewcol->pack_start (*render_pb, false);
   viewcol->pack_start (*render_txt, true);
   viewcol->add_attribute (*render_pb, "pixbuf", exifColumns.icon);
   viewcol->add_attribute (*render_txt, "markup", exifColumns.field);
 
   render_pb->property_ypad() = 0;
   render_txt->property_ypad() = 0;
   render_pb->property_yalign() = 0;
   render_txt->property_yalign() = 0;
  
   exifTree->append_column (*viewcol); 
   
   Gtk::TreeView::Column *viewcolv = Gtk::manage(new Gtk::TreeView::Column ("Value"));
   Gtk::CellRendererText *render_txtv = Gtk::manage(new Gtk::CellRendererText());
   viewcolv->pack_start (*render_txtv, true);
   viewcolv->add_attribute (*render_txtv, "markup", exifColumns.value);
   
   render_txtv->property_ypad() = 0;
   
   exifTree->append_column (*viewcolv); 
   
   pack_start (*scrolledWindow);
   
   show_all ();
}

ExifPanel::~ExifPanel () {
}

void ExifPanel::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    setImageData (idata);
    
    enableListener ();
}

void ExifPanel::write (ProcParams* pp, ParamsEdited* pedited) {


}

void ExifPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {


}
        
void ExifPanel::setImageData (const ImageMetaData* id) {
    
    idata = id; 
    exifTreeModel->clear ();

    if( !idata )
        return;

    std::vector<ExifPair>  list=idata->getExifData ();

    // sort elements for grouping into different folders
    std::sort( list.begin(),list.end());
    Glib::ustring currentGroup="";
    Gtk::TreeModel::Children chldr=exifTreeModel->children();
    for( std::vector<ExifPair>::iterator iter= list.begin(); iter != list.end(); iter++ ) {
        if( iter->group != currentGroup ){
            chldr = addGroup (exifTreeModel->children(), iter->group );
            currentGroup = iter->group;
        }
        if (iter->value.length() > 80 )
            iter->value = iter->value.substr(0, 80) + " [...]";
        addTag (chldr, iter->field, iter->value );
    }

}

Gtk::TreeModel::Children ExifPanel::addGroup (const Gtk::TreeModel::Children& root, Glib::ustring group )
{
    Gtk::TreeModel::Row row = *(exifTreeModel->append(root));

	row[exifColumns.field] = Glib::ustring("<b>") + group + "</b>";

	return row.children();
}

Gtk::TreeModel::Children ExifPanel::addTag (const Gtk::TreeModel::Children& root, Glib::ustring field, Glib::ustring value )
{
    Gtk::TreeModel::Row row = *(exifTreeModel->append(root));

    row[exifColumns.orig_value] = value;

    row[exifColumns.field] = field;
    row[exifColumns.value] = value;
    row[exifColumns.value] = Glib::ustring("<i>") + value +  Glib::ustring("</i>");
    
    return row.children();
}


void ExifPanel::notifyListener () {
    
    if (listener)
        listener->panelChanged (EvExif, M("HISTORY_CHANGED"));
}
