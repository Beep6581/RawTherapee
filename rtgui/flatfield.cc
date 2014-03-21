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
#include "flatfield.h"
#include "options.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

FlatField::FlatField () : Gtk::VBox(), FoldableToolPanel(this)
{
	set_border_width(4);

	hbff = Gtk::manage(new Gtk::HBox());
	flatFieldFile = Gtk::manage(new MyFileChooserButton(M("TP_FLATFIELD_LABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
	flatFieldFilePersister.reset(new FileChooserLastFolderPersister(flatFieldFile, options.lastFlatfieldDir));
	ffLabel = Gtk::manage(new Gtk::Label(M("GENERAL_FILE")));
	flatFieldFileReset = Gtk::manage(new Gtk::Button());
	flatFieldFileReset->set_image (*Gtk::manage(new RTImage ("gtk-cancel.png")));
	hbff->pack_start(*ffLabel, Gtk::PACK_SHRINK, 4);
	hbff->pack_start(*flatFieldFile);
	hbff->pack_start(*flatFieldFileReset, Gtk::PACK_SHRINK, 4);
	flatFieldAutoSelect = Gtk::manage(new Gtk::CheckButton((M("TP_FLATFIELD_AUTOSELECT"))));
	ffInfo = Gtk::manage(new Gtk::Label(""));
	ffInfo->set_alignment(0,0); //left align
	flatFieldBlurRadius = Gtk::manage(new Adjuster (M("TP_FLATFIELD_BLURRADIUS"),0,200,2,32));
	flatFieldBlurRadius->setAdjusterListener (this);
	if (flatFieldBlurRadius->delay < 1000) flatFieldBlurRadius->delay = 1000;
	flatFieldBlurRadius->show();

	Gtk::HBox* hbffbt = Gtk::manage (new Gtk::HBox ());
	hbffbt->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_FLATFIELD_BLURTYPE") +": ")));
	flatFieldBlurType = Gtk::manage (new MyComboBoxText ());
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_AREA"));
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_VERTICAL"));
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_HORIZONTAL"));
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_VERTHORIZ"));
	flatFieldBlurType->set_active(0);
	hbffbt->pack_end (*flatFieldBlurType);

	pack_start( *hbff, Gtk::PACK_SHRINK, 4);
	pack_start( *flatFieldAutoSelect, Gtk::PACK_SHRINK, 4);
	pack_start( *ffInfo, Gtk::PACK_SHRINK, 4);
	pack_start( *hbffbt, Gtk::PACK_SHRINK, 4);
	pack_start( *flatFieldBlurRadius, Gtk::PACK_SHRINK, 4);

	flatFieldFileconn = flatFieldFile->signal_file_set().connect ( sigc::mem_fun(*this, &FlatField::flatFieldFileChanged), true);
	flatFieldFileReset->signal_clicked().connect( sigc::mem_fun(*this, &FlatField::flatFieldFile_Reset), true );
	flatFieldAutoSelectconn = flatFieldAutoSelect->signal_toggled().connect ( sigc::mem_fun(*this, &FlatField::flatFieldAutoSelectChanged), true);
	flatFieldBlurTypeconn = flatFieldBlurType->signal_changed().connect( sigc::mem_fun(*this, &FlatField::flatFieldBlurTypeChanged) );
	lastShortcutPath = "";
}

void FlatField::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	flatFieldAutoSelectconn.block (true);
	flatFieldBlurTypeconn.block (true);

	if(pedited ){
		flatFieldAutoSelect->set_inconsistent (!pedited->raw.ff_AutoSelect);
		flatFieldBlurRadius->setEditedState( pedited->raw.ff_BlurRadius ? Edited : UnEdited );
		if( !pedited->raw.ff_BlurType )
			flatFieldBlurType->set_active(procparams::RAWParams::numFlatFileBlurTypes); // No name
	}
	if (safe_file_test (pp->raw.ff_file, Glib::FILE_TEST_EXISTS))
		flatFieldFile->set_filename (pp->raw.ff_file);
	else
		flatFieldFile_Reset();
	hbff->set_sensitive( !pp->raw.ff_AutoSelect );

	lastFFAutoSelect = pp->raw.ff_AutoSelect;
	if( pp->raw.ff_AutoSelect  && ffp && !batchMode){
		// retrieve the auto-selected ff filename
		rtengine::RawImage *img = ffp->getFF();
		if( img ){
			ffInfo->set_text( Glib::ustring::compose("%1: f/%2", Glib::path_get_basename(img->get_filename()), img->get_aperture()) ); // !!! need to add focallength in mm and format aperture to ##.#
		}else{
			ffInfo->set_text(Glib::ustring(M("TP_PREPROCESS_NO_FOUND")));
		}
	}
	else ffInfo->set_text("");

	flatFieldAutoSelect ->set_active(pp->raw.ff_AutoSelect);
	flatFieldBlurRadius->setValue (pp->raw.ff_BlurRadius);
	flatFieldBlurType->set_active(procparams::RAWParams::numFlatFileBlurTypes);

	//flatFieldBlurType
	for( size_t i=0; i< procparams::RAWParams::numFlatFileBlurTypes;i++)
		if( pp->raw.ff_BlurType == procparams::RAWParams::ff_BlurTypestring[i]){
			flatFieldBlurType->set_active(i);
			break;
		}

	flatFieldAutoSelect->set_active (pp->raw.ff_AutoSelect);
	flatFieldBlurRadius->setValue (pp->raw.ff_BlurRadius);

	ffChanged = false;

	flatFieldAutoSelectconn.block (false);
	flatFieldBlurTypeconn.block (false);
	enableListener ();
}

void FlatField::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.ff_file = flatFieldFile->get_filename();
	pp->raw.ff_AutoSelect = flatFieldAutoSelect->get_active();
	pp->raw.ff_BlurRadius = flatFieldBlurRadius->getIntValue();

	int currentRow = flatFieldBlurType->get_active_row_number();
	if( currentRow>=0 && currentRow < procparams::RAWParams::numFlatFileBlurTypes)
		pp->raw.ff_BlurType = procparams::RAWParams::ff_BlurTypestring[currentRow];

	if (pedited) {
		pedited->raw.ff_file = ffChanged;
		pedited->raw.ff_AutoSelect = !flatFieldAutoSelect->get_inconsistent();
		pedited->raw.ff_BlurRadius = flatFieldBlurRadius->getEditedState ();
		pedited->raw.ff_BlurType = flatFieldBlurType->get_active_row_number() != procparams::RAWParams::numFlatFileBlurTypes;
	}

}

void FlatField::adjusterChanged (Adjuster* a, double newval)
{
	if (listener) {

		Glib::ustring value = a->getTextValue();

		listener->panelChanged (EvFlatFieldBlurRadius,  value );
	}
}

void FlatField::setBatchMode(bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	flatFieldBlurRadius->showEditedCB ();
}

void FlatField::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	flatFieldBlurRadius->setDefault( defParams->raw.ff_BlurRadius);

	if (pedited) {
		flatFieldBlurRadius->setDefaultEditedState( pedited->raw.ff_BlurRadius ? Edited : UnEdited);
	} else {
		flatFieldBlurRadius->setDefaultEditedState( Irrelevant );
	}
}

void FlatField::flatFieldFileChanged()
{
	ffChanged=true;
    if (listener)
        listener->panelChanged (EvFlatFieldFile, Glib::path_get_basename(flatFieldFile->get_filename()));
}

void FlatField::flatFieldFile_Reset()
{
	ffChanged=true;

// caution: I had to make this hack, because set_current_folder() doesn't work correctly!
//          Because szeva doesn't exist since he was committed to happy hunting ground in Issue 316
//          we can use him now for this hack
	flatFieldFile->set_filename (options.lastFlatfieldDir + "/szeva");
// end of the hack

	if (!options.lastFlatfieldDir.empty())
		flatFieldFile->set_current_folder(options.lastFlatfieldDir);

	ffInfo->set_text("");
    if (listener)
        listener->panelChanged (EvFlatFieldFile, M("GENERAL_NONE") );
}

void FlatField::flatFieldBlurTypeChanged ()
{
	int  curSelection = flatFieldBlurType->get_active_row_number();

	Glib::ustring s="";
	if( curSelection>=0 && curSelection < procparams::RAWParams::numFlatFileBlurTypes)
	    s = flatFieldBlurType->get_active_text();

    if (listener)
        listener->panelChanged (EvFlatFieldBlurType, s);
}

void FlatField::flatFieldAutoSelectChanged()
{
    if (batchMode) {
        if (flatFieldAutoSelect->get_inconsistent()) {
        	flatFieldAutoSelect->set_inconsistent (false);
        	flatFieldAutoSelectconn.block (true);
        	flatFieldAutoSelect->set_active (false);
        	flatFieldAutoSelectconn.block (false);
        }
        else if (lastFFAutoSelect)
        	flatFieldAutoSelect->set_inconsistent (true);

        lastFFAutoSelect = flatFieldAutoSelect->get_active ();
    }
	hbff->set_sensitive( !flatFieldAutoSelect->get_active() );

    if( flatFieldAutoSelect->get_active()  && ffp && !batchMode){
 	  // retrieve the auto-selected ff filename
       rtengine::RawImage *img = ffp->getFF();
      if( img ){
        ffInfo->set_text( Glib::ustring::compose("%1: f/%2", Glib::path_get_basename(img->get_filename()), img->get_aperture()) ); // !!! need to add focallength in mm and format aperture to ##.#
      }else{
    	  ffInfo->set_text(Glib::ustring(M("TP_PREPROCESS_NO_FOUND")));
      }
    }
    else{ffInfo->set_text("");}

    if (listener)
        listener->panelChanged (EvFlatFieldAutoSelect, flatFieldAutoSelect->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));

}

void FlatField::setShortcutPath(Glib::ustring path)
{
	if (path == "") return;
#ifdef WIN32
	// Dirty workaround, waiting for a clean solution by using exceptions!
	if (!safe_is_shortcut_dir(path))
#endif
	{
		if (lastShortcutPath != "") {
			try {
				flatFieldFile->remove_shortcut_folder(lastShortcutPath);
			}
			catch (Glib::Error &err) {}
		}
		lastShortcutPath = path;
		try {
			flatFieldFile->add_shortcut_folder(path);
		}
		catch (Glib::Error &err) {}
	}
}
