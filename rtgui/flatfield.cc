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

FlatField::FlatField () : FoldableToolPanel(this, "flatfield", M("TP_FLATFIELD_LABEL"))
{
	hbff = Gtk::manage(new Gtk::HBox());
	hbff->set_spacing(2);
	flatFieldFile = Gtk::manage(new MyFileChooserButton(M("TP_FLATFIELD_LABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
	flatFieldFilePersister.reset(new FileChooserLastFolderPersister(flatFieldFile, options.lastFlatfieldDir));
	ffLabel = Gtk::manage(new Gtk::Label(M("GENERAL_FILE")));
	flatFieldFileReset = Gtk::manage(new Gtk::Button());
	flatFieldFileReset->set_image (*Gtk::manage(new RTImage ("gtk-cancel.png")));
	hbff->pack_start(*ffLabel, Gtk::PACK_SHRINK, 0);
	hbff->pack_start(*flatFieldFile);
	hbff->pack_start(*flatFieldFileReset, Gtk::PACK_SHRINK, 0);
	flatFieldAutoSelect = Gtk::manage(new Gtk::CheckButton((M("TP_FLATFIELD_AUTOSELECT"))));
	ffInfo = Gtk::manage(new Gtk::Label(""));
	ffInfo->set_alignment(0,0); //left align
	flatFieldBlurRadius = Gtk::manage(new Adjuster (M("TP_FLATFIELD_BLURRADIUS"),0,200,2,32));
	flatFieldBlurRadius->setAdjusterListener (this);
	if (flatFieldBlurRadius->delay < 1000) flatFieldBlurRadius->delay = 1000;
	flatFieldBlurRadius->show();

	Gtk::HBox* hbffbt = Gtk::manage (new Gtk::HBox ());
	hbffbt->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_FLATFIELD_BLURTYPE") +":")));
	hbffbt->set_spacing(4);
	flatFieldBlurType = Gtk::manage (new MyComboBoxText ());
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_AREA"));
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_VERTICAL"));
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_HORIZONTAL"));
	flatFieldBlurType->append_text(M("TP_FLATFIELD_BT_VERTHORIZ"));
	flatFieldBlurType->set_active(0);
	hbffbt->pack_end (*flatFieldBlurType);

	flatFieldClipControl = Gtk::manage (new Adjuster(M("TP_FLATFIELD_CLIPCONTROL"), 0., 100., 1., 0.));
	flatFieldClipControl->setAdjusterListener(this);
	flatFieldClipControl->addAutoButton("");
	if (flatFieldClipControl->delay < 1000) flatFieldClipControl->delay = 1000;
	flatFieldClipControl->show();
	flatFieldClipControl->set_tooltip_markup (M("TP_FLATFIELD_CLIPCONTROL_TOOLTIP"));

	pack_start( *hbff, Gtk::PACK_SHRINK, 0);
	pack_start( *flatFieldAutoSelect, Gtk::PACK_SHRINK, 0);
	pack_start( *ffInfo, Gtk::PACK_SHRINK, 0);
	pack_start( *hbffbt, Gtk::PACK_SHRINK, 0);
	pack_start( *flatFieldBlurRadius, Gtk::PACK_SHRINK, 0);
	pack_start( *flatFieldClipControl, Gtk::PACK_SHRINK, 0);

	flatFieldFileconn = flatFieldFile->signal_file_set().connect ( sigc::mem_fun(*this, &FlatField::flatFieldFileChanged), true);
	flatFieldFileReset->signal_clicked().connect( sigc::mem_fun(*this, &FlatField::flatFieldFile_Reset), true );
	flatFieldAutoSelectconn = flatFieldAutoSelect->signal_toggled().connect ( sigc::mem_fun(*this, &FlatField::flatFieldAutoSelectChanged), true);
	flatFieldBlurTypeconn = flatFieldBlurType->signal_changed().connect( sigc::mem_fun(*this, &FlatField::flatFieldBlurTypeChanged) );
	lastShortcutPath = "";
	
	// Set filename filters
	b_filter_asCurrent = false;
	Gtk::FileFilter *filter_any = Gtk::manage(new Gtk::FileFilter);
    filter_any->add_pattern("*");
    filter_any->set_name(M("TP_FLATFIELD_FILEDLGFILTERANY"));
    flatFieldFile->add_filter (*filter_any);
    
    // filters for all supported non-raw extensions       
    for (size_t i=0; i<options.parseExtensions.size(); i++) {
        if (options.parseExtensionsEnabled[i] && options.parseExtensions[i].uppercase()!="JPG" && options.parseExtensions[i].uppercase()!="JPEG" && options.parseExtensions[i].uppercase()!="PNG" && options.parseExtensions[i].uppercase()!="TIF" && options.parseExtensions[i].uppercase()!="TIFF"  )
        {    
	        Gtk::FileFilter *filter_ff = Gtk::manage(new Gtk::FileFilter);
	        filter_ff->add_pattern("*." + options.parseExtensions[i]);
   	        filter_ff->add_pattern("*." + options.parseExtensions[i].uppercase());
	        filter_ff->set_name(options.parseExtensions[i].uppercase());
	        flatFieldFile->add_filter (*filter_ff);
	        //printf("adding filter %s \n",options.parseExtensions[i].uppercase().c_str());
        }
    }
}

void FlatField::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	flatFieldAutoSelectconn.block (true);
	flatFieldBlurTypeconn.block (true);

	//flatFieldBlurType
	for( size_t i=0; i< procparams::RAWParams::numFlatFileBlurTypes;i++)
		if( pp->raw.ff_BlurType == procparams::RAWParams::ff_BlurTypestring[i]){
			flatFieldBlurType->set_active(i);
			break;
		}

	if (multiImage || pp->raw.ff_BlurType == procparams::RAWParams::ff_BlurTypestring[procparams::RAWParams::area_ff])
		flatFieldClipControl->show();
	else
		flatFieldClipControl->hide();

	flatFieldAutoSelect->set_active (pp->raw.ff_AutoSelect);
	flatFieldBlurRadius->setValue (pp->raw.ff_BlurRadius);
	flatFieldClipControl->setValue (pp->raw.ff_clipControl);
	flatFieldClipControl->setAutoValue (pp->raw.ff_AutoClipControl);

	if(pedited ){
		flatFieldAutoSelect->set_inconsistent (!pedited->raw.ff_AutoSelect);
		flatFieldBlurRadius->setEditedState( pedited->raw.ff_BlurRadius ? Edited : UnEdited );
		flatFieldClipControl->setEditedState( pedited->raw.ff_clipControl ? Edited : UnEdited );
		flatFieldClipControl->setAutoInconsistent(multiImage && !pedited->raw.ff_AutoClipControl);
		if( !pedited->raw.ff_BlurType )
			flatFieldBlurType->set_active(procparams::RAWParams::numFlatFileBlurTypes); // No name
	}
	if (safe_file_test (pp->raw.ff_file, Glib::FILE_TEST_EXISTS))
		flatFieldFile->set_filename (pp->raw.ff_file);
	else
		flatFieldFile_Reset();
	hbff->set_sensitive( !pp->raw.ff_AutoSelect );

	lastFFAutoSelect = pp->raw.ff_AutoSelect;
	lastFFAutoClipCtrl = pp->raw.ff_AutoClipControl;

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

	ffChanged = false;

	flatFieldAutoSelectconn.block (false);
	flatFieldBlurTypeconn.block (false);
	enableListener ();
	
	// Add filter with the current file extension if the current file is raw
	if (ffp && !batchMode){
	
		if (b_filter_asCurrent){
			//First, remove last filter_asCurrent if it was set for a raw file
			std::vector<const Gtk::FileFilter*> filters = flatFieldFile->list_filters();
			flatFieldFile->remove_filter(**(filters.end()-1));
			b_filter_asCurrent = false;
		}

		Glib::ustring fname = Glib::path_get_basename(ffp->GetCurrentImageFilePath());
		Glib::ustring filetype;

		if (fname!=""){
			// get image filetype, set filter to the same as current image's filetype
			std::string::size_type idx;
			idx = fname.rfind('.');
			if(idx != std::string::npos){
				filetype = fname.substr(idx+1);
				//exclude non-raw
				israw = filetype.uppercase()!="JPG" && filetype.uppercase()!="JPEG" && filetype.uppercase()!="PNG" && filetype.uppercase()!="TIF" && filetype.uppercase()!="TIFF";

				if (israw)
				{
					b_filter_asCurrent = true; //prevent re-adding this filter on every pp3 file read
					Gtk::FileFilter *filter_asCurrent = Gtk::manage(new Gtk::FileFilter);
					filter_asCurrent->add_pattern("*." + filetype);
					filter_asCurrent->set_name(M("TP_FLATFIELD_FILEDLGFILTERFF") + " (" + filetype + ")");
					flatFieldFile->add_filter (*filter_asCurrent);
					flatFieldFile->set_filter (*filter_asCurrent);
				}
			}
		}
	}

}

void FlatField::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.ff_file = flatFieldFile->get_filename();
	pp->raw.ff_AutoSelect = flatFieldAutoSelect->get_active();
	pp->raw.ff_BlurRadius = flatFieldBlurRadius->getIntValue();
	pp->raw.ff_clipControl = flatFieldClipControl->getIntValue();
	pp->raw.ff_AutoClipControl = flatFieldClipControl->getAutoValue();

	int currentRow = flatFieldBlurType->get_active_row_number();
	if( currentRow>=0 && currentRow < procparams::RAWParams::numFlatFileBlurTypes)
		pp->raw.ff_BlurType = procparams::RAWParams::ff_BlurTypestring[currentRow];

	if (pedited) {
		pedited->raw.ff_file = ffChanged;
		pedited->raw.ff_AutoSelect = !flatFieldAutoSelect->get_inconsistent();
		pedited->raw.ff_BlurRadius = flatFieldBlurRadius->getEditedState ();
		pedited->raw.ff_clipControl = flatFieldClipControl->getEditedState ();
		pedited->raw.ff_AutoClipControl = !flatFieldClipControl->getAutoInconsistent();
		pedited->raw.ff_BlurType = flatFieldBlurType->get_active_row_number() != procparams::RAWParams::numFlatFileBlurTypes;
	}

}

void FlatField::adjusterChanged (Adjuster* a, double newval)
{
	if (listener) {

		Glib::ustring value = a->getTextValue();

		if (a==flatFieldBlurRadius)
			listener->panelChanged (EvFlatFieldBlurRadius,  value);
		else if (a==flatFieldClipControl)
			listener->panelChanged (EvFlatFieldClipControl,  value);
	}
}

void FlatField::adjusterAutoToggled (Adjuster* a, bool newval) {

	if (multiImage) {
		if (flatFieldClipControl->getAutoInconsistent()) {
			flatFieldClipControl->setAutoInconsistent(false);
			flatFieldClipControl->setAutoValue(false);
		}
		else if (lastFFAutoClipCtrl)
			flatFieldClipControl->setAutoInconsistent(true);

		lastFFAutoClipCtrl = flatFieldClipControl->getAutoValue();

	}

	if (listener) {
		if(a==flatFieldClipControl) {
			if (flatFieldClipControl->getAutoInconsistent())
				listener->panelChanged (EvFlatFieldAutoClipControl, M("GENERAL_UNCHANGED"));
			else if (flatFieldClipControl->getAutoValue())
				listener->panelChanged (EvFlatFieldAutoClipControl, M("GENERAL_ENABLED"));
			else
				listener->panelChanged (EvFlatFieldAutoClipControl, M("GENERAL_DISABLED"));
		}
	}
}

void FlatField::setBatchMode(bool batchMode)
{
	ToolPanel::setBatchMode (batchMode);
	flatFieldBlurRadius->showEditedCB ();
	flatFieldClipControl->showEditedCB ();
}

void FlatField::setAdjusterBehavior (bool clipctrladd) {
	flatFieldClipControl->setAddMode(clipctrladd);
}

void FlatField::trimValues (rtengine::procparams::ProcParams* pp) {
	flatFieldClipControl->trimValue(pp->raw.ff_clipControl);
}

void FlatField::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	flatFieldBlurRadius->setDefault( defParams->raw.ff_BlurRadius);
	flatFieldClipControl->setDefault( defParams->raw.ff_clipControl);

	if (pedited) {
		flatFieldBlurRadius->setDefaultEditedState( pedited->raw.ff_BlurRadius ? Edited : UnEdited);
		flatFieldClipControl->setDefaultEditedState( pedited->raw.ff_clipControl ? Edited : UnEdited);
	} else {
		flatFieldBlurRadius->setDefaultEditedState( Irrelevant );
		flatFieldClipControl->setDefaultEditedState( Irrelevant );
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

	if (multiImage || curSelection == procparams::RAWParams::area_ff)
		flatFieldClipControl->show();
	else
		flatFieldClipControl->hide();

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
