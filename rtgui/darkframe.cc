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
#include "darkframe.h"
#include "options.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include <sstream>
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

DarkFrame::DarkFrame () : Gtk::VBox(), FoldableToolPanel(this)
{
	hbdf = Gtk::manage(new Gtk::HBox());
	darkFrameFile = Gtk::manage(new MyFileChooserButton(M("TP_DARKFRAME_LABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
	darkFrameFilePersister.reset(new FileChooserLastFolderPersister(darkFrameFile, options.lastDarkframeDir));
	dfLabel = Gtk::manage(new Gtk::Label(M("GENERAL_FILE")));
	btnReset = Gtk::manage(new Gtk::Button());
	btnReset->set_image (*Gtk::manage(new RTImage ("gtk-cancel.png")));
	hbdf->pack_start(*dfLabel, Gtk::PACK_SHRINK, 4);
	hbdf->pack_start(*darkFrameFile);
	hbdf->pack_start(*btnReset, Gtk::PACK_SHRINK, 4);
	dfAuto = Gtk::manage(new Gtk::CheckButton((M("TP_DARKFRAME_AUTOSELECT"))));
	dfInfo = Gtk::manage(new Gtk::Label(""));
	dfInfo->set_alignment(0,0); //left align

	pack_start( *hbdf, Gtk::PACK_SHRINK, 4);
	pack_start( *dfAuto, Gtk::PACK_SHRINK, 4);
	pack_start( *dfInfo, Gtk::PACK_SHRINK, 4);

	dfautoconn = dfAuto->signal_toggled().connect ( sigc::mem_fun(*this, &DarkFrame::dfAutoChanged), true);
	dfFile = darkFrameFile->signal_file_set().connect ( sigc::mem_fun(*this, &DarkFrame::darkFrameChanged), true);
	btnReset->signal_clicked().connect( sigc::mem_fun(*this, &DarkFrame::darkFrameReset), true );
}

void DarkFrame::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
	disableListener ();
	dfautoconn.block(true);

	if(pedited ){
		dfAuto->set_inconsistent(!pedited->raw.dfAuto );
	}
	if (safe_file_test (pp->raw.dark_frame, Glib::FILE_TEST_EXISTS))
		darkFrameFile->set_filename (pp->raw.dark_frame);
	hbdf->set_sensitive( !pp->raw.df_autoselect );

	lastDFauto = pp->raw.df_autoselect;

	if( pp->raw.df_autoselect  && dfp && !batchMode){
		// retrieve the auto-selected df filename
		rtengine::RawImage *img = dfp->getDF();
		if( img ){
			dfInfo->set_text( Glib::ustring::compose("%1: %2ISO %3s", Glib::path_get_basename(img->get_filename()), img->get_ISOspeed(), img->get_shutter()) );
		}else{
			dfInfo->set_text(Glib::ustring(M("TP_PREPROCESS_NO_FOUND")));
		}
	}
	else dfInfo->set_text("");

	dfAuto->set_active( pp->raw.df_autoselect );
	dfChanged = false;

	dfautoconn.block(false);
	enableListener ();
}

void DarkFrame::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.dark_frame = darkFrameFile->get_filename();
	pp->raw.df_autoselect = dfAuto->get_active();

	if (pedited) {
		pedited->raw.darkFrame = dfChanged;
		pedited->raw.dfAuto = !dfAuto->get_inconsistent();
	}

}

void DarkFrame::dfAutoChanged()
{
    if (batchMode) {
        if (dfAuto->get_inconsistent()) {
        	dfAuto->set_inconsistent (false);
        	dfautoconn.block (true);
        	dfAuto->set_active (false);
        	dfautoconn.block (false);
        }
        else if (lastDFauto)
        	dfAuto->set_inconsistent (true);

        lastDFauto = dfAuto->get_active ();
    }

    if(dfAuto->get_active() && dfp && !batchMode){
 	 // retrieve the auto-selected df filename
      rtengine::RawImage *img = dfp->getDF();
      if( img ){
        dfInfo->set_text( Glib::ustring::compose("%1: %2ISO %3s", Glib::path_get_basename(img->get_filename()), img->get_ISOspeed(), img->get_shutter()) );
      }else{
		dfInfo->set_text(Glib::ustring(M("TP_PREPROCESS_NO_FOUND")));
	  }
    }
    else{dfInfo->set_text("");}

	hbdf->set_sensitive( !dfAuto->get_active() );
    if (listener)
        listener->panelChanged (EvPreProcessAutoDF, dfAuto->get_active()?M("GENERAL_ENABLED"):M("GENERAL_DISABLED"));
}

void DarkFrame::darkFrameChanged()
{
	dfChanged=true;
    if (listener)
        listener->panelChanged (EvPreProcessDFFile, Glib::path_get_basename(darkFrameFile->get_filename()));
}

void DarkFrame::darkFrameReset()
{
	dfChanged=true;
	//darkFrameFile->set_current_name("");
	darkFrameFile->set_filename ("");

	if (!options.lastDarkframeDir.empty())
		darkFrameFile->set_current_folder(options.lastDarkframeDir);

	dfInfo->set_text("");
    if (listener)
        listener->panelChanged (EvPreProcessDFFile, M("GENERAL_NONE"));

}
