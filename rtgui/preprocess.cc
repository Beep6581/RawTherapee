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
#include <preprocess.h>
#include <options.h>
#include <guiutils.h>
using namespace rtengine;
using namespace rtengine::procparams;

PreProcess::PreProcess (): Gtk::VBox(), FoldableToolPanel(this)
{
	hbdf = Gtk::manage(new Gtk::HBox());
	darkFrameFile = Gtk::manage(new Gtk::FileChooserButton(M("TP_PREPROCESS_DARKFRAME"), Gtk::FILE_CHOOSER_ACTION_OPEN));
	dfLabel = Gtk::manage(new Gtk::Label(M("TP_PREPROCESS_DARKFRAME")));
	btnReset = Gtk::manage(new Gtk::Button());
	btnReset->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-cancel"), Gtk::ICON_SIZE_BUTTON)));
	hbdf->pack_start(*dfLabel, Gtk::PACK_SHRINK, 4);
	hbdf->pack_start(*darkFrameFile);
	hbdf->pack_start(*btnReset, Gtk::PACK_SHRINK, 4);
	dfAuto = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_DFAUTOSELECT"))));
	
	caAutocorrect = Gtk::manage(new Gtk::CheckButton((M("PREFERENCES_CACORRECTION"))));
	caRed = Gtk::manage(new Adjuster (M("PREFERENCES_CARED"),-4.0,4.0,0.1,0));
	caRed->setAdjusterListener (this);
	caRed->show();
	caBlue = Gtk::manage(new Adjuster (M("PREFERENCES_CABLUE"),-4.0,4.0,0.1,0));
	caBlue->setAdjusterListener (this);
	caBlue->show();
//exposi
	/*exPos = Gtk::manage(new Adjuster (M("PREFERENCES_EXPOS"),0.2,4.0,0.1,1));
	exPos->setAdjusterListener (this);
	exPos->show();
	exPreser = Gtk::manage(new Adjuster (M("PREFERENCES_PRESER"),0,2.5,0.1,0));
	exPreser->setAdjusterListener (this);
	exPreser->show();*/
	PexPos = Gtk::manage(new Adjuster (M("PREFERENCES_EXPOS"),0.2,4.0,0.1,1));
	PexPos->setAdjusterListener (this);
	PexPos->show();
	PexPreser = Gtk::manage(new Adjuster (M("PREFERENCES_PRESER"),0,2.5,0.1,0));
	PexPreser->setAdjusterListener (this);
	PexPreser->show();

	
	hotDeadPixel = Gtk::manage(new Gtk::CheckButton((M("PREFERENCES_HOTDEADPIXFILT"))));

	lineDenoise = Gtk::manage(new Adjuster (M("PREFERENCES_LINEDENOISE"),0,1000,1,0));
	lineDenoise->setAdjusterListener (this);
	lineDenoise->show();

	greenEqThreshold = Gtk::manage(new Adjuster (M("PREFERENCES_GREENEQUIL"),0,100,1,0));
	greenEqThreshold->setAdjusterListener (this);
	greenEqThreshold->show();

    pack_start( *hbdf, Gtk::PACK_SHRINK, 4);
    pack_start( *dfAuto, Gtk::PACK_SHRINK, 4);
    pack_start( *Gtk::manage (new  Gtk::HSeparator()));
    pack_start( *hotDeadPixel, Gtk::PACK_SHRINK, 4);
    pack_start( *Gtk::manage (new  Gtk::HSeparator()));
    pack_start( *caAutocorrect, Gtk::PACK_SHRINK, 4);
	pack_start( *caRed, Gtk::PACK_SHRINK, 4);
    pack_start( *caBlue, Gtk::PACK_SHRINK, 4);
	pack_start( *PexPos, Gtk::PACK_SHRINK, 4);//exposi
    pack_start( *PexPreser, Gtk::PACK_SHRINK, 4);
	
    pack_start( *Gtk::manage (new  Gtk::HSeparator()));
    pack_start( *lineDenoise, Gtk::PACK_SHRINK, 4);
    pack_start( *Gtk::manage (new  Gtk::HSeparator()));
    pack_start( *greenEqThreshold, Gtk::PACK_SHRINK, 4);

   caacsconn = caAutocorrect->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::caCorrectionChanged), true);
   dfautoconn = dfAuto->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::dfAutoChanged), true);
   hdpixelconn = hotDeadPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::hotDeadPixelChanged), true);
   dfFile = darkFrameFile->signal_file_set().connect ( sigc::mem_fun(*this, &PreProcess::darkFrameChanged), true);
   btnReset->signal_clicked().connect( sigc::mem_fun(*this, &PreProcess::darkFrameReset), true );
}


void PreProcess::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
   disableListener ();
   caacsconn.block (true);
   dfautoconn.block(true);
   hdpixelconn.block (true);

   if(pedited ){
	   dfAuto->set_inconsistent(!pedited->raw.dfAuto );
	   caAutocorrect->set_inconsistent(!pedited->raw.caCorrection);
	   caRed->setEditedState( pedited->raw.caRed ? Edited : UnEdited );
	   caBlue->setEditedState( pedited->raw.caBlue ? Edited : UnEdited );
	   PexPos->setEditedState( pedited->raw.exPos ? Edited : UnEdited );
	   PexPreser->setEditedState( pedited->raw.exPreser ? Edited : UnEdited );
	   
	   //exposure
	   hotDeadPixel->set_inconsistent (!pedited->raw.hotDeadPixel);
	   lineDenoise->setEditedState( pedited->raw.linenoise ? Edited : UnEdited );
	   greenEqThreshold->setEditedState( pedited->raw.greenEq ? Edited : UnEdited );
   }

   if (Glib::file_test (pp->raw.dark_frame, Glib::FILE_TEST_EXISTS))
      darkFrameFile->set_filename (pp->raw.dark_frame);
   else if( !options.rtSettings.darkFramesPath.empty() )
	   darkFrameFile->set_current_folder( options.rtSettings.darkFramesPath );

   lastCA  = pp->raw.ca_autocorrect;
   lastHot = pp->raw.hotdeadpix_filt;
   lastDFauto = pp->raw.df_autoselect;

   dfAuto->set_active( pp->raw.df_autoselect );
   caAutocorrect->set_active(pp->raw.ca_autocorrect);
	caRed->setValue (pp->raw.cared);
	caBlue->setValue (pp->raw.cablue);
	PexPos->setValue (pp->raw.expos);
	PexPreser->setValue (pp->raw.preser);//exposi
	
   hotDeadPixel->set_active (pp->raw.hotdeadpix_filt);
   lineDenoise->setValue (pp->raw.linenoise);
   greenEqThreshold->setValue (pp->raw.greenthresh);

   dfChanged = false;


   caacsconn.block (false);
   dfautoconn.block(false);
   hdpixelconn.block (false);

   enableListener ();
}

void PreProcess::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
	pp->raw.dark_frame = darkFrameFile->get_filename();
	pp->raw.df_autoselect = dfAuto->get_active();
	pp->raw.ca_autocorrect = caAutocorrect->get_active();
	pp->raw.cared = (double)caRed->getValue();
	pp->raw.cablue = (double)caBlue->getValue();
	pp->raw.expos = (double)PexPos->getValue();
	pp->raw.preser = (double)PexPreser->getValue();//exposi
	
	pp->raw.hotdeadpix_filt = hotDeadPixel->get_active();
	pp->raw.linenoise = (int)lineDenoise->getValue();
	pp->raw.greenthresh = (int)greenEqThreshold->getValue();

	if (pedited) {
		pedited->raw.darkFrame = dfChanged;
		pedited->raw.dfAuto = !dfAuto->get_inconsistent();
		pedited->raw.linenoise = lineDenoise->getEditedState ();
		pedited->raw.greenEq= greenEqThreshold->getEditedState ();
		pedited->raw.caCorrection = !caAutocorrect->get_inconsistent();
		pedited->raw.caRed = caRed->getEditedState ();
		pedited->raw.caBlue = caBlue->getEditedState ();
		pedited->raw.exPos = PexPos->getEditedState ();
		pedited->raw.exPreser = PexPreser->getEditedState ();//exposi
		
		pedited->raw.hotDeadPixel = !hotDeadPixel->get_inconsistent();
	}
}

void PreProcess::adjusterChanged (Adjuster* a, double newval)
{
    if (listener)
        listener->panelChanged (EvPreProcess,  Glib::ustring("params") );
}

void PreProcess::setBatchMode(bool batchMode)
{
   ToolPanel::setBatchMode (batchMode);
	caRed->showEditedCB ();
	caBlue->showEditedCB ();
	PexPos->showEditedCB ();
	PexPreser->showEditedCB ();//exposi
	
   lineDenoise->showEditedCB ();
   greenEqThreshold->showEditedCB ();
}

void PreProcess::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
	lineDenoise->setDefault( defParams->raw.linenoise);
	caRed->setDefault( defParams->raw.cared);
	caBlue->setDefault( defParams->raw.cablue);
	PexPos->setDefault( defParams->raw.expos);
	PexPreser->setDefault( defParams->raw.preser);
	
	greenEqThreshold->setDefault (defParams->raw.greenthresh);
	if (pedited) {
		lineDenoise->setDefaultEditedState( pedited->raw.linenoise ? Edited : UnEdited);
		caRed->setDefaultEditedState( pedited->raw.caRed ? Edited : UnEdited);
		caBlue->setDefaultEditedState( pedited->raw.caBlue ? Edited : UnEdited);
		PexPos->setDefaultEditedState( pedited->raw.exPos ? Edited : UnEdited);
		PexPreser->setDefaultEditedState( pedited->raw.exPreser ? Edited : UnEdited);
		
		greenEqThreshold->setDefaultEditedState(pedited->raw.greenEq ? Edited : UnEdited);
	}else{
		lineDenoise->setDefaultEditedState( Irrelevant );
		caRed->setDefaultEditedState( Irrelevant );
		caBlue->setDefaultEditedState( Irrelevant );
		PexPos->setDefaultEditedState( Irrelevant );
		PexPreser->setDefaultEditedState( Irrelevant );
		
		greenEqThreshold->setDefaultEditedState(Irrelevant );
	}
}

void PreProcess::caCorrectionChanged()
{
    if (batchMode) {
        if (caAutocorrect->get_inconsistent()) {
        	caAutocorrect->set_inconsistent (false);
        	caacsconn.block (true);
            caAutocorrect->set_active (false);
            caacsconn.block (false);
        }
        else if (lastCA)
        	caAutocorrect->set_inconsistent (true);

        lastCA = caAutocorrect->get_active ();
    }
    if (listener)
        listener->panelChanged (EvPreProcess, Glib::ustring(M("PREFERENCES_CACORRECTION"))+"="+(caAutocorrect->get_active()?"ON":"OFF") );
}

void PreProcess::dfAutoChanged()
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
	hbdf->set_sensitive( !dfAuto->get_active() );
    if (listener)
        listener->panelChanged (EvPreProcess, Glib::ustring(M("TP_PREPROCESS_DFAUTOSELECT"))+"="+(dfAuto->get_active()?"ON":"OFF") );
}

void PreProcess::hotDeadPixelChanged ()
{
    if (batchMode) {
        if (hotDeadPixel->get_inconsistent()) {
        	hotDeadPixel->set_inconsistent (false);
        	hdpixelconn.block (true);
        	hotDeadPixel->set_active (false);
        	hdpixelconn.block (false);
        }
        else if (lastHot)
        	hotDeadPixel->set_inconsistent (true);

        lastHot = hotDeadPixel->get_active ();
    }
    if (listener)
        listener->panelChanged (EvPreProcess, Glib::ustring(M("PREFERENCES_HOTDEADPIXFILT"))+"="+(hotDeadPixel->get_active()?"ON":"OFF") );
}
void PreProcess::darkFrameChanged()
{
	dfChanged=true;
    if (listener)
        listener->panelChanged (EvPreProcess, Glib::ustring(M("TP_PREPROCESS_DARKFRAME"))+"="+darkFrameFile->get_filename());
}

void PreProcess::darkFrameReset()
{
	dfChanged=true;
	darkFrameFile->set_current_name("");
	darkFrameFile->set_filename ("");
    if (listener)
        listener->panelChanged (EvPreProcess, Glib::ustring(M("TP_PREPROCESS_DARKFRAME"))+"=0" );

}
