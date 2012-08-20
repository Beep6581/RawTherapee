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
#include <iomanip>
#include "icmpanel.h"
#include "options.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include "../rtengine/iccstore.h"
#include "../rtengine/dcp.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

extern Options options;

ICMPanel::ICMPanel () : Gtk::VBox(), FoldableToolPanel(this), iunchanged(NULL), icmplistener(NULL), lastRefFilename("") {

//    set_border_width (4);

    ipDialog = Gtk::manage (new MyFileChooserButton (M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    ipDialog->set_tooltip_text (M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    ipDialogPersister.reset(new FileChooserLastFolderPersister(ipDialog, options.lastIccDir));

    Gtk::Label* ilab = Gtk::manage (new Gtk::Label ());
    ilab->set_alignment (0.0, 0.5);
    ilab->set_markup (Glib::ustring("<b>") + M("TP_ICM_INPUTPROFILE") + "</b>");
    pack_start (*ilab, Gtk::PACK_SHRINK, 4);
    
    inone = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTNONE")));
    inone->set_tooltip_text (M("TP_ICM_INPUTNONE_TOOLTIP"));
    pack_start (*inone, Gtk::PACK_SHRINK, 4);

    iembedded = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTEMBEDDED")));
    iembedded->set_tooltip_text (M("TP_ICM_INPUTEMBEDDED_TOOLTIP"));
    pack_start (*iembedded, Gtk::PACK_SHRINK, 4);

    icamera = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCAMERA")));
    icamera->set_tooltip_text (M("TP_ICM_INPUTCAMERA_TOOLTIP"));
    pack_start (*icamera, Gtk::PACK_SHRINK, 4);

    icameraICC = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCAMERAICC")));
    icameraICC->set_tooltip_text (M("TP_ICM_INPUTCAMERAICC_TOOLTIP"));
    pack_start (*icameraICC, Gtk::PACK_SHRINK, 4);

    ifromfile = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCUSTOM")+":"));
    Gtk::HBox* ffbox = Gtk::manage (new Gtk::HBox ());
    ifromfile->set_tooltip_text (M("TP_ICM_INPUTCUSTOM_TOOLTIP"));
    ffbox->pack_start (*ifromfile, Gtk::PACK_SHRINK);
    ffbox->pack_start (*ipDialog);

    pack_start (*ffbox, Gtk::PACK_SHRINK, 4);

    opts = icamera->get_group();
    icameraICC->set_group (opts);
    iembedded->set_group (opts);
    ifromfile->set_group (opts);
    inone->set_group (opts);

    Gtk::HBox* hb = Gtk::manage (new Gtk::HBox ());
    hb->show ();
    Gtk::Label* ppl = Gtk::manage (new Gtk::Label (M("TP_ICM_PREFERREDPROFILE")+":"));
    ppl->show ();
    prefprof = Gtk::manage (new MyComboBoxText ());
    prefprof->append_text (M("TP_ICM_PREFERREDPROFILE_1"));
    prefprof->append_text (M("TP_ICM_PREFERREDPROFILE_2"));
    prefprof->append_text (M("TP_ICM_PREFERREDPROFILE_3"));
    prefprof->append_text (M("TP_ICM_PREFERREDPROFILE_4"));
    prefprof->show ();
    hb->pack_start(*ppl, Gtk::PACK_SHRINK, 4);
    hb->pack_start(*prefprof);  
    pack_start (*hb, Gtk::PACK_SHRINK, 4);

    ckbToneCurve = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_TONECURVE")));
    ckbToneCurve->set_sensitive (false);
    ckbToneCurve->set_tooltip_text (M("TP_ICM_TONECURVE_TOOLTIP"));
    pack_start (*ckbToneCurve, Gtk::PACK_SHRINK, 4);

    ckbBlendCMSMatrix = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_BLENDCMSMATRIX")));
    ckbBlendCMSMatrix->set_sensitive (false);
    ckbBlendCMSMatrix->set_tooltip_text (M("TP_ICM_BLENDCMSMATRIX_TOOLTIP"));
    pack_start (*ckbBlendCMSMatrix, Gtk::PACK_SHRINK, 4);

    saveRef = Gtk::manage (new Gtk::Button (M("TP_ICM_SAVEREFERENCE")));
    saveRef->set_image (*Gtk::manage (new RTImage ("gtk-save-large.png")));
    pack_start (*saveRef, Gtk::PACK_SHRINK, 4);


    Gtk::HSeparator* hsep1 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep1, Gtk::PACK_SHRINK, 2);

    Gtk::Label* wlab = Gtk::manage (new Gtk::Label ());
    wlab->set_alignment (0.0, 0.5);
    wlab->set_markup (Glib::ustring("<b>") + M("TP_ICM_WORKINGPROFILE") + "</b>");
    
    pack_start (*wlab, Gtk::PACK_SHRINK, 4);
    wnames = Gtk::manage (new MyComboBoxText ());
    pack_start (*wnames, Gtk::PACK_SHRINK, 4);

    Gtk::HSeparator* hsep2 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep2, Gtk::PACK_SHRINK, 2);

    Gtk::Label* olab = Gtk::manage (new Gtk::Label ());
    olab->set_alignment (0.0, 0.5);
    olab->set_markup (Glib::ustring("<b>") + M("TP_ICM_OUTPUTPROFILE") + "</b>");

    pack_start (*olab, Gtk::PACK_SHRINK, 4);
    onames = Gtk::manage (new MyComboBoxText ());
    pack_start (*onames, Gtk::PACK_SHRINK, 4);

    std::vector<std::string> wpnames = rtengine::getWorkingProfiles ();
    for (size_t i=0; i<wpnames.size(); i++)
        wnames->append_text (wpnames[i]);
  
 
	Gtk::HSeparator* hsep22 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep22, Gtk::PACK_SHRINK, 2);
	
    Gtk::Label* galab = Gtk::manage (new Gtk::Label ());
    galab->set_alignment (0.0, 0.5);
    galab->set_markup (Glib::ustring("<b>") + M("TP_GAMMA_OUTPUT") + "</b>");
	
    pack_start (*galab, Gtk::PACK_SHRINK, 4);
    wgamma = Gtk::manage (new MyComboBoxText ());
    pack_start (*wgamma, Gtk::PACK_SHRINK, 4);
		
	Gtk::HSeparator* hsep23 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep23, Gtk::PACK_SHRINK, 2);
	
	freegamma = Gtk::manage(new Gtk::CheckButton((M("TP_GAMMA_FREE"))));	
	freegamma->set_active (false);
	pack_start( *freegamma);
	
	gampos = Gtk::manage(new Adjuster (M("TP_GAMMA_CURV"),1,3.5,0.01,2.22));
	gampos->setAdjusterListener (this);
	if (gampos->delay < 1000) gampos->delay = 1000;
	gampos->show();
	slpos = Gtk::manage(new Adjuster (M("TP_GAMMA_SLOP"),0,15,0.01,4.5));
	slpos->setAdjusterListener (this); 
	if (slpos->delay < 1000) slpos->delay = 1000;
	slpos->show();
	pack_start( *gampos, Gtk::PACK_SHRINK, 4);//gamma
	pack_start( *slpos, Gtk::PACK_SHRINK, 4);//slope
		
	gamcsconn = freegamma->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::GamChanged));

		
		
    std::vector<std::string> wpgamma = rtengine::getGamma ();
    for (size_t i=0; i<wpgamma.size(); i++)
        wgamma->append_text (wpgamma[i]);

    onames->append_text (M("TP_ICM_NOICM"));
    onames->set_active (0);

    std::vector<std::string> opnames = iccStore->getOutputProfiles ();
    for (size_t i=0; i<opnames.size(); i++)
        onames->append_text (opnames[i]);

    wnames->set_active (0);
    onames->set_active (0);
	wgamma->set_active (0);

    Gtk::FileFilter filter_icc;
    filter_icc.set_name(M("TP_ICM_FILEDLGFILTERICM"));
    filter_icc.add_pattern("*.dcp");
    filter_icc.add_pattern("*.DCP");
    filter_icc.add_pattern("*.dng");
    filter_icc.add_pattern("*.DNG");
    filter_icc.add_pattern("*.icc");
    filter_icc.add_pattern("*.icm");
    filter_icc.add_pattern("*.ICC");
    filter_icc.add_pattern("*.ICM");
    Gtk::FileFilter filter_any;
    filter_any.set_name(M("TP_ICM_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");

    ipDialog->add_filter (filter_icc);
    ipDialog->add_filter (filter_any);

    oldip = "";

    wnames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::wpChanged) );
    onames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::opChanged) );
    wgamma->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::gpChanged) );
    prefprof->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::prefProfChanged) );
	
    icamera->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    icameraICC->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    iembedded->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    ifromfile->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    ckbBlendCMSMatrix->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::iccTogglesChanged) );
    ckbToneCurve->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::toneCurveChanged) );

    ipc = ipDialog->signal_selection_changed().connect( sigc::mem_fun(*this, &ICMPanel::ipSelectionChanged) );
    saveRef->signal_pressed().connect( sigc::mem_fun(*this, &ICMPanel::saveReferencePressed) );

    show_all ();
}

void ICMPanel::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    ipc.block (true);
    if (pp->icm.input == "(none)" && icamera->get_state()!=Gtk::STATE_INSENSITIVE) {
        inone->set_active (true); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
        ckbBlendCMSMatrix->set_sensitive (false);
    }
    else if (pp->icm.input == "(embedded)" || ((pp->icm.input == "(camera)" || pp->icm.input=="") && icamera->get_state()==Gtk::STATE_INSENSITIVE)) {
        iembedded->set_active (true); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
        ckbBlendCMSMatrix->set_sensitive (false);
    }
    else if ((pp->icm.input == "(cameraICC)") && icameraICC->get_state()!=Gtk::STATE_INSENSITIVE) {
        icameraICC->set_active (true); prefprof->set_sensitive (true); ckbToneCurve->set_sensitive (true);
        ckbBlendCMSMatrix->set_sensitive (true);
    }
    else if ((pp->icm.input == "(cameraICC)") && icameraICC->get_state()==Gtk::STATE_INSENSITIVE) {
    	// this is the case when (cameraICC) is instructed by packaged profiles, but ICC file is not found
    	// therefore falling back UI to explicitly reflect the (camera) option
    	icamera->set_active (true);
        prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);   // RT's own are always single-illuminant and tone curve disabled
    	ckbBlendCMSMatrix->set_sensitive (false);
    }
    else if ((pp->icm.input == "(camera)" || pp->icm.input=="") && icamera->get_state()!=Gtk::STATE_INSENSITIVE) {
        icamera->set_active (true);
        ckbBlendCMSMatrix->set_sensitive (false); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
    }
    else {
        ifromfile->set_active (true);
        oldip = pp->icm.input.substr(5);  // cut of "file:"
        ipDialog->set_filename (pp->icm.input.substr(5));
        ckbBlendCMSMatrix->set_sensitive (true);  prefprof->set_sensitive (true); ckbToneCurve->set_sensitive (true);
    }

    wnames->set_active_text (pp->icm.working);   
    wgamma->set_active_text (pp->icm.gamma);    
	
    if (pp->icm.output==ColorManagementParams::NoICMString)
        onames->set_active_text (M("TP_ICM_NOICM"));
    else
        onames->set_active_text (pp->icm.output);

    if (onames->get_active_row_number()==-1)
        onames->set_active_text (M("TP_ICM_NOICM"));

    prefprof->set_active(pp->icm.preferredProfile-1);

    ckbToneCurve->set_active (pp->icm.toneCurve);
    ckbBlendCMSMatrix->set_active (pp->icm.blendCMSMatrix);
	onames->set_sensitive(wgamma->get_active_row_number()==0 || freegamma->get_active()); //"default"
	wgamma->set_sensitive(!freegamma->get_active());
	
    if (pedited) {
        iunchanged->set_active (!pedited->icm.input);
        ckbToneCurve->set_sensitive (false);
        ckbBlendCMSMatrix->set_sensitive (false);
        if (!pedited->icm.working)
            wnames->set_active_text(M("GENERAL_UNCHANGED"));
        if (!pedited->icm.output)
            onames->set_active_text(M("GENERAL_UNCHANGED"));
        if (!pedited->icm.preferredProfile)
            prefprof->set_active_text(M("GENERAL_UNCHANGED"));
        if (!pedited->icm.gamma){
            wgamma->set_active_text(M("GENERAL_UNCHANGED"));
            wgamma->set_active_text(M("GENERAL_UNCHANGED"));
			}
	    gampos->setEditedState      (pedited->icm.gampos ? Edited : UnEdited);
        slpos->setEditedState  	 (pedited->icm.slpos ? Edited : UnEdited);
		
    }
	
	gamcsconn.block (true);
    freegamma->set_active (pp->icm.freegamma);
    gamcsconn.block (false);

	lastgamfree = pp->icm.freegamma;
	gampos->setValue (pp->icm.gampos);
	slpos->setValue (pp->icm.slpos);
       
       
    ipc.block (false);

    enableListener ();
}

void ICMPanel::write (ProcParams* pp, ParamsEdited* pedited) {

    if (inone->get_active())
        pp->icm.input = "(none)";
    else if (iembedded->get_active ())
        pp->icm.input = "(embedded)";
    else if (icamera->get_active ())
        pp->icm.input = "(camera)";
    else if (icameraICC->get_active ())
        pp->icm.input = "(cameraICC)";
    else {
        if (safe_file_test (ipDialog->get_filename (), Glib::FILE_TEST_EXISTS) && !safe_file_test (ipDialog->get_filename (), Glib::FILE_TEST_IS_DIR))
        pp->icm.input = "file:"+ipDialog->get_filename ();
        else
            pp->icm.input = "";  // just a directory

        Glib::ustring p=Glib::path_get_dirname(ipDialog->get_filename ());
     }

    pp->icm.working = wnames->get_active_text ();
    pp->icm.gamma = wgamma->get_active_text ();
    pp->icm.preferredProfile = prefprof->get_active_row_number()+1;
   
    if (onames->get_active_text()==M("TP_ICM_NOICM"))
        pp->icm.output  = ColorManagementParams::NoICMString;
    else
        pp->icm.output  = onames->get_active_text();
		pp->icm.freegamma = freegamma->get_active();
    pp->icm.toneCurve = ckbToneCurve->get_active ();
    pp->icm.blendCMSMatrix = ckbBlendCMSMatrix->get_active ();
	pp->icm.gampos =(double) gampos->getValue();
	pp->icm.slpos =(double) slpos->getValue();
	
    if (pedited) {
        pedited->icm.input = !iunchanged->get_active ();
        pedited->icm.working = wnames->get_active_text()!=M("GENERAL_UNCHANGED");
        pedited->icm.output = onames->get_active_text()!=M("GENERAL_UNCHANGED");
        pedited->icm.preferredProfile = prefprof->get_active_text()!=M("GENERAL_UNCHANGED");
        pedited->icm.toneCurve = !ckbToneCurve->get_inconsistent ();
        pedited->icm.blendCMSMatrix = !ckbBlendCMSMatrix->get_inconsistent ();
        pedited->icm.gamma = wgamma->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->icm.freegamma =!freegamma->get_inconsistent();
        pedited->icm.gampos          = gampos->getEditedState ();
        pedited->icm.slpos   		 = slpos->getEditedState ();
		
    }
}
void ICMPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {
  gampos->setDefault (defParams->icm.gampos);
  slpos->setDefault (defParams->icm.slpos);
 
  if (pedited) {
          gampos->setDefaultEditedState (pedited->icm.gampos ? Edited : UnEdited);
          slpos->setDefaultEditedState (pedited->icm.slpos ? Edited : UnEdited);

  }
  else {
          gampos->setDefaultEditedState (Irrelevant);
          slpos->setDefaultEditedState (Irrelevant);

  }
  }
  
 void ICMPanel::setAdjusterBehavior (bool gammaadd, bool slopeadd) {
	gampos->setAddMode (gammaadd);
	slpos->setAddMode (slopeadd);
}
 
  
void ICMPanel::adjusterChanged (Adjuster* a, double newval) {

    if (listener && freegamma->get_active()) {

        Glib::ustring costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), newval);

        if (a==gampos) 
            listener->panelChanged (EvGAMPOS, costr);
		else if (a==slpos)
            listener->panelChanged (EvSLPOS, costr);
		
    }
}

void ICMPanel::wpChanged () {

    if (listener)
        listener->panelChanged (EvWProfile, wnames->get_active_text ());
}

void ICMPanel::gpChanged () {

    if (listener)
        {listener->panelChanged (EvGAMMA, wgamma->get_active_text ());
		onames->set_sensitive(wgamma->get_active_row_number()==0); //"default"
		 }
}

void ICMPanel::prefProfChanged() {
    if (listener)
        listener->panelChanged (EvPrefProfile, prefprof->get_active_text ());
}

void ICMPanel::toneCurveChanged() {
    if (listener)
        listener->panelChanged (EvDCPToneCurve, "");
}

void ICMPanel::ipChanged () {

    std::string profname;
    if (inone->get_active()) {
        profname = "(none)";
        ckbBlendCMSMatrix->set_sensitive(false); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
    }
    else if (iembedded->get_active ()) {
        profname = "(embedded)";
        ckbBlendCMSMatrix->set_sensitive(false); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
    }
    else if (icamera->get_active ()) {
        profname = "(camera)";
        ckbBlendCMSMatrix->set_sensitive(false); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
    }
    else if (icameraICC->get_active ()) {
        profname = "(cameraICC)";
        ckbBlendCMSMatrix->set_sensitive(true); prefprof->set_sensitive (false); ckbToneCurve->set_sensitive (false);
    }
    else {
        profname = ipDialog->get_filename ();
        ckbBlendCMSMatrix->set_sensitive(true); prefprof->set_sensitive (true); ckbToneCurve->set_sensitive (true);
    }
    
    if (listener && profname!=oldip)
        listener->panelChanged (EvIProfile, profname);

    oldip = profname;
}

void ICMPanel::iccTogglesChanged() {
    if (listener) listener->panelChanged (EvIProfile, "");
}

void ICMPanel::GamChanged() {
    if (batchMode) {
        if (freegamma->get_inconsistent()) {
            freegamma->set_inconsistent (false);
            gamcsconn.block (true);
            freegamma->set_active (false);
            gamcsconn.block (false);
        }
        else if (lastgamfree)
            freegamma->set_inconsistent (true);

        lastgamfree = freegamma->get_active ();
    }
    
    if (listener) {
        if (freegamma->get_active()){
            listener->panelChanged (EvGAMFREE, M("GENERAL_ENABLED"));
			onames->set_sensitive(!freegamma->get_active());//disabled choice
			wgamma->set_sensitive(!freegamma->get_active());
		}
        else {
            listener->panelChanged (EvGAMFREE, M("GENERAL_DISABLED"));
		   	onames->set_sensitive(!freegamma->get_active() && wgamma->get_active_row_number()==0);
			wgamma->set_sensitive(!freegamma->get_active()); 
		}
    }
}

void ICMPanel::opChanged () {

    if (listener)
        listener->panelChanged (EvOProfile, onames->get_active_text());
}

void ICMPanel::setRawMeta (bool raw, const rtengine::ImageData* pMeta) {

    disableListener ();

    icamera->set_active (raw);
    iembedded->set_active (!raw);
    icamera->set_sensitive (raw);
    icameraICC->set_sensitive (raw && (iccStore->getStdProfile(pMeta->getCamera()) != NULL || dcpStore->getStdProfile(pMeta->getCamera()) != NULL));
    iembedded->set_sensitive (!raw);  

    enableListener ();
}

void ICMPanel::ipSelectionChanged() {

	if (ipDialog->get_filename() == "")
		return;

	ipChanged();
}

void ICMPanel::saveReferencePressed () {

    if (!icmplistener)
        return;
    Gtk::FileChooserDialog dialog(M("TP_ICM_SAVEREFERENCEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    FileChooserLastFolderPersister persister(&dialog, options.lastProfilingReferenceDir);
    dialog.set_current_name (lastRefFilename);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_tif;
    filter_tif.set_name(M("SAVEDLG_TIFFFILTER"));
    filter_tif.add_pattern("*.tif");
    filter_tif.add_pattern("*.tiff");
    dialog.add_filter(filter_tif);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("TP_ICM_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

    //dialog.set_do_overwrite_confirmation (true);

    bool done = false;
    do {
        int result = dialog.run();
        if (result != Gtk::RESPONSE_OK) {
            done = true;
        } else {
            std::string fname = dialog.get_filename();
            Glib::ustring ext = getExtension(fname);
            if (ext != "tif" && ext != "tiff")
                fname += ".tif";
            if (confirmOverwrite(dialog, fname)) {
                icmplistener->saveInputICCReference (fname);
                lastRefFilename = Glib::path_get_basename (fname);
                done = true;
            }
        }
    } while (!done);
    return;
}

void ICMPanel::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    iunchanged = Gtk::manage (new Gtk::RadioButton (M("GENERAL_UNCHANGED")));
    iunchanged->set_group (opts);
    pack_start (*iunchanged, Gtk::PACK_SHRINK, 4);
    reorder_child (*iunchanged, 5);
    removeIfThere (this, saveRef);
    onames->append_text (M("GENERAL_UNCHANGED"));
    wnames->append_text (M("GENERAL_UNCHANGED"));
	wgamma->append_text (M("GENERAL_UNCHANGED"));
    prefprof->append_text (M("GENERAL_UNCHANGED"));
	gampos->showEditedCB ();
	slpos->showEditedCB ();


}

