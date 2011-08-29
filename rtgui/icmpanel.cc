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
#include <icmpanel.h>
#include <options.h>
#include <guiutils.h>
#include <safegtk.h>
#include <iccstore.h>

using namespace rtengine;
using namespace rtengine::procparams;

extern Options options;

Glib::ustring ICMPanel::lastICCWorkDir;

ICMPanel::ICMPanel () : Gtk::VBox(), FoldableToolPanel(this), iunchanged(NULL), icmplistener(NULL) {

//    set_border_width (4);

    ipDialog = Gtk::manage (new Gtk::FileChooserButton (M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));

    Gtk::Label* ilab = Gtk::manage (new Gtk::Label ());
    ilab->set_alignment (0.0, 0.5);
    ilab->set_markup (Glib::ustring("<b>") + M("TP_ICM_INPUTPROFILE") + "</b>");
    pack_start (*ilab, Gtk::PACK_SHRINK, 4);
    
    inone = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTNONE")));
    pack_start (*inone, Gtk::PACK_SHRINK, 4);

    iembedded = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTEMBEDDED")));
    pack_start (*iembedded, Gtk::PACK_SHRINK, 4);

    icamera = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCAMERA")));
    pack_start (*icamera, Gtk::PACK_SHRINK, 4);

    ifromfile = Gtk::manage (new Gtk::RadioButton (M("TP_ICM_INPUTCUSTOM")+":"));
    Gtk::HBox* ffbox = Gtk::manage (new Gtk::HBox ());
    ffbox->pack_start (*ifromfile, Gtk::PACK_SHRINK);
    ffbox->pack_start (*ipDialog);

    pack_start (*ffbox, Gtk::PACK_SHRINK, 4);

    opts = icamera->get_group();
    iembedded->set_group (opts);
    ifromfile->set_group (opts);
    inone->set_group (opts);

    igamma = Gtk::manage (new Gtk::CheckButton (M("TP_ICM_GAMMABEFOREINPUT")));
    igamma->set_sensitive (false);
    pack_start (*igamma, Gtk::PACK_SHRINK, 4);

    saveRef = Gtk::manage (new Gtk::Button (M("TP_ICM_SAVEREFERENCE")));
    saveRef->set_image (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
    pack_start (*saveRef, Gtk::PACK_SHRINK, 4);


    Gtk::HSeparator* hsep1 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep1, Gtk::PACK_SHRINK, 2);

    Gtk::Label* wlab = Gtk::manage (new Gtk::Label ());
    wlab->set_alignment (0.0, 0.5);
    wlab->set_markup (Glib::ustring("<b>") + M("TP_ICM_WORKINGPROFILE") + "</b>");
    
    pack_start (*wlab, Gtk::PACK_SHRINK, 4);
    wnames = Gtk::manage (new Gtk::ComboBoxText ());    
    pack_start (*wnames, Gtk::PACK_SHRINK, 4);

    Gtk::HSeparator* hsep2 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep2, Gtk::PACK_SHRINK, 2);

    Gtk::Label* olab = Gtk::manage (new Gtk::Label ());
    olab->set_alignment (0.0, 0.5);
    olab->set_markup (Glib::ustring("<b>") + M("TP_ICM_OUTPUTPROFILE") + "</b>");

    pack_start (*olab, Gtk::PACK_SHRINK, 4);
    onames = Gtk::manage (new Gtk::ComboBoxText ());    
    pack_start (*onames, Gtk::PACK_SHRINK, 4);

    std::vector<std::string> wpnames = rtengine::getWorkingProfiles ();
    for (int i=0; i<wpnames.size(); i++)
        wnames->append_text (wpnames[i]);
  
 
	Gtk::HSeparator* hsep22 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep22, Gtk::PACK_SHRINK, 2);
	
    Gtk::Label* galab = Gtk::manage (new Gtk::Label ());
    galab->set_alignment (0.0, 0.5);
    galab->set_markup (Glib::ustring("<b>") + M("TP_GAMMA_OUTPUT") + "</b>");
	
    pack_start (*galab, Gtk::PACK_SHRINK, 4);
    wgamma = Gtk::manage (new Gtk::ComboBoxText ());    
    pack_start (*wgamma, Gtk::PACK_SHRINK, 4);
		
	Gtk::HSeparator* hsep23 = Gtk::manage (new Gtk::HSeparator ());
    pack_start (*hsep23, Gtk::PACK_SHRINK, 2);
	
	freegamma = Gtk::manage(new Gtk::CheckButton((M("TP_GAMMA_FREE"))));	
	freegamma->set_active (false);
	pack_start( *freegamma);
	
	g_ampos = Gtk::manage(new Adjuster (M("TP_GAMMA_CURV"),1,3.5,0.01,2.22));
	g_ampos->setAdjusterListener (this);
	if (g_ampos->delay < 1000) g_ampos->delay = 1000;
	g_ampos->show();
	s_lpos = Gtk::manage(new Adjuster (M("TP_GAMMA_SLOP"),0,15,0.01,4.5));
	s_lpos->setAdjusterListener (this); 
	if (s_lpos->delay < 1000) s_lpos->delay = 1000;
	s_lpos->show();
	pack_start( *g_ampos, Gtk::PACK_SHRINK, 4);//gamma
	pack_start( *s_lpos, Gtk::PACK_SHRINK, 4);//slope
		
	gamcsconn = freegamma->signal_toggled().connect ( sigc::mem_fun(*this, &ICMPanel::GamChanged));

		
		
    std::vector<std::string> wpgamma = rtengine::getGamma ();
    for (int i=0; i<wpgamma.size(); i++)
        wgamma->append_text (wpgamma[i]);

    onames->append_text (M("TP_ICM_NOICM"));
    onames->set_active (0);

    std::vector<std::string> opnames = iccStore->getOutputProfiles ();
    for (int i=0; i<opnames.size(); i++)
        onames->append_text (opnames[i]);

    wnames->set_active (0);
    onames->set_active (0);
	wgamma->set_active (0);

    Gtk::FileFilter filter_icc;
    filter_icc.set_name(M("TP_ICM_FILEDLGFILTERICM"));
    filter_icc.add_pattern("*.icc");
    filter_icc.add_pattern("*.icm");
    filter_icc.add_pattern("*.ICC");
    filter_icc.add_pattern("*.ICM");
    Gtk::FileFilter filter_any;
    filter_any.set_name(M("TP_ICM_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");

    ipDialog->add_filter (filter_icc);
    ipDialog->add_filter (filter_any);

    ipDialog->set_current_folder ( lastICCWorkDir.empty() ? options.rtSettings.iccDirectory : lastICCWorkDir);

    oldip = "";

    wnames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::wpChanged) );
    onames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::opChanged) );
    wgamma->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::gpChanged) );
	
    icamera->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    iembedded->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    ifromfile->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    igamma->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::profAppGammaChanged) );
    ipc = ipDialog->signal_selection_changed().connect( sigc::mem_fun(*this, &ICMPanel::ipSelectionChanged) );
    saveRef->signal_pressed().connect( sigc::mem_fun(*this, &ICMPanel::saveReferencePressed) );

    show_all ();
}

void ICMPanel::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    ipc.block (true);
    if (pp->icm.input == "(none)" && icamera->get_state()!=Gtk::STATE_INSENSITIVE) {
        inone->set_active (true);
        igamma->set_sensitive (false);
    }
    else if (pp->icm.input == "(embedded)" || ((pp->icm.input == "(camera)" || pp->icm.input=="") && icamera->get_state()==Gtk::STATE_INSENSITIVE)) {
        iembedded->set_active (true);
        igamma->set_sensitive (false);
    }
    else if ((pp->icm.input == "(camera)" || pp->icm.input=="") && icamera->get_state()!=Gtk::STATE_INSENSITIVE) {
        icamera->set_active (true);
        igamma->set_sensitive (false);
    }
    else {
        ifromfile->set_active (true);
        oldip = pp->icm.input.substr(5);
        ipDialog->set_filename (pp->icm.input.substr(5));
        igamma->set_sensitive (true);
    }

    wnames->set_active_text (pp->icm.working);   
    wgamma->set_active_text (pp->icm.gamma);    
	
    if (pp->icm.output==ColorManagementParams::NoICMString)
        onames->set_active_text (M("TP_ICM_NOICM"));
    else
        onames->set_active_text (pp->icm.output);

    if (onames->get_active_row_number()==-1)
        onames->set_active_text (M("TP_ICM_NOICM"));

    igamma->set_active (pp->icm.gammaOnInput);
	onames->set_sensitive(wgamma->get_active_row_number()==0 || freegamma->get_active()); //"default"
	wgamma->set_sensitive(!freegamma->get_active());
	
    if (pedited) {
        iunchanged->set_active (!pedited->icm.input);
        igamma->set_sensitive (false);
        if (!pedited->icm.working)
            wnames->set_active_text(M("GENERAL_UNCHANGED"));
        if (!pedited->icm.output)
            onames->set_active_text(M("GENERAL_UNCHANGED"));
        if (!pedited->icm.gamma){
            wgamma->set_active_text(M("GENERAL_UNCHANGED"));
            wgamma->set_active_text(M("GENERAL_UNCHANGED"));
			}
	    g_ampos->setEditedState      (pedited->icm.gampos ? Edited : UnEdited);
        s_lpos->setEditedState  	 (pedited->icm.slpos ? Edited : UnEdited);
		
    }
	
	gamcsconn.block (true);
    freegamma->set_active (pp->icm.freegamma);
    gamcsconn.block (false);

	lastgamfree = pp->icm.freegamma;
	g_ampos->setValue (pp->icm.gampos);
	s_lpos->setValue (pp->icm.slpos);
       
       
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
    else {
        pp->icm.input = "file:"+ipDialog->get_filename ();

        Glib::ustring p=Glib::path_get_dirname(ipDialog->get_filename ());
        if (p!=options.rtSettings.iccDirectory) {
            lastICCWorkDir=p;
        }
     }

    pp->icm.working = wnames->get_active_text ();
    pp->icm.gamma = wgamma->get_active_text ();
   
    if (onames->get_active_text()==M("TP_ICM_NOICM"))
        pp->icm.output  = ColorManagementParams::NoICMString;
    else
        pp->icm.output  = onames->get_active_text();
		pp->icm.freegamma = freegamma->get_active();
    pp->icm.gammaOnInput = igamma->get_active ();
	pp->icm.gampos =(double) g_ampos->getValue();
	pp->icm.slpos =(double) s_lpos->getValue();
	
    if (pedited) {
        pedited->icm.input = !iunchanged->get_active ();
        pedited->icm.working = wnames->get_active_text()!=M("GENERAL_UNCHANGED");
        pedited->icm.output = onames->get_active_text()!=M("GENERAL_UNCHANGED");
        pedited->icm.gammaOnInput = !ifromfile->get_active ();
        pedited->icm.gamma = wgamma->get_active_text()!=M("GENERAL_UNCHANGED");
		pedited->icm.freegamma =!freegamma->get_inconsistent();
        pedited->icm.gampos          = g_ampos->getEditedState ();
        pedited->icm.slpos   		 = s_lpos->getEditedState ();
		
    }
}
void ICMPanel::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {
  g_ampos->setDefault (defParams->icm.gampos);
  s_lpos->setDefault (defParams->icm.slpos);
 
  if (pedited) {
          g_ampos->setDefaultEditedState (pedited->icm.gampos ? Edited : UnEdited);
          s_lpos->setDefaultEditedState (pedited->icm.slpos ? Edited : UnEdited);

  }
  else {
          g_ampos->setDefaultEditedState (Irrelevant);
          s_lpos->setDefaultEditedState (Irrelevant);

  }
  }
void ICMPanel::adjusterChanged (Adjuster* a, double newval) {

    if (listener && freegamma->get_active()) {

        Glib::ustring costr = Glib::ustring::format ((int)a->getValue());

        if (a==g_ampos) 
            listener->panelChanged (EvGAMPOS, costr);
		else if (a==s_lpos)
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

void ICMPanel::ipChanged () {

    std::string profname;
    if (inone->get_active()) {
        profname = "(none)";
        igamma->set_sensitive (false);
    }
    else if (iembedded->get_active ()) {
        profname = "(embedded)";
        igamma->set_sensitive (false);
    }
    else if (icamera->get_active ()) {
        profname = "(camera)";
        igamma->set_sensitive (false);
    }
    else {
        profname = ipDialog->get_filename ();
        igamma->set_sensitive (true);
    }
    
    if (listener && profname!=oldip)
        listener->panelChanged (EvIProfile, profname);

    oldip = profname;
}

void ICMPanel::profAppGammaChanged() {
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

void ICMPanel::setRaw (bool raw) {

    disableListener ();

    icamera->set_active (raw);
    iembedded->set_active (!raw);
    icamera->set_sensitive (raw);
    iembedded->set_sensitive (!raw);  

    enableListener ();
}

void ICMPanel::ipSelectionChanged () {

    if (ipDialog->get_filename () == "")
        return;

        ipChanged ();
}

void ICMPanel::saveReferencePressed () {

    if (!icmplistener)
        return;
    Gtk::FileChooserDialog dialog(M("TP_ICM_SAVEREFERENCEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_tif;
    filter_tif.set_name(M("SAVEDLG_TIFFFILTER"));
    filter_tif.add_pattern("*.tif");
    dialog.add_filter(filter_tif);

    dialog.set_do_overwrite_confirmation (true);

    if (dialog.run()==Gtk::RESPONSE_OK) 
        icmplistener->saveInputICCReference (dialog.get_filename());
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
	g_ampos->showEditedCB ();
	s_lpos->showEditedCB ();


}

