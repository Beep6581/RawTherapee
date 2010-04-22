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

using namespace rtengine;
using namespace rtengine::procparams;

extern Options options;

ICMPanel::ICMPanel () : ToolPanel(), icmplistener(NULL), iunchanged(NULL) {

//    set_border_width (4);

    ipDialog = Gtk::manage (new Gtk::FileChooserButton (M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    opDialog = Gtk::manage (new Gtk::FileChooserButton (M("TP_ICM_INPUTDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN));

    Gtk::Label* ilab = Gtk::manage (new Gtk::Label ());
    ilab->set_alignment (0.0, 0.5);
    ilab->set_markup (Glib::ustring("<b>") + M("TP_ICM_INPUTPROFILE") + "</b>");
    pack_start (*ilab, Gtk::PACK_SHRINK, 4);
    
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

    onames->append_text (M("TP_ICM_NOICM"));
    onames->set_active (0);

    std::vector<std::string> opnames = rtengine::getOutputProfiles ();
    for (int i=0; i<opnames.size(); i++)
        onames->append_text (opnames[i]);

    wnames->set_active (0);
    onames->set_active (0);

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
    opDialog->add_filter (filter_icc);
    opDialog->add_filter (filter_any);

    if (Glib::file_test (options.rtSettings.iccDirectory, Glib::FILE_TEST_IS_DIR)) {
        ipDialog->set_current_folder (options.rtSettings.iccDirectory);
        opDialog->set_current_folder (options.rtSettings.iccDirectory);
    }    

    oldip = "";

    wnames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::wpChanged) );
    onames->signal_changed().connect( sigc::mem_fun(*this, &ICMPanel::opChanged) );
    icamera->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    iembedded->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    ifromfile->signal_toggled().connect( sigc::mem_fun(*this, &ICMPanel::ipChanged) );
    ipc = ipDialog->signal_selection_changed().connect( sigc::mem_fun(*this, &ICMPanel::ipSelectionChanged) );
    saveRef->signal_pressed().connect( sigc::mem_fun(*this, &ICMPanel::saveReferencePressed) );

    show_all ();
}

void ICMPanel::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();

    ipc.block (true);
    if (pp->icm.input == "(embedded)" || ((pp->icm.input == "(camera)" || pp->icm.input=="") && icamera->get_state()==Gtk::STATE_INSENSITIVE)) {
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
    if (pp->icm.output=="No ICM: sRGB output")
        onames->set_active_text (M("TP_ICM_NOICM"));
    else
        onames->set_active_text (pp->icm.output);

    if (onames->get_active_row_number()==-1)
        onames->set_active_text (M("TP_ICM_NOICM"));

    igamma->set_active (pp->icm.gammaOnInput);

    if (pedited) {
        iunchanged->set_active (!pedited->icm.input);
        igamma->set_sensitive (false);
        if (!pedited->icm.working)
            wnames->set_active_text("(Unchanged)");
        if (!pedited->icm.output)
            onames->set_active_text("(Unchanged)");
    }
        
    ipc.block (false);

    enableListener ();
}

void ICMPanel::write (ProcParams* pp, ParamsEdited* pedited) {

    if (iembedded->get_active ())
        pp->icm.input = "(embedded)";
    else if (icamera->get_active ())
        pp->icm.input = "(camera)";
    else
        pp->icm.input = "file:"+ipDialog->get_filename ();

    pp->icm.working = wnames->get_active_text ();
    
    if (onames->get_active_text()==M("TP_ICM_NOICM"))
        pp->icm.output  = "No ICM: sRGB output";
    else
        pp->icm.output  = onames->get_active_text();
    pp->icm.gammaOnInput = igamma->get_active ();
    
    if (pedited) {
        pedited->icm.input = !iunchanged->get_active ();
        pedited->icm.working = wnames->get_active_text()!="(Unchanged)";
        pedited->icm.output = onames->get_active_text()!="(Unchanged)";
        pedited->icm.gammaOnInput = !ifromfile->get_active ();
    }
}

void ICMPanel::wpChanged () {

    if (listener)
        listener->panelChanged (EvWProfile, wnames->get_active_text ());
}

void ICMPanel::ipChanged () {

    std::string profname;
    if (iembedded->get_active ()) {
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
    else
        ipChanged ();

}
void ICMPanel::saveReferencePressed () {

    if (!icmplistener)
        return;
    Gtk::FileChooserDialog dialog(M("TP_ICM_SAVEREFERENCEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_jpg;
    filter_jpg.set_name(M("SAVEDLG_JPGFILTER"));
    filter_jpg.add_pattern("*.jpg");
    dialog.add_filter(filter_jpg);

    dialog.set_do_overwrite_confirmation (true);

    if (dialog.run()==Gtk::RESPONSE_OK) 
        icmplistener->saveInputICCReference (dialog.get_filename());

}
void ICMPanel::setBatchMode (bool batchMode) {

    ToolPanel::setBatchMode (batchMode);
    iunchanged = Gtk::manage (new Gtk::RadioButton ("(Unchanged)"));
    iunchanged->set_group (opts);
    pack_start (*iunchanged, Gtk::PACK_SHRINK, 4);
    reorder_child (*iunchanged, 5);
    removeIfThere (this, saveRef);
    onames->append_text ("(Unchanged)");
    wnames->append_text ("(Unchanged)");
}

