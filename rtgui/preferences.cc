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
#include <sigc++/class_slot.h>
#include "preferences.h"
#include <multilangmgr.h>
#include <splash.h>
#include <cachemanager.h>
#include <addsetids.h>
#include <dfmanager.h>
#include <sstream>

extern Options options;
extern Glib::ustring argv0;

Preferences::Preferences  (RTWindow *rtwindow):parent(rtwindow)  {
  
    set_title (M("MAIN_BUTTON_PREFERENCES"));

    moptions.copyFrom (&options);

    set_size_request (650, 650);
    set_border_width (4);

    Gtk::VBox* mainvb = get_vbox ();
    set_has_separator (false);

    Gtk::Notebook* nb = Gtk::manage (new Gtk::Notebook ());
    mainvb->pack_start (*nb);

    Gtk::HSeparator* hsep1 = Gtk::manage (new Gtk::HSeparator ());
    mainvb->pack_start (*hsep1, Gtk::PACK_SHRINK, 2);

    Gtk::HBox* buttonpanel = Gtk::manage (new Gtk::HBox ());
    mainvb->pack_end (*buttonpanel, Gtk::PACK_SHRINK, 2);

//    Gtk::Button* load   = Gtk::manage (new Gtk::Button (M("GENERAL_LOAD")));
//    Gtk::Button* save   = Gtk::manage (new Gtk::Button (M("GENERAL_SAVE")));
    Gtk::Button* about  = Gtk::manage (new Gtk::Button (M("GENERAL_ABOUT")));
    Gtk::Button* ok     = Gtk::manage (new Gtk::Button (M("GENERAL_OK")));
    Gtk::Button* cancel = Gtk::manage (new Gtk::Button (M("GENERAL_CANCEL")));

//    save->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
//    load->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));
    about->set_image (*Gtk::manage(new Gtk::Image (argv0+"/images/logoicon16.png")));
    ok->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-ok"), Gtk::ICON_SIZE_BUTTON)));
    cancel->set_image (*Gtk::manage(new Gtk::Image (Gtk::StockID("gtk-cancel"), Gtk::ICON_SIZE_BUTTON)));


//    load->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::loadPressed) );
//    save->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::savePressed) );
    about->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::aboutPressed) );
    ok->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::okPressed) );
    cancel->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::cancelPressed) );

//    buttonpanel->pack_start (*load, Gtk::PACK_SHRINK, 4);
//    buttonpanel->pack_start (*save, Gtk::PACK_SHRINK, 4);
    buttonpanel->pack_start (*about, Gtk::PACK_SHRINK, 4);
    buttonpanel->pack_end (*ok, Gtk::PACK_SHRINK, 4);
    buttonpanel->pack_end (*cancel, Gtk::PACK_SHRINK, 4);

    nb->append_page (*getGeneralPanel(),        M("PREFERENCES_TAB_GENERAL"));
    nb->append_page (*getProcParamsPanel(),     M("PREFERENCES_TAB_IMPROC"));
    nb->append_page (*getFileBrowserPanel(),    M("PREFERENCES_TAB_BROWSER"));
    nb->append_page (*getColorManagementPanel(),M("PREFERENCES_TAB_COLORMGR"));
    nb->append_page (*getBatchProcPanel(),      M("PREFERENCES_BATCH_PROCESSING"));
    nb->set_current_page (0);

    fillPreferences ();

    show_all_children ();
    set_modal (true);
}

Gtk::Widget* Preferences::getBatchProcPanel () {

    Gtk::VBox* mvbpp = Gtk::manage (new Gtk::VBox ());

    Gtk::ScrolledWindow* behscrollw = Gtk::manage (new Gtk::ScrolledWindow ());
    behscrollw->set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    Gtk::Frame* behFrame = Gtk::manage (new Gtk::Frame (M("PREFERENCES_BEHAVIOR")));
    behFrame->add (*behscrollw);
    mvbpp->pack_start (*behFrame);
//    mvbpp->pack_start (*behFrame, Gtk::PACK_SHRINK, 2);
    Gtk::TreeView* behTreeView = Gtk::manage (new Gtk::TreeView ());
    behscrollw->add (*behTreeView);

    behModel = Gtk::TreeStore::create (behavColumns);
    behTreeView->set_model (behModel);
    
    behTreeView->append_column (M("PREFERENCES_PROPERTY"), behavColumns.label); 
    behTreeView->append_column_editable (M("PREFERENCES_ADD"), behavColumns.badd); 
    behTreeView->append_column_editable (M("PREFERENCES_SET"), behavColumns.bset); 
    
    Gtk::CellRendererToggle* cr_add = static_cast<Gtk::CellRendererToggle*> (behTreeView->get_column (1)->get_first_cell_renderer());
    Gtk::CellRendererToggle* cr_set = static_cast<Gtk::CellRendererToggle*> (behTreeView->get_column (2)->get_first_cell_renderer());

    cr_add->set_radio (true);
    cr_add->set_property("xalign", 0.0f);
    sigc::connection addc = cr_add->signal_toggled().connect (sigc::mem_fun (*this, &Preferences::behAddRadioToggled));
    cr_set->set_radio (true);
    cr_set->set_property("xalign", 0.0f);
    sigc::connection setc = cr_set->signal_toggled().connect (sigc::mem_fun (*this, &Preferences::behSetRadioToggled));

    behTreeView->get_column (1)->add_attribute (*cr_add, "visible", behavColumns.visible);
    behTreeView->get_column (1)->set_sizing(Gtk::TREE_VIEW_COLUMN_FIXED);
    behTreeView->get_column (1)->set_fixed_width (50);
    behTreeView->get_column (2)->add_attribute (*cr_set, "visible", behavColumns.visible);
    behTreeView->get_column (2)->set_sizing(Gtk::TREE_VIEW_COLUMN_FIXED);
    behTreeView->get_column (2)->set_fixed_width (50);

    // fill model
    Gtk::TreeModel::iterator mi, ci;

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_EXPOSURE_LABEL"));
    appendBehavList (mi, M("TP_EXPOSURE_EXPCOMP"), ADDSET_TC_EXPCOMP, false);
    appendBehavList (mi, M("TP_EXPOSURE_BRIGHTNESS"), ADDSET_TC_BRIGHTNESS, false);
    appendBehavList (mi, M("TP_EXPOSURE_BLACKLEVEL"), ADDSET_TC_BLACKLEVEL, false);
    appendBehavList (mi, M("TP_EXPOSURE_CONTRAST"), ADDSET_TC_CONTRAST, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_SHADOWSHLIGHTS_LABEL"));
    appendBehavList (mi, M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), ADDSET_SH_HIGHLIGHTS, false);
    appendBehavList (mi, M("TP_SHADOWSHLIGHTS_SHADOWS"), ADDSET_SH_SHADOWS, false);
    appendBehavList (mi, M("TP_SHADOWSHLIGHTS_LOCALCONTR"), ADDSET_SH_LOCALCONTRAST, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_LABCURVE_LABEL"));
    appendBehavList (mi, M("TP_LABCURVE_BRIGHTNESS"), ADDSET_LC_BRIGHTNESS, false);
    appendBehavList (mi, M("TP_LABCURVE_CONTRAST"), ADDSET_LC_CONTRAST, false);
	appendBehavList (mi, M("TP_LABCURVE_SATURATION"), ADDSET_LC_SATURATION, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_SHARPENING_LABEL"));
    appendBehavList (mi, M("TP_SHARPENING_AMOUNT"), ADDSET_SHARP_AMOUNT, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_LUMADENOISE_LABEL"));
    appendBehavList (mi, M("TP_LUMADENOISE_EDGETOLERANCE"), ADDSET_LD_EDGETOLERANCE, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_WBALANCE_LABEL"));
    appendBehavList (mi, M("TP_WBALANCE_TEMPERATURE"), ADDSET_WB_TEMPERATURE, true);
    appendBehavList (mi, M("TP_WBALANCE_GREEN"), ADDSET_WB_GREEN, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_COLORBOOST_LABEL"));
    appendBehavList (mi, M("TP_COLORBOOST_AMOUNT"), ADDSET_CBOOST_AMOUNT, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_COLORSHIFT_LABEL"));
    appendBehavList (mi, M("TP_COLORSHIFT_BLUEYELLOW"), ADDSET_CS_BLUEYELLOW, false);
    appendBehavList (mi, M("TP_COLORSHIFT_GREENMAGENTA"), ADDSET_CS_GREENMAGENTA, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_ROTATE_LABEL"));
    appendBehavList (mi, M("TP_ROTATE_DEGREE"), ADDSET_ROTATE_DEGREE, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_DISTORTION_LABEL"));
    appendBehavList (mi, M("TP_DISTORTION_AMOUNT"), ADDSET_DIST_AMOUNT, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_PERSPECTIVE_LABEL"));
    appendBehavList (mi, M("TP_PERSPECTIVE_HORIZONTAL")+", "+M("TP_PERSPECTIVE_VERTICAL"), ADDSET_PERSPECTIVE, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_CACORRECTION_LABEL"));
    appendBehavList (mi, M("TP_CACORRECTION_BLUE")+", "+M("TP_CACORRECTION_RED"), ADDSET_CA, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_VIGNETTING_LABEL"));
    appendBehavList (mi, M("TP_VIGNETTING_AMOUNT"), ADDSET_VIGN_AMOUNT, false);

    behTreeView->expand_all ();

    return mvbpp;
}

void Preferences::appendBehavList (Gtk::TreeModel::iterator& parent, Glib::ustring label, int id, bool set) {

    Gtk::TreeModel::iterator ci = behModel->append (parent->children());
    ci->set_value (behavColumns.label, label);
    ci->set_value (behavColumns.visible, true);
    ci->set_value (behavColumns.badd, !set);
    ci->set_value (behavColumns.bset, set);
    ci->set_value (behavColumns.addsetid, id);
}

void Preferences::behAddRadioToggled (const Glib::ustring& path) {

    Gtk::TreeModel::iterator iter = behModel->get_iter (path);
    bool set = iter->get_value (behavColumns.bset);
    iter->set_value (behavColumns.bset, false);
    iter->set_value (behavColumns.badd, true);
}

void Preferences::behSetRadioToggled (const Glib::ustring& path) {

    Gtk::TreeModel::iterator iter = behModel->get_iter (path);
    bool add = iter->get_value (behavColumns.badd);
    iter->set_value (behavColumns.bset, true);
    iter->set_value (behavColumns.badd, false);
}

Gtk::Widget* Preferences::getProcParamsPanel () {

    Gtk::VBox* mvbpp = Gtk::manage (new Gtk::VBox ());

    Gtk::Frame* fpp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_IMPROCPARAMS")));
    Gtk::Label* drlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_FORRAW")+":"));
    rprofiles = Gtk::manage (new Gtk::ComboBoxText ());
    Gtk::Label* drimg = Gtk::manage (new Gtk::Label (M("PREFERENCES_FORIMAGE")+":"));
    iprofiles = Gtk::manage (new Gtk::ComboBoxText ());  
    Gtk::Table* defpt = Gtk::manage (new Gtk::Table (2, 2));
    defpt->attach (*drlab, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    defpt->attach (*rprofiles, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach (*drimg, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    defpt->attach (*iprofiles, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    fpp->add (*defpt);

    mvbpp->pack_start (*fpp, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_PROFILEHANDLING")));
    Gtk::VBox* vbdp = Gtk::manage (new Gtk::VBox ());
    saveParamsFile = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_PROFILESAVEINPUT")));
    vbdp->pack_start (*saveParamsFile, Gtk::PACK_SHRINK, 4);
    saveParamsCache = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_PROFILESAVECACHE")));
    vbdp->pack_start (*saveParamsCache, Gtk::PACK_SHRINK, 4); 
    Gtk::Label* lplab = Gtk::manage (new Gtk::Label (M("PREFERENCES_PROFILELOADPR")+":"));
    loadParamsPreference = Gtk::manage (new Gtk::ComboBoxText ());
    loadParamsPreference->append_text (M("PREFERENCES_PROFILEPRCACHE"));
    loadParamsPreference->append_text (M("PREFERENCES_PROFILEPRFILE"));
    Gtk::HBox* hb41 = Gtk::manage (new Gtk::HBox ());
    hb41->pack_start (*lplab, Gtk::PACK_SHRINK, 4);
    hb41->pack_start (*loadParamsPreference);
    vbdp->pack_start (*hb41, Gtk::PACK_SHRINK, 4);
    fdp->add (*vbdp);
    mvbpp->pack_start (*fdp, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdf = Gtk::manage (new Gtk::Frame (M("PREFERENCES_DARKFRAME")) );
    Gtk::HBox* hb42 = Gtk::manage (new Gtk::HBox ());
    darkFrameDir = Gtk::manage(new Gtk::FileChooserButton(M("PREFERENCES_DIRDARKFRAMES"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label *dfLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_DIRDARKFRAMES")));
    hb42->pack_start(*dfLab , Gtk::PACK_SHRINK, 4 );
    hb42->pack_start(*darkFrameDir);
    dfLabel = Gtk::manage(new Gtk::Label("Found:"));
    Gtk::VBox* vbdf = Gtk::manage (new Gtk::VBox ());
    vbdf->pack_start( *hb42, Gtk::PACK_SHRINK, 4);
    vbdf->pack_start( *dfLabel, Gtk::PACK_SHRINK, 4 );
    fdf->add( *vbdf );
    mvbpp->pack_start ( *fdf , Gtk::PACK_SHRINK, 4);
    mvbpp->set_border_width (4);

    std::vector<Glib::ustring> pnames;
    if (options.multiUser)
        parseDir (Options::rtdir + "/" + options.profilePath, pnames, paramFileExtension);
    parseDir (argv0 + "/" + options.profilePath, pnames, paramFileExtension);
    for (int i=0; i<pnames.size(); i++) {
        rprofiles->append_text (pnames[i]);
        iprofiles->append_text (pnames[i]);
    }

    dfconn = darkFrameDir->signal_file_set().connect ( sigc::mem_fun(*this, &Preferences::darkFrameChanged), true);

    return mvbpp;
}

Gtk::Widget* Preferences::getColorManagementPanel () {

    Gtk::VBox* mvbcm = Gtk::manage (new Gtk::VBox ());
    mvbcm->set_border_width (4);

    Gtk::Label* intlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_CMETRICINTENT")+":"));
    intent = Gtk::manage (new Gtk::ComboBoxText ());
    intent->append_text (M("PREFERENCES_INTENT_PERCEPTUAL"));
    intent->append_text (M("PREFERENCES_INTENT_RELATIVE"));
    intent->append_text (M("PREFERENCES_INTENT_SATURATION"));
    intent->append_text (M("PREFERENCES_INTENT_ABSOLUTE"));

    iccDir = Gtk::manage (new Gtk::FileChooserButton (M("PREFERENCES_ICCDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label* pdlabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_ICCDIR")+":"));

    monProfile = Gtk::manage (new Gtk::FileChooserButton (M("PREFERENCES_MONITORICC"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    Gtk::Label* mplabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_MONITORICC")+":"));

    Gtk::Table* colt = Gtk::manage (new Gtk::Table (3, 2));
    colt->attach (*intlab, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colt->attach (*intent, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    colt->attach (*pdlabel, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colt->attach (*iccDir, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);
    colt->attach (*mplabel, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colt->attach (*monProfile, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    mvbcm->pack_start (*colt, Gtk::PACK_SHRINK, 4);

    return mvbcm;
}

Gtk::Widget* Preferences::getGeneralPanel () {

    Gtk::VBox* mvbsd = new Gtk::VBox ();

    Gtk::Frame* fworklflow = new Gtk::Frame (M("PREFERENCES_WORKFLOW"));
    Gtk::HBox* hbworkflow = new Gtk::HBox ();
    hbworkflow->set_border_width (4);
    Gtk::Label* flayoutlab = new Gtk::Label (M("PREFERENCES_EDITORLAYOUT")+":");
    editorLayout = new Gtk::ComboBoxText ();

    editorLayout->append_text (M("PREFERENCES_SINGLETAB"));
    editorLayout->append_text (M("PREFERENCES_MULTITAB"));
    editorLayout->set_active (1);

    hbworkflow->pack_start (*flayoutlab, Gtk::PACK_SHRINK, 4);
    hbworkflow->pack_start (*editorLayout);
    fworklflow->add (*hbworkflow);
    mvbsd->pack_start (*fworklflow, Gtk::PACK_SHRINK, 4);
     
    Gtk::Frame* flang = new Gtk::Frame (M("PREFERENCES_DEFAULTLANG"));
    Gtk::HBox* hblang = new Gtk::HBox ();
    hblang->set_border_width (4);
    Gtk::Label* langlab = new Gtk::Label (M("PREFERENCES_SELECTLANG")+":");
    languages = new Gtk::ComboBoxText ();

    std::vector<Glib::ustring> langs;
    parseDir (argv0 + "/languages", langs, "");
    for (int i=0; i<langs.size(); i++) {
	if ("default" != langs[i] && "README" != langs[i] && "LICENSE" != langs[i]) {
  	    languages->append_text (langs[i]);
	}
    }

    Gtk::Label* langw = new Gtk::Label (Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")");
    hblang->pack_start (*langlab, Gtk::PACK_SHRINK, 4);
    hblang->pack_start (*languages);
    hblang->pack_end (*langw, Gtk::PACK_SHRINK, 4);
    flang->add (*hblang);
    mvbsd->pack_start (*flang, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* ftheme = new Gtk::Frame (M("PREFERENCES_DEFAULTTHEME"));
    Gtk::VBox* vbftheme = new Gtk::VBox ();
    Gtk::HBox* hbUseSystemTheme = new Gtk::HBox ();
    hbUseSystemTheme->set_border_width (4);
    Gtk::Label* useNextStart = new Gtk::Label (Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")");
    chUseSystemTheme =  new Gtk::CheckButton (M("PREFERENCES_USESYSTEMTHEME"));
    hbUseSystemTheme->pack_start(*chUseSystemTheme, Gtk::PACK_SHRINK);
    hbUseSystemTheme->pack_start (*useNextStart, Gtk::PACK_SHRINK, 4);
    vbftheme->pack_start(*hbUseSystemTheme, Gtk::PACK_SHRINK, 4);

    hbtheme = new Gtk::HBox ();
    hbtheme->set_border_width (4);
    Gtk::Label* themelab = new Gtk::Label (M("PREFERENCES_SELECTTHEME")+":");
    theme = new Gtk::ComboBoxText ();

    theme->append_text (Glib::ustring("(")+M("PREFERENCES_GTKTHEME")+")");
    theme->set_active (0);
    std::vector<Glib::ustring> themes;
    parseDir (argv0 + "/themes", themes, "");
    for (int i=0; i<themes.size(); i++) 
        theme->append_text (themes[i]);

    Gtk::Label* fontlab = new Gtk::Label (M("PREFERENCES_SELECTFONT")+":");
    fontbutton = new Gtk::FontButton ();
    fontbutton->set_use_size(true);
    fontbutton->set_font_name(options.font);

    hbtheme->pack_start (*themelab, Gtk::PACK_SHRINK, 4);
    hbtheme->pack_start (*theme);
    hbtheme->pack_start (*fontlab, Gtk::PACK_SHRINK, 4);
    hbtheme->pack_start (*fontbutton);
    vbftheme->pack_end(*hbtheme, Gtk::PACK_SHRINK, 4);
    ftheme->add (*vbftheme);
    mvbsd->pack_start (*ftheme, Gtk::PACK_SHRINK, 4);
  
//-----

    Gtk::Frame* frl = new Gtk::Frame (M("PREFERENCES_CLIPPINGIND"));
    blinkClipped = new Gtk::CheckButton (M("PREFERENCES_BLINKCLIPPED"));
    Gtk::VBox* vbrl = new Gtk::VBox ();
    vbrl->set_border_width (4);
    vbrl->pack_start (*blinkClipped, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* vbhl = new Gtk::HBox ();
    Gtk::Label* hll = new Gtk::Label (M("PREFERENCES_HLTHRESHOLD")+": ");
    hlThresh = new Gtk::SpinButton ();
    hlThresh->set_digits (0);
    hlThresh->set_increments (1, 10);
    hlThresh->set_range (0, 255);
    vbhl->pack_start (*hll, Gtk::PACK_SHRINK, 8);
    vbhl->pack_start (*hlThresh, Gtk::PACK_SHRINK, 8);

    vbrl->pack_start (*vbhl, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* vbsh = new Gtk::HBox ();
    Gtk::Label* shl = new Gtk::Label (M("PREFERENCES_SHTHRESHOLD")+": ");
    shThresh = new Gtk::SpinButton ();
    shThresh->show ();
    shThresh->set_digits (0);
    shThresh->set_increments (1, 10);
    shThresh->set_range (0, 255);
    vbsh->pack_start (*shl, Gtk::PACK_SHRINK, 8);
    vbsh->pack_start (*shThresh, Gtk::PACK_SHRINK, 8);
    vbrl->pack_start (*vbsh, Gtk::PACK_SHRINK, 4);

    frl->add (*vbrl);  
    mvbsd->pack_start (*frl, Gtk::PACK_SHRINK, 4);

//-----
    Gtk::Frame* fdf = new Gtk::Frame (M("PREFERENCES_DATEFORMAT"));

    Gtk::HBox* hb6 = new Gtk::HBox ();
    Gtk::VBox* dfvb = new Gtk::VBox ();
    Gtk::Label* dflab = new Gtk::Label (M("PREFERENCES_DATEFORMAT")+":");
    hb6->pack_start (*dflab, Gtk::PACK_SHRINK,4);
    dateformat = new Gtk::Entry ();
    dateformat->set_tooltip_markup (M("PREFERENCES_DATEFORMATHINT"));
    dflab->set_tooltip_markup (M("PREFERENCES_DATEFORMATHINT"));
    hb6->pack_start (*dateformat);
    dfvb->pack_start (*hb6, Gtk::PACK_SHRINK, 4);
    fdf->add (*dfvb);
    dfvb->set_border_width (4);

    mvbsd->pack_start (*fdf, Gtk::PACK_SHRINK, 4);

  //-----
    Gtk::Frame* fdg = new Gtk::Frame (M("PREFERENCES_EXTERNALEDITOR"));
    Gtk::VBox* dgvb = new Gtk::VBox ();

    Gtk::HBox* hb7c = new Gtk::HBox ();
    edOther = new Gtk::RadioButton (M("PREFERENCES_EDITORCMDLINE")+":");
    hb7c->pack_start (*edOther, Gtk::PACK_SHRINK,4);
    editorToSendTo = new Gtk::Entry ();
    hb7c->pack_start (*editorToSendTo);
    Gtk::RadioButton::Group ge = edOther->get_group();
  
#ifdef __APPLE__
  Gtk::HBox* hb7 = new Gtk::HBox ();
  edGimp = new Gtk::RadioButton ("GIMP");
  hb7->pack_start (*edGimp, Gtk::PACK_SHRINK,4);
  dgvb->pack_start (*hb7, Gtk::PACK_SHRINK, 4);
  edGimp->set_group (ge);
 
  Gtk::HBox* hb7b = new Gtk::HBox ();
  edPS = new Gtk::RadioButton (M("PREFERENCES_PSPATH")+":");
  hb7b->pack_start (*edPS, Gtk::PACK_SHRINK,4);
  psDir = new Gtk::FileChooserButton (M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER);
  hb7b->pack_start (*psDir);
  dgvb->pack_start (*hb7b, Gtk::PACK_SHRINK, 4);
  edPS->set_group (ge);
#elif defined _WIN32
  Gtk::HBox* hb7 = new Gtk::HBox ();
  edGimp = new Gtk::RadioButton (M("PREFERENCES_GIMPPATH")+":");
  hb7->pack_start (*edGimp, Gtk::PACK_SHRINK,4);
  gimpDir = new Gtk::FileChooserButton (M("PREFERENCES_GIMPPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER);
  hb7->pack_start (*gimpDir);
  dgvb->pack_start (*hb7, Gtk::PACK_SHRINK, 4);
  edGimp->set_group (ge);
 
  Gtk::HBox* hb7b = new Gtk::HBox ();
  edPS = new Gtk::RadioButton (M("PREFERENCES_PSPATH")+":");
  hb7b->pack_start (*edPS, Gtk::PACK_SHRINK,4);
  psDir = new Gtk::FileChooserButton (M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER);
  hb7b->pack_start (*psDir);
  dgvb->pack_start (*hb7b, Gtk::PACK_SHRINK, 4);
  edPS->set_group (ge);
#else
    Gtk::HBox* hb7 = new Gtk::HBox ();
    edGimp = new Gtk::RadioButton ("GIMP");
    hb7->pack_start (*edGimp, Gtk::PACK_SHRINK,4);
    dgvb->pack_start (*hb7, Gtk::PACK_SHRINK, 4);
    edGimp->set_group (ge);
#endif

    dgvb->pack_start (*hb7c, Gtk::PACK_SHRINK, 4);
    dgvb->set_border_width (4);
    fdg->add (*dgvb);
    mvbsd->pack_start (*fdg, Gtk::PACK_SHRINK, 4);

    mvbsd->set_border_width (4);

    tconn = theme->signal_changed().connect( sigc::mem_fun(*this, &Preferences::themeChanged) );
    fconn = fontbutton->signal_font_set().connect( sigc::mem_fun(*this, &Preferences::fontChanged) );
    usethcon = chUseSystemTheme->signal_clicked ().connect( sigc::mem_fun(*this, &Preferences::useThemeChanged) );
    
    return mvbsd;
}

Gtk::Widget* Preferences::getFileBrowserPanel () {

    Gtk::VBox* mvbfb = new Gtk::VBox ();
    mvbfb->set_border_width (4);

    Gtk::Frame* fsd = new Gtk::Frame (M("PREFERENCES_STARTUPIMDIR"));

    sdcurrent = new Gtk::RadioButton (M("PREFERENCES_DIRSOFTWARE"));
    sdlast    = new Gtk::RadioButton (M("PREFERENCES_DIRLAST"));
    sdhome    = new Gtk::RadioButton (M("PREFERENCES_DIRHOME"));
    sdother   = new Gtk::RadioButton (M("PREFERENCES_DIROTHER")+": ");
    startupdir = new Gtk::Entry ();

    Gtk::Button* sdselect = new Gtk::Button ("");
    sdselect->set_image (*(new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

    Gtk::RadioButton::Group opts = sdcurrent->get_group();
    sdlast->set_group (opts);
    sdhome->set_group (opts);
    sdother->set_group (opts);

    Gtk::VBox* vbsd = new Gtk::VBox ();
    vbsd->pack_start (*sdcurrent, Gtk::PACK_SHRINK,0);
    vbsd->pack_start (*sdlast, Gtk::PACK_SHRINK,0);
    vbsd->pack_start (*sdhome, Gtk::PACK_SHRINK,0);
    Gtk::HBox* otherbox = new Gtk::HBox ();
    otherbox->pack_start (*sdother, Gtk::PACK_SHRINK);
    otherbox->pack_start (*startupdir);
    otherbox->pack_end (*sdselect, Gtk::PACK_SHRINK, 4);
    vbsd->pack_start (*otherbox, Gtk::PACK_SHRINK, 0);
    vbsd->set_border_width (4);

    fsd->add (*vbsd);
    mvbfb->pack_start (*fsd, Gtk::PACK_SHRINK, 4);

    sdselect->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::selectStartupDir) );

//---


    Gtk::Frame* fro = new Gtk::Frame (M("PREFERENCES_FBROWSEROPTS"));
    showDateTime = new Gtk::CheckButton (M("PREFERENCES_SHOWDATETIME"));
    showBasicExif = new Gtk::CheckButton (M("PREFERENCES_SHOWBASICEXIF"));
    Gtk::VBox* vbro = new Gtk::VBox ();
    overlayedFileNames = new Gtk::CheckButton (M("PREFERENCES_OVERLAY_FILENAMES"));
    vbro->set_border_width (4);
    vbro->pack_start (*showDateTime, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*showBasicExif, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*overlayedFileNames, Gtk::PACK_SHRINK, 4); 

    fro->add (*vbro);  

    Gtk::Frame* fre = new Gtk::Frame (M("PREFERENCES_PARSEDEXT"));
    Gtk::VBox* vbre = new Gtk::VBox ();
    vbre->set_border_width (4);
    Gtk::HBox* hb0 = new Gtk::HBox ();
    Gtk::Label* elab = new Gtk::Label (M("PREFERENCES_PARSEDEXTADD")+":");
    hb0->pack_start (*elab, Gtk::PACK_SHRINK, 4);
    extension = new Gtk::Entry ();
    hb0->pack_start (*extension);
    addExt = new Gtk::Button ();
    delExt = new Gtk::Button ();
    addExt->set_tooltip_text (M("PREFERENCES_PARSEDEXTADDHINT"));
    delExt->set_tooltip_text (M("PREFERENCES_PARSEDEXTDELHINT"));
    Gtk::Image* addExtImg = new Gtk::Image (argv0+"/images/list-add12.png");
    Gtk::Image* delExtImg = new Gtk::Image (argv0+"/images/list-remove12r.png");
    addExt->add (*addExtImg);
    delExt->add (*delExtImg);
    hb0->pack_end (*delExt, Gtk::PACK_SHRINK, 4);
    hb0->pack_end (*addExt, Gtk::PACK_SHRINK, 4);
    extensions = new Gtk::TreeView ();
    Gtk::ScrolledWindow* hscrollw = new Gtk::ScrolledWindow ();
    hscrollw->set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    hscrollw->add (*extensions);
    extensionModel = Gtk::ListStore::create (extensionColumns);
    extensions->set_model (extensionModel);
    extensions->append_column_editable("Enabled", extensionColumns.enabled);
    extensions->append_column("Extension", extensionColumns.ext);
    extensions->set_headers_visible (false);  
    vbre->pack_start (*hscrollw);
    vbre->pack_start (*hb0, Gtk::PACK_SHRINK, 0);

    fre->add (*vbre);  

    Gtk::Frame* frc = new Gtk::Frame (M("PREFERENCES_CACHEOPTS"));
    Gtk::VBox* vbc = new Gtk::VBox ();
    frc->add (*vbc);  
    vbc->set_border_width (4);

    Gtk::Label* cflab = new Gtk::Label (M("PREFERENCES_CACHETHUMBFORM")+":");
    cformat = new Gtk::ComboBoxText ();
    cformat->append_text (M("PREFERENCES_CACHEFORMAT1"));
    cformat->append_text (M("PREFERENCES_CACHEFORMAT2"));
    cformat->append_text (M("PREFERENCES_CACHEFORMAT1")+", 16 bit");
    Gtk::HBox* hb2 = new Gtk::HBox ();
    hb2->pack_start (*cflab, Gtk::PACK_SHRINK, 4);
    hb2->pack_start (*cformat);
    vbc->pack_start (*hb2, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hb3 = new Gtk::HBox ();
    Gtk::Label* chlab = new Gtk::Label (M("PREFERENCES_CACHETHUMBHEIGHT")+":");
    maxThumbSize = new Gtk::SpinButton ();
    hb3->pack_start (*chlab, Gtk::PACK_SHRINK, 8);
    hb3->pack_start (*maxThumbSize, Gtk::PACK_SHRINK, 8);

    maxThumbSize->set_digits (0);
    maxThumbSize->set_increments (1, 10);
    maxThumbSize->set_range (40, 400);
    vbc->pack_start (*hb3, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hb4 = new Gtk::HBox ();
    Gtk::Label* celab = new Gtk::Label (M("PREFERENCES_CACHEMAXENTRIES")+":");
    maxCacheEntries = new Gtk::SpinButton ();
    hb4->pack_start (*celab, Gtk::PACK_SHRINK, 8);
    hb4->pack_start (*maxCacheEntries, Gtk::PACK_SHRINK, 8);

    maxCacheEntries->set_digits (0);
    maxCacheEntries->set_increments (1, 10);
    maxCacheEntries->set_range (10, 100000);
    vbc->pack_start (*hb4, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hb5 = new Gtk::HBox ();
    clearThumbnails = new Gtk::Button (M("PREFERENCES_CACHECLEARTHUMBS"));
    clearProfiles = new Gtk::Button (M("PREFERENCES_CACHECLEARPROFILES"));
    clearAll = new Gtk::Button (M("PREFERENCES_CACHECLEARALL"));
    hb5->pack_start (*clearThumbnails, Gtk::PACK_SHRINK, 8);
    hb5->pack_start (*clearProfiles, Gtk::PACK_SHRINK, 8);
    hb5->pack_start (*clearAll, Gtk::PACK_SHRINK, 8);
    vbc->pack_start (*hb5, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hb6 = new Gtk::HBox ();
    Gtk::VBox* vb6 = new Gtk::VBox ();

    vb6->pack_start (*fro);
    vb6->pack_end (*frc);
    hb6->pack_start (*vb6);
    hb6->pack_start (*fre);

    mvbfb->pack_start (*hb6, Gtk::PACK_SHRINK, 4);
  
//  mvbfb->pack_start (*fro, Gtk::PACK_SHRINK, 4);
//  mvbfb->pack_start (*fre);
//  mvbfb->pack_start (*frc, Gtk::PACK_SHRINK, 4);

    addExt->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::addExtPressed) );
    delExt->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::delExtPressed) );
    extension->signal_activate().connect( sigc::mem_fun(*this, &Preferences::addExtPressed) );
    clearThumbnails->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::clearThumbImagesPressed) );
    clearProfiles->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::clearProfilesPressed) );
    clearAll->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::clearAllPressed) );

    return mvbfb;
}

void Preferences::parseDir (Glib::ustring dirname, std::vector<Glib::ustring>& items, Glib::ustring ext) {

    // process directory
    Glib::Dir* dir = NULL;
    try {
        dir = new Glib::Dir (dirname);
    }
    catch (const Glib::FileError& fe) {
        return;
    }
    dirname = dirname + "/";
    for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
      Glib::ustring fname = dirname + *i;
      Glib::ustring sname = *i;
      // ignore directories
      if (!Glib::file_test (fname, Glib::FILE_TEST_IS_DIR) && sname.size() >= ext.size() && sname.substr (sname.size()-ext.size(), ext.size()).casefold() == ext) 
            items.push_back (sname.substr(0,sname.size()-ext.size()));
    }
    delete dir;
}

void Preferences::storePreferences () {

    moptions.defProfRaw          = rprofiles->get_active_text();
    moptions.defProfImg          = iprofiles->get_active_text();
    moptions.dateFormat          = dateformat->get_text();
    moptions.fbShowDateTime  = showDateTime->get_active ();
    moptions.fbShowBasicExif = showBasicExif->get_active ();
    moptions.blinkClipped    = blinkClipped->get_active ();
    moptions.highlightThreshold = (int)hlThresh->get_value ();
    moptions.shadowThreshold = (int)shThresh->get_value ();
    moptions.language        = languages->get_active_text ();
    moptions.theme           = theme->get_active_text ();
    moptions.useSystemTheme  = chUseSystemTheme->get_active ();
    moptions.font            = fontbutton->get_font_name();
#ifdef _WIN32    
    moptions.gimpDir        = gimpDir->get_filename ();
    moptions.psDir          = psDir->get_filename ();
#elif defined __APPLE__
    moptions.psDir          = psDir->get_filename (); 
#endif	
    moptions.customEditorProg = editorToSendTo->get_text ();
    if (edGimp->get_active ())
        moptions.editorToSendTo = 1;
#ifdef _WIN32    
    else if (edPS->get_active ())
        moptions.editorToSendTo = 2;
#elif defined __APPLE__   
    else if (edPS->get_active ())
      moptions.editorToSendTo = 2; 
#endif	
    else if (edOther->get_active ())
        moptions.editorToSendTo = 3;


    moptions.rtSettings.monitorProfile      = monProfile->get_filename ();
	moptions.rtSettings.iccDirectory        = iccDir->get_filename ();
	moptions.rtSettings.colorimetricIntent  = intent->get_active_row_number ();

    if (dmethodBatch->get_active_row_number()==0)
        moptions.rtSettings.demosaicMethodBatch = "eahd";
    else if (dmethodBatch->get_active_row_number()==1)
        moptions.rtSettings.demosaicMethodBatch = "hphd";
    else if (dmethodBatch->get_active_row_number()==2)
        moptions.rtSettings.demosaicMethodBatch = "vng4";
    else if (dmethodBatch->get_active_row_number()==3)
        moptions.rtSettings.demosaicMethodBatch = "amaze";
    else if (dmethodBatch->get_active_row_number()==4)
        moptions.rtSettings.demosaicMethodBatch = "dcb";
    else if (dmethodBatch->get_active_row_number()==5)
        moptions.rtSettings.demosaicMethodBatch = "ahd";
    else if (dmethodBatch->get_active_row_number()==6)
        moptions.rtSettings.demosaicMethodBatch = "bilinear";

    if (sdcurrent->get_active ()) 
        moptions.startupDir = STARTUPDIR_CURRENT;
    else if (sdhome->get_active ()) 
        moptions.startupDir = STARTUPDIR_HOME;
    else if (sdlast->get_active ()) 
        moptions.startupDir = STARTUPDIR_LAST;
    else if (sdother->get_active ()) {
        moptions.startupDir = STARTUPDIR_CUSTOM;
        moptions.startupPath = startupdir->get_text();
    }

    moptions.parseExtensions.clear ();
    moptions.parseExtensionsEnabled.clear ();
    Gtk::TreeNodeChildren c = extensionModel->children ();
    for (int i=0; i<c.size(); i++) {
        moptions.parseExtensions.push_back (c[i][extensionColumns.ext]);
        moptions.parseExtensionsEnabled.push_back (c[i][extensionColumns.enabled]);
    }
    
    if (cformat->get_active_row_number() == 0)
        moptions.thumbnailFormat = FT_Custom;
    else if (cformat->get_active_row_number() == 1)
        moptions.thumbnailFormat = FT_Jpeg;
    else if (cformat->get_active_row_number() == 2)
        moptions.thumbnailFormat = FT_Custom16;

    moptions.maxThumbnailHeight = (int)maxThumbSize->get_value ();
    moptions.maxCacheEntries = (int)maxCacheEntries->get_value ();
    moptions.overlayedFileNames = overlayedFileNames->get_active ();
    
    moptions.saveParamsFile = saveParamsFile->get_active ();
    moptions.saveParamsCache = saveParamsCache->get_active ();
    moptions.paramsLoadLocation = (PPLoadLocation)loadParamsPreference->get_active_row_number ();

    moptions.rtSettings.darkFramesPath =   darkFrameDir->get_filename();

    int i = 0;
    moptions.baBehav.resize (ADDSET_PARAM_NUM);
    for (Gtk::TreeIter sections=behModel->children().begin();  sections!=behModel->children().end(); sections++)
        for (Gtk::TreeIter adjs=sections->children().begin();  adjs!=sections->children().end(); adjs++) 
            moptions.baBehav[adjs->get_value (behavColumns.addsetid)] = adjs->get_value (behavColumns.badd);

    moptions.tabbedUI = (bool)editorLayout->get_active_row_number();
}

void Preferences::fillPreferences () {

    tconn.block (true);
    dfconn.block (true);

    rprofiles->set_active_text (moptions.defProfRaw);
    iprofiles->set_active_text (moptions.defProfImg);
    dateformat->set_text (moptions.dateFormat);
    if (Glib::file_test (moptions.rtSettings.monitorProfile, Glib::FILE_TEST_EXISTS)) 
        monProfile->set_filename (moptions.rtSettings.monitorProfile);
    if (Glib::file_test (moptions.rtSettings.iccDirectory, Glib::FILE_TEST_IS_DIR)) 
        iccDir->set_filename (moptions.rtSettings.iccDirectory);
	intent->set_active (moptions.rtSettings.colorimetricIntent);
    languages->set_active_text (moptions.language);
    theme->set_active_text (moptions.theme);
    chUseSystemTheme->set_active(moptions.useSystemTheme);
    fontbutton->set_font_name(moptions.font);
    showDateTime->set_active (moptions.fbShowDateTime);
    showBasicExif->set_active (moptions.fbShowBasicExif);
    blinkClipped->set_active (moptions.blinkClipped);
    hlThresh->set_value (moptions.highlightThreshold);
    shThresh->set_value (moptions.shadowThreshold);

    edGimp->set_active (moptions.editorToSendTo==1);
    edOther->set_active (moptions.editorToSendTo==3);
#ifdef _WIN32    
    edPS->set_active (moptions.editorToSendTo==2);
    if (Glib::file_test (moptions.gimpDir, Glib::FILE_TEST_IS_DIR)) 
        gimpDir->set_filename (moptions.gimpDir);
    if (Glib::file_test (moptions.psDir, Glib::FILE_TEST_IS_DIR)) 
        psDir->set_filename (moptions.psDir);
#elif defined __APPLE__
  edPS->set_active (moptions.editorToSendTo==2);
  if (Glib::file_test (moptions.psDir, Glib::FILE_TEST_IS_DIR))
    psDir->set_filename (moptions.psDir); 
#endif	
    editorToSendTo->set_text (moptions.customEditorProg);


    if (moptions.rtSettings.demosaicMethodBatch=="eahd")
        dmethodBatch->set_active (0);
    else if (moptions.rtSettings.demosaicMethodBatch=="hphd")
        dmethodBatch->set_active (1);
    else if (moptions.rtSettings.demosaicMethodBatch=="vng4")
        dmethodBatch->set_active (2);
	else if (moptions.rtSettings.demosaicMethodBatch=="amaze")//Emil's code for AMaZE
        dmethodBatch->set_active (3);
    else if (moptions.rtSettings.demosaicMethodBatch=="dcb")
        dmethodBatch->set_active (4);
    else if (moptions.rtSettings.demosaicMethodBatch=="ahd")
        dmethodBatch->set_active (5);
    else if (moptions.rtSettings.demosaicMethodBatch=="bilinear")
         dmethodBatch->set_active (6);

    if (moptions.startupDir==STARTUPDIR_CURRENT) 
        sdcurrent->set_active ();
    else if (moptions.startupDir==STARTUPDIR_LAST) 
        sdlast->set_active ();
    else if (moptions.startupDir==STARTUPDIR_HOME) 
        sdhome->set_active ();
    else if (moptions.startupDir==STARTUPDIR_CUSTOM) {
        sdother->set_active ();
        startupdir->set_text (moptions.startupPath);
    }
    
    extensionModel->clear ();
    for (int i=0; i<moptions.parseExtensions.size(); i++) {
        Gtk::TreeRow row = *(extensionModel->append());
        row[extensionColumns.enabled] = moptions.parseExtensionsEnabled[i];
        row[extensionColumns.ext]     = moptions.parseExtensions[i];
    }
       
    if (moptions.thumbnailFormat == FT_Custom)
        cformat->set_active (0);
    else if (moptions.thumbnailFormat == FT_Jpeg)
        cformat->set_active (1);
    else if (moptions.thumbnailFormat == FT_Custom16)
        cformat->set_active (2);
    
    maxThumbSize->set_value (moptions.maxThumbnailHeight);
    maxCacheEntries->set_value (moptions.maxCacheEntries);
    overlayedFileNames->set_active (moptions.overlayedFileNames);
    
    saveParamsFile->set_active (moptions.saveParamsFile);
    saveParamsCache->set_active (moptions.saveParamsCache);
    loadParamsPreference->set_active (moptions.paramsLoadLocation);    

    editorLayout->set_active(moptions.tabbedUI);
    darkFrameDir->set_filename( moptions.rtSettings.darkFramesPath );
    updateDFinfos();

    addc.block (true);
    setc.block (true);
    if (moptions.baBehav.size() == ADDSET_PARAM_NUM) {
		for (int i=0; i<moptions.baBehav.size(); i++) 
			for (Gtk::TreeIter sections=behModel->children().begin();  sections!=behModel->children().end(); sections++) 
				for (Gtk::TreeIter adjs=sections->children().begin();  adjs!=sections->children().end(); adjs++) 
					if (adjs->get_value (behavColumns.addsetid) == i) {
						adjs->set_value (behavColumns.badd, moptions.baBehav[i]==1);
						adjs->set_value (behavColumns.bset, moptions.baBehav[i]!=1);
						break;
					}
    }
    addc.block (false);
    setc.block (false);
    tconn.block (false);
    dfconn.block (false);
}

void Preferences::loadPressed () {

    moptions.copyFrom (&options);
    fillPreferences ();
}

void Preferences::savePressed () {

    storePreferences ();
    options.copyFrom (&moptions);
    Options::save ();
}

void Preferences::okPressed () {

    storePreferences ();
    workflowUpdate();
    options.copyFrom (&moptions);   
    hide ();
}

void Preferences::cancelPressed () {

	// set the initial theme back
	if (theme->get_active_text () != options.theme)
		switchThemeTo(options.theme);

	// set the initial font back
	if (fontbutton->get_font_name() != options.font)
		switchFontTo(options.font);
    hide ();
}

void Preferences::selectStartupDir () {

    Gtk::FileChooserDialog dialog(M("PREFERENCES_DIRSELECTDLG"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER);
//    dialog.set_transient_for(*this);

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) 
        startupdir->set_text (dialog.get_filename());
}

void Preferences::aboutPressed () {

    Splash* splash = new Splash (-1);
    splash->set_transient_for (*this);
    splash->set_modal (true);   
    splash->show ();
}

void Preferences::themeChanged () {

	switchThemeTo(theme->get_active_text ());
}

void Preferences::fontChanged () {

	switchFontTo(fontbutton->get_font_name());
}

void Preferences::switchThemeTo(Glib::ustring newTheme) {

	std::vector<Glib::ustring> files;
	files.push_back (argv0+"/themes/"+newTheme);
	Gtk::RC::set_default_files (files);
	Gtk::RC::reparse_all (Gtk::Settings::get_default());
	GdkEventClient event = { GDK_CLIENT_EVENT, NULL, TRUE, gdk_atom_intern("_GTK_READ_RCFILES", FALSE), 8 };
	gdk_event_send_clientmessage_toall ((GdkEvent*)&event);
}

void Preferences::workflowUpdate (){

    if(moptions.tabbedUI != options.tabbedUI)
    {
        parent->MoveFileBrowserToMain();
        parent->SetMainCurrent();
        if(moptions.tabbedUI)
            parent->epanel->hide_all();
        else
           parent->epanel->show_all();
    }
}

void Preferences::switchFontTo(Glib::ustring newFont) {

	Gtk::RC::parse_string (Glib::ustring::compose(
			"style \"clearlooks-default\" { font_name = \"%1\" }", newFont));
	Gtk::RC::reparse_all (Gtk::Settings::get_default());
	GdkEventClient event = { GDK_CLIENT_EVENT, NULL, TRUE, gdk_atom_intern("_GTK_READ_RCFILES", FALSE), 8 };
	gdk_event_send_clientmessage_toall ((GdkEvent*)&event);
}

void Preferences::useThemeChanged(){

    if(!chUseSystemTheme->get_active()){
        hbtheme->set_sensitive(true);
        fontbutton->set_sensitive(true);  
    }
    else{
        hbtheme->set_sensitive(false);
        fontbutton->set_sensitive(false);  
    }
}

void Preferences::addExtPressed () {

  Gtk::TreeNodeChildren c = extensionModel->children ();
  for (int i=0; i<c.size(); i++)
    if (c[i][extensionColumns.ext] == extension->get_text ())
        return;

  Gtk::TreeRow row = *(extensionModel->append());

  row[extensionColumns.enabled] = true;
  row[extensionColumns.ext]     = extension->get_text ();
}

void Preferences::delExtPressed () {

    extensionModel->erase (extensions->get_selection()->get_selected ()); 
}

void Preferences::clearProfilesPressed () {

    cacheMgr->clearProfiles ();
}

void Preferences::clearThumbImagesPressed () {

    cacheMgr->clearThumbImages ();
}

void Preferences::clearAllPressed () {

    cacheMgr->clearAll ();
}

void Preferences::darkFrameChanged ()
{
	Glib::ustring s(darkFrameDir->get_filename());
	if( s.compare( rtengine::dfm.getPathname()) !=0 ){
	   rtengine::dfm.init( s );
	   updateDFinfos();
	}
}
void Preferences::updateDFinfos()
{
    int t1,t2;
    rtengine::dfm.getStat(t1,t2);
    std::ostringstream s;
    s << "Found: "<< t1 << " shots, " << t2 << " templates";
    dfLabel->set_text(s.str());
}
