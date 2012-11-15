/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>, Oliver Duis <www.oliverduis.de>
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
#include "multilangmgr.h"
#include "splash.h"
#include "cachemanager.h"
#include "addsetids.h"
#include "../rtengine/dfmanager.h"
#include "../rtengine/ffmanager.h"
#include <sstream>
#include "../rtengine/safegtk.h"
#include "rtimage.h"

extern Options options;
extern Glib::ustring argv0;

Preferences::Preferences  (RTWindow *rtwindow):parent(rtwindow)  {
  
	splash = NULL;

    set_title (M("MAIN_BUTTON_PREFERENCES"));

    moptions.copyFrom (&options);

    /*
     * Do not increase height, since it's not visible on e.g. smaller netbook screens
     * Default height is about 620 pixels actually, that's why we do not set the height anymore
     * Netbook users will most certainly set a smaller font and use the "slimUI" mode,
     * so they'll be able to shrink the pref window and close it.
     */
    set_size_request (650, -1);
    set_default_size (options.preferencesWidth, options.preferencesHeight);
    set_border_width(4);

    Gtk::VBox* mainvb = get_vbox ();
    mainvb->set_spacing(8);
    set_has_separator (false);

    Gtk::Notebook* nb = Gtk::manage (new Gtk::Notebook ());
    mainvb->pack_start (*nb);

    Gtk::HBox* buttonpanel = Gtk::manage (new Gtk::HBox ());
    buttonpanel->set_spacing(8);
    mainvb->pack_start (*buttonpanel, Gtk::PACK_SHRINK, 0);

    Gtk::Button* about  = Gtk::manage (new Gtk::Button (M("GENERAL_ABOUT")));
    Gtk::Button* ok     = Gtk::manage (new Gtk::Button (M("GENERAL_OK")));
    Gtk::Button* cancel = Gtk::manage (new Gtk::Button (M("GENERAL_CANCEL")));

    about->set_image (*Gtk::manage(new RTImage ("rt-logo.png")));
    ok->set_image (*Gtk::manage(new RTImage ("gtk-apply.png")));
    cancel->set_image (*Gtk::manage(new RTImage ("gtk-cancel.png")));


    about->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::aboutPressed) );
    ok->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::okPressed) );
    cancel->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::cancelPressed) );

    buttonpanel->pack_start (*about, Gtk::PACK_SHRINK, 0);
    buttonpanel->pack_end (*ok, Gtk::PACK_SHRINK, 0);
    buttonpanel->pack_end (*cancel, Gtk::PACK_SHRINK, 0);

    nb->append_page (*getGeneralPanel(),        M("PREFERENCES_TAB_GENERAL"));
    nb->append_page (*getProcParamsPanel(),     M("PREFERENCES_TAB_IMPROC"));
    nb->append_page (*getFileBrowserPanel(),    M("PREFERENCES_TAB_BROWSER"));
    nb->append_page (*getColorManagementPanel(),M("PREFERENCES_TAB_COLORMGR"));
    nb->append_page (*getBatchProcPanel(),      M("PREFERENCES_BATCH_PROCESSING"));
    nb->append_page (*getSoundPanel(),          M("PREFERENCES_TAB_SOUND"));
    nb->set_current_page (0);

    fillPreferences ();

    show_all_children ();
    set_modal (true);
}


Preferences::~Preferences () {

    options.preferencesWidth = get_width();
    options.preferencesHeight = get_height();
}

Gtk::Widget* Preferences::getBatchProcPanel () {

    Gtk::VBox* mvbpp = Gtk::manage (new Gtk::VBox ());
    mvbpp->set_border_width(4);

    Gtk::ScrolledWindow* behscrollw = Gtk::manage (new Gtk::ScrolledWindow ());
    behscrollw->set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    behscrollw->set_border_width(8);
    behscrollw->set_size_request(-1, 60);
    Gtk::Frame* behFrame = Gtk::manage (new Gtk::Frame (M("PREFERENCES_BEHAVIOR")));
    behFrame->add (*behscrollw);
    //mvbpp->pack_start (*behFrame);
    mvbpp->pack_start (*behFrame, Gtk::PACK_EXPAND_WIDGET, 4);
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

    /*
     *   The TRUE/FALSE values of appendBehavList are replaced by the one defined in options.cc,
     */
    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_EXPOSURE_LABEL"));
    appendBehavList (mi, M("TP_EXPOSURE_EXPCOMP"), ADDSET_TC_EXPCOMP, false);
    appendBehavList (mi, M("TP_EXPOSURE_COMPRHIGHLIGHTS"), ADDSET_TC_HLCOMPAMOUNT, false);
    appendBehavList (mi, M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), ADDSET_TC_HLCOMPTHRESH, false);
    appendBehavList (mi, M("TP_EXPOSURE_BLACKLEVEL"), ADDSET_TC_BLACKLEVEL, false);
    appendBehavList (mi, M("TP_EXPOSURE_COMPRSHADOWS"), ADDSET_TC_SHCOMP, false);
    appendBehavList (mi, M("TP_EXPOSURE_BRIGHTNESS"), ADDSET_TC_BRIGHTNESS, false);
    appendBehavList (mi, M("TP_EXPOSURE_CONTRAST"), ADDSET_TC_CONTRAST, false);
    appendBehavList (mi, M("TP_EXPOSURE_SATURATION"), ADDSET_TC_SATURATION, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_SHADOWSHLIGHTS_LABEL"));
    appendBehavList (mi, M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), ADDSET_SH_HIGHLIGHTS, false);
    appendBehavList (mi, M("TP_SHADOWSHLIGHTS_SHADOWS"), ADDSET_SH_SHADOWS, false);
    appendBehavList (mi, M("TP_SHADOWSHLIGHTS_LOCALCONTR"), ADDSET_SH_LOCALCONTRAST, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_LABCURVE_LABEL"));
    appendBehavList (mi, M("TP_LABCURVE_BRIGHTNESS"), ADDSET_LC_BRIGHTNESS, false);
    appendBehavList (mi, M("TP_LABCURVE_CONTRAST"), ADDSET_LC_CONTRAST, false);
	appendBehavList (mi, M("TP_LABCURVE_CHROMATICITY"), ADDSET_LC_CHROMATICITY, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_SHARPENING_LABEL"));
    appendBehavList (mi, M("TP_SHARPENING_AMOUNT"), ADDSET_SHARP_AMOUNT, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_SHARPENEDGE_LABEL"));
    appendBehavList (mi, M("TP_SHARPENEDGE_PASSES"), ADDSET_SHARPENEDGE_PASS, false);
    appendBehavList (mi, M("TP_SHARPENEDGE_AMOUNT"), ADDSET_SHARPENEDGE_AMOUNT, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_SHARPENMICRO_LABEL"));
    appendBehavList (mi, M("TP_SHARPENMICRO_AMOUNT"), ADDSET_SHARPENMICRO_AMOUNT, false);
    appendBehavList (mi, M("TP_SHARPENMICRO_UNIFORMITY"), ADDSET_SHARPENMICRO_UNIFORMITY, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_DIRPYRDENOISE_LABEL"));
    appendBehavList (mi, M("TP_DIRPYRDENOISE_LUMA")+", "+M("TP_DIRPYRDENOISE_CHROMA"), ADDSET_DIRPYRDN_CHLUM, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_GAMMA"), ADDSET_DIRPYRDN_GAMMA, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_WBALANCE_LABEL"));
    appendBehavList (mi, M("TP_WBALANCE_TEMPERATURE"), ADDSET_WB_TEMPERATURE, true);
    appendBehavList (mi, M("TP_WBALANCE_GREEN"), ADDSET_WB_GREEN, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_COLORAPP_LABEL"));
    appendBehavList (mi, M("TP_COLORAPP_CIECAT_DEGREE"),ADDSET_CAT_DEGREE, true);
    appendBehavList (mi, M("TP_COLORAPP_ADAPTSCENE"),ADDSET_CAT_ADAPTSCENE, true);
    appendBehavList (mi, M("TP_COLORAPP_ADAPTVIEWING"),ADDSET_CAT_ADAPTVIEWING, true);
    appendBehavList (mi, M("TP_COLORAPP_LIGHT"),ADDSET_CAT_LIGHT, true);
    appendBehavList (mi, M("TP_COLORAPP_BRIGHT"),ADDSET_CAT_BRIGHT, true);
    appendBehavList (mi, M("TP_COLORAPP_CHROMA"),ADDSET_CAT_CHROMA, true);
    appendBehavList (mi, M("TP_COLORAPP_RSTPRO"),ADDSET_CAT_RSTPRO, true);
    appendBehavList (mi, M("TP_COLORAPP_CONTRAST"),ADDSET_CAT_CONTRAST, true);
    appendBehavList (mi, M("TP_COLORAPP_CONTRAST_Q"),ADDSET_CAT_CONTRAST_Q, true);
    appendBehavList (mi, M("TP_COLORAPP_CHROMA_S"),ADDSET_CAT_CHROMA_S, true);
    appendBehavList (mi, M("TP_COLORAPP_CHROMA_M"),ADDSET_CAT_CHROMA_M, true);
    appendBehavList (mi, M("TP_COLORAPP_HUE"),ADDSET_CAT_HUE, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_VIBRANCE_LABEL"));
    appendBehavList (mi, M("TP_VIBRANCE_PASTELS"), ADDSET_VIBRANCE_PASTELS, false);
    appendBehavList (mi, M("TP_VIBRANCE_SATURATED"), ADDSET_VIBRANCE_SATURATED, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_GAMMA_OUTPUT"));
    appendBehavList (mi, M("TP_GAMMA_CURV"), ADDSET_FREE_OUPUT_GAMMA, false);
    appendBehavList (mi, M("TP_GAMMA_SLOP"), ADDSET_FREE_OUTPUT_SLOPE, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_CHMIXER_LABEL"));
    appendBehavList (mi, M("TP_CHMIXER_RED")+", "+M("TP_CHMIXER_GREEN")+", "+M("TP_CHMIXER_BLUE"), ADDSET_CHMIXER, false);

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
    appendBehavList (mi, M("TP_VIGNETTING_RADIUS"), ADDSET_VIGN_RADIUS, false);
    appendBehavList (mi, M("TP_VIGNETTING_STRENGTH"), ADDSET_VIGN_STRENGTH, false);
    appendBehavList (mi, M("TP_VIGNETTING_CENTER_X")+", "+M("TP_VIGNETTING_CENTER_Y"), ADDSET_VIGN_CENTER, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_DIRPYREQUALIZER_LABEL"));
    appendBehavList (mi, M("TP_EXPOSURE_CONTRAST")+", "+M("TP_DIRPYREQUALIZER_THRESHOLD"), ADDSET_DIRPYREQ, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_PREPROCESS_LABEL"));
    appendBehavList (mi, M("TP_PREPROCESS_GREENEQUIL"), ADDSET_PREPROCESS_GREENEQUIL, false);
    appendBehavList (mi, M("TP_PREPROCESS_LINEDENOISE"), ADDSET_PREPROCESS_LINEDENOISE, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_EXPOSCORR_LABEL"));
    appendBehavList (mi, M("TP_RAWEXPOS_LINEAR"), ADDSET_RAWEXPOS_LINEAR, false);
    appendBehavList (mi, M("TP_RAWEXPOS_PRESER"), ADDSET_RAWEXPOS_PRESER, false);
    appendBehavList (mi, M("TP_RAWEXPOS_BLACKS"), ADDSET_RAWEXPOS_BLACKS, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_CHROMATABERR_LABEL"));
    appendBehavList (mi, M("TP_RAWCACORR_CARED")+", "+M("TP_RAWCACORR_CABLUE"), ADDSET_RAWCACORR, true);

    behTreeView->expand_all ();

    chOverwriteOutputFile =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_OVERWRITEOUTPUTFILE")) );
    mvbpp->pack_start(*chOverwriteOutputFile, Gtk::PACK_SHRINK, 4);

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
    //bool set = iter->get_value (behavColumns.bset);
    iter->set_value (behavColumns.bset, false);
    iter->set_value (behavColumns.badd, true);
}

void Preferences::behSetRadioToggled (const Glib::ustring& path) {

    Gtk::TreeModel::iterator iter = behModel->get_iter (path);
    //bool add = iter->get_value (behavColumns.badd);
    iter->set_value (behavColumns.bset, true);
    iter->set_value (behavColumns.badd, false);
}

Gtk::Widget* Preferences::getProcParamsPanel () {

    Gtk::VBox* mvbpp = Gtk::manage (new Gtk::VBox ());

    Gtk::Frame* fpp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_IMPROCPARAMS")));
    Gtk::Label* drlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_FORRAW")+":", Gtk::ALIGN_LEFT));
    rprofiles = Gtk::manage (new Gtk::ComboBoxText ());
    rprofiles->set_size_request(50, -1);
    rprofiles->signal_changed().connect( sigc::mem_fun(*this, &Preferences::forRAWComboChanged) );
    forRAWComboChanged(); // update the tooltip
    Gtk::Label* drimg = Gtk::manage (new Gtk::Label (M("PREFERENCES_FORIMAGE")+":", Gtk::ALIGN_LEFT));
    iprofiles = Gtk::manage (new Gtk::ComboBoxText ());
    iprofiles->set_size_request(50, -1);
    iprofiles->signal_changed().connect( sigc::mem_fun(*this, &Preferences::forImageComboChanged) );
    forImageComboChanged(); // update the tooltip
    Gtk::Table* defpt = Gtk::manage (new Gtk::Table (2, 2));
    defpt->attach (*drlab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach (*rprofiles, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    defpt->attach (*drimg, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach (*iprofiles, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    fpp->add (*defpt);

    mvbpp->pack_start (*fpp, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_PROFILEHANDLING")));
    Gtk::VBox* vbdp = Gtk::manage (new Gtk::VBox ());
    vbdp->set_border_width (4);
    saveParamsFile = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_PROFILESAVEINPUT")));
    vbdp->pack_start (*saveParamsFile, Gtk::PACK_SHRINK, 4);
    saveParamsCache = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_PROFILESAVECACHE")));
    vbdp->pack_start (*saveParamsCache, Gtk::PACK_SHRINK, 4);
    Gtk::Label* lplab = Gtk::manage (new Gtk::Label (M("PREFERENCES_PROFILELOADPR")+":"));
    loadParamsPreference = Gtk::manage (new Gtk::ComboBoxText ());
    loadParamsPreference->append_text (M("PREFERENCES_PROFILEPRCACHE"));
    loadParamsPreference->append_text (M("PREFERENCES_PROFILEPRFILE"));
    Gtk::HBox* hb41 = Gtk::manage (new Gtk::HBox ());
    hb41->pack_start (*lplab, Gtk::PACK_SHRINK, 0);
    hb41->pack_start (*loadParamsPreference, Gtk::PACK_EXPAND_WIDGET, 0);
    hb41->set_spacing(4);
    vbdp->pack_start (*hb41, Gtk::PACK_EXPAND_WIDGET, 4);
    fdp->add (*vbdp);
    mvbpp->pack_start (*fdp, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdf = Gtk::manage (new Gtk::Frame (M("PREFERENCES_DARKFRAME")) );
    Gtk::HBox* hb42 = Gtk::manage (new Gtk::HBox ());
    darkFrameDir = Gtk::manage(new Gtk::FileChooserButton(M("PREFERENCES_DIRDARKFRAMES"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label *dfLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_DIRDARKFRAMES")+":"));
    hb42->pack_start(*dfLab , Gtk::PACK_SHRINK, 4 );
    hb42->pack_start(*darkFrameDir, Gtk::PACK_EXPAND_WIDGET, 4);
    dfLabel = Gtk::manage(new Gtk::Label("Found:"));
    Gtk::VBox* vbdf = Gtk::manage (new Gtk::VBox ());
    vbdf->pack_start( *hb42, Gtk::PACK_SHRINK, 4);
    vbdf->pack_start( *dfLabel, Gtk::PACK_SHRINK, 4 );
    fdf->add( *vbdf );
    mvbpp->pack_start ( *fdf , Gtk::PACK_SHRINK, 4);
    mvbpp->set_border_width (4);

    //dfconn = darkFrameDir->signal_file_set().connect ( sigc::mem_fun(*this, &Preferences::darkFrameChanged), true);
    dfconn = darkFrameDir->signal_current_folder_changed().connect ( sigc::mem_fun(*this, &Preferences::darkFrameChanged), true);

    // FLATFIELD
    Gtk::Frame* fff = Gtk::manage (new Gtk::Frame (M("PREFERENCES_FLATFIELD")) );
    Gtk::HBox* hb43 = Gtk::manage (new Gtk::HBox ());
    flatFieldDir = Gtk::manage(new Gtk::FileChooserButton(M("PREFERENCES_FLATFIELDSDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label *ffLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_FLATFIELDSDIR")));
    hb43->pack_start(*ffLab , Gtk::PACK_SHRINK, 4 );
    hb43->pack_start(*flatFieldDir);
    ffLabel = Gtk::manage(new Gtk::Label("Found:"));
    Gtk::VBox* vbff = Gtk::manage (new Gtk::VBox ());
    vbff->pack_start( *hb43, Gtk::PACK_SHRINK, 4);
    vbff->pack_start( *ffLabel, Gtk::PACK_SHRINK, 4 );
    fff->add( *vbff );
    mvbpp->pack_start ( *fff , Gtk::PACK_SHRINK, 4);
    mvbpp->set_border_width (4);

    //ffconn = flatFieldDir->signal_file_set().connect ( sigc::mem_fun(*this, &Preferences::flatFieldChanged), true);
    ffconn = flatFieldDir->signal_current_folder_changed().connect ( sigc::mem_fun(*this, &Preferences::flatFieldChanged), true);
	
    std::vector<Glib::ustring> pnames;
    parseDir (options.getUserProfilePath(), pnames, paramFileExtension);
    parseDir (options.getGlobalProfilePath(), pnames, paramFileExtension);
    for (size_t i=0; i<pnames.size(); i++) {
        rprofiles->append_text (pnames[i]);
        iprofiles->append_text (pnames[i]);
    }

    Gtk::Frame* fmd = Gtk::manage (new Gtk::Frame (M("PREFERENCES_METADATA")));
    Gtk::VBox* vbmd = Gtk::manage (new Gtk::VBox ());
    ckbTunnelMetaData = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_TUNNELMETADATA")));
    vbmd->pack_start (*ckbTunnelMetaData, Gtk::PACK_SHRINK, 4);
    fmd->add (*vbmd);
    mvbpp->pack_start (*fmd, Gtk::PACK_SHRINK, 4);

    return mvbpp;
}

Gtk::Widget* Preferences::getColorManagementPanel () {

    Gtk::VBox* mvbcm = Gtk::manage (new Gtk::VBox ());
    mvbcm->set_border_width (4);

    Gtk::Label* intlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_CMETRICINTENT")+":", Gtk::ALIGN_LEFT));
    intent = Gtk::manage (new Gtk::ComboBoxText ());
    intent->append_text (M("PREFERENCES_INTENT_PERCEPTUAL"));
    intent->append_text (M("PREFERENCES_INTENT_RELATIVE"));
    intent->append_text (M("PREFERENCES_INTENT_SATURATION"));
    intent->append_text (M("PREFERENCES_INTENT_ABSOLUTE"));

    iccDir = Gtk::manage (new Gtk::FileChooserButton (M("PREFERENCES_ICCDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label* pdlabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_ICCDIR")+":", Gtk::ALIGN_LEFT));

    monProfile = Gtk::manage (new Gtk::FileChooserButton (M("PREFERENCES_MONITORICC"), Gtk::FILE_CHOOSER_ACTION_OPEN));
    Gtk::Label* mplabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_MONITORICC")+":", Gtk::ALIGN_LEFT));

	cbAutoMonProfile = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_AUTOMONPROFILE")));
	autoMonProfileConn  = cbAutoMonProfile->signal_toggled().connect (sigc::mem_fun(*this, &Preferences::autoMonProfileToggled));

    Gtk::Table* colt = Gtk::manage (new Gtk::Table (3, 2));
    colt->attach (*intlab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colt->attach (*intent, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colt->attach (*pdlabel, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colt->attach (*iccDir, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colt->attach (*mplabel, 0, 1, 2, 3, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colt->attach (*monProfile, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colt->attach (*cbAutoMonProfile, 1, 2, 3, 4, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    mvbcm->pack_start (*colt, Gtk::PACK_SHRINK, 4);

    autoMonProfileToggled();
    //Gtk::Frame* fdp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_OUTPUTDEVICE")));
    Gtk::VBox* vbdp = Gtk::manage (new Gtk::VBox ());
    vbdp->set_border_width (4);
    //Gtk::Label* viewlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_VIEW")+":"));
    Gtk::Label* viewlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_VIEW")+":", Gtk::ALIGN_LEFT));

    view = Gtk::manage (new Gtk::ComboBoxText ());
    view->append_text (M("PREFERENCES_D50"));
    view->append_text (M("PREFERENCES_D55"));
    view->append_text (M("PREFERENCES_D60"));
    view->append_text (M("PREFERENCES_D65"));
    view->append_text (M("PREFERENCES_BLACKBODY"));
    view->append_text (M("PREFERENCES_FLUOF2"));
    view->append_text (M("PREFERENCES_FLUOF7"));
    view->append_text (M("PREFERENCES_FLUOF11"));

    Gtk::Label* greylab = Gtk::manage (new Gtk::Label (M("PREFERENCES_GREY")+":", Gtk::ALIGN_LEFT));
    grey = Gtk::manage (new Gtk::ComboBoxText ());
    grey->append_text (M("PREFERENCES_GREY05"));
    grey->append_text (M("PREFERENCES_GREY10"));
    grey->append_text (M("PREFERENCES_GREY15"));
    grey->append_text (M("PREFERENCES_GREY18"));
    grey->append_text (M("PREFERENCES_GREY23"));
    grey->append_text (M("PREFERENCES_GREY30"));
    grey->append_text (M("PREFERENCES_GREY40"));

    Gtk::Label* restartNeeded1 = Gtk::manage( new Gtk::Label (Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    Gtk::Label* restartNeeded2 = Gtk::manage( new Gtk::Label (Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );

    Gtk::Table* colo = Gtk::manage (new Gtk::Table (2, 3));
    colo->attach (*viewlab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colo->attach (*view, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colo->attach (*restartNeeded1, 2, 3, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colo->attach (*greylab, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colo->attach (*grey, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colo->attach (*restartNeeded2, 2, 3, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);

    mvbcm->pack_start (*colo, Gtk::PACK_SHRINK, 4);

    return mvbcm;
}

Gtk::Widget* Preferences::getGeneralPanel () {

    Gtk::VBox* mvbsd = Gtk::manage( new Gtk::VBox () );

    Gtk::Frame* fworklflow = Gtk::manage(  new Gtk::Frame (M("PREFERENCES_WORKFLOW")) );
    Gtk::HBox* hbworkflow = Gtk::manage( new Gtk::HBox () );
    hbworkflow->set_border_width (4);
    Gtk::Label* flayoutlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_EDITORLAYOUT")+":") );
    editorLayout = Gtk::manage( new Gtk::ComboBoxText () );
    editorLayout->set_size_request(45, -1);

    editorLayout->append_text (M("PREFERENCES_SINGLETAB"));
    editorLayout->append_text (M("PREFERENCES_SINGLETABVERTAB"));
    editorLayout->append_text (M("PREFERENCES_MULTITAB"));
    editorLayout->append_text (M("PREFERENCES_MULTITABDUALMON"));
    editorLayout->set_active (2);
    editorLayout->signal_changed().connect( sigc::mem_fun(*this, &Preferences::layoutComboChanged) );
    layoutComboChanged(); // update the tooltip

    hbworkflow->pack_start (*flayoutlab, Gtk::PACK_SHRINK, 4);
    hbworkflow->pack_start (*editorLayout);
    Gtk::Label* lNextStart = Gtk::manage( new Gtk::Label (Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    hbworkflow->pack_end (*lNextStart, Gtk::PACK_SHRINK, 4);

    Gtk::VBox* vbworkflow = Gtk::manage( new Gtk::VBox () );
    vbworkflow->pack_start (*hbworkflow, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hbworkflow2 = Gtk::manage( new Gtk::HBox () );
    ckbHistogramPositionLeft =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_HISTOGRAMPOSITIONLEFT")) );
    hbworkflow2->pack_start (*ckbHistogramPositionLeft, Gtk::PACK_SHRINK, 4);
    ckbShowProfileSelector =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SHOWPROFILESELECTOR")) );
    hbworkflow2->pack_start (*ckbShowProfileSelector, Gtk::PACK_SHRINK, 4);
    ckbSquareDetailWindow =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SQUAREDETAILWINDOW")) );
    hbworkflow2->pack_start (*ckbSquareDetailWindow, Gtk::PACK_SHRINK, 4);
    vbworkflow->pack_start (*hbworkflow2, Gtk::PACK_SHRINK, 4);
    
    Gtk::HBox* hbworkflow3 = Gtk::manage( new Gtk::HBox () );
    ckbFileBrowserToolbarSingleRow =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_FILEBROWSERTOOLBARSINGLEROW")) );

    hbworkflow3->pack_start (*ckbFileBrowserToolbarSingleRow, Gtk::PACK_SHRINK, 4);
    vbworkflow->pack_start (*hbworkflow3, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hbworkflow4 = Gtk::manage( new Gtk::HBox () );

    Gtk::Label* hb4label =  Gtk::manage( new Gtk::Label (M("PREFERENCES_TP_LABEL")) );
    hbworkflow4->pack_start (*hb4label, Gtk::PACK_SHRINK, 4);
    ckbHideTPVScrollbar =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_TP_VSCROLLBAR")) );
    hbworkflow4->pack_start (*ckbHideTPVScrollbar, Gtk::PACK_SHRINK, 4);

    ckbUseIconNoText =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_TP_USEICONORTEXT")) );
    hbworkflow4->pack_start (*ckbUseIconNoText, Gtk::PACK_SHRINK, 4);

    vbworkflow->pack_start (*hbworkflow4, Gtk::PACK_SHRINK, 4);

    fworklflow->add (*vbworkflow);
    mvbsd->pack_start (*fworklflow, Gtk::PACK_SHRINK, 4);
     
    Gtk::Frame* flang = Gtk::manage( new Gtk::Frame (M("PREFERENCES_DEFAULTLANG")) );
    Gtk::HBox* hblang = Gtk::manage( new Gtk::HBox () );
    hblang->set_border_width (4);

    ckbLangAutoDetect =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_LANGAUTODETECT")) );

    Gtk::Label* langlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTLANG")+":") );
    languages = Gtk::manage( new Gtk::ComboBoxText () );

    std::vector<Glib::ustring> langs;
    parseDir (argv0 + "/languages", langs, "");
    for (size_t i=0; i<langs.size(); i++) {
        if ("default" != langs[i] && "README" != langs[i] && "LICENSE" != langs[i]) {
            languages->append_text (langs[i]);
        }
    }

    Gtk::Label* langw = Gtk::manage( new Gtk::Label (Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    hblang->pack_start (*ckbLangAutoDetect, Gtk::PACK_SHRINK, 4);
    hblang->pack_start (*langlab, Gtk::PACK_SHRINK, 8);
    hblang->pack_start (*languages);
    hblang->pack_end (*langw, Gtk::PACK_SHRINK, 4);
    flang->add (*hblang);
    mvbsd->pack_start (*flang, Gtk::PACK_SHRINK, 4);

    langAutoDetectConn  = ckbLangAutoDetect->signal_toggled().connect (sigc::mem_fun(*this, &Preferences::langAutoDetectToggled));

    Gtk::Frame* ftheme = Gtk::manage( new Gtk::Frame (M("PREFERENCES_DEFAULTTHEME")) );
    Gtk::VBox* vbftheme = Gtk::manage( new Gtk::VBox () );
    vbftheme->set_border_width(4);
    vbftheme->set_spacing(4);
    Gtk::HBox* hbUseSystemTheme = Gtk::manage( new Gtk::HBox () );
    hbUseSystemTheme->set_spacing(4);
    chUseSystemTheme =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_USESYSTEMTHEME")) );
    Gtk::Label* useNextStart = Gtk::manage( new Gtk::Label (Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );

    Gtk::Label* cutOverlayLabel = Gtk::manage( new Gtk::Label (M("PREFERENCES_CUTOVERLAYBRUSH") + ":") );
    butCropCol= Gtk::manage( new Gtk::ColorButton() );
    butCropCol->set_use_alpha(true);

    hbUseSystemTheme->pack_start(*chUseSystemTheme, Gtk::PACK_SHRINK);
    hbUseSystemTheme->pack_start (*useNextStart, Gtk::PACK_SHRINK, 0);
    hbUseSystemTheme->pack_end (*butCropCol, Gtk::PACK_SHRINK, 0);
    hbUseSystemTheme->pack_end (*cutOverlayLabel, Gtk::PACK_SHRINK, 0);
    vbftheme->pack_start(*hbUseSystemTheme, Gtk::PACK_SHRINK, 0);

    hbtheme = Gtk::manage( new Gtk::HBox () );
    hbtheme->set_spacing (4);
    Gtk::Label* themelab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTTHEME")+":") );
    theme = Gtk::manage( new Gtk::ComboBoxText () );

    theme->append_text (Glib::ustring("(")+M("PREFERENCES_GTKTHEME")+")");
    theme->set_active (0);
    std::vector<Glib::ustring> themes;
    parseDir (argv0 + "/themes", themes, ".gtkrc");
    for (size_t i=0; i<themes.size(); i++)
        theme->append_text (themes[i]);

    Gtk::Label* fontlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTFONT")+":") );
    fontbutton = Gtk::manage( new Gtk::FontButton ());
    fontbutton->set_use_size(true);
    fontbutton->set_font_name(options.font);

    hbtheme->pack_start (*themelab, Gtk::PACK_SHRINK, 0);
    hbtheme->pack_start (*theme);
    hbtheme->pack_start (*fontlab, Gtk::PACK_SHRINK, 0);
    hbtheme->pack_start (*fontbutton);
    vbftheme->pack_start(*hbtheme, Gtk::PACK_SHRINK, 0);

    slimUI = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SLIMUI")) );
    vbftheme->pack_start(*slimUI, Gtk::PACK_SHRINK, 0);

    ftheme->add (*vbftheme);
    mvbsd->pack_start (*ftheme, Gtk::PACK_SHRINK, 0);

//-----

    Gtk::HBox* hbcd = Gtk::manage( new Gtk::HBox () );
    hbcd->set_spacing(4);

    Gtk::Frame* frl = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CLIPPINGIND")));
    blinkClipped = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_BLINKCLIPPED")) );
    Gtk::VBox* vbrl = Gtk::manage( new Gtk::VBox () );
    vbrl->set_border_width(4);
    vbrl->set_spacing (4);
    vbrl->pack_start (*blinkClipped, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* vbhl = Gtk::manage( new Gtk::HBox () );
    vbhl->set_spacing(4);
    Gtk::Label* hll = Gtk::manage( new Gtk::Label (M("PREFERENCES_HLTHRESHOLD")+": "));
    hlThresh = Gtk::manage( new Gtk::SpinButton () );
    hlThresh->set_digits (0);
    hlThresh->set_increments (1, 10);
    hlThresh->set_range (0, 255);
    vbhl->pack_start (*hll, Gtk::PACK_SHRINK, 0);
    vbhl->pack_end (*hlThresh, Gtk::PACK_SHRINK, 0);

    vbrl->pack_start (*vbhl, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* vbsh = Gtk::manage( new Gtk::HBox () );
    vbsh->set_spacing (4);
    Gtk::Label* shl = Gtk::manage( new Gtk::Label (M("PREFERENCES_SHTHRESHOLD")+": ") );
    shThresh = Gtk::manage( new Gtk::SpinButton () );
    shThresh->show ();
    shThresh->set_digits (0);
    shThresh->set_increments (1, 10);
    shThresh->set_range (0, 255);
    vbsh->pack_start (*shl, Gtk::PACK_SHRINK, 0);
    vbsh->pack_end (*shThresh, Gtk::PACK_SHRINK, 0);
    vbrl->pack_start (*vbsh, Gtk::PACK_SHRINK, 0);

    frl->add (*vbrl);  
    hbcd->pack_start (*frl, true, true, 0);

//-----
    Gtk::VBox* dfpfvb = Gtk::manage( new Gtk::VBox () );
    dfpfvb->set_spacing (4);

    Gtk::Frame* fdf = Gtk::manage( new Gtk::Frame (M("PREFERENCES_DATEFORMATFRAME")) );

    Gtk::HBox* hb6 = Gtk::manage( new Gtk::HBox () );
    hb6->set_border_width (4);
    hb6->set_spacing (4);
    Gtk::Label* dflab = Gtk::manage( new Gtk::Label (M("PREFERENCES_DATEFORMAT")+":", Gtk::ALIGN_LEFT));
    hb6->pack_start (*dflab, Gtk::PACK_SHRINK,4);
    dateformat = Gtk::manage( new Gtk::Entry () );
    dateformat->set_tooltip_markup (M("PREFERENCES_DATEFORMATHINT"));
    dflab->set_tooltip_markup (M("PREFERENCES_DATEFORMATHINT"));
    hb6->pack_start (*dflab, Gtk::PACK_SHRINK, 0);
    hb6->pack_end (*dateformat, Gtk::PACK_SHRINK, 0);
    fdf->add (*hb6);

    dfpfvb->pack_start (*fdf, true, true, 0);

//-----
    Gtk::Frame* pff = Gtk::manage( new Gtk::Frame (M("PREFERENCES_PANFACTORFRAME")) );

    Gtk::HBox* pfhb = Gtk::manage( new Gtk::HBox () );
    pfhb->set_border_width(4);
    pfhb->set_spacing(4);
    Gtk::Label* pfl = Gtk::manage( new Gtk::Label (M("PREFERENCES_PANFACTORLABEL") + ":", Gtk::ALIGN_LEFT));
    panFactor = Gtk::manage( new Gtk::SpinButton () );
    panFactor->set_digits (0);
    panFactor->set_increments (1, 5);
    panFactor->set_range (1, 10);
    pfhb->pack_start (*pfl, Gtk::PACK_SHRINK, 0);
    pfhb->pack_end (*panFactor, Gtk::PACK_SHRINK, 0);
    pff->add (*pfhb);

    dfpfvb->pack_start (*pff, true, true, 0);
    hbcd->pack_start (*dfpfvb, Gtk::PACK_SHRINK, 4);
    mvbsd->pack_start (*hbcd, Gtk::PACK_SHRINK, 4);

  //-----
    Gtk::Frame* fdg = Gtk::manage( new Gtk::Frame (M("PREFERENCES_EXTERNALEDITOR")) );
    Gtk::VBox* dgvb = Gtk::manage( new Gtk::VBox () );

    Gtk::HBox* hb7c = Gtk::manage( new Gtk::HBox () );
    edOther = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_EDITORCMDLINE")+":"));
    hb7c->pack_start (*edOther, Gtk::PACK_SHRINK,4);
    editorToSendTo = Gtk::manage( new Gtk::Entry () );
    hb7c->pack_start (*editorToSendTo);
    Gtk::RadioButton::Group ge = edOther->get_group();
  
#ifdef __APPLE__
  Gtk::HBox* hb7 = Gtk::manage( new Gtk::HBox () );
  edGimp = Gtk::manage( new Gtk::RadioButton ("GIMP") );
  hb7->pack_start (*edGimp, Gtk::PACK_SHRINK,4);
  dgvb->pack_start (*hb7, Gtk::PACK_SHRINK, 4);
  edGimp->set_group (ge);
 
  Gtk::HBox* hb7b = Gtk::manage( new Gtk::HBox () );
  edPS = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_PSPATH")+":"));
  hb7b->pack_start (*edPS, Gtk::PACK_SHRINK,4);
  psDir = Gtk::manage( new Gtk::FileChooserButton (M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) );
  hb7b->pack_start (*psDir);
  dgvb->pack_start (*hb7b, Gtk::PACK_SHRINK, 4);
  edPS->set_group (ge);
#elif defined WIN32
  Gtk::HBox* hb7 = Gtk::manage( new Gtk::HBox () );
  edGimp = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_GIMPPATH")+":") );
  hb7->pack_start (*edGimp, Gtk::PACK_SHRINK,4);
  gimpDir = Gtk::manage( new Gtk::FileChooserButton (M("PREFERENCES_GIMPPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) );
  hb7->pack_start (*gimpDir);
  dgvb->pack_start (*hb7, Gtk::PACK_SHRINK, 4);
  edGimp->set_group (ge);
 
  Gtk::HBox* hb7b = Gtk::manage( new Gtk::HBox ());
  edPS = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_PSPATH")+":") );
  hb7b->pack_start (*edPS, Gtk::PACK_SHRINK,4);
  psDir = Gtk::manage( new Gtk::FileChooserButton (M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) );
  hb7b->pack_start (*psDir);
  dgvb->pack_start (*hb7b, Gtk::PACK_SHRINK, 4);
  edPS->set_group (ge);
#else
    Gtk::HBox* hb7 = Gtk::manage( new Gtk::HBox () );
    edGimp = Gtk::manage( new Gtk::RadioButton ("GIMP") );
    hb7->pack_start (*edGimp, Gtk::PACK_SHRINK,4);
    dgvb->pack_start (*hb7, Gtk::PACK_SHRINK, 4);
    edGimp->set_group (ge);
#endif

    dgvb->pack_start (*hb7c, Gtk::PACK_SHRINK, 4);
    dgvb->set_border_width (4);
    fdg->add (*dgvb);
    mvbsd->pack_start (*fdg, Gtk::PACK_SHRINK, 4);


    // Custom profile builder box
    Gtk::Frame* cpfrm = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CUSTPROFBUILD")) );

    Gtk::HBox* cphb = Gtk::manage( new Gtk::HBox () );
    cphb->set_border_width (4);
    cphb->set_spacing (4);

    Gtk::Label* cplab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CUSTPROFBUILDPATH")+":") );
    cphb->pack_start (*cplab, Gtk::PACK_SHRINK,4);

    txtCustProfBuilderPath = Gtk::manage( new Gtk::Entry () );
    txtCustProfBuilderPath->set_tooltip_markup (M("PREFERENCES_CUSTPROFBUILDHINT"));
    cphb->set_tooltip_markup (M("PREFERENCES_CUSTPROFBUILDHINT"));
    cphb->pack_start (*txtCustProfBuilderPath);
    
    cpfrm->add (*cphb);

    mvbsd->pack_start (*cpfrm, Gtk::PACK_SHRINK, 4);


    mvbsd->set_border_width (4);

    tconn = theme->signal_changed().connect( sigc::mem_fun(*this, &Preferences::themeChanged) );
    sconn = slimUI->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::themeChanged) );
    fconn = fontbutton->signal_font_set().connect( sigc::mem_fun(*this, &Preferences::fontChanged) );
    usethcon = chUseSystemTheme->signal_clicked ().connect( sigc::mem_fun(*this, &Preferences::useThemeChanged) );
    
    return mvbsd;
}

Gtk::Widget* Preferences::getFileBrowserPanel () {

    Gtk::VBox* mvbfb = Gtk::manage( new Gtk::VBox () );
    mvbfb->set_border_width (4);

    Gtk::Frame* fsd = Gtk::manage( new Gtk::Frame (M("PREFERENCES_STARTUPIMDIR")) );

    sdcurrent = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIRSOFTWARE")) );
    sdlast    = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIRLAST")) );
    sdhome    = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIRHOME")) );
    sdother   = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIROTHER")+": ") );
    startupdir = Gtk::manage( new Gtk::Entry () );

    Gtk::Button* sdselect = Gtk::manage( new Gtk::Button ("") );
    sdselect->set_image (*Gtk::manage(new RTImage ("gtk-open.png")));

    Gtk::RadioButton::Group opts = sdcurrent->get_group();
    sdlast->set_group (opts);
    sdhome->set_group (opts);
    sdother->set_group (opts);

    Gtk::VBox* vbsd = Gtk::manage( new Gtk::VBox () );
    vbsd->pack_start (*sdcurrent, Gtk::PACK_SHRINK,0);
    vbsd->pack_start (*sdlast, Gtk::PACK_SHRINK,0);
    vbsd->pack_start (*sdhome, Gtk::PACK_SHRINK,0);
    Gtk::HBox* otherbox = Gtk::manage( new Gtk::HBox () );
    otherbox->pack_start (*sdother, Gtk::PACK_SHRINK);
    otherbox->pack_start (*startupdir);
    otherbox->pack_end (*sdselect, Gtk::PACK_SHRINK, 4);
    vbsd->pack_start (*otherbox, Gtk::PACK_SHRINK, 0);
    vbsd->set_border_width (4);

    fsd->add (*vbsd);
    mvbfb->pack_start (*fsd, Gtk::PACK_SHRINK, 4);

    sdselect->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::selectStartupDir) );

//---


    Gtk::Frame* fro = Gtk::manage( new Gtk::Frame (M("PREFERENCES_FBROWSEROPTS")) );
    showDateTime = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SHOWDATETIME")) );
    showBasicExif = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SHOWBASICEXIF")) );
    showExpComp = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SHOWEXPOSURECOMPENSATION")) );
    Gtk::VBox* vbro = Gtk::manage( new Gtk::VBox () );
    Gtk::HBox* hbro1 = Gtk::manage( new Gtk::HBox () );
    overlayedFileNames = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_OVERLAY_FILENAMES")) );
	ckbInternalThumbIfUntouched = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_INTERNALTHUMBIFUNTOUCHED")));

    vbro->set_border_width (4);
    vbro->pack_start (*showDateTime, Gtk::PACK_SHRINK, 0);
    hbro1->pack_start (*showBasicExif, Gtk::PACK_SHRINK, 0);
    hbro1->pack_start (*showExpComp, Gtk::PACK_SHRINK, 4);
    vbro->pack_start (*hbro1, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*overlayedFileNames, Gtk::PACK_SHRINK, 0);
	vbro->pack_start (*ckbInternalThumbIfUntouched, Gtk::PACK_SHRINK, 0);

    fro->add (*vbro);  


    Gtk::Frame* frmnu = Gtk::manage( new Gtk::Frame (M("PREFERENCES_MENUOPTIONS")) );
    ckbmenuGroupRank = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPRANK")) );
    ckbmenuGroupLabel = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPLABEL")) );
    ckbmenuGroupFileOperations = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPFILEOPERATIONS")) );
    ckbmenuGroupProfileOperations = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPPROFILEOPERATIONS")) );
    ckbmenuGroupExtProg = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPEXTPROGS")) );
    Gtk::VBox* vbmnu = Gtk::manage( new Gtk::VBox () );

    vbmnu->set_border_width (4);
    vbmnu->pack_start (*ckbmenuGroupRank, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupLabel, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupFileOperations, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupProfileOperations, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupExtProg, Gtk::PACK_SHRINK, 0);

	frmnu->add (*vbmnu);


    Gtk::Frame* fre = Gtk::manage( new Gtk::Frame (M("PREFERENCES_PARSEDEXT")) );
    Gtk::VBox* vbre = Gtk::manage( new Gtk::VBox () );
    vbre->set_border_width (4);
    Gtk::HBox* hb0 = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* elab = Gtk::manage( new Gtk::Label (M("PREFERENCES_PARSEDEXTADD")+":") );
    hb0->pack_start (*elab, Gtk::PACK_SHRINK, 4);
    extension = Gtk::manage( new Gtk::Entry () );
    extension->set_width_chars(5);
    hb0->pack_start (*extension);
    addExt = Gtk::manage( new Gtk::Button () );
    delExt = Gtk::manage( new Gtk::Button () );
    addExt->set_tooltip_text (M("PREFERENCES_PARSEDEXTADDHINT"));
    delExt->set_tooltip_text (M("PREFERENCES_PARSEDEXTDELHINT"));
    Gtk::Image* addExtImg = Gtk::manage( new RTImage ("list-add-small.png") );
    Gtk::Image* delExtImg = Gtk::manage( new RTImage ("list-remove-red-small.png") );
    addExt->add (*addExtImg);
    delExt->add (*delExtImg);
    hb0->pack_end (*delExt, Gtk::PACK_SHRINK, 4);
    hb0->pack_end (*addExt, Gtk::PACK_SHRINK, 4);
    extensions = Gtk::manage( new Gtk::TreeView () );
    Gtk::ScrolledWindow* hscrollw = Gtk::manage( new Gtk::ScrolledWindow () );
    hscrollw->set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    hscrollw->add (*extensions);
    extensionModel = Gtk::ListStore::create (extensionColumns);
    extensions->set_model (extensionModel);
    extensions->append_column_editable("Enabled", extensionColumns.enabled);
    extensions->append_column("Extension", extensionColumns.ext);
    extensions->set_headers_visible (false);  
    vbre->pack_start (*hscrollw);
    vbre->pack_start (*hb0, Gtk::PACK_SHRINK, 4);

    fre->add (*vbre);  

    Gtk::Frame* frc = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CACHEOPTS")) );
    Gtk::VBox* vbc = Gtk::manage( new Gtk::VBox () );
    frc->add (*vbc);  
    vbc->set_border_width (4);

    Gtk::Label* cflab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CACHETHUMBFORM")+":") );
    cflab->set_alignment(Gtk::ALIGN_LEFT, Gtk::ALIGN_TOP);
    cformat = Gtk::manage( new Gtk::ComboBoxText () );
    cformat->set_size_request(50, -1);
    cformat->append_text (M("PREFERENCES_CACHEFORMAT1"));
    cformat->append_text (M("PREFERENCES_CACHEFORMAT2"));
    cformat->append_text (M("PREFERENCES_CACHEFORMAT1")+", 16 bit");
    cformat->signal_changed().connect( sigc::mem_fun(*this, &Preferences::cacheFormatComboChanged) );
    cacheFormatComboChanged(); // update the tooltip
    vbc->pack_start (*cflab, Gtk::PACK_SHRINK, 2);
    vbc->pack_start (*cformat);

    Gtk::HBox* hb3 = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* chlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CACHETHUMBHEIGHT")+":") );
    maxThumbSize = Gtk::manage( new Gtk::SpinButton () );
    hb3->pack_start (*chlab, Gtk::PACK_SHRINK, 4);
    hb3->pack_start (*maxThumbSize, Gtk::PACK_SHRINK, 4);

    maxThumbSize->set_digits (0);
    maxThumbSize->set_increments (1, 10);
    maxThumbSize->set_range (40, 800);
    vbc->pack_start (*hb3, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hb4 = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* celab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CACHEMAXENTRIES")+":") );
    maxCacheEntries = Gtk::manage( new Gtk::SpinButton () );
    hb4->pack_start (*celab, Gtk::PACK_SHRINK, 4);
    hb4->pack_start (*maxCacheEntries, Gtk::PACK_SHRINK, 4);

    maxCacheEntries->set_digits (0);
    maxCacheEntries->set_increments (1, 10);
    maxCacheEntries->set_range (10, 100000);
    vbc->pack_start (*hb4, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hb5 = Gtk::manage( new Gtk::HBox () );
    clearThumbnails = Gtk::manage( new Gtk::Button (M("PREFERENCES_CACHECLEARTHUMBS")) );
    clearProfiles = Gtk::manage( new Gtk::Button (M("PREFERENCES_CACHECLEARPROFILES")) );
    clearAll = Gtk::manage( new Gtk::Button (M("PREFERENCES_CACHECLEARALL")) );
    hb5->pack_start (*clearThumbnails, Gtk::PACK_SHRINK, 4);
    hb5->pack_start (*clearProfiles, Gtk::PACK_SHRINK, 4);
    hb5->pack_start (*clearAll, Gtk::PACK_SHRINK, 4);
    vbc->pack_start (*hb5, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hb6 = Gtk::manage( new Gtk::HBox () );
    Gtk::VBox* vb6 = Gtk::manage( new Gtk::VBox () );

    vb6->pack_start (*fro);
    vb6->pack_start (*frmnu);
    vb6->pack_end (*frc);
    hb6->pack_start (*vb6);
    hb6->pack_start (*fre);
    hb6->set_spacing(4);

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

Gtk::Widget* Preferences::getSoundPanel () {
    Gtk::VBox* pSnd = new Gtk::VBox ();

    ckbSndEnable = Gtk::manage( new Gtk::CheckButton (M("GENERAL_ENABLE")));
    sndEnableConn  = ckbSndEnable->signal_toggled().connect (sigc::mem_fun(*this, &Preferences::sndEnableToggled));

    pSnd->pack_start (*ckbSndEnable, Gtk::PACK_SHRINK, 8);

    Gtk::Label* lSndHelp = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_HELP")));
    pSnd->pack_start (*lSndHelp, Gtk::PACK_SHRINK, 4);

    // BatchQueueDone
    Gtk::HBox* pBatchQueueDone = Gtk::manage( new Gtk::HBox() );

    Gtk::Label* lSndBatchQueueDone = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_BATCHQUEUEDONE") + Glib::ustring(":")));
    pBatchQueueDone->pack_start (*lSndBatchQueueDone, Gtk::PACK_SHRINK, 12);
    
    txtSndBatchQueueDone =  Gtk::manage (new Gtk::Entry());
    pBatchQueueDone->pack_end (*txtSndBatchQueueDone, Gtk::PACK_EXPAND_WIDGET, 4);
    
    pSnd->pack_start (*pBatchQueueDone, Gtk::PACK_SHRINK, 4);

    // LngEditProcDone
    Gtk::HBox* pSndLngEditProcDone = Gtk::manage( new Gtk::HBox() );

    Gtk::Label* lSndLngEditProcDone = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_LNGEDITPROCDONE") + Glib::ustring(":")));
    pSndLngEditProcDone->pack_start (*lSndLngEditProcDone, Gtk::PACK_SHRINK, 12);
    
    txtSndLngEditProcDone =  Gtk::manage (new Gtk::Entry());
    pSndLngEditProcDone->pack_start (*txtSndLngEditProcDone, Gtk::PACK_EXPAND_WIDGET, 4);

    Gtk::Label* lSndLngEditProcDoneSecs = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_TRESHOLDSECS") + Glib::ustring(":")));
    pSndLngEditProcDone->pack_start (*lSndLngEditProcDoneSecs, Gtk::PACK_SHRINK, 12);
 
    spbSndLngEditProcDoneSecs = Gtk::manage( new Gtk::SpinButton () );
    spbSndLngEditProcDoneSecs->set_digits (1);
    spbSndLngEditProcDoneSecs->set_increments (0.5, 1);
    spbSndLngEditProcDoneSecs->set_range (0, 10);
    pSndLngEditProcDone->pack_end (*spbSndLngEditProcDoneSecs, Gtk::PACK_SHRINK, 4);

    pSnd->pack_start (*pSndLngEditProcDone, Gtk::PACK_SHRINK, 4);

    pSnd->set_border_width (4);

    sndEnableToggled();

    return pSnd;
}

void Preferences::parseDir (Glib::ustring dirname, std::vector<Glib::ustring>& items, Glib::ustring ext) {

    if (dirname.empty())
        return;

    // process directory
    Glib::Dir* dir = NULL;
    try {
        dir = new Glib::Dir (dirname);
    }
    catch (const Glib::FileError& fe) {
        return;
    }
    for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
      Glib::ustring fname = Glib::build_filename(dirname, *i);
      Glib::ustring sname = *i;
      // ignore directories
      if (!safe_file_test (fname, Glib::FILE_TEST_IS_DIR) && sname.size() >= ext.size() && sname.substr (sname.size()-ext.size(), ext.size()).casefold() == ext) 
            items.push_back (sname.substr(0,sname.size()-ext.size()));
    }
    delete dir;
}

void Preferences::storePreferences () {

    moptions.defProfRaw          = rprofiles->get_active_text();
    if (moptions.defProfRaw.empty()) moptions.defProfRaw = DEFPROFILE_RAW;
    moptions.defProfImg          = iprofiles->get_active_text();
    if (moptions.defProfImg.empty()) moptions.defProfImg = DEFPROFILE_IMG;
    moptions.dateFormat          = dateformat->get_text();
    moptions.panAccelFactor      = (int)panFactor->get_value();
    moptions.fbShowDateTime  = showDateTime->get_active ();
    moptions.fbShowBasicExif = showBasicExif->get_active ();
    moptions.fbShowExpComp   = showExpComp->get_active ();
    moptions.menuGroupRank              = ckbmenuGroupRank->get_active();
    moptions.menuGroupLabel             = ckbmenuGroupLabel->get_active();
    moptions.menuGroupFileOperations    = ckbmenuGroupFileOperations->get_active();
    moptions.menuGroupProfileOperations = ckbmenuGroupProfileOperations->get_active();
    moptions.menuGroupExtProg           = ckbmenuGroupExtProg->get_active();
    moptions.blinkClipped    = blinkClipped->get_active ();
    moptions.highlightThreshold = (int)hlThresh->get_value ();
    moptions.shadowThreshold = (int)shThresh->get_value ();
    moptions.language        = languages->get_active_text ();
    moptions.languageAutoDetect = ckbLangAutoDetect->get_active ();
    moptions.theme           = theme->get_active_text ();
    moptions.slimUI          = slimUI->get_active ();
    moptions.useSystemTheme  = chUseSystemTheme->get_active ();
     
    Gdk::Color cropCol=butCropCol->get_color();
    moptions.cutOverlayBrush[0]=cropCol.get_red_p();
    moptions.cutOverlayBrush[1]=cropCol.get_green_p();
    moptions.cutOverlayBrush[2]=cropCol.get_blue_p();
    moptions.cutOverlayBrush[3]=butCropCol->get_alpha()/65535.0;

    moptions.font            = fontbutton->get_font_name();
#ifdef WIN32    
    moptions.gimpDir        = gimpDir->get_filename ();
    moptions.psDir          = psDir->get_filename ();
#elif defined __APPLE__
    moptions.psDir          = psDir->get_filename (); 
#endif
    moptions.customEditorProg = editorToSendTo->get_text ();
    if (edGimp->get_active ())
        moptions.editorToSendTo = 1;
#ifdef WIN32    
    else if (edPS->get_active ())
        moptions.editorToSendTo = 2;
#elif defined __APPLE__   
    else if (edPS->get_active ())
      moptions.editorToSendTo = 2; 
#endif
    else if (edOther->get_active ())
        moptions.editorToSendTo = 3;

    moptions.customProfileBuilder = txtCustProfBuilderPath->get_text();

    moptions.rtSettings.monitorProfile      = monProfile->get_filename ();
    moptions.rtSettings.autoMonitorProfile  = cbAutoMonProfile->get_active ();
    moptions.rtSettings.iccDirectory        = iccDir->get_filename ();
    moptions.rtSettings.colorimetricIntent  = intent->get_active_row_number ();
    moptions.rtSettings.viewingdevice       = view->get_active_row_number ();
    moptions.rtSettings.viewingdevicegrey   = grey->get_active_row_number ();

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
    for (size_t i=0; i<c.size(); i++) {
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
    moptions.internalThumbIfUntouched = ckbInternalThumbIfUntouched->get_active ();
    
    moptions.saveParamsFile = saveParamsFile->get_active ();
    moptions.saveParamsCache = saveParamsCache->get_active ();
    moptions.paramsLoadLocation = (PPLoadLocation)loadParamsPreference->get_active_row_number ();

    moptions.tunnelMetaData = ckbTunnelMetaData->get_active ();

    moptions.rtSettings.darkFramesPath =   darkFrameDir->get_filename();
    moptions.rtSettings.flatFieldsPath =   flatFieldDir->get_filename();

    moptions.baBehav.resize (ADDSET_PARAM_NUM);
    for (Gtk::TreeIter sections=behModel->children().begin();  sections!=behModel->children().end(); sections++)
        for (Gtk::TreeIter adjs=sections->children().begin();  adjs!=sections->children().end(); adjs++) 
            moptions.baBehav[adjs->get_value (behavColumns.addsetid)] = adjs->get_value (behavColumns.badd);

    int editorMode=editorLayout->get_active_row_number();
    moptions.tabbedUI = (editorMode>1);
    moptions.multiDisplayMode = editorMode==3 ? 1:0;
    moptions.mainNBVertical = editorMode==1;

    moptions.histogramPosition = ckbHistogramPositionLeft->get_active() ? 1 : 2;
    moptions.showProfileSelector = ckbShowProfileSelector->get_active();
    moptions.squareDetailWindow = ckbSquareDetailWindow->get_active();
    moptions.FileBrowserToolbarSingleRow = ckbFileBrowserToolbarSingleRow->get_active();
    moptions.hideTPVScrollbar = ckbHideTPVScrollbar->get_active();
    moptions.overwriteOutputFile = chOverwriteOutputFile->get_active ();
    moptions.UseIconNoText = ckbUseIconNoText->get_active();

    // Sounds
    moptions.sndEnable = ckbSndEnable->get_active ();
    moptions.sndBatchQueueDone = txtSndBatchQueueDone->get_text ();
    moptions.sndLngEditProcDone     = txtSndLngEditProcDone->get_text ();
    moptions.sndLngEditProcDoneSecs = spbSndLngEditProcDoneSecs->get_value ();
}

void Preferences::fillPreferences () {

    tconn.block (true);
    sconn.block (true);
    dfconn.block (true);
    ffconn.block (true);

    rprofiles->set_active_text (moptions.defProfRaw);
    iprofiles->set_active_text (moptions.defProfImg);
    dateformat->set_text (moptions.dateFormat);
    panFactor->set_value(moptions.panAccelFactor);
    if (safe_file_test (moptions.rtSettings.monitorProfile, Glib::FILE_TEST_EXISTS)) 
        monProfile->set_filename (moptions.rtSettings.monitorProfile);
    if (moptions.rtSettings.monitorProfile.empty())
        monProfile->set_current_folder (moptions.rtSettings.iccDirectory);
    cbAutoMonProfile->set_active(moptions.rtSettings.autoMonitorProfile);

    if (Glib::file_test (moptions.rtSettings.iccDirectory, Glib::FILE_TEST_IS_DIR)) 
        iccDir->set_current_folder (moptions.rtSettings.iccDirectory);
    intent->set_active (moptions.rtSettings.colorimetricIntent);
    view->set_active (moptions.rtSettings.viewingdevice);
    grey->set_active (moptions.rtSettings.viewingdevicegrey);

    languages->set_active_text (moptions.language);
    ckbLangAutoDetect->set_active (moptions.languageAutoDetect);
    theme->set_active_text (moptions.theme);
    slimUI->set_active(moptions.slimUI);
    chUseSystemTheme->set_active(moptions.useSystemTheme);

    Gdk::Color cropCol;
    cropCol.set_rgb_p(moptions.cutOverlayBrush[0],moptions.cutOverlayBrush[1],moptions.cutOverlayBrush[2]);
    butCropCol->set_color(cropCol);
    butCropCol->set_alpha ( (unsigned short)(moptions.cutOverlayBrush[3]*65535.0));

    fontbutton->set_font_name(moptions.font);
    showDateTime->set_active (moptions.fbShowDateTime);
    showBasicExif->set_active (moptions.fbShowBasicExif);
    showExpComp->set_active (moptions.fbShowExpComp);
    ckbmenuGroupRank->set_active(moptions.menuGroupRank);
    ckbmenuGroupLabel->set_active(moptions.menuGroupLabel);
    ckbmenuGroupFileOperations->set_active(moptions.menuGroupFileOperations);
    ckbmenuGroupProfileOperations->set_active(moptions.menuGroupProfileOperations);
    ckbmenuGroupExtProg->set_active(moptions.menuGroupExtProg);

    blinkClipped->set_active (moptions.blinkClipped);
    hlThresh->set_value (moptions.highlightThreshold);
    shThresh->set_value (moptions.shadowThreshold);

    edGimp->set_active (moptions.editorToSendTo==1);
    edOther->set_active (moptions.editorToSendTo==3);
#ifdef WIN32    
    edPS->set_active (moptions.editorToSendTo==2);
    if (safe_file_test (moptions.gimpDir, Glib::FILE_TEST_IS_DIR)) 
        gimpDir->set_filename (moptions.gimpDir);
    if (safe_file_test (moptions.psDir, Glib::FILE_TEST_IS_DIR)) 
        psDir->set_filename (moptions.psDir);
#elif defined __APPLE__
  edPS->set_active (moptions.editorToSendTo==2);
  if (safe_file_test (moptions.psDir, Glib::FILE_TEST_IS_DIR))
    psDir->set_filename (moptions.psDir); 
#endif
    editorToSendTo->set_text (moptions.customEditorProg);

    txtCustProfBuilderPath->set_text(moptions.customProfileBuilder);

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
    for (size_t i=0; i<moptions.parseExtensions.size(); i++) {
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
    ckbInternalThumbIfUntouched->set_active(moptions.internalThumbIfUntouched);
    
    saveParamsFile->set_active (moptions.saveParamsFile);
    saveParamsCache->set_active (moptions.saveParamsCache);
    loadParamsPreference->set_active (moptions.paramsLoadLocation);    

    ckbTunnelMetaData->set_active (moptions.tunnelMetaData); 

    if (!moptions.tabbedUI)
        editorLayout->set_active(moptions.mainNBVertical ? 1 : 0);
    else 
        editorLayout->set_active(moptions.multiDisplayMode ? 3 : 2);

    ckbHistogramPositionLeft->set_active(moptions.histogramPosition==1);
    ckbShowProfileSelector->set_active(moptions.showProfileSelector);
    ckbSquareDetailWindow->set_active(moptions.squareDetailWindow);
    ckbFileBrowserToolbarSingleRow->set_active(moptions.FileBrowserToolbarSingleRow);
    ckbHideTPVScrollbar->set_active(moptions.hideTPVScrollbar);
    ckbUseIconNoText->set_active(moptions.UseIconNoText);

    //darkFrameDir->set_filename( moptions.rtSettings.darkFramesPath );
    //updateDFinfos();
    darkFrameDir->set_current_folder( moptions.rtSettings.darkFramesPath );
    darkFrameChanged ();
    
    //flatFieldDir->set_filename( moptions.rtSettings.flatFieldsPath );
    //updateFFinfos();
    flatFieldDir->set_current_folder( moptions.rtSettings.flatFieldsPath );
    flatFieldChanged ();

    addc.block (true);
    setc.block (true);
    if (moptions.baBehav.size() == ADDSET_PARAM_NUM) {
        for (size_t i=0; i<moptions.baBehav.size(); i++)
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
    sconn.block (false);
    dfconn.block (false);
    ffconn.block (false);

    chOverwriteOutputFile->set_active (moptions.overwriteOutputFile);

    // Sounds
    ckbSndEnable->set_active (moptions.sndEnable);
    txtSndBatchQueueDone->set_text (moptions.sndBatchQueueDone);
    txtSndLngEditProcDone->set_text (moptions.sndLngEditProcDone);
    spbSndLngEditProcDoneSecs->set_value (moptions.sndLngEditProcDoneSecs);
}

/*
void Preferences::loadPressed () {

    moptions.copyFrom (&options);
    fillPreferences ();
}

void Preferences::savePressed () {

    storePreferences ();
    options.copyFrom (&moptions);
    Options::save ();
}
*/

void Preferences::autoMonProfileToggled () {
	monProfile->set_sensitive(!cbAutoMonProfile->get_active());
}

void Preferences::sndEnableToggled () {
	txtSndBatchQueueDone->set_sensitive(ckbSndEnable->get_active());
	txtSndLngEditProcDone->set_sensitive(ckbSndEnable->get_active());
	spbSndLngEditProcDoneSecs->set_sensitive(ckbSndEnable->get_active());
}

void Preferences::langAutoDetectToggled () {
	languages->set_sensitive(!ckbLangAutoDetect->get_active());
}

void Preferences::okPressed () {

    storePreferences ();
    workflowUpdate();
    options.copyFrom (&moptions);   
    options.filterOutParsedExtensions();
    Options::save ();
    hide ();
}

void Preferences::cancelPressed () {

	// set the initial theme back
	if (theme->get_active_text () != options.theme) {
		RTImage::setPaths(options);
		RTImage::updateImages();
		switchThemeTo(options.theme, options.slimUI);
	}

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

    splash = new Splash (*this);
    splash->set_transient_for (*this);
    splash->signal_delete_event().connect( sigc::mem_fun(*this, &Preferences::splashClosed) );
    splash->show ();
}

void Preferences::themeChanged () {

	moptions.theme = theme->get_active_text ();
	moptions.useSystemTheme = chUseSystemTheme->get_active ();
	RTImage::setPaths(moptions);
	RTImage::updateImages();
	switchThemeTo(theme->get_active_text (), slimUI->get_active());
}

void Preferences::forRAWComboChanged () {
	rprofiles->set_tooltip_text(rprofiles->get_active_text());
}

void Preferences::forImageComboChanged () {
	iprofiles->set_tooltip_text(iprofiles->get_active_text());
}

void Preferences::layoutComboChanged () {
	editorLayout->set_tooltip_text(editorLayout->get_active_text());
}

void Preferences::cacheFormatComboChanged () {
	cformat->set_tooltip_text(cformat->get_active_text());
}


void Preferences::fontChanged () {

	switchFontTo(fontbutton->get_font_name());
}

void Preferences::switchThemeTo(Glib::ustring newTheme, bool slimInterface) {

	std::vector<Glib::ustring> files;
	files.push_back (argv0+"/themes/"+newTheme+".gtkrc");
	if (slimInterface)
		files.push_back (argv0+"/themes/slim");
	Gtk::RC::set_default_files (files);

#ifndef WIN32
   // For an unknown reason, gtkmm 2.22 don't know the gtk-button-images property, while it exists in the documentation...
   // Anyway, the problem was Linux only
   static Glib::RefPtr<Gtk::Settings> settings = Gtk::Settings::get_default();
   if (settings)
      settings->property_gtk_button_images().set_value(true);
   else
      printf("Error: no default settings to update!\n");
#endif

   Gtk::RC::reparse_all (Gtk::Settings::get_default());
   GdkEventClient event = { GDK_CLIENT_EVENT, NULL, TRUE, gdk_atom_intern("_GTK_READ_RCFILES", FALSE), 8 };
   gdk_event_send_clientmessage_toall ((GdkEvent*)&event);
}

void Preferences::workflowUpdate (){

    if(moptions.tabbedUI != options.tabbedUI) {
        parent->MoveFileBrowserToMain();
        parent->SetMainCurrent();
        if(moptions.tabbedUI)
            parent->epanel->hide_all();
        else
           parent->epanel->show_all();
    }
    if(moptions.hideTPVScrollbar != options.hideTPVScrollbar) {
    	// Update the tool panels
   		parent->updateTPVScrollbar (moptions.hideTPVScrollbar);
    }
    if(moptions.UseIconNoText != options.UseIconNoText) {
    	// Update the tool's tab titles
    	parent->updateTabsUsesIcons(moptions.UseIconNoText);
    }
    if(moptions.FileBrowserToolbarSingleRow != options.FileBrowserToolbarSingleRow) {
    	// Update the position of the Query toolbar
    	parent->updateFBQueryTB(moptions.FileBrowserToolbarSingleRow);
    }
    if(moptions.histogramPosition != options.histogramPosition) {
    	// Update the position of the Histogram
    	parent->updateHistogramPosition(options.histogramPosition, moptions.histogramPosition);
    }
    if(moptions.showProfileSelector != options.showProfileSelector) {
    	// Update the position of the Profile selector
    	parent->updateTPProfileSelector(moptions.showProfileSelector);
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
  for (size_t i=0; i<c.size(); i++)
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
	//Glib::ustring s(darkFrameDir->get_filename());
	Glib::ustring s(darkFrameDir->get_current_folder());
	//if( s.compare( rtengine::dfm.getPathname()) !=0 ){
	   rtengine::dfm.init( s );
	   updateDFinfos();
	//}
}

void Preferences::flatFieldChanged ()
{
	//Glib::ustring s(flatFieldDir->get_filename());
	Glib::ustring s(flatFieldDir->get_current_folder());
	//if( s.compare( rtengine::ffm.getPathname()) !=0 ){
	   rtengine::ffm.init( s );
	   updateFFinfos();
	//}
}

void Preferences::updateDFinfos()
{
    int t1,t2;
    rtengine::dfm.getStat(t1,t2);
    Glib::ustring s = Glib::ustring::compose("%1: %2 %3, %4 %5", M("PREFERENCES_DARKFRAMEFOUND"), t1, M("PREFERENCES_DARKFRAMESHOTS"), t2, M("PREFERENCES_DARKFRAMETEMPLATES"));
    dfLabel->set_text(s);
}

void Preferences::updateFFinfos()
{
    int t1,t2;
    rtengine::ffm.getStat(t1,t2);
    Glib::ustring s = Glib::ustring::compose("%1: %2 %3, %4 %5", M("PREFERENCES_FLATFIELDFOUND"), t1, M("PREFERENCES_FLATFIELDSHOTS"), t2, M("PREFERENCES_FLATFIELDTEMPLATES"));
    ffLabel->set_text(s);
}

bool Preferences::splashClosed(GdkEventAny* event) {
	delete splash;
	splash = NULL;
	return true;
}
