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
#include <sigc++/slot.h>
#include "preferences.h"
#include "multilangmgr.h"
#include "splash.h"
#include "cachemanager.h"
#include "addsetids.h"
#include "../rtengine/dfmanager.h"
#include "../rtengine/ffmanager.h"
#include <sstream>
#include "rtimage.h"
#ifdef _OPENMP
#include <omp.h>
#endif

extern Options options;
extern Glib::ustring argv0;
Glib::RefPtr<Gtk::CssProvider> themecss;
Glib::RefPtr<Gtk::CssProvider> fontcss;

Preferences::Preferences  (RTWindow *rtwindow)
    : Gtk::Dialog (M("MAIN_BUTTON_PREFERENCES"), *rtwindow, true)
    , splash (nullptr)
    , rprofiles (nullptr)
    , iprofiles (nullptr)
    , parent (rtwindow)
{
    regex = Glib::Regex::create(THEMEREGEXSTR, Glib::RegexCompileFlags::REGEX_CASELESS);

    moptions.copyFrom (&options);

    /*
     * Do not increase height, since it's not visible on e.g. smaller netbook
     * screens. The default height is about 620 pixels currently, that's why
     * we do not set the height anymore. Netbook users will most certainly set
     * a smaller font, so they'll be able to shrink the Preferences window and
     * close it.
     */
    set_size_request (650, -1);
    set_default_size (options.preferencesWidth, options.preferencesHeight);

    Gtk::Box* mainBox = get_content_area ();
//GTK318
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
    mainBox->set_spacing(8);
#endif
//GTK318
    //set_has_separator (false);

    Gtk::Notebook* nb = Gtk::manage (new Gtk::Notebook ());
    nb->set_name ("PrefNotebook");
    mainBox->pack_start (*nb);

    Gtk::Button* about  = Gtk::manage (new Gtk::Button (M("GENERAL_ABOUT")));
    Gtk::Button* ok     = Gtk::manage (new Gtk::Button (M("GENERAL_OK")));
    Gtk::Button* cancel = Gtk::manage (new Gtk::Button (M("GENERAL_CANCEL")));

    about->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::aboutPressed) );
    ok->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::okPressed) );
    cancel->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::cancelPressed) );

    get_action_area()->pack_start (*about);
    get_action_area()->pack_end (*ok);
    get_action_area()->pack_end (*cancel);

    nb->append_page (*getGeneralPanel(),        M("PREFERENCES_TAB_GENERAL"));
    nb->append_page (*getProcParamsPanel(),     M("PREFERENCES_TAB_IMPROC"));
    nb->append_page (*getFileBrowserPanel(),    M("PREFERENCES_TAB_BROWSER"));
    nb->append_page (*getColorManagementPanel(), M("PREFERENCES_TAB_COLORMGR"));
    nb->append_page (*getBatchProcPanel(),      M("PREFERENCES_BATCH_PROCESSING"));
    nb->append_page (*getPerformancePanel(),    M("PREFERENCES_TAB_PERFORMANCE"));
    // Sounds only on Windows and Linux
#if defined(WIN32) || defined(__linux__)
    nb->append_page (*getSoundPanel(),          M("PREFERENCES_TAB_SOUND"));
#endif
    nb->set_current_page (0);

    profileStore.addListener(this);

    fillPreferences ();

    show_all_children ();
}


Preferences::~Preferences ()
{

    profileStore.removeListener(this);
    get_size(options.preferencesWidth, options.preferencesHeight);
}

int Preferences::getThemeRowNumber(Glib::ustring& longThemeFName)
{

    if (regex->match(longThemeFName + ".css", matchInfo)) {
        for (size_t i=0 ; i<themeFNames.size(); ++i) {
            if (themeFNames.at(i).longFName == longThemeFName) {
                return (int)i;
            }
        }
    }
    return -1;
}

Gtk::Widget* Preferences::getBatchProcPanel ()
{

    Gtk::VBox* mvbpp = Gtk::manage (new Gtk::VBox ());

    Gtk::ScrolledWindow* behscrollw = Gtk::manage (new Gtk::ScrolledWindow ());
    behscrollw->set_policy (Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    behscrollw->set_size_request(-1, 60);
    Gtk::VBox* vbbeh = Gtk::manage( new Gtk::VBox () );
    vbbeh->pack_start (*behscrollw, Gtk::PACK_EXPAND_WIDGET);
    Gtk::Frame* behFrame = Gtk::manage (new Gtk::Frame (M("PREFERENCES_BEHAVIOR")));
    behFrame->add (*vbbeh);
    //mvbpp->pack_start (*behFrame);
    mvbpp->pack_start (*behFrame, Gtk::PACK_EXPAND_WIDGET, 4);
    Gtk::TreeView* behTreeView = Gtk::manage (new Gtk::TreeView ());
    behscrollw->add (*behTreeView);

    behModel = Gtk::TreeStore::create (behavColumns);
    behTreeView->set_model (behModel);

    behTreeView->append_column (M("PREFERENCES_PROPERTY"), behavColumns.label);
    behTreeView->append_column_editable (M("PREFERENCES_ADD"), behavColumns.badd);
    behTreeView->append_column_editable (M("PREFERENCES_SET"), behavColumns.bset);

    Gtk::CellRendererToggle* cr_add = static_cast<Gtk::CellRendererToggle*> (behTreeView->get_column (1)->get_first_cell());
    Gtk::CellRendererToggle* cr_set = static_cast<Gtk::CellRendererToggle*> (behTreeView->get_column (2)->get_first_cell());

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
    mi->set_value (behavColumns.label, M("TP_RETINEX_LABEL"));
    appendBehavList (mi, M("TP_RETINEX_STRENGTH"), ADDSET_RETI_STR, false);
    appendBehavList (mi, M("TP_RETINEX_NEIGHBOR"), ADDSET_RETI_NEIGH, false);
    appendBehavList (mi, M("TP_RETINEX_VARIANCE"), ADDSET_RETI_VART, false);
    appendBehavList (mi, M("TP_RETINEX_GAMMA"), ADDSET_RETI_GAM, false);
    appendBehavList (mi, M("TP_RETINEX_SLOPE"), ADDSET_RETI_SLO, false);
    appendBehavList (mi, M("TP_RETINEX_GAIN"), ADDSET_RETI_GAIN, false);
    appendBehavList (mi, M("TP_RETINEX_OFFSET"), ADDSET_RETI_OFFS, false);
    appendBehavList (mi, M("TP_RETINEX_THRESHOLD"), ADDSET_RETI_LIMD, false);

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
    //  appendBehavList (mi, M("TP_DIRPYRDENOISE_LUMA")+", "+M("TP_DIRPYRDENOISE_CHROMA"), ADDSET_DIRPYRDN_CHLUM, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_LUMA"), ADDSET_DIRPYRDN_LUMA, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_LDETAIL"), ADDSET_DIRPYRDN_LUMDET, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_CHROMA"), ADDSET_DIRPYRDN_CHROMA, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_RED"), ADDSET_DIRPYRDN_CHROMARED, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_BLUE"), ADDSET_DIRPYRDN_CHROMABLUE, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_GAMMA"), ADDSET_DIRPYRDN_GAMMA, true);
    appendBehavList (mi, M("TP_DIRPYRDENOISE_PASSES"), ADDSET_DIRPYRDN_PASSES, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_WBALANCE_LABEL"));
    appendBehavList (mi, M("TP_WBALANCE_TEMPERATURE"), ADDSET_WB_TEMPERATURE, true);
    appendBehavList (mi, M("TP_WBALANCE_GREEN"), ADDSET_WB_GREEN, true);
    appendBehavList (mi, M("TP_WBALANCE_EQBLUERED"), ADDSET_WB_EQUAL, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_COLORAPP_LABEL"));
    appendBehavList (mi, M("TP_COLORAPP_CIECAT_DEGREE"), ADDSET_CAT_DEGREE, true);
    appendBehavList (mi, M("TP_COLORAPP_ADAPTSCENE"), ADDSET_CAT_ADAPTSCENE, true);
    appendBehavList (mi, M("TP_COLORAPP_LIGHT"), ADDSET_CAT_LIGHT, true);
    appendBehavList (mi, M("TP_COLORAPP_BRIGHT"), ADDSET_CAT_BRIGHT, true);
    appendBehavList (mi, M("TP_COLORAPP_CHROMA"), ADDSET_CAT_CHROMA, true);
    appendBehavList (mi, M("TP_COLORAPP_RSTPRO"), ADDSET_CAT_RSTPRO, true);
    appendBehavList (mi, M("TP_COLORAPP_CONTRAST"), ADDSET_CAT_CONTRAST, true);
    appendBehavList (mi, M("TP_COLORAPP_CONTRAST_Q"), ADDSET_CAT_CONTRAST_Q, true);
    appendBehavList (mi, M("TP_COLORAPP_CHROMA_S"), ADDSET_CAT_CHROMA_S, true);
    appendBehavList (mi, M("TP_COLORAPP_CHROMA_M"), ADDSET_CAT_CHROMA_M, true);
    appendBehavList (mi, M("TP_COLORAPP_HUE"), ADDSET_CAT_HUE, true);
    appendBehavList (mi, M("TP_COLORAPP_ADAPTVIEWING"), ADDSET_CAT_ADAPTVIEWING, true);
    appendBehavList (mi, M("TP_COLORAPP_BADPIXSL"), ADDSET_CAT_BADPIX, true);

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
    appendBehavList (mi, M("TP_CHMIXER_RED") + ", " + M("TP_CHMIXER_GREEN") + ", " + M("TP_CHMIXER_BLUE"), ADDSET_CHMIXER, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_BWMIX_LABEL"));
    appendBehavList (mi, M("TP_BWMIX_MIXC"), ADDSET_BLACKWHITE_HUES, false);
    appendBehavList (mi, M("TP_BWMIX_GAMMA"), ADDSET_BLACKWHITE_GAMMA, false);

    mi = behModel->append ();
    mi->set_value( behavColumns.label, M("TP_FILMSIMULATION_LABEL") );
    appendBehavList( mi, M( "TP_FILMSIMULATION_STRENGTH" ), ADDSET_FILMSIMULATION_STRENGTH, true );

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_COLORTONING_LABEL"));
    appendBehavList (mi, M("TP_COLORTONING_SPLITCOCO"), ADDSET_COLORTONING_SPLIT , true);
    appendBehavList (mi, M("TP_COLORTONING_SATURATIONTHRESHOLD"), ADDSET_COLORTONING_SATTHRESHOLD , true);
    appendBehavList (mi, M("TP_COLORTONING_SATURATEDOPACITY"), ADDSET_COLORTONING_SATOPACITY , true);
    appendBehavList (mi, M("TP_COLORTONING_BALANCE"), ADDSET_COLORTONING_BALANCE , true);
    appendBehavList (mi, M("TP_COLORTONING_STRENGTH"), ADDSET_COLORTONING_STRENGTH , true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_ROTATE_LABEL"));
    appendBehavList (mi, M("TP_ROTATE_DEGREE"), ADDSET_ROTATE_DEGREE, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_DISTORTION_LABEL"));
    appendBehavList (mi, M("TP_DISTORTION_AMOUNT"), ADDSET_DIST_AMOUNT, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_PERSPECTIVE_LABEL"));
    appendBehavList (mi, M("TP_PERSPECTIVE_HORIZONTAL") + ", " + M("TP_PERSPECTIVE_VERTICAL"), ADDSET_PERSPECTIVE, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_GRADIENT_LABEL"));
    appendBehavList (mi, M("TP_GRADIENT_DEGREE"), ADDSET_GRADIENT_DEGREE, false);
    appendBehavList (mi, M("TP_GRADIENT_FEATHER"), ADDSET_GRADIENT_FEATHER, false);
    appendBehavList (mi, M("TP_GRADIENT_STRENGTH"), ADDSET_GRADIENT_STRENGTH, false);
    appendBehavList (mi, M("TP_GRADIENT_CENTER_X") + ", " + M("TP_GRADIENT_CENTER_Y"), ADDSET_GRADIENT_CENTER, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_PCVIGNETTE_LABEL"));
    appendBehavList (mi, M("TP_PCVIGNETTE_STRENGTH"), ADDSET_PCVIGNETTE_STRENGTH, false);
    appendBehavList (mi, M("TP_PCVIGNETTE_FEATHER"), ADDSET_PCVIGNETTE_FEATHER, false);
    appendBehavList (mi, M("TP_PCVIGNETTE_ROUNDNESS"), ADDSET_PCVIGNETTE_ROUNDNESS, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_CACORRECTION_LABEL"));
    appendBehavList (mi, M("TP_CACORRECTION_BLUE") + ", " + M("TP_CACORRECTION_RED"), ADDSET_CA, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_VIGNETTING_LABEL"));
    appendBehavList (mi, M("TP_VIGNETTING_AMOUNT"), ADDSET_VIGN_AMOUNT, false);
    appendBehavList (mi, M("TP_VIGNETTING_RADIUS"), ADDSET_VIGN_RADIUS, false);
    appendBehavList (mi, M("TP_VIGNETTING_STRENGTH"), ADDSET_VIGN_STRENGTH, false);
    appendBehavList (mi, M("TP_VIGNETTING_CENTER_X") + ", " + M("TP_VIGNETTING_CENTER_Y"), ADDSET_VIGN_CENTER, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_DIRPYREQUALIZER_LABEL"));
    appendBehavList (mi, M("TP_EXPOSURE_CONTRAST"), ADDSET_DIRPYREQ, true);
    appendBehavList (mi, M("TP_DIRPYREQUALIZER_THRESHOLD"), ADDSET_DIRPYREQ_THRESHOLD, true);
    appendBehavList (mi, M("TP_DIRPYREQUALIZER_SKIN"), ADDSET_DIRPYREQ_SKINPROTECT, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_WAVELET_LABEL"));
    appendBehavList (mi, M("TP_WAVELET_LEVELS"), ADDSET_WA_THRES, true);
    //  appendBehavList (mi, M("TP_WAVELET_CONTRAST"), ADDSET_WA, true);
    appendBehavList (mi, M("TP_WAVELET_THRESHOLD"), ADDSET_WA_THRESHOLD, true);
    appendBehavList (mi, M("TP_WAVELET_THRESHOLD2"), ADDSET_WA_THRESHOLD2, true);
    appendBehavList (mi, M("TP_WAVELET_CHRO"), ADDSET_WA_CHRO, true);
    appendBehavList (mi, M("TP_WAVELET_CHR"), ADDSET_WA_CHROMA, true);
    appendBehavList (mi, M("TP_WAVELET_SKIN"), ADDSET_WA_SKINPROTECT, true);
    appendBehavList (mi, M("TP_WAVELET_EDRAD"), ADDSET_WA_EDGRAD, true);
    appendBehavList (mi, M("TP_WAVELET_EDVAL"), ADDSET_WA_EDGVAL, true);
    appendBehavList (mi, M("TP_WAVELET_RESCON"), ADDSET_WA_RESCON, true);
    appendBehavList (mi, M("TP_WAVELET_THR"), ADDSET_WA_THRR, true);
    appendBehavList (mi, M("TP_WAVELET_RESCONH"), ADDSET_WA_RESCONH, true);
    appendBehavList (mi, M("TP_WAVELET_THRH"), ADDSET_WA_THRRH, true);
    appendBehavList (mi, M("TP_WAVELET_RESCHRO"), ADDSET_WA_RESCHRO, true);
    appendBehavList (mi, M("TP_WAVELET_TMSTRENGTH"), ADDSET_WA_TMRS, true);
    appendBehavList (mi, M("TP_WAVELET_SKY"), ADDSET_WA_SKYPROTECT, true);
    appendBehavList (mi, M("TP_WAVELET_CONTRA"), ADDSET_WA_CONTRAST, true);
    appendBehavList (mi, M("TP_WAVELET_STRENGTH"), ADDSET_WA_STRENGTH, true);
    appendBehavList (mi, M("TP_WAVELET_COMPGAMMA"), ADDSET_WA_GAMMA, true);
    appendBehavList (mi, M("TP_WAVELET_EDGEDETECT"), ADDSET_WA_EDGEDETECT, true);
    appendBehavList (mi, M("TP_WAVELET_EDGEDETECTTHR"), ADDSET_WA_EDGEDETECTTHR, true);
    appendBehavList (mi, M("TP_WAVELET_EDGEDETECTTHR2"), ADDSET_WA_EDGEDETECTTHR2, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_PREPROCESS_LABEL"));
    appendBehavList (mi, M("TP_PREPROCESS_GREENEQUIL"), ADDSET_PREPROCESS_GREENEQUIL, false);
    appendBehavList (mi, M("TP_PREPROCESS_LINEDENOISE"), ADDSET_PREPROCESS_LINEDENOISE, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_EXPOS_WHITEPOINT_LABEL"));
    appendBehavList (mi, M("TP_RAWEXPOS_LINEAR"), ADDSET_RAWEXPOS_LINEAR, false);
    appendBehavList (mi, M("TP_RAWEXPOS_PRESER"), ADDSET_RAWEXPOS_PRESER, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_RAWEXPOS_BLACKS"));
    appendBehavList (mi, M("TP_RAWEXPOS_RGB"), ADDSET_RAWEXPOS_BLACKS, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_FLATFIELD_LABEL"));
    appendBehavList (mi, M("TP_FLATFIELD_CLIPCONTROL"), ADDSET_RAWFFCLIPCONTROL, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("TP_CHROMATABERR_LABEL"));
    appendBehavList (mi, M("TP_RAWCACORR_CARED") + ", " + M("TP_RAWCACORR_CABLUE"), ADDSET_RAWCACORR, true);

    behTreeView->expand_all ();

    behAddAll = Gtk::manage( new Gtk::Button (M("PREFERENCES_BEHADDALL")) );
    behSetAll = Gtk::manage( new Gtk::Button (M("PREFERENCES_BEHSETALL")) );
    behAddAll->set_tooltip_markup (M("PREFERENCES_BEHADDALLHINT"));
    behSetAll->set_tooltip_markup (M("PREFERENCES_BEHSETALLHINT"));

    behAddAll->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::behAddAllPressed) );
    behSetAll->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::behSetAllPressed) );

    Gtk::HBox* buttonpanel1 = Gtk::manage (new Gtk::HBox ());
    //buttonpanel1->set_spacing(8);
    buttonpanel1->pack_end (*behSetAll, Gtk::PACK_SHRINK, 4);
    buttonpanel1->pack_end (*behAddAll, Gtk::PACK_SHRINK, 4);
    vbbeh->pack_start (*buttonpanel1, Gtk::PACK_SHRINK, 4);

    chOverwriteOutputFile =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_OVERWRITEOUTPUTFILE")) );
    mvbpp->pack_start(*chOverwriteOutputFile, Gtk::PACK_SHRINK, 4);

    return mvbpp;
}

void Preferences::appendBehavList (Gtk::TreeModel::iterator& parent, Glib::ustring label, int id, bool set)
{

    Gtk::TreeModel::iterator ci = behModel->append (parent->children());
    ci->set_value (behavColumns.label, label);
    ci->set_value (behavColumns.visible, true);
    ci->set_value (behavColumns.badd, !set);
    ci->set_value (behavColumns.bset, set);
    ci->set_value (behavColumns.addsetid, id);
}

void Preferences::behAddRadioToggled (const Glib::ustring& path)
{

    Gtk::TreeModel::iterator iter = behModel->get_iter (path);
    //bool set = iter->get_value (behavColumns.bset);
    iter->set_value (behavColumns.bset, false);
    iter->set_value (behavColumns.badd, true);
}

void Preferences::behSetRadioToggled (const Glib::ustring& path)
{

    Gtk::TreeModel::iterator iter = behModel->get_iter (path);
    //bool add = iter->get_value (behavColumns.badd);
    iter->set_value (behavColumns.bset, true);
    iter->set_value (behavColumns.badd, false);
}

Gtk::Widget* Preferences::getProcParamsPanel ()
{

    Gtk::VBox* mvbpp = Gtk::manage (new Gtk::VBox ());

    Gtk::Frame* fpp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_IMPROCPARAMS")));
    Gtk::VBox* vbpp = Gtk::manage (new Gtk::VBox ());
    Gtk::Label* drlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_FORRAW") + ":", Gtk::ALIGN_START));
    rprofiles = Gtk::manage (new ProfileStoreComboBox ());
    setExpandAlignProperties(rprofiles, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    rprofiles->set_size_request(50, -1);
    rpconn = rprofiles->signal_changed().connect( sigc::mem_fun(*this, &Preferences::forRAWComboChanged) );
    Gtk::Label* drimg = Gtk::manage (new Gtk::Label (M("PREFERENCES_FORIMAGE") + ":", Gtk::ALIGN_START));
    iprofiles = Gtk::manage (new ProfileStoreComboBox ());
    iprofiles->set_size_request(50, -1);
    setExpandAlignProperties(iprofiles, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    ipconn = iprofiles->signal_changed().connect( sigc::mem_fun(*this, &Preferences::forImageComboChanged) );
    Gtk::Table* defpt = Gtk::manage (new Gtk::Table (2, 2));
    defpt->attach (*drlab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach (*rprofiles, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    defpt->attach (*drimg, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach (*iprofiles, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    vbpp->pack_start (*defpt, Gtk::PACK_SHRINK, 4);
    useBundledProfiles = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_USEBUNDLEDPROFILES")));
    bpconn = useBundledProfiles->signal_clicked().connect ( sigc::mem_fun(*this, &Preferences::bundledProfilesChanged) );
    vbpp->pack_start (*useBundledProfiles, Gtk::PACK_SHRINK, 4);
    fpp->add (*vbpp);
    mvbpp->pack_start (*fpp, Gtk::PACK_SHRINK, 4);

    // Custom profile builder box
    Gtk::Frame* cpfrm = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CUSTPROFBUILD")) );
    Gtk::Label* cplab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CUSTPROFBUILDPATH") + ":", Gtk::ALIGN_START) );
    txtCustProfBuilderPath = Gtk::manage( new Gtk::Entry () );
    txtCustProfBuilderPath->set_tooltip_markup (M("PREFERENCES_CUSTPROFBUILDHINT"));
    Gtk::Label* cpltypelab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CUSTPROFBUILDKEYFORMAT") + ":", Gtk::ALIGN_START) );
    custProfBuilderLabelType = Gtk::manage (new Gtk::ComboBoxText ());
    custProfBuilderLabelType->append (M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_TID"));
    custProfBuilderLabelType->append (M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_NAME"));
    custProfBuilderLabelType->append (M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_TID") + "_" + M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_NAME"));
    Gtk::Table* cpbt = Gtk::manage (new Gtk::Table (2, 2));
    cpbt->attach (*cplab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    cpbt->attach (*txtCustProfBuilderPath, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    cpbt->attach (*cpltypelab, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    cpbt->attach (*custProfBuilderLabelType, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    cpfrm->add (*cpbt);
    mvbpp->pack_start (*cpfrm, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdp = Gtk::manage (new Gtk::Frame (M("PREFERENCES_PROFILEHANDLING")));
    Gtk::VBox* vbdp = Gtk::manage (new Gtk::VBox ());
    saveParamsFile = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_PROFILESAVEINPUT")));
    vbdp->pack_start (*saveParamsFile, Gtk::PACK_SHRINK, 4);
    saveParamsCache = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_PROFILESAVECACHE")));
    vbdp->pack_start (*saveParamsCache, Gtk::PACK_SHRINK, 4);
    Gtk::Label* lplab = Gtk::manage (new Gtk::Label (M("PREFERENCES_PROFILELOADPR") + ":"));
    loadParamsPreference = Gtk::manage (new Gtk::ComboBoxText ());
    loadParamsPreference->append (M("PREFERENCES_PROFILEPRCACHE"));
    loadParamsPreference->append (M("PREFERENCES_PROFILEPRFILE"));
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
    Gtk::Label *dfLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_DIRDARKFRAMES") + ":"));
    hb42->pack_start(*dfLab , Gtk::PACK_SHRINK, 4 );
    hb42->pack_start(*darkFrameDir, Gtk::PACK_EXPAND_WIDGET, 4);
    dfLabel = Gtk::manage(new Gtk::Label("Found:"));
    Gtk::VBox* vbdf = Gtk::manage (new Gtk::VBox ());
    vbdf->pack_start( *hb42, Gtk::PACK_SHRINK, 4);
    vbdf->pack_start( *dfLabel, Gtk::PACK_SHRINK, 4 );
    fdf->add( *vbdf );
    mvbpp->pack_start ( *fdf , Gtk::PACK_SHRINK, 4);

    //dfconn = darkFrameDir->signal_file_set().connect ( sigc::mem_fun(*this, &Preferences::darkFrameChanged), true);
    dfconn = darkFrameDir->signal_current_folder_changed().connect ( sigc::mem_fun(*this, &Preferences::darkFrameChanged), true);

    // FLATFIELD
    Gtk::Frame* fff = Gtk::manage (new Gtk::Frame (M("PREFERENCES_FLATFIELD")) );
    Gtk::HBox* hb43 = Gtk::manage (new Gtk::HBox ());
    flatFieldDir = Gtk::manage(new Gtk::FileChooserButton(M("PREFERENCES_FLATFIELDSDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label *ffLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_FLATFIELDSDIR") + ":"));
    hb43->pack_start(*ffLab , Gtk::PACK_SHRINK, 4 );
    hb43->pack_start(*flatFieldDir);
    ffLabel = Gtk::manage(new Gtk::Label("Found:"));
    Gtk::VBox* vbff = Gtk::manage (new Gtk::VBox ());
    vbff->pack_start( *hb43, Gtk::PACK_SHRINK, 4);
    vbff->pack_start( *ffLabel, Gtk::PACK_SHRINK, 4 );
    fff->add( *vbff );
    mvbpp->pack_start ( *fff , Gtk::PACK_SHRINK, 4);

    //ffconn = flatFieldDir->signal_file_set().connect ( sigc::mem_fun(*this, &Preferences::flatFieldChanged), true);
    ffconn = flatFieldDir->signal_current_folder_changed().connect ( sigc::mem_fun(*this, &Preferences::flatFieldChanged), true);

    //Cluts Dir
    Gtk::Frame* clutsDirFrame = Gtk::manage (new Gtk::Frame (M("PREFERENCES_FILMSIMULATION")) );
    Gtk::HBox* clutsDirBox = Gtk::manage (new Gtk::HBox ());
    clutsDir = Gtk::manage(new Gtk::FileChooserButton(M("PREFERENCES_CLUTSDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    Gtk::Label *clutsDirLabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_CLUTSDIR") + ":"));
    Gtk::Label* clutsRestartNeeded = Gtk::manage( new Gtk::Label (Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    clutsDirBox->pack_start( *clutsDirLabel, Gtk::PACK_SHRINK, 4 );
    clutsDirBox->pack_start( *clutsDir );
    clutsDirBox->pack_start( *clutsRestartNeeded, Gtk::PACK_SHRINK, 4 );
    clutsDirFrame->add( *clutsDirBox );
    mvbpp->pack_start( *clutsDirFrame, Gtk::PACK_SHRINK, 4 );

    Gtk::Frame* fmd = Gtk::manage (new Gtk::Frame (M("PREFERENCES_METADATA")));
    Gtk::VBox* vbmd = Gtk::manage (new Gtk::VBox ());
    ckbTunnelMetaData = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_TUNNELMETADATA")));
    vbmd->pack_start (*ckbTunnelMetaData, Gtk::PACK_SHRINK, 4);
    fmd->add (*vbmd);
    mvbpp->pack_start (*fmd, Gtk::PACK_SHRINK, 4);

    return mvbpp;
}

Gtk::Widget* Preferences::getPerformancePanel ()
{
    Gtk::VBox* mainContainer = Gtk::manage( new Gtk::VBox () );
    mainContainer->set_spacing(4);

    Gtk::Frame* fprevdemo = Gtk::manage (new Gtk::Frame (M("PREFERENCES_PREVDEMO")));
    Gtk::HBox* hbprevdemo = Gtk::manage (new Gtk::HBox (false, 4));
    Gtk::Label* lprevdemo = Gtk::manage (new Gtk::Label (M("PREFERENCES_PREVDEMO_LABEL")));
    cprevdemo = Gtk::manage (new Gtk::ComboBoxText ());
    cprevdemo->append (M("PREFERENCES_PREVDEMO_FAST"));
    cprevdemo->append (M("PREFERENCES_PREVDEMO_SIDECAR"));
    cprevdemo->set_active (1);
    hbprevdemo->pack_start (*lprevdemo, Gtk::PACK_SHRINK);
    hbprevdemo->pack_start (*cprevdemo);
    fprevdemo->add (*hbprevdemo);
    mainContainer->pack_start (*fprevdemo, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* ftiffserialize = Gtk::manage (new Gtk::Frame (M("PREFERENCES_SERIALIZE_TIFF_READ")));
    Gtk::HBox* htiffserialize = Gtk::manage (new Gtk::HBox (false, 4));
    ctiffserialize = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SERIALIZE_TIFF_READ_LABEL")) );
    ctiffserialize->set_tooltip_text(M("PREFERENCES_SERIALIZE_TIFF_READ_TOOLTIP"));
    htiffserialize->pack_start (*ctiffserialize);
    ftiffserialize->add (*htiffserialize);
    mainContainer->pack_start (*ftiffserialize, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fclut = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CLUTSCACHE")) );
    Gtk::HBox* clutCacheSizeHB = Gtk::manage( new Gtk::HBox () );
    clutCacheSizeHB->set_spacing(4);
    Gtk::Label* CLUTLl = Gtk::manage( new Gtk::Label (M("PREFERENCES_CLUTSCACHE_LABEL") + ":", Gtk::ALIGN_START));
    clutCacheSizeSB = Gtk::manage( new Gtk::SpinButton () );
    clutCacheSizeSB->set_digits (0);
    clutCacheSizeSB->set_increments (1, 5);
    clutCacheSizeSB->set_max_length(2);  // Will this be sufficient? :)
#ifdef _OPENMP
    clutCacheSizeSB->set_range (1, 3 * omp_get_num_procs());
#else
    clutCacheSizeSB->set_range (1, 12);
#endif
    clutCacheSizeHB->pack_start (*CLUTLl, Gtk::PACK_SHRINK, 0);
    clutCacheSizeHB->pack_end (*clutCacheSizeSB, Gtk::PACK_SHRINK, 0);
    fclut->add (*clutCacheSizeHB);
    mainContainer->pack_start (*fclut, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* finspect = Gtk::manage(  new Gtk::Frame (M("PREFERENCES_INSPECT_LABEL")) );
    Gtk::HBox* maxIBuffersHB = Gtk::manage( new Gtk::HBox () );
    maxIBuffersHB->set_spacing(4);
    maxIBuffersHB->set_tooltip_text(M("PREFERENCES_INSPECT_MAXBUFFERS_TOOLTIP"));
    Gtk::Label* maxIBufferLbl = Gtk::manage( new Gtk::Label (M("PREFERENCES_INSPECT_MAXBUFFERS_LABEL") + ":", Gtk::ALIGN_START));
    maxInspectorBuffersSB = Gtk::manage( new Gtk::SpinButton () );
    maxInspectorBuffersSB->set_digits (0);
    maxInspectorBuffersSB->set_increments (1, 5);
    maxInspectorBuffersSB->set_max_length(2);
    maxInspectorBuffersSB->set_range (1, 12);  // ... we have to set a limit, 12 seem to be enough even for systems with tons of RAM
    maxIBuffersHB->pack_start (*maxIBufferLbl, Gtk::PACK_SHRINK, 0);
    maxIBuffersHB->pack_end (*maxInspectorBuffersSB, Gtk::PACK_SHRINK, 0);
    finspect->add(*maxIBuffersHB);
    mainContainer->pack_start(*finspect, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdenoise = Gtk::manage( new Gtk::Frame (M("PREFERENCES_NOISE")) );
    Gtk::VBox* vbdenoise = Gtk::manage( new Gtk::VBox (Gtk::PACK_SHRINK, 4) );

    Gtk::Label* lreloadneeded2 = Gtk::manage (new Gtk::Label (M("PREFERENCES_IMG_RELOAD_NEEDED"), Gtk::ALIGN_START));
    Gtk::HBox* threadLimitHB = Gtk::manage (new Gtk::HBox (Gtk::PACK_SHRINK, 4));
    threadLimitHB->set_tooltip_text(M("PREFERENCES_RGBDTL_TOOLTIP"));
    Gtk::Label* RGBDTLl = Gtk::manage( new Gtk::Label (M("PREFERENCES_RGBDTL_LABEL") + ":", Gtk::ALIGN_START));
    rgbDenoiseTreadLimitSB = Gtk::manage( new Gtk::SpinButton () );
    rgbDenoiseTreadLimitSB->set_digits (0);
    rgbDenoiseTreadLimitSB->set_increments (1, 5);
    rgbDenoiseTreadLimitSB->set_max_length(2);  // Will this be sufficient? :)
#ifdef _OPENMP
    int maxThreadNumber = omp_get_max_threads();
#else
    int maxThreadNumber = 10;
#endif
    rgbDenoiseTreadLimitSB->set_range (0, maxThreadNumber);
    threadLimitHB->pack_start (*RGBDTLl, Gtk::PACK_SHRINK, 2);
    threadLimitHB->pack_end (*rgbDenoiseTreadLimitSB, Gtk::PACK_SHRINK, 2);

    Gtk::Label* dnlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_LEVDN") + ":", Gtk::ALIGN_START));
    Gtk::Label* dnautlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_LEVAUTDN") + ":", Gtk::ALIGN_START));
    Gtk::Label* dnautsimpllab = Gtk::manage (new Gtk::Label (M("PREFERENCES_SIMPLAUT") + ":", Gtk::ALIGN_START));
    Gtk::Label* dntilab = Gtk::manage (new Gtk::Label (M("PREFERENCES_TINB") + ":", Gtk::ALIGN_START));
    Gtk::Label* dnwavlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_WAVLEV") + ":", Gtk::ALIGN_START));
    Gtk::Label* dnlisslab = Gtk::manage (new Gtk::Label (M("PREFERENCES_LISS") + ":", Gtk::ALIGN_START));

    dnv = Gtk::manage (new Gtk::ComboBoxText ());
    dnv->append (M("PREFERENCES_MIN"));
    dnv->append (M("PREFERENCES_SMA"));
    dnv->append (M("PREFERENCES_MED"));
    dnv->append (M("PREFERENCES_MAX"));
    dnaut = Gtk::manage (new Gtk::ComboBoxText ());
    dnaut->append (M("PREFERENCES_AUTLOW"));
    dnaut->append (M("PREFERENCES_AUTSTD"));

    dnautsimpl = Gtk::manage (new Gtk::ComboBoxText ());
    dnautsimpl->append (M("PREFERENCES_STDAUT"));
    dnautsimpl->append (M("PREFERENCES_EXPAUT"));

    dnliss = Gtk::manage (new Gtk::ComboBoxText ());
    dnliss->append (M("PREFERENCES_AUTLISVLOW"));//very low
    dnliss->append (M("PREFERENCES_AUTLISLOW"));//low
    dnliss->append (M("PREFERENCES_AUTLISSTD"));//med
    dnliss->append (M("PREFERENCES_AUTLISMAX"));//max

    dnti = Gtk::manage (new Gtk::ComboBoxText ());
    dnti->append (M("PREFERENCES_TISTD"));
    dnti->append (M("PREFERENCES_TIMAX"));

    dnwavlev = Gtk::manage (new Gtk::ComboBoxText ());
    dnwavlev->append (M("PREFERENCES_WLZER"));
    dnwavlev->append (M("PREFERENCES_WLONE"));
    dnwavlev->append (M("PREFERENCES_WLTWO"));

    Gtk::Table* colon = Gtk::manage (new Gtk::Table (6, 2));
    colon->attach (*dnlab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colon->attach (*dnv, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colon->attach (*dnautlab, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colon->attach (*dnaut, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colon->attach (*dnautsimpllab, 0, 1, 2, 3, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colon->attach (*dnautsimpl, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colon->attach (*dnlisslab, 0, 1, 3, 4, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colon->attach (*dnliss, 1, 2, 3, 4, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colon->attach (*dntilab, 0, 1, 4, 5, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colon->attach (*dnti, 1, 2, 4, 5, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    colon->attach (*dnwavlab, 0, 1, 5, 6, Gtk::FILL, Gtk::SHRINK, 2, 2);
    colon->attach (*dnwavlev, 1, 2, 5, 6, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);

    vbdenoise->pack_start (*lreloadneeded2, Gtk::PACK_SHRINK);
    vbdenoise->pack_start (*colon, Gtk::PACK_SHRINK);
    vbdenoise->pack_start(*threadLimitHB, Gtk::PACK_SHRINK);
    // <--- To be hard-coded and removed once tested
    cbdaubech = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_DAUB_LABEL"), Gtk::ALIGN_START));
    cbdaubech->set_tooltip_markup (M("PREFERENCES_DAUB_TOOLTIP"));
//   vbdenoise->pack_start (*cbdaubech, Gtk::PACK_SHRINK);
    // --->
    fdenoise->add (*vbdenoise);
    mainContainer->pack_start (*fdenoise, Gtk::PACK_SHRINK, 4);

    return mainContainer;
}

Gtk::Widget* Preferences::getColorManagementPanel ()
{

    Gtk::VBox* mvbcm = Gtk::manage (new Gtk::VBox ());
    mvbcm->set_spacing (4);

    iccDir = Gtk::manage (new Gtk::FileChooserButton (M("PREFERENCES_ICCDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(iccDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* pdlabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_ICCDIR") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(pdlabel, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Grid* iccdgrid = Gtk::manage (new Gtk::Grid ());
    setExpandAlignProperties (iccdgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    iccdgrid->set_column_spacing (4);

    iccdgrid->attach (*pdlabel, 0, 0, 1, 1);
    iccdgrid->attach (*iccDir, 1, 0, 1, 1);

    iccDir->signal_selection_changed ().connect (sigc::mem_fun (this, &Preferences::iccDirChanged));

    mvbcm->pack_start(*iccdgrid, Gtk::PACK_SHRINK);

    //-------------------------  MONITOR ----------------------

    Gtk::Frame* fmonitor = Gtk::manage( new Gtk::Frame (M("PREFERENCES_MONITOR")) );
    Gtk::Grid* gmonitor = Gtk::manage( new Gtk::Grid () );
    gmonitor->set_column_spacing (4);

    monProfile = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(monProfile, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* mplabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_MONPROFILE") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(mplabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    monIntent = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(monIntent, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* milabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_MONINTENT")+":", Gtk::ALIGN_START));
    setExpandAlignProperties(milabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    monProfile->append (M("PREFERENCES_PROFILE_NONE"));
    monProfile->set_active (0);

    const std::vector<Glib::ustring> profiles = rtengine::ICCStore::getInstance ()->getProfiles (rtengine::ICCStore::ProfileType::MONITOR);
    for (const auto profile : profiles) {
        monProfile->append (profile);
    }

    // same order as the enum
    monIntent->append (M("PREFERENCES_INTENT_PERCEPTUAL"));
    monIntent->append (M("PREFERENCES_INTENT_RELATIVE"));
    monIntent->append (M("PREFERENCES_INTENT_ABSOLUTE"));
    monIntent->set_active (1);
    monIntent->set_size_request(120, -1);

    monBPC = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_CMMBPC")));
    setExpandAlignProperties(monBPC, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    monBPC->set_active (true);

#if defined(WIN32) // Auto-detection not implemented for Linux, see issue 851
    cbAutoMonProfile = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_AUTOMONPROFILE")));
    setExpandAlignProperties(cbAutoMonProfile, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    autoMonProfileConn  = cbAutoMonProfile->signal_toggled().connect (sigc::mem_fun(*this, &Preferences::autoMonProfileToggled));
#endif

    int row = 0;
    gmonitor->attach (*mplabel, 0, row, 1, 1);
#if defined(__APPLE__) // monitor profile not supported on apple
    Gtk::Label *osxwarn = Gtk::manage (new Gtk::Label (M("PREFERENCES_MONPROFILE_WARNOSX"), Gtk::ALIGN_LEFT));
    setExpandAlignProperties(osxwarn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    gmonitor->attach (*osxwarn, 1, row, 1, 1);
#else
    gmonitor->attach (*monProfile, 1, row, 1, 1);
#endif
    ++row;
#if defined(WIN32)
    gmonitor->attach (*cbAutoMonProfile, 1, row, 1, 1);
    ++row;
#endif
    gmonitor->attach (*milabel, 0, row, 1, 1);
    gmonitor->attach (*monIntent, 1, row, 1, 1);
    ++row;
    gmonitor->attach (*monBPC, 0, row, 2, 1);

#if defined(WIN32)
    autoMonProfileToggled();
#endif

    fmonitor->add(*gmonitor);

    mvbcm->pack_start(*fmonitor, Gtk::PACK_SHRINK);

    //-------------------------  PRINTER ----------------------

    Gtk::Frame* fprinter = Gtk::manage( new Gtk::Frame (M("PREFERENCES_PRINTER")) );
    Gtk::Grid* gprinter = Gtk::manage( new Gtk::Grid () );
    gprinter->set_column_spacing (4);
    prtProfile = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(prtProfile, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* pplabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_PRTPROFILE") + ":"));
    setExpandAlignProperties(pplabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    prtIntent = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(prtIntent, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* pilabel = Gtk::manage (new Gtk::Label (M("PREFERENCES_PRTINTENT")+":"));
    setExpandAlignProperties(pilabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    prtProfile->append (M("PREFERENCES_PROFILE_NONE"));
    prtProfile->set_active (0);

    const std::vector<Glib::ustring> prtprofiles = rtengine::ICCStore::getInstance ()->getProfiles (rtengine::ICCStore::ProfileType::PRINTER);
    for (const auto prtprofile : prtprofiles)
        prtProfile->append (prtprofile);

    // same order as the enum
    prtIntent->append (M("PREFERENCES_INTENT_PERCEPTUAL"));
    prtIntent->append (M("PREFERENCES_INTENT_RELATIVE"));
    prtIntent->append (M("PREFERENCES_INTENT_ABSOLUTE"));
    prtIntent->set_active (1);

    prtBPC = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_CMMBPC")));
    setExpandAlignProperties(prtBPC, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    prtBPC->set_active (true);

    row = 0;
    gprinter->attach (*pplabel, 0, row, 1, 1);
    gprinter->attach (*prtProfile, 1, row, 1, 1);
    ++row;
    gprinter->attach (*pilabel, 0, row, 1, 1);
    gprinter->attach (*prtIntent, 1, row, 1, 1);
    ++row;
    gprinter->attach (*prtBPC, 0, row, 2, 1);

#if defined(WIN32)
    autoMonProfileToggled();
#endif

    fprinter->add(*gprinter);

    mvbcm->pack_start(*fprinter, Gtk::PACK_SHRINK);

    //-------------------------  CIECAM ----------------------

    Gtk::Label* viewlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_VIEW") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(viewlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    view = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(view, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    view->append (M("PREFERENCES_D50"));
    view->append (M("PREFERENCES_D55"));
    view->append (M("PREFERENCES_D60"));
    view->append (M("PREFERENCES_D65"));
    view->append (M("PREFERENCES_BLACKBODY"));
    view->append (M("PREFERENCES_FLUOF2"));
    view->append (M("PREFERENCES_FLUOF7"));
    view->append (M("PREFERENCES_FLUOF11"));

    Gtk::Label* greylab = Gtk::manage (new Gtk::Label (M("PREFERENCES_GREY") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(greylab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    grey = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(grey, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    grey->append (M("PREFERENCES_GREY05"));
    grey->append (M("PREFERENCES_GREY10"));
    grey->append (M("PREFERENCES_GREY15"));
    grey->append (M("PREFERENCES_GREY18"));
    grey->append (M("PREFERENCES_GREY23"));
    grey->append (M("PREFERENCES_GREY30"));
    grey->append (M("PREFERENCES_GREY40"));

    Gtk::Label* greySclab = Gtk::manage (new Gtk::Label (M("PREFERENCES_GREYSC") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(greySclab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    greySc = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(greySc, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    greySc->append (M("PREFERENCES_GREYSCA"));
    greySc->append (M("PREFERENCES_GREYSC18"));

    Gtk::Frame* fcielab = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CIEART_FRAME")) );
    setExpandAlignProperties(fcielab, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    Gtk::Grid* colo = Gtk::manage (new Gtk::Grid ());
    setExpandAlignProperties(colo, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Label* lreloadneeded1 = Gtk::manage (new Gtk::Label (M("PREFERENCES_IMG_RELOAD_NEEDED"), Gtk::ALIGN_START));
    setExpandAlignProperties(lreloadneeded1, true, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    colo->attach (*lreloadneeded1, 0, 0, 2, 1);
    colo->attach (*viewlab, 0, 1, 1, 1);
    colo->attach (*view, 1, 1, 1, 1);
    colo->attach (*greylab, 0, 2, 1, 1);
    colo->attach (*grey, 1, 2, 1, 1);
    colo->attach (*greySclab, 0, 3, 1, 1);
    colo->attach (*greySc, 1, 3, 1, 1);
    cbciecamfloat = Gtk::manage (new Gtk::CheckButton (M("PREFERENCES_CIEART_LABEL")));
    setExpandAlignProperties(cbciecamfloat, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    colo->attach (*cbciecamfloat, 0, 4, 2, 1);
    cbciecamfloat->set_tooltip_markup (M("PREFERENCES_CIEART_TOOLTIP"));
    fcielab->add (*colo);

    mvbcm->pack_start (*fcielab, Gtk::PACK_SHRINK, 4);

    return mvbcm;
}

Gtk::Widget* Preferences::getGeneralPanel ()
{

    Gtk::Grid* mvbsd = Gtk::manage( new Gtk::Grid () );
    mvbsd->set_column_spacing(4);
    mvbsd->set_row_spacing(4);

    Gtk::Frame* fworklflow = Gtk::manage (new Gtk::Frame (M("PREFERENCES_WORKFLOW")));
    setExpandAlignProperties(fworklflow, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    Gtk::Grid* workflowGrid = Gtk::manage (new Gtk::Grid());
    workflowGrid->set_column_spacing(4);
    workflowGrid->set_row_spacing(4);
    setExpandAlignProperties(workflowGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    Gtk::Label* flayoutlab = Gtk::manage (new Gtk::Label (M("PREFERENCES_EDITORLAYOUT") + ":"));
    setExpandAlignProperties(flayoutlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    editorLayout = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(editorLayout, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    editorLayout->append (M("PREFERENCES_SINGLETAB"));
    editorLayout->append (M("PREFERENCES_SINGLETABVERTAB"));
    editorLayout->append (M("PREFERENCES_MULTITAB"));
    editorLayout->append (M("PREFERENCES_MULTITABDUALMON"));
    editorLayout->set_active (2);
    Gtk::CellRendererText* cellRenderer = dynamic_cast<Gtk::CellRendererText*>(editorLayout->get_first_cell());
    cellRenderer->property_ellipsize() = Pango::ELLIPSIZE_MIDDLE;
    cellRenderer->property_ellipsize_set() = true;
    editorLayout->signal_changed().connect (sigc::mem_fun(*this, &Preferences::layoutComboChanged));
    layoutComboChanged(); // update the tooltip
    Gtk::Label* lNextStart = Gtk::manage( new Gtk::Label (Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    setExpandAlignProperties(lNextStart, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*flayoutlab, Gtk::POS_LEFT, 1, 1);
    workflowGrid->attach_next_to(*editorLayout, *flayoutlab, Gtk::POS_RIGHT, 1, 1);
    workflowGrid->attach_next_to(*lNextStart, *editorLayout, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* curveBBoxPosL = Gtk::manage (new Gtk::Label (M("PREFERENCES_CURVEBBOXPOS") + ":"));
    setExpandAlignProperties(curveBBoxPosL, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    curveBBoxPosC = Gtk::manage (new Gtk::ComboBoxText ());
    setExpandAlignProperties(curveBBoxPosC, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    curveBBoxPosC->append (M("PREFERENCES_CURVEBBOXPOS_ABOVE"));
    curveBBoxPosC->append (M("PREFERENCES_CURVEBBOXPOS_RIGHT"));
    curveBBoxPosC->append (M("PREFERENCES_CURVEBBOXPOS_BELOW"));
    curveBBoxPosC->append (M("PREFERENCES_CURVEBBOXPOS_LEFT"));
    curveBBoxPosC->set_active (1);
    Gtk::Label* curveBBoxPosRestartL = Gtk::manage (new Gtk::Label (Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")"));
    setExpandAlignProperties(curveBBoxPosRestartL, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*curveBBoxPosL, *flayoutlab, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*curveBBoxPosC, *editorLayout, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*curveBBoxPosRestartL, *lNextStart, Gtk::POS_BOTTOM, 1, 1);

    ckbHistogramPositionLeft =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_HISTOGRAMPOSITIONLEFT")) );
    setExpandAlignProperties(ckbHistogramPositionLeft, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    ckbHistogramWorking =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_HISTOGRAMWORKING")) );
    setExpandAlignProperties(ckbHistogramWorking, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    ckbHistogramWorking->set_tooltip_markup (M("PREFERENCES_HISTOGRAM_TOOLTIP"));
    workflowGrid->attach_next_to(*ckbHistogramPositionLeft, *curveBBoxPosL, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*ckbHistogramWorking, *curveBBoxPosC, Gtk::POS_BOTTOM, 2, 1);

    ckbFileBrowserToolbarSingleRow =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_FILEBROWSERTOOLBARSINGLEROW")) );
    setExpandAlignProperties(ckbFileBrowserToolbarSingleRow, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    ckbShowFilmStripToolBar =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_SHOWFILMSTRIPTOOLBAR")) );
    setExpandAlignProperties(ckbShowFilmStripToolBar, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    workflowGrid->attach_next_to(*ckbFileBrowserToolbarSingleRow, *ckbHistogramPositionLeft, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*ckbShowFilmStripToolBar, *ckbHistogramWorking, Gtk::POS_BOTTOM, 2, 1);

    Gtk::Label* hb4label =  Gtk::manage( new Gtk::Label (M("PREFERENCES_TP_LABEL")) );
    setExpandAlignProperties(hb4label, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    ckbHideTPVScrollbar =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_TP_VSCROLLBAR")) );
    setExpandAlignProperties(ckbHideTPVScrollbar, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    ckbUseIconNoText =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_TP_USEICONORTEXT")) );
    setExpandAlignProperties(ckbUseIconNoText, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*hb4label, *ckbFileBrowserToolbarSingleRow, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*ckbHideTPVScrollbar, *hb4label, Gtk::POS_RIGHT, 1, 1);
    workflowGrid->attach_next_to(*ckbUseIconNoText, *ckbHideTPVScrollbar, Gtk::POS_RIGHT, 1, 1);

    fworklflow->add (*workflowGrid);
    mvbsd->attach_next_to(*fworklflow, Gtk::POS_TOP, 2, 1);

    // ---------------------------------------------

    Gtk::Frame* flang = Gtk::manage( new Gtk::Frame (M("PREFERENCES_DEFAULTLANG")) );
    setExpandAlignProperties(flang, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    Gtk::Grid* langGrid = Gtk::manage( new Gtk::Grid() );
    langGrid->set_column_spacing(4);
    langGrid->set_row_spacing(4);
    setExpandAlignProperties(langGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);

    ckbLangAutoDetect =  Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_LANGAUTODETECT")) );
    setExpandAlignProperties(ckbLangAutoDetect, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    Gtk::Label* langlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTLANG") + ":") );
    setExpandAlignProperties(langlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    languages = Gtk::manage( new Gtk::ComboBoxText () );
    setExpandAlignProperties(languages, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    std::vector<Glib::ustring> langs;
    parseDir (argv0 + "/languages", langs, "");

    for (size_t i = 0; i < langs.size(); i++) {
        if ("default" != langs[i] && "README" != langs[i] && "LICENSE" != langs[i]) {
            languages->append (langs[i]);
        }
    }

    Gtk::Label* langw = Gtk::manage( new Gtk::Label (Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    setExpandAlignProperties(langw, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    langGrid->attach_next_to(*ckbLangAutoDetect, Gtk::POS_LEFT, 3, 1);
    langGrid->attach_next_to(*langlab, *ckbLangAutoDetect, Gtk::POS_BOTTOM, 1, 1);
    langGrid->attach_next_to(*languages, *langlab, Gtk::POS_RIGHT, 1, 1);
    langGrid->attach_next_to(*langw, *languages, Gtk::POS_RIGHT, 1, 1);
    flang->add (*langGrid);
    mvbsd->attach_next_to(*flang, *fworklflow, Gtk::POS_BOTTOM, 2, 1);

    // ---------------------------------------------

    Gtk::Frame* ftheme = Gtk::manage( new Gtk::Frame (M("PREFERENCES_DEFAULTTHEME")) );
    setExpandAlignProperties(ftheme, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    Gtk::Grid* themeGrid = Gtk::manage( new Gtk::Grid() );
    themeGrid->set_column_spacing(4);
    themeGrid->set_row_spacing(4);
    setExpandAlignProperties(themeGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    Gtk::Label* themelab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTTHEME") + ":") );
    setExpandAlignProperties(themelab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    theme = Gtk::manage( new Gtk::ComboBoxText () );
    setExpandAlignProperties(theme, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    theme->set_active (0);
    parseThemeDir (Glib::build_filename(argv0, "themes"));

    for (size_t i = 0; i < themeFNames.size(); i++) {
        theme->append (themeFNames.at(i).shortFName);
    }

    themeGrid->attach_next_to(*themelab, Gtk::POS_LEFT, 1, 1);
    themeGrid->attach_next_to(*theme, *themelab, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* fontlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTFONT")) );
    setExpandAlignProperties(fontlab, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    fontButton = Gtk::manage( new Gtk::FontButton ());
    setExpandAlignProperties(fontButton, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    fontButton->set_use_size(true);
    fontButton->set_font_name(Glib::ustring::compose("%1 %2", options.fontFamily == "default" ? "sans" : options.fontFamily, options.fontSize));

    themeGrid->attach_next_to(*fontlab, *theme, Gtk::POS_RIGHT, 1, 1);
    themeGrid->attach_next_to(*fontButton, *fontlab, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* cpfontlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_SELECTFONT_COLPICKER") + ":") );
    setExpandAlignProperties(cpfontlab, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    colorPickerFontButton = Gtk::manage( new Gtk::FontButton ());
    setExpandAlignProperties(fontButton, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    colorPickerFontButton->set_use_size(true);
    colorPickerFontButton->set_font_name(Glib::ustring::compose("%1 %2", options.CPFontFamily == "default" ? "sans" : options.CPFontFamily, options.CPFontSize));

    themeGrid->attach_next_to(*cpfontlab, *fontButton, Gtk::POS_RIGHT, 1, 1);
    themeGrid->attach_next_to(*colorPickerFontButton, *cpfontlab, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* cutOverlayLabel = Gtk::manage( new Gtk::Label (M("PREFERENCES_CUTOVERLAYBRUSH") + ":") );
    setExpandAlignProperties(cutOverlayLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    butCropCol = Gtk::manage( new Gtk::ColorButton() );
    setExpandAlignProperties(butCropCol, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    butCropCol->set_use_alpha(true);
    themeGrid->attach_next_to(*cutOverlayLabel, *themelab, Gtk::POS_BOTTOM, 1, 1);
    themeGrid->attach_next_to(*butCropCol, *cutOverlayLabel, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* navGuideLabel = Gtk::manage( new Gtk::Label (M("PREFERENCES_NAVGUIDEBRUSH") + ":") );
    setExpandAlignProperties(navGuideLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    butNavGuideCol = Gtk::manage( new Gtk::ColorButton() );
    setExpandAlignProperties(butNavGuideCol, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    butNavGuideCol->set_use_alpha(true);
    themeGrid->attach_next_to(*navGuideLabel, *butCropCol, Gtk::POS_RIGHT, 2, 1);
    themeGrid->attach_next_to(*butNavGuideCol, *navGuideLabel, Gtk::POS_RIGHT, 1, 1);

    ftheme->add (*themeGrid);
    mvbsd->attach_next_to(*ftheme, *flang, Gtk::POS_BOTTOM, 2, 1);

    // ---------------------------------------------

    Gtk::Frame* fclip = Gtk::manage( new Gtk::Frame (M("PREFERENCES_CLIPPINGIND")));
    setExpandAlignProperties(fclip, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Grid* clipGrid = Gtk::manage( new Gtk::Grid() );
    clipGrid->set_column_spacing(4);
    clipGrid->set_row_spacing(4);
    setExpandAlignProperties(clipGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    Gtk::Label* hll = Gtk::manage( new Gtk::Label (M("PREFERENCES_HLTHRESHOLD") + ": "));
    setExpandAlignProperties(hll, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    hlThresh = Gtk::manage( new Gtk::SpinButton () );
    setExpandAlignProperties(hlThresh, false, false, Gtk::ALIGN_END, Gtk::ALIGN_BASELINE);
    hlThresh->set_digits (0);
    hlThresh->set_increments (1, 10);
    hlThresh->set_range (0, 255);
    clipGrid->attach_next_to(*hll, Gtk::POS_LEFT, 1, 1);
    clipGrid->attach_next_to(*hlThresh, *hll, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* shl = Gtk::manage( new Gtk::Label (M("PREFERENCES_SHTHRESHOLD") + ": ") );
    setExpandAlignProperties(shl, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    shThresh = Gtk::manage( new Gtk::SpinButton () );
    setExpandAlignProperties(shThresh, false, false, Gtk::ALIGN_END, Gtk::ALIGN_BASELINE);
    shThresh->show ();
    shThresh->set_digits (0);
    shThresh->set_increments (1, 10);
    shThresh->set_range (0, 255);
    clipGrid->attach_next_to(*shl, *hll, Gtk::POS_BOTTOM, 1, 1);
    clipGrid->attach_next_to(*shThresh, *shl, Gtk::POS_RIGHT, 1, 1);

    fclip->add (*clipGrid);
    mvbsd->attach_next_to(*fclip, *ftheme, Gtk::POS_BOTTOM, 1, 1);

    // ---------------------------------------------

    Gtk::Frame* fnav = Gtk::manage( new Gtk::Frame (M("PREFERENCES_NAVIGATIONFRAME")) );
    setExpandAlignProperties(fclip, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Grid* navigationGrid = Gtk::manage( new Gtk::Grid() );
    navigationGrid->set_column_spacing(4);
    navigationGrid->set_row_spacing(4);
    setExpandAlignProperties(fclip, false, false, Gtk::ALIGN_START, Gtk::ALIGN_FILL);

    Gtk::Label* panFactorLabel = Gtk::manage( new Gtk::Label (M("PREFERENCES_PANFACTORLABEL") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(panFactorLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    panFactor = Gtk::manage( new Gtk::SpinButton () );
    setExpandAlignProperties(panFactor, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    panFactor->set_digits (0);
    panFactor->set_increments (1, 5);
    panFactor->set_range (1, 10);
    navigationGrid->attach_next_to(*panFactorLabel, Gtk::POS_LEFT, 1, 1);
    navigationGrid->attach_next_to(*panFactor, *panFactorLabel, Gtk::POS_RIGHT, 1, 1);

    rememberZoomPanCheckbutton = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_REMEMBERZOOMPAN")) );
    setExpandAlignProperties(rememberZoomPanCheckbutton, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    rememberZoomPanCheckbutton->set_tooltip_text(M("PREFERENCES_REMEMBERZOOMPAN_TOOLTIP"));

    navigationGrid->attach_next_to(*rememberZoomPanCheckbutton, *panFactorLabel, Gtk::POS_BOTTOM, 2, 1);

    fnav->add (*navigationGrid);
    mvbsd->attach_next_to(*fnav, *fclip, Gtk::POS_RIGHT, 1, 1);

    // ---------------------------------------------

    Gtk::Frame* fdg = Gtk::manage( new Gtk::Frame (M("PREFERENCES_EXTERNALEDITOR")) );
    setExpandAlignProperties(fdg, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Grid* externaleditorGrid = Gtk::manage( new Gtk::Grid() );
    externaleditorGrid->set_column_spacing(4);
    externaleditorGrid->set_row_spacing(4);
    setExpandAlignProperties(externaleditorGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    edOther = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_EDITORCMDLINE") + ":"));
    setExpandAlignProperties(edOther, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    editorToSendTo = Gtk::manage( new Gtk::Entry () );
    setExpandAlignProperties(editorToSendTo, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    Gtk::RadioButton::Group ge = edOther->get_group();

#ifdef __APPLE__
    edGimp = Gtk::manage( new Gtk::RadioButton ("GIMP") );
    setExpandAlignProperties(edGimp, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    edGimp->set_group (ge);
    externaleditorGrid->attach_next_to(*edGimp, Gtk::POS_TOP, 2, 1);

    edPS = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_PSPATH") + ":"));
    setExpandAlignProperties(edPS, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    psDir = Gtk::manage( new Gtk::FileChooserButton (M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) );
    setExpandAlignProperties(psDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edPS, *edGimp, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*psDir, *edPS, Gtk::POS_RIGHT, 1, 1);
    edPS->set_group (ge);

    externaleditorGrid->attach_next_to(*edOther, *edPS, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*editorToSendTo, *edOther, Gtk::POS_RIGHT, 1, 1);
#elif defined WIN32
    edGimp = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_GIMPPATH") + ":") );
    setExpandAlignProperties(edGimp, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    gimpDir = Gtk::manage( new Gtk::FileChooserButton (M("PREFERENCES_GIMPPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) );
    setExpandAlignProperties(gimpDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edGimp, Gtk::POS_TOP, 1, 1);
    externaleditorGrid->attach_next_to(*gimpDir, *edGimp, Gtk::POS_RIGHT, 1, 1);
    edGimp->set_group (ge);

    edPS = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_PSPATH") + ":") );
    setExpandAlignProperties(edPS, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    psDir = Gtk::manage( new Gtk::FileChooserButton (M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) );
    setExpandAlignProperties(psDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edPS, *edGimp, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*psDir, *edPS, Gtk::POS_RIGHT, 1, 1);
    edPS->set_group (ge);

    externaleditorGrid->attach_next_to(*edOther, *edPS, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*editorToSendTo, *edOther, Gtk::POS_RIGHT, 1, 1);
#else
    edGimp = Gtk::manage( new Gtk::RadioButton ("GIMP") );
    setExpandAlignProperties(edGimp, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edGimp, Gtk::POS_TOP, 2, 1);
    edGimp->set_group (ge);

    externaleditorGrid->attach_next_to(*edOther, *edGimp, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*editorToSendTo, *edOther, Gtk::POS_RIGHT, 1, 1);
#endif

    fdg->add (*externaleditorGrid);
    mvbsd->attach_next_to(*fdg, *fclip, Gtk::POS_BOTTOM, 2, 1);

    langAutoDetectConn = ckbLangAutoDetect->signal_toggled().connect (sigc::mem_fun(*this, &Preferences::langAutoDetectToggled));
    tconn = theme->signal_changed().connect( sigc::mem_fun(*this, &Preferences::themeChanged) );
    fconn = fontButton->signal_font_set().connect( sigc::mem_fun(*this, &Preferences::fontChanged) );

    return mvbsd;
}

Gtk::Widget* Preferences::getFileBrowserPanel ()
{

    Gtk::VBox* mvbfb = Gtk::manage( new Gtk::VBox () );

    Gtk::Frame* fsd = Gtk::manage( new Gtk::Frame (M("PREFERENCES_STARTUPIMDIR")) );

    sdcurrent = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIRSOFTWARE")) );
    sdlast    = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIRLAST")) );
    sdhome    = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIRHOME")) );
    sdother   = Gtk::manage( new Gtk::RadioButton (M("PREFERENCES_DIROTHER") + ": ") );
    startupdir = Gtk::manage( new Gtk::Entry () );

    Gtk::Button* sdselect = Gtk::manage( new Gtk::Button () );
    sdselect->set_image (*Gtk::manage(new RTImage ("gtk-open.png")));

    Gtk::RadioButton::Group opts = sdcurrent->get_group();
    sdlast->set_group (opts);
    sdhome->set_group (opts);
    sdother->set_group (opts);

    Gtk::VBox* vbsd = Gtk::manage( new Gtk::VBox () );
    vbsd->pack_start (*sdcurrent, Gtk::PACK_SHRINK, 0);
    vbsd->pack_start (*sdlast, Gtk::PACK_SHRINK, 0);
    vbsd->pack_start (*sdhome, Gtk::PACK_SHRINK, 0);
    Gtk::HBox* otherbox = Gtk::manage( new Gtk::HBox () );
    otherbox->pack_start (*sdother, Gtk::PACK_SHRINK);
    otherbox->pack_start (*startupdir);
    otherbox->pack_end (*sdselect, Gtk::PACK_SHRINK, 4);
    vbsd->pack_start (*otherbox, Gtk::PACK_SHRINK, 0);

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
    Gtk::HBox* hbro0 = Gtk::manage( new Gtk::HBox () );
    overlayedFileNames = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_OVERLAY_FILENAMES")) );
    filmStripOverlayedFileNames = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_OVERLAY_FILENAMES_FILMSTRIP")) );
    sameThumbSize = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_FSTRIP_SAME_THUMB_HEIGHT")) );
    sameThumbSize->set_tooltip_text(M("PREFERENCES_FSTRIP_SAME_THUMB_HEIGHT_HINT"));
    ckbInternalThumbIfUntouched = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_INTERNALTHUMBIFUNTOUCHED")));

    vbro->pack_start (*showDateTime, Gtk::PACK_SHRINK, 0);
    Gtk::Label* dflab = Gtk::manage( new Gtk::Label (M("PREFERENCES_DATEFORMAT") + ":", Gtk::ALIGN_START));
    dateformat = Gtk::manage( new Gtk::Entry () );
    dateformat->set_tooltip_markup (M("PREFERENCES_DATEFORMATHINT"));
    dflab->set_tooltip_markup (M("PREFERENCES_DATEFORMATHINT"));
    hbro0->pack_start (*dflab, Gtk::PACK_SHRINK, 4);
    hbro0->pack_start (*dateformat, Gtk::PACK_SHRINK, 0);

    vbro->pack_start (*hbro0, Gtk::PACK_SHRINK, 0);
    hbro1->pack_start (*showBasicExif, Gtk::PACK_SHRINK, 0);
    hbro1->pack_start (*showExpComp, Gtk::PACK_SHRINK, 4);
    vbro->pack_start (*hbro1, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*overlayedFileNames, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*filmStripOverlayedFileNames, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*sameThumbSize, Gtk::PACK_SHRINK, 0);
    vbro->pack_start (*ckbInternalThumbIfUntouched, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hbrecent = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* labrecent = Gtk::manage( new Gtk::Label (M("PREFERENCES_MAXRECENTFOLDERS") + ":") );
    maxRecentFolders = Gtk::manage( new Gtk::SpinButton () );
    hbrecent->pack_start (*labrecent, Gtk::PACK_SHRINK, 4);
    hbrecent->pack_start (*maxRecentFolders, Gtk::PACK_SHRINK, 4);
    maxRecentFolders->set_digits (0);
    maxRecentFolders->set_increments (1, 5);
    maxRecentFolders->set_range (1, 25);
    vbro->pack_start (*hbrecent, Gtk::PACK_SHRINK, 4);

    fro->add (*vbro);


    Gtk::Frame* frmnu = Gtk::manage( new Gtk::Frame (M("PREFERENCES_MENUOPTIONS")) );
    ckbmenuGroupRank = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPRANK")) );
    ckbmenuGroupLabel = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPLABEL")) );
    ckbmenuGroupFileOperations = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPFILEOPERATIONS")) );
    ckbmenuGroupProfileOperations = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPPROFILEOPERATIONS")) );
    ckbmenuGroupExtProg = Gtk::manage( new Gtk::CheckButton (M("PREFERENCES_MENUGROUPEXTPROGS")) );
    Gtk::VBox* vbmnu = Gtk::manage( new Gtk::VBox () );

    vbmnu->pack_start (*ckbmenuGroupRank, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupLabel, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupFileOperations, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupProfileOperations, Gtk::PACK_SHRINK, 0);
    vbmnu->pack_start (*ckbmenuGroupExtProg, Gtk::PACK_SHRINK, 0);

    frmnu->add (*vbmnu);


    Gtk::Frame* fre = Gtk::manage( new Gtk::Frame (M("PREFERENCES_PARSEDEXT")) );
    Gtk::VBox* vbre = Gtk::manage( new Gtk::VBox () );
    Gtk::HBox* hb0 = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* elab = Gtk::manage( new Gtk::Label (M("PREFERENCES_PARSEDEXTADD") + ":") );
    hb0->pack_start (*elab, Gtk::PACK_SHRINK, 4);
    extension = Gtk::manage( new Gtk::Entry () );
    extension->set_width_chars(5);
    extension->set_max_width_chars(5);
    hb0->pack_start (*extension);
    addExt = Gtk::manage( new Gtk::Button () );
    delExt = Gtk::manage( new Gtk::Button () );
    moveExtUp = Gtk::manage( new Gtk::Button () );
    moveExtDown = Gtk::manage( new Gtk::Button () );
    addExt->set_tooltip_text (M("PREFERENCES_PARSEDEXTADDHINT"));
    delExt->set_tooltip_text (M("PREFERENCES_PARSEDEXTDELHINT"));
    moveExtUp->set_tooltip_text (M("PREFERENCES_PARSEDEXTUPHINT"));
    moveExtDown->set_tooltip_text (M("PREFERENCES_PARSEDEXTDOWNHINT"));
    Gtk::Image* addExtImg = Gtk::manage( new RTImage ("list-add-small.png") );
    Gtk::Image* delExtImg = Gtk::manage( new RTImage ("list-remove-red-small.png") );
    Gtk::Image* moveExtUpImg = Gtk::manage( new RTImage ("arrow-up-small.png") );
    Gtk::Image* moveExtDownImg = Gtk::manage( new RTImage ("arrow-down-small.png") );
    addExt->add (*addExtImg);
    delExt->add (*delExtImg);
    moveExtUp->set_image (*moveExtUpImg);
    moveExtDown->set_image (*moveExtDownImg);
    hb0->pack_end (*moveExtDown, Gtk::PACK_SHRINK, 4);
    hb0->pack_end (*moveExtUp, Gtk::PACK_SHRINK, 4);
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

    Gtk::HBox* hb3 = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* chlab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CACHETHUMBHEIGHT") + ":") );
    maxThumbSize = Gtk::manage( new Gtk::SpinButton () );
    hb3->pack_start (*chlab, Gtk::PACK_SHRINK, 4);
    hb3->pack_start (*maxThumbSize, Gtk::PACK_SHRINK, 4);

    maxThumbSize->set_digits (0);
    maxThumbSize->set_increments (1, 10);
    maxThumbSize->set_range (40, 800);
    vbc->pack_start (*hb3, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hb4 = Gtk::manage( new Gtk::HBox () );
    Gtk::Label* celab = Gtk::manage( new Gtk::Label (M("PREFERENCES_CACHEMAXENTRIES") + ":") );
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
    moveExtUp->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::moveExtUpPressed) );
    moveExtDown->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::moveExtDownPressed) );
    extension->signal_activate().connect( sigc::mem_fun(*this, &Preferences::addExtPressed) );
    clearThumbnails->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::clearThumbImagesPressed) );
    clearProfiles->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::clearProfilesPressed) );
    clearAll->signal_clicked().connect( sigc::mem_fun(*this, &Preferences::clearAllPressed) );

    return mvbfb;
}

Gtk::Widget* Preferences::getSoundPanel ()
{
    Gtk::VBox* pSnd = new Gtk::VBox ();

    ckbSndEnable = Gtk::manage( new Gtk::CheckButton (M("GENERAL_ENABLE")));
    sndEnableConn  = ckbSndEnable->signal_toggled().connect (sigc::mem_fun(*this, &Preferences::sndEnableToggled));

    pSnd->pack_start (*ckbSndEnable, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hblSndHelp = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* lSndHelp = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_HELP")));
    hblSndHelp->pack_start (*lSndHelp, Gtk::PACK_SHRINK, 4);
    pSnd->pack_start (*hblSndHelp, Gtk::PACK_SHRINK, 4);

    // BatchQueueDone
    Gtk::HBox* pBatchQueueDone = Gtk::manage( new Gtk::HBox() );

    Gtk::Label* lSndBatchQueueDone = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_BATCHQUEUEDONE") + Glib::ustring(":")));
    pBatchQueueDone->pack_start (*lSndBatchQueueDone, Gtk::PACK_SHRINK, 4);

    txtSndBatchQueueDone =  Gtk::manage (new Gtk::Entry());
    pBatchQueueDone->pack_end (*txtSndBatchQueueDone, Gtk::PACK_EXPAND_WIDGET, 4);

    pSnd->pack_start (*pBatchQueueDone, Gtk::PACK_SHRINK, 4);

    // LngEditProcDone
    Gtk::HBox* pSndLngEditProcDone = Gtk::manage( new Gtk::HBox() );

    Gtk::Label* lSndLngEditProcDone = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_LNGEDITPROCDONE") + Glib::ustring(":")));
    pSndLngEditProcDone->pack_start (*lSndLngEditProcDone, Gtk::PACK_SHRINK, 4);

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

    sndEnableToggled();

    return pSnd;
}

void Preferences::parseDir (Glib::ustring dirname, std::vector<Glib::ustring>& items, Glib::ustring ext)
{

    if (dirname.empty()) {
        return;
    }

    // process directory
    Glib::Dir* dir = nullptr;

    try {
        dir = new Glib::Dir (dirname);
    } catch (const Glib::Error& e) {
        return;
    }

    for (Glib::DirIterator i = dir->begin(); i != dir->end(); ++i) {
        Glib::ustring fname = Glib::build_filename(dirname, *i);
        Glib::ustring sname = *i;

        // ignore directories
        if (!Glib::file_test (fname, Glib::FILE_TEST_IS_DIR) && sname.size() >= ext.size() && sname.substr (sname.size() - ext.size(), ext.size()).casefold() == ext) {
            items.push_back (sname.substr(0, sname.size() - ext.size()));
        }
    }

    std::sort(items.begin(), items.end());
    delete dir;
}

void Preferences::parseThemeDir (Glib::ustring dirname)
{

    if (dirname.empty()) {
        return;
    }

    // process directory
    Glib::Dir* dir = nullptr;

    try {
        dir = new Glib::Dir (dirname);
    } catch (const Glib::Error& e) {
        return;
    }

    for (Glib::DirIterator i = dir->begin(); i != dir->end(); ++i) {
        Glib::ustring fname = Glib::build_filename(dirname, *i);
        Glib::ustring sname = *i;

        bool keepIt = false;

        // ignore directories and filter out unsupported theme
        if (regex->match(sname, matchInfo) && !Glib::file_test (fname, Glib::FILE_TEST_IS_DIR) && sname.size() >= 4) {
            Glib::ustring fname2 = matchInfo.fetch(1);
            Glib::ustring minMinor = matchInfo.fetch(2);
            Glib::ustring maxMinor = matchInfo.fetch(3);

            if (!minMinor.empty()) {
                guint64 minMinorVal = g_ascii_strtoll(minMinor.c_str(), 0, 0);
                if ((guint64)GTK_MINOR_VERSION >= minMinorVal) {
                    keepIt = true;
                }
            }
            if (!maxMinor.empty()) {
                guint64 maxMinorVal = g_ascii_strtoll(maxMinor.c_str(), 0, 0);
                if ((guint64)GTK_MINOR_VERSION <= maxMinorVal) {
                    keepIt = true;
                }
            }
            if (keepIt) {
                themeFNames.push_back(ThemeFilename(matchInfo.fetch(1), sname.substr(0, sname.size() - 4)));
            }
        }
    }
    std::sort(themeFNames.begin(), themeFNames.end(), [] (const ThemeFilename& firstDir, const ThemeFilename& secondDir)
            {
                return firstDir.longFName < secondDir.longFName;
            });

    delete dir;
}

void Preferences::storePreferences ()
{

    // With the new mechanism, we can't be sure of the availability of the DEFPROFILE_RAW & DEFPROFILE_IMG profiles,
    // because useBundledProfiles may be false. We're now using DEFPROFILE_INTERNAL instead, which is always available.
    moptions.defProfRaw          = rprofiles->getFullPathFromActiveRow();

    if (moptions.defProfRaw.empty()) {
        moptions.defProfRaw = DEFPROFILE_INTERNAL;
    }

    moptions.defProfImg          = iprofiles->getFullPathFromActiveRow();

    if (moptions.defProfImg.empty()) {
        moptions.defProfImg = DEFPROFILE_INTERNAL;
    }

    moptions.dateFormat          = dateformat->get_text();
    moptions.panAccelFactor      = (int)panFactor->get_value();
    moptions.rememberZoomAndPan = rememberZoomPanCheckbutton->get_active();
    moptions.fbShowDateTime  = showDateTime->get_active ();
    moptions.fbShowBasicExif = showBasicExif->get_active ();
    moptions.fbShowExpComp   = showExpComp->get_active ();
    moptions.menuGroupRank              = ckbmenuGroupRank->get_active();
    moptions.menuGroupLabel             = ckbmenuGroupLabel->get_active();
    moptions.menuGroupFileOperations    = ckbmenuGroupFileOperations->get_active();
    moptions.menuGroupProfileOperations = ckbmenuGroupProfileOperations->get_active();
    moptions.menuGroupExtProg           = ckbmenuGroupExtProg->get_active();
    moptions.highlightThreshold = (int)hlThresh->get_value ();
    moptions.shadowThreshold = (int)shThresh->get_value ();
    moptions.language        = languages->get_active_text ();
    moptions.languageAutoDetect = ckbLangAutoDetect->get_active ();
    moptions.theme           = themeFNames.at(theme->get_active_row_number ()).longFName;

    Gdk::RGBA cropCol = butCropCol->get_rgba();
    moptions.cutOverlayBrush[0] = cropCol.get_red();
    moptions.cutOverlayBrush[1] = cropCol.get_green();
    moptions.cutOverlayBrush[2] = cropCol.get_blue();
    moptions.cutOverlayBrush[3] = butCropCol->get_alpha() / 65535.0;

    Gdk::RGBA NavGuideCol = butNavGuideCol->get_rgba();
    moptions.navGuideBrush[0] = NavGuideCol.get_red();
    moptions.navGuideBrush[1] = NavGuideCol.get_green();
    moptions.navGuideBrush[2] = NavGuideCol.get_blue();
    moptions.navGuideBrush[3] = butNavGuideCol->get_alpha() / 65535.0;

    Pango::FontDescription fd(fontButton->get_font_name());
    moptions.fontFamily      = fd.get_family();
    moptions.fontSize        = fd.get_size() / Pango::SCALE;

    Pango::FontDescription cpfd(colorPickerFontButton->get_font_name());
    moptions.CPFontFamily    = cpfd.get_family();
    moptions.CPFontSize      = cpfd.get_size() / Pango::SCALE;

#ifdef WIN32
    moptions.gimpDir        = gimpDir->get_filename ();
    moptions.psDir          = psDir->get_filename ();
#elif defined __APPLE__
    moptions.psDir          = psDir->get_filename ();
#endif
    moptions.customEditorProg = editorToSendTo->get_text ();

    if (edGimp->get_active ()) {
        moptions.editorToSendTo = 1;
    }

#ifdef WIN32
    else if (edPS->get_active ()) {
        moptions.editorToSendTo = 2;
    }

#elif defined __APPLE__
    else if (edPS->get_active ()) {
        moptions.editorToSendTo = 2;
    }

#endif
    else if (edOther->get_active ()) {
        moptions.editorToSendTo = 3;
    }

    moptions.CPBPath = txtCustProfBuilderPath->get_text();
    moptions.CPBKeys = CPBKeyType(custProfBuilderLabelType->get_active_row_number());

    if (!prtProfile->get_active_row_number()) {
        moptions.rtSettings.printerProfile = "";
    } else {
        moptions.rtSettings.printerProfile = prtProfile->get_active_text ();
    }
    switch (prtIntent->get_active_row_number ()) {
    default:
    case 0:
        moptions.rtSettings.printerIntent = rtengine::RI_PERCEPTUAL;
        break;
    case 1:
        moptions.rtSettings.printerIntent = rtengine::RI_RELATIVE;
        break;
    case 2:
        moptions.rtSettings.printerIntent = rtengine::RI_ABSOLUTE;
        break;
    }
    moptions.rtSettings.printerBPC = prtBPC->get_active ();

#if !defined(__APPLE__) // monitor profile not supported on apple
    if (!monProfile->get_active_row_number()) {
        moptions.rtSettings.monitorProfile = "";
    } else {
        moptions.rtSettings.monitorProfile = monProfile->get_active_text ();
    }
    switch (monIntent->get_active_row_number ()) {
    default:
    case 0:
        moptions.rtSettings.monitorIntent = rtengine::RI_PERCEPTUAL;
        break;
    case 1:
        moptions.rtSettings.monitorIntent = rtengine::RI_RELATIVE;
        break;
    case 2:
        moptions.rtSettings.monitorIntent = rtengine::RI_ABSOLUTE;
        break;
    }
    moptions.rtSettings.monitorBPC = monBPC->get_active ();
#if defined(WIN32)
    moptions.rtSettings.autoMonitorProfile  = cbAutoMonProfile->get_active ();
#endif
#endif

    moptions.rtSettings.iccDirectory        = iccDir->get_filename ();
    moptions.rtSettings.viewingdevice       = view->get_active_row_number ();
    moptions.rtSettings.viewingdevicegrey   = grey->get_active_row_number ();
    moptions.rtSettings.viewinggreySc   = greySc->get_active_row_number ();
    //  moptions.rtSettings.autocielab            = cbAutocielab->get_active ();
    moptions.rtSettings.ciecamfloat             = cbciecamfloat->get_active ();
    moptions.rtSettings.HistogramWorking            = ckbHistogramWorking->get_active ();
    moptions.rtSettings.leveldnv   = dnv->get_active_row_number ();
    moptions.rtSettings.leveldnti   = dnti->get_active_row_number ();
    moptions.rtSettings.leveldnliss   = dnliss->get_active_row_number ();
    moptions.rtSettings.leveldnaut   = dnaut->get_active_row_number ();
    moptions.rtSettings.nrwavlevel   = dnwavlev->get_active_row_number ();
    moptions.rtSettings.leveldnautsimpl   = dnautsimpl->get_active_row_number ();
    moptions.rtSettings.daubech             = cbdaubech->get_active ();

    moptions.prevdemo = (prevdemo_t)cprevdemo->get_active_row_number ();
    moptions.serializeTiffRead = ctiffserialize->get_active();

    if (sdcurrent->get_active ()) {
        moptions.startupDir = STARTUPDIR_CURRENT;
    } else if (sdhome->get_active ()) {
        moptions.startupDir = STARTUPDIR_HOME;
    } else if (sdlast->get_active ()) {
        moptions.startupDir = STARTUPDIR_LAST;
    } else if (sdother->get_active ()) {
        moptions.startupDir = STARTUPDIR_CUSTOM;
        moptions.startupPath = startupdir->get_text();
    }

    moptions.parseExtensions.clear ();
    moptions.parseExtensionsEnabled.clear ();
    Gtk::TreeNodeChildren c = extensionModel->children ();

    for (size_t i = 0; i < c.size(); i++) {
        moptions.parseExtensions.push_back (c[i][extensionColumns.ext]);
        moptions.parseExtensionsEnabled.push_back (c[i][extensionColumns.enabled]);
    }

    moptions.maxRecentFolders = (int)maxRecentFolders->get_value();
    moptions.maxThumbnailHeight = (int)maxThumbSize->get_value ();
    moptions.maxCacheEntries = (int)maxCacheEntries->get_value ();
    moptions.overlayedFileNames = overlayedFileNames->get_active ();
    moptions.filmStripOverlayedFileNames = filmStripOverlayedFileNames->get_active();
    moptions.sameThumbSize = sameThumbSize->get_active();
    moptions.internalThumbIfUntouched = ckbInternalThumbIfUntouched->get_active ();

    moptions.saveParamsFile = saveParamsFile->get_active ();
    moptions.saveParamsCache = saveParamsCache->get_active ();
    moptions.paramsLoadLocation = (PPLoadLocation)loadParamsPreference->get_active_row_number ();
    moptions.useBundledProfiles = useBundledProfiles->get_active ();

    moptions.tunnelMetaData = ckbTunnelMetaData->get_active ();

    moptions.rtSettings.darkFramesPath =   darkFrameDir->get_filename();
    moptions.rtSettings.flatFieldsPath =   flatFieldDir->get_filename();

    moptions.clutsDir = clutsDir->get_filename();

    moptions.baBehav.resize (ADDSET_PARAM_NUM);

    for (Gtk::TreeIter sections = behModel->children().begin();  sections != behModel->children().end(); sections++)
        for (Gtk::TreeIter adjs = sections->children().begin();  adjs != sections->children().end(); adjs++) {
            moptions.baBehav[adjs->get_value (behavColumns.addsetid)] = adjs->get_value (behavColumns.badd);
        }

    int editorMode = editorLayout->get_active_row_number();
    moptions.tabbedUI = (editorMode > 1);
    moptions.multiDisplayMode = editorMode == 3 ? 1 : 0;
    moptions.mainNBVertical = editorMode == 1;

    moptions.curvebboxpos = curveBBoxPosC->get_active_row_number();
    moptions.histogramPosition = ckbHistogramPositionLeft->get_active() ? 1 : 2;
    moptions.FileBrowserToolbarSingleRow = ckbFileBrowserToolbarSingleRow->get_active();
    moptions.showFilmStripToolBar = ckbShowFilmStripToolBar->get_active();
    moptions.hideTPVScrollbar = ckbHideTPVScrollbar->get_active();
    moptions.overwriteOutputFile = chOverwriteOutputFile->get_active ();
    moptions.UseIconNoText = ckbUseIconNoText->get_active();

    moptions.rgbDenoiseThreadLimit = rgbDenoiseTreadLimitSB->get_value_as_int();
    moptions.clutCacheSize = clutCacheSizeSB->get_value_as_int();
    moptions.maxInspectorBuffers = maxInspectorBuffersSB->get_value_as_int();

    // Sounds only on Windows and Linux
#if defined(WIN32) || defined(__linux__)
    moptions.sndEnable = ckbSndEnable->get_active ();
    moptions.sndBatchQueueDone = txtSndBatchQueueDone->get_text ();
    moptions.sndLngEditProcDone     = txtSndLngEditProcDone->get_text ();
    moptions.sndLngEditProcDoneSecs = spbSndLngEditProcDoneSecs->get_value ();
#endif
}

void Preferences::fillPreferences ()
{

    tconn.block (true);
    fconn.block (true);
    sconn.block (true);
    dfconn.block (true);
    ffconn.block (true);
    rpconn.block(true);
    ipconn.block(true);
    bpconn.block(true);

    rprofiles->setActiveRowFromFullPath (moptions.defProfRaw);
    forRAWComboChanged(); // update the tooltip
    iprofiles->setActiveRowFromFullPath (moptions.defProfImg);
    forImageComboChanged(); // update the tooltip
    dateformat->set_text (moptions.dateFormat);
    panFactor->set_value (moptions.panAccelFactor);
    rememberZoomPanCheckbutton->set_active (moptions.rememberZoomAndPan);
    ctiffserialize->set_active(moptions.serializeTiffRead);

    setActiveTextOrIndex (*prtProfile, moptions.rtSettings.printerProfile, 0);
    switch (moptions.rtSettings.printerIntent) {
    default:
    case rtengine::RI_PERCEPTUAL:
        prtIntent->set_active (0);
        break;
    case rtengine::RI_RELATIVE:
        prtIntent->set_active (1);
        break;
    case rtengine::RI_ABSOLUTE:
        prtIntent->set_active (2);
        break;
    }
    prtBPC->set_active (moptions.rtSettings.printerBPC);

#if !defined(__APPLE__) // monitor profile not supported on apple
    setActiveTextOrIndex (*monProfile, moptions.rtSettings.monitorProfile, 0);
    switch (moptions.rtSettings.monitorIntent) {
    default:
    case rtengine::RI_PERCEPTUAL:
        monIntent->set_active (0);
        break;
    case rtengine::RI_RELATIVE:
        monIntent->set_active (1);
        break;
    case rtengine::RI_ABSOLUTE:
        monIntent->set_active (2);
        break;
    }
    monBPC->set_active (moptions.rtSettings.monitorBPC);
#if defined(WIN32)
    cbAutoMonProfile->set_active(moptions.rtSettings.autoMonitorProfile);
#endif
#endif

    if (Glib::file_test (moptions.rtSettings.iccDirectory, Glib::FILE_TEST_IS_DIR)) {
        iccDir->set_current_folder (moptions.rtSettings.iccDirectory);
    }

    view->set_active (moptions.rtSettings.viewingdevice);
    grey->set_active (moptions.rtSettings.viewingdevicegrey);
    greySc->set_active (moptions.rtSettings.viewinggreySc);
    dnv->set_active (moptions.rtSettings.leveldnv);
    dnti->set_active (moptions.rtSettings.leveldnti);
    dnliss->set_active (moptions.rtSettings.leveldnliss);
    dnaut->set_active (moptions.rtSettings.leveldnaut);
    dnautsimpl->set_active (moptions.rtSettings.leveldnautsimpl);
    dnwavlev->set_active (moptions.rtSettings.nrwavlevel);
    cprevdemo->set_active (moptions.prevdemo);
    cbdaubech->set_active (moptions.rtSettings.daubech);

//  cbAutocielab->set_active (moptions.rtSettings.autocielab);
    cbciecamfloat->set_active (moptions.rtSettings.ciecamfloat);
    ckbHistogramWorking->set_active (moptions.rtSettings.HistogramWorking);
    languages->set_active_text (moptions.language);
    ckbLangAutoDetect->set_active (moptions.languageAutoDetect);
    int themeNbr = getThemeRowNumber(moptions.theme);
    theme->set_active (themeNbr==-1 ? 0 : themeNbr);

    Gdk::RGBA cropCol;
    cropCol.set_rgba(moptions.cutOverlayBrush[0], moptions.cutOverlayBrush[1], moptions.cutOverlayBrush[2]);
    butCropCol->set_rgba(cropCol);
    butCropCol->set_alpha ( (unsigned short)(moptions.cutOverlayBrush[3] * 65535.0));

    Gdk::RGBA NavGuideCol;
    NavGuideCol.set_rgba(moptions.navGuideBrush[0], moptions.navGuideBrush[1], moptions.navGuideBrush[2]);
    butNavGuideCol->set_rgba(NavGuideCol);
    butNavGuideCol->set_alpha ( (unsigned short)(moptions.navGuideBrush[3] * 65535.0));

    fontButton->set_font_name(Glib::ustring::compose("%1 %2", options.fontFamily == "default" ? "sans" : options.fontFamily, options.fontSize));
    colorPickerFontButton->set_font_name(Glib::ustring::compose("%1 %2", options.CPFontFamily == "default" ? "sans" : options.CPFontFamily, options.CPFontSize));

    showDateTime->set_active (moptions.fbShowDateTime);
    showBasicExif->set_active (moptions.fbShowBasicExif);
    showExpComp->set_active (moptions.fbShowExpComp);
    ckbmenuGroupRank->set_active(moptions.menuGroupRank);
    ckbmenuGroupLabel->set_active(moptions.menuGroupLabel);
    ckbmenuGroupFileOperations->set_active(moptions.menuGroupFileOperations);
    ckbmenuGroupProfileOperations->set_active(moptions.menuGroupProfileOperations);
    ckbmenuGroupExtProg->set_active(moptions.menuGroupExtProg);

    hlThresh->set_value (moptions.highlightThreshold);
    shThresh->set_value (moptions.shadowThreshold);

    edGimp->set_active (moptions.editorToSendTo == 1);
    edOther->set_active (moptions.editorToSendTo == 3);
#ifdef WIN32
    edPS->set_active (moptions.editorToSendTo == 2);

    if (Glib::file_test (moptions.gimpDir, Glib::FILE_TEST_IS_DIR)) {
        gimpDir->set_current_folder (moptions.gimpDir);
    }

    if (Glib::file_test (moptions.psDir, Glib::FILE_TEST_IS_DIR)) {
        psDir->set_current_folder (moptions.psDir);
    }

#elif defined __APPLE__
    edPS->set_active (moptions.editorToSendTo == 2);

    if (Glib::file_test (moptions.psDir, Glib::FILE_TEST_IS_DIR)) {
        psDir->set_current_folder (moptions.psDir);
    }

#endif
    editorToSendTo->set_text (moptions.customEditorProg);

    txtCustProfBuilderPath->set_text(moptions.CPBPath);
    custProfBuilderLabelType->set_active(moptions.CPBKeys);


    if (moptions.startupDir == STARTUPDIR_CURRENT) {
        sdcurrent->set_active ();
    } else if (moptions.startupDir == STARTUPDIR_LAST) {
        sdlast->set_active ();
    } else if (moptions.startupDir == STARTUPDIR_HOME) {
        sdhome->set_active ();
    } else if (moptions.startupDir == STARTUPDIR_CUSTOM) {
        sdother->set_active ();
        startupdir->set_text (moptions.startupPath);
    }

    extensionModel->clear ();

    for (size_t i = 0; i < moptions.parseExtensions.size(); i++) {
        Gtk::TreeRow row = *(extensionModel->append());
        row[extensionColumns.enabled] = moptions.parseExtensionsEnabled[i];
        row[extensionColumns.ext]     = moptions.parseExtensions[i];
    }

    maxThumbSize->set_value (moptions.maxThumbnailHeight);
    maxRecentFolders->set_value(moptions.maxRecentFolders);
    maxCacheEntries->set_value (moptions.maxCacheEntries);
    overlayedFileNames->set_active (moptions.overlayedFileNames);
    filmStripOverlayedFileNames->set_active(moptions.filmStripOverlayedFileNames);
    sameThumbSize->set_active(moptions.sameThumbSize);
    ckbInternalThumbIfUntouched->set_active(moptions.internalThumbIfUntouched);

    saveParamsFile->set_active (moptions.saveParamsFile);
    saveParamsCache->set_active (moptions.saveParamsCache);
    loadParamsPreference->set_active (moptions.paramsLoadLocation);
    useBundledProfiles->set_active (moptions.useBundledProfiles);

    ckbTunnelMetaData->set_active (moptions.tunnelMetaData);

    if (!moptions.tabbedUI) {
        editorLayout->set_active(moptions.mainNBVertical ? 1 : 0);
    } else {
        editorLayout->set_active(moptions.multiDisplayMode ? 3 : 2);
    }

    curveBBoxPosC->set_active(moptions.curvebboxpos);
    ckbHistogramPositionLeft->set_active(moptions.histogramPosition == 1);
//   ckbHistogramWorking->set_active(moptions.histogramWorking==1);
    ckbFileBrowserToolbarSingleRow->set_active(moptions.FileBrowserToolbarSingleRow);
    ckbShowFilmStripToolBar->set_active(moptions.showFilmStripToolBar);
    ckbHideTPVScrollbar->set_active(moptions.hideTPVScrollbar);
    ckbUseIconNoText->set_active(moptions.UseIconNoText);

    rgbDenoiseTreadLimitSB->set_value(moptions.rgbDenoiseThreadLimit);
    clutCacheSizeSB->set_value(moptions.clutCacheSize);
    maxInspectorBuffersSB->set_value(moptions.maxInspectorBuffers);

    darkFrameDir->set_current_folder( moptions.rtSettings.darkFramesPath );
    darkFrameChanged ();

    flatFieldDir->set_current_folder( moptions.rtSettings.flatFieldsPath );
    flatFieldChanged ();

    clutsDir->set_current_folder( moptions.clutsDir );

    addc.block (true);
    setc.block (true);

    if (moptions.baBehav.size() == ADDSET_PARAM_NUM) {
        for (size_t i = 0; i < moptions.baBehav.size(); i++)
            for (Gtk::TreeIter sections = behModel->children().begin();  sections != behModel->children().end(); sections++)
                for (Gtk::TreeIter adjs = sections->children().begin();  adjs != sections->children().end(); adjs++)
                    if (adjs->get_value (behavColumns.addsetid) == (int)i) {
                        adjs->set_value (behavColumns.badd, moptions.baBehav[i] == 1);
                        adjs->set_value (behavColumns.bset, moptions.baBehav[i] != 1);
                        break;
                    }
    }

    addc.block (false);
    setc.block (false);
    fconn.block (false);
    tconn.block (false);
    sconn.block (false);
    dfconn.block (false);
    ffconn.block (false);
    rpconn.block(true);
    ipconn.block(true);
    bpconn.block(false);

    chOverwriteOutputFile->set_active (moptions.overwriteOutputFile);

    // Sounds only on Windows and Linux
#if defined(WIN32) || defined(__linux__)
    ckbSndEnable->set_active (moptions.sndEnable);
    txtSndBatchQueueDone->set_text (moptions.sndBatchQueueDone);
    txtSndLngEditProcDone->set_text (moptions.sndLngEditProcDone);
    spbSndLngEditProcDoneSecs->set_value (moptions.sndLngEditProcDoneSecs);
#endif
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

#if defined(WIN32)
void Preferences::autoMonProfileToggled ()
{
    monProfile->set_sensitive(!cbAutoMonProfile->get_active());
}
#endif
/*
void Preferences::autocielabToggled () {
//  cbAutocielab->set_sensitive(cbAutocielab->get_active());
}
*/
void Preferences::sndEnableToggled ()
{
    txtSndBatchQueueDone->set_sensitive(ckbSndEnable->get_active());
    txtSndLngEditProcDone->set_sensitive(ckbSndEnable->get_active());
    spbSndLngEditProcDoneSecs->set_sensitive(ckbSndEnable->get_active());
}

void Preferences::langAutoDetectToggled ()
{
    languages->set_sensitive(!ckbLangAutoDetect->get_active());
}

void Preferences::okPressed ()
{

    storePreferences ();
    workflowUpdate();
    options.copyFrom (&moptions);
    options.filterOutParsedExtensions();
    Options::save ();
    hide ();
}

void Preferences::cancelPressed ()
{
    // set the initial theme back
    if (themeFNames.at(theme->get_active_row_number ()).longFName != options.theme) {
        RTImage::setPaths(options);
        RTImage::updateImages();
        switchThemeTo(options.theme);
    }

    // set the initial font back
    Pango::FontDescription fd(fontButton->get_font_name());
    if (fd.get_family() != options.fontFamily && (fd.get_size() / Pango::SCALE) != options.fontSize) {
        switchFontTo(options.fontFamily == "default" ? "sans" : options.fontFamily, options.fontSize);
    }

    // update the profileStore
    if (useBundledProfiles->get_active () != options.useBundledProfiles) {
        // we have to rescan with the old value;
        bpconn.block(true);
        useBundledProfiles->set_active (false);
        bundledProfilesChanged();
        bpconn.block(false);
    }

    hide ();
}

void Preferences::selectStartupDir ()
{

    Gtk::FileChooserDialog dialog (getToplevelWindow (this), M("PREFERENCES_DIRSELECTDLG"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER);
//    dialog.set_transient_for(*this);

    //Add response buttons the the dialog:
    dialog.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(M("GENERAL_OPEN"), Gtk::RESPONSE_OK);

    int result = dialog.run();

    if (result == Gtk::RESPONSE_OK) {
        startupdir->set_text (dialog.get_filename());
    }
}

void Preferences::aboutPressed ()
{

    splash = new Splash (*this);
    splash->set_transient_for (*this);
    splash->signal_delete_event().connect( sigc::mem_fun(*this, &Preferences::splashClosed) );
    splash->show ();
}

void Preferences::themeChanged ()
{

    moptions.theme = themeFNames.at(theme->get_active_row_number ()).longFName;
    RTImage::setPaths(moptions);
    RTImage::updateImages();
    switchThemeTo(moptions.theme);
}

void Preferences::forRAWComboChanged ()
{
    if (!rprofiles) {
        return;
    }

    const ProfileStoreEntry *selectedEntry = rprofiles->getSelectedEntry();

    if (!selectedEntry) {
        return;
    }

    if (selectedEntry->type == PSET_FOLDER) {
        rpconn.block(true);
        rprofiles->set_active(currRawRow);
        rpconn.block(false);
    } else {
        currRawRow = rprofiles->get_active();
    }

    rprofiles->set_tooltip_text(selectedEntry->label);
}

void Preferences::forImageComboChanged ()
{
    if (!iprofiles) {
        return;
    }

    const ProfileStoreEntry *selectedEntry = iprofiles->getSelectedEntry();

    if (!selectedEntry) {
        return;
    }

    if (selectedEntry->type == PSET_FOLDER) {
        ipconn.block(true);
        iprofiles->set_active(currImgRow);
        ipconn.block(false);
    } else {
        currImgRow = rprofiles->get_active();
    }

    iprofiles->set_tooltip_text(iprofiles->getSelectedEntry()->label);
}

void Preferences::layoutComboChanged ()
{
    editorLayout->set_tooltip_text(editorLayout->get_active_text());
}

void Preferences::bundledProfilesChanged ()
{
    rpconn.block (true);
    ipconn.block (true);

    // parseProfiles does use options.useBundledProfiles, so we temporarily change its value
    bool currValue = options.useBundledProfiles;
    options.useBundledProfiles = useBundledProfiles->get_active ();

    // rescan the file's tree
    profileStore.parseProfiles(); // This will call Preferences::updateProfileList in return

    // restoring back the old value
    options.useBundledProfiles = currValue;

    ipconn.block (false);
    rpconn.block (false);
}

void Preferences::iccDirChanged ()
{
    const auto currentSelection = monProfile->get_active_text ();
    const auto profiles = rtengine::ICCStore::getInstance ()->getProfilesFromDir (iccDir->get_filename ());

    monProfile->remove_all();

    monProfile->append (M("PREFERENCES_PROFILE_NONE"));

    for (const auto& profile : profiles)
        monProfile->append (profile);

    setActiveTextOrIndex(*monProfile, currentSelection, 0);
}

void Preferences::storeCurrentValue()
{
    // TODO: Find a way to get and restore the current selection; the following line can't work anymore
    storedValueRaw = rprofiles->getFullPathFromActiveRow();
    storedValueImg = iprofiles->getFullPathFromActiveRow();
}

void Preferences::updateProfileList()
{
    rprofiles->updateProfileList();
    iprofiles->updateProfileList();
}

void Preferences::restoreValue()
{
    if (!rprofiles->setActiveRowFromFullPath(storedValueRaw)) {
        moptions.defProfRaw = DEFPROFILE_INTERNAL;
        rpconn.block(true);
        rprofiles->setInternalEntry();
        rpconn.block(false);
    }

    currRawRow = rprofiles->get_active();

    if (!iprofiles->setActiveRowFromFullPath(storedValueImg)) {
        moptions.defProfImg = DEFPROFILE_INTERNAL;
        ipconn.block(true);
        iprofiles->setInternalEntry();
        ipconn.block(false);
    }

    currImgRow = iprofiles->get_active();

    storedValueRaw = "";
    storedValueImg = "";
}

void Preferences::switchThemeTo(Glib::ustring newTheme)
{

    Glib::ustring filename(Glib::build_filename(argv0, "themes", newTheme + ".css"));

    if (!themecss) {
        themecss = Gtk::CssProvider::create();
        Glib::RefPtr<Gdk::Screen> screen = Gdk::Screen::get_default();
        Gtk::StyleContext::add_provider_for_screen(screen, themecss, GTK_STYLE_PROVIDER_PRIORITY_USER);
    }

    try {
        themecss->load_from_path (filename);
    } catch (Glib::Error &err) {
        printf("Error: Can't load css file \"%s\"\nMessage: %s\n", filename.c_str(), err.what().c_str());
    } catch (...) {
        printf("Error: Can't load css file \"%s\"\n", filename.c_str());
    }
}

void Preferences::fontChanged ()
{

    Pango::FontDescription fd(fontButton->get_font_name());
    switchFontTo(fd.get_family(), fd.get_size() / Pango::SCALE);
}

void Preferences::switchFontTo(const Glib::ustring &newFontFamily, const int newFontSize)
{

    if (newFontFamily != "default") {
        if (!fontcss) {
            fontcss = Gtk::CssProvider::create();
            Glib::RefPtr<Gdk::Screen> screen = Gdk::Screen::get_default();
            Gtk::StyleContext::add_provider_for_screen(screen, fontcss, GTK_STYLE_PROVIDER_PRIORITY_USER);
        }

        try {
            //GTK318
            #if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
        	fontcss->load_from_data (Glib::ustring::compose("* { font-family: %1; font-size: %2px }", newFontFamily, newFontSize));
            #else
            fontcss->load_from_data (Glib::ustring::compose("* { font-family: %1; font-size: %2pt }", newFontFamily, newFontSize));
            #endif
            //GTK318
        } catch (Glib::Error &err) {
            printf("Error: \"%s\"\n", err.what().c_str());
        } catch (...) {
            printf("Error: Can't find the font named \"%s\"\n", newFontFamily.c_str());
        }
    }
}

void Preferences::workflowUpdate ()
{

    if(moptions.tabbedUI != options.tabbedUI) {
        parent->MoveFileBrowserToMain();
        parent->CloseOpenEditors();
        parent->SetMainCurrent();

        if(moptions.tabbedUI) {
            parent->epanel->hide();
            parent->set_title_decorated("");
        } else {
            parent->epanel->show_all();
            parent->set_title_decorated(parent->epanel->getFileName());
        }
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

    if(moptions.showFilmStripToolBar != options.showFilmStripToolBar) {
        // Update the visibility of FB toolbar
        parent->updateFBToolBarVisibility(moptions.showFilmStripToolBar);
    }

    if(moptions.histogramPosition != options.histogramPosition) {
        // Update the position of the Histogram
        parent->updateHistogramPosition(options.histogramPosition, moptions.histogramPosition);
    }

    if(  moptions.rtSettings.printerProfile != options.rtSettings.printerProfile
       ||moptions.rtSettings.printerBPC     != options.rtSettings.printerBPC
       ||moptions.rtSettings.printerIntent  != options.rtSettings.printerIntent)
    {
        // Update the position of the Histogram
        parent->updateProfiles(moptions.rtSettings.printerProfile, moptions.rtSettings.printerIntent, moptions.rtSettings.printerBPC);
    }

}

void Preferences::addExtPressed ()
{

    Gtk::TreeNodeChildren c = extensionModel->children ();

    for (size_t i = 0; i < c.size(); i++)
        if (c[i][extensionColumns.ext] == extension->get_text ()) {
            return;
        }

    Gtk::TreeRow row = *(extensionModel->append());

    row[extensionColumns.enabled] = true;
    row[extensionColumns.ext]     = extension->get_text ();
}

void Preferences::delExtPressed ()
{

    extensionModel->erase (extensions->get_selection()->get_selected ());
}

void Preferences::moveExtUpPressed ()
{
    const Glib::RefPtr<Gtk::TreeSelection> selection = extensions->get_selection ();
    if (!selection)
        return;

    const Gtk::TreeModel::iterator selected = selection->get_selected ();
    if (!selected || selected == extensionModel->children ().begin ())
        return;

    Gtk::TreeModel::iterator previous = selected;
    --previous;
    extensionModel->iter_swap (selected, previous);
}

void Preferences::moveExtDownPressed ()
{
    const Glib::RefPtr<Gtk::TreeSelection> selection = extensions->get_selection ();
    if (!selection)
        return;

    const Gtk::TreeModel::iterator selected = selection->get_selected ();
    if (!selected)
        return;

    Gtk::TreeModel::iterator next = selected;
    if (++next)
        extensionModel->iter_swap (selected, next);
}

void Preferences::clearProfilesPressed ()
{

    cacheMgr->clearProfiles ();
}

void Preferences::clearThumbImagesPressed ()
{

    cacheMgr->clearImages ();
}

void Preferences::clearAllPressed ()
{

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
    int t1, t2;
    rtengine::dfm.getStat(t1, t2);
    Glib::ustring s = Glib::ustring::compose("%1: %2 %3, %4 %5", M("PREFERENCES_DARKFRAMEFOUND"), t1, M("PREFERENCES_DARKFRAMESHOTS"), t2, M("PREFERENCES_DARKFRAMETEMPLATES"));
    dfLabel->set_text(s);
}

void Preferences::updateFFinfos()
{
    int t1, t2;
    rtengine::ffm.getStat(t1, t2);
    Glib::ustring s = Glib::ustring::compose("%1: %2 %3, %4 %5", M("PREFERENCES_FLATFIELDFOUND"), t1, M("PREFERENCES_FLATFIELDSHOTS"), t2, M("PREFERENCES_FLATFIELDTEMPLATES"));
    ffLabel->set_text(s);
}

bool Preferences::splashClosed(GdkEventAny* event)
{
    delete splash;
    splash = nullptr;
    return true;
}

void Preferences::behAddAllPressed ()
{

    if (moptions.baBehav.size() == ADDSET_PARAM_NUM) {
        for (size_t i = 0; i < moptions.baBehav.size(); i++)
            for (Gtk::TreeIter sections = behModel->children().begin();  sections != behModel->children().end(); sections++)
                for (Gtk::TreeIter adjs = sections->children().begin();  adjs != sections->children().end(); adjs++)
                    if (adjs->get_value (behavColumns.addsetid) == (int)i) {
                        adjs->set_value (behavColumns.badd, true);
                        adjs->set_value (behavColumns.bset, false);
                        break;
                    }
    }
}

void Preferences::behSetAllPressed ()
{

    if (moptions.baBehav.size() == ADDSET_PARAM_NUM) {
        for (size_t i = 0; i < moptions.baBehav.size(); i++)
            for (Gtk::TreeIter sections = behModel->children().begin();  sections != behModel->children().end(); sections++)
                for (Gtk::TreeIter adjs = sections->children().begin();  adjs != sections->children().end(); adjs++)
                    if (adjs->get_value (behavColumns.addsetid) == (int)i) {
                        adjs->set_value (behavColumns.badd, false);
                        adjs->set_value (behavColumns.bset, true);
                        break;
                    }
    }
}
