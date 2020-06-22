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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <sigc++/slot.h>
#include "preferences.h"
#include "multilangmgr.h"
#include "splash.h"
#include "cachemanager.h"
#include "addsetids.h"
#include "../rtengine/dfmanager.h"
#include "../rtengine/ffmanager.h"
#include "../rtengine/iccstore.h"
#include "../rtengine/procparams.h"
#include <sstream>
#include "rtimage.h"
#include "rtwindow.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
void placeSpinBox(Gtk::Container* where, Gtk::SpinButton* &spin, const std::string &labelText, int digits, int inc0, int inc1, int maxLength, int range0, int range1, const std::string &toolTip = "") {
    Gtk::HBox* HB = Gtk::manage ( new Gtk::HBox () );
    HB->set_spacing (4);
    if (!toolTip.empty()) {
        HB->set_tooltip_text (M (toolTip));
    }
    Gtk::Label* label = Gtk::manage ( new Gtk::Label (M (labelText) + ":", Gtk::ALIGN_START));
    spin = Gtk::manage ( new Gtk::SpinButton () );
    spin->set_digits (digits);
    spin->set_increments (inc0, inc1);
    spin->set_max_length (maxLength); // Will this be sufficient? :)
    spin->set_range (range0, range1);
    HB->pack_start (*label, Gtk::PACK_SHRINK, 0);
    HB->pack_end (*spin, Gtk::PACK_SHRINK, 0);
    where->add(*HB);
}
}

extern Glib::ustring argv0;
Glib::RefPtr<Gtk::CssProvider> themecss;
Glib::RefPtr<Gtk::CssProvider> fontcss;

Preferences::Preferences(RTWindow *rtwindow)
    : Gtk::Dialog(M("MAIN_BUTTON_PREFERENCES"), *rtwindow, true)
    , regex(Glib::Regex::create (THEMEREGEXSTR, Glib::RegexCompileFlags::REGEX_CASELESS))
    , splash(nullptr)
    , rprofiles(nullptr)
    , iprofiles(nullptr)
    , parent(rtwindow)
    , newFont(false)
    , newCPFont(false)
{

    moptions.copyFrom(&options);

    set_size_request(650, -1);
    set_default_size(options.preferencesWidth, options.preferencesHeight);

    Pango::FontDescription defaultFont = get_style_context()->get_font();
    initialFontFamily = defaultFont.get_family();
    initialFontSize = defaultFont.get_size() / Pango::SCALE;

    Gtk::Box* mainBox = get_content_area();
//GTK318
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
    mainBox->set_spacing(8);
#endif
//GTK318

    Gtk::Notebook* nb = Gtk::manage(new Gtk::Notebook());
    nb->set_scrollable(true);
    nb->set_name("PrefNotebook");
    mainBox->pack_start(*nb);

    Gtk::Button* about  = Gtk::manage(new Gtk::Button(M("GENERAL_ABOUT")));
    Gtk::Button* ok     = Gtk::manage(new Gtk::Button(M("GENERAL_OK")));
    Gtk::Button* cancel = Gtk::manage(new Gtk::Button(M("GENERAL_CANCEL")));

    about->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::aboutPressed));
    ok->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::okPressed));
    cancel->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::cancelPressed));

    get_action_area()->pack_start(*about);
    get_action_area()->pack_end(*ok);
    get_action_area()->pack_end(*cancel);

    nb->append_page(*getGeneralPanel(), M("PREFERENCES_TAB_GENERAL"));
    nb->append_page(*getImageProcessingPanel(), M("PREFERENCES_TAB_IMPROC"));
    nb->append_page(*getDynamicProfilePanel(), M("PREFERENCES_TAB_DYNAMICPROFILE"));
    nb->append_page(*getFileBrowserPanel(), M("PREFERENCES_TAB_BROWSER"));
    nb->append_page(*getColorManPanel(), M("PREFERENCES_TAB_COLORMGR"));
    nb->append_page(*getBatchProcPanel(), M("PREFERENCES_BATCH_PROCESSING"));
    nb->append_page(*getPerformancePanel(), M("PREFERENCES_TAB_PERFORMANCE"));
    // Sounds only on Windows and Linux
#if defined(WIN32) || defined(__linux__)
    nb->append_page(*getSoundsPanel(), M("PREFERENCES_TAB_SOUND"));
#endif
    nb->set_current_page(0);

    ProfileStore::getInstance()->addListener(this);

    fillPreferences();

    show_all_children();
}


Preferences::~Preferences()
{

    ProfileStore::getInstance()->removeListener(this);
    get_size(options.preferencesWidth, options.preferencesHeight);
}

int Preferences::getThemeRowNumber (const Glib::ustring& longThemeFName)
{

    if (regex->match(longThemeFName + ".css", matchInfo)) {
        for (size_t i = 0 ; i < themeFNames.size(); ++i) {
            if (themeFNames.at(i).longFName == longThemeFName) {
                return (int)i;
            }
        }
    }

    return -1;
}

Gtk::Widget* Preferences::getBatchProcPanel()
{
    swBatchProc = Gtk::manage(new Gtk::ScrolledWindow());
    swBatchProc->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::VBox* vbBatchProc = Gtk::manage (new Gtk::VBox ());

    Gtk::ScrolledWindow* behscrollw = Gtk::manage(new Gtk::ScrolledWindow());
    behscrollw->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    behscrollw->set_size_request(-1, 60);
    Gtk::VBox* vbbeh = Gtk::manage(new Gtk::VBox());
    vbbeh->pack_start(*behscrollw, Gtk::PACK_EXPAND_WIDGET);
    Gtk::Frame* behFrame = Gtk::manage(new Gtk::Frame(M("PREFERENCES_BEHAVIOR")));
    behFrame->add(*vbbeh);
    vbBatchProc->pack_start (*behFrame, Gtk::PACK_EXPAND_WIDGET, 4);
    Gtk::TreeView* behTreeView = Gtk::manage(new Gtk::TreeView());
    behscrollw->add(*behTreeView);

    behModel = Gtk::TreeStore::create(behavColumns);
    behTreeView->set_model(behModel);

    behTreeView->append_column(M("PREFERENCES_PROPERTY"), behavColumns.label);
    behTreeView->append_column_editable(M("PREFERENCES_ADD"), behavColumns.badd);
    behTreeView->append_column_editable(M("PREFERENCES_SET"), behavColumns.bset);

    Gtk::CellRendererToggle* cr_add = static_cast<Gtk::CellRendererToggle*>(behTreeView->get_column(1)->get_first_cell());
    Gtk::CellRendererToggle* cr_set = static_cast<Gtk::CellRendererToggle*>(behTreeView->get_column(2)->get_first_cell());

    cr_add->set_radio(true);
    cr_add->set_property("xalign", 0.0f);
    sigc::connection addc = cr_add->signal_toggled().connect(sigc::mem_fun(*this, &Preferences::behAddRadioToggled));
    cr_set->set_radio(true);
    cr_set->set_property("xalign", 0.0f);
    sigc::connection setc = cr_set->signal_toggled().connect(sigc::mem_fun(*this, &Preferences::behSetRadioToggled));

    behTreeView->get_column(1)->add_attribute(*cr_add, "visible", behavColumns.visible);
    behTreeView->get_column(2)->add_attribute(*cr_set, "visible", behavColumns.visible);

    // fill model
    Gtk::TreeModel::iterator mi, ci;

    /*
     * The TRUE/FALSE values of appendBehavList are replaced by the one defined in options.cc,
     */
    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_EXPOSURE_LABEL"));
    appendBehavList(mi, M("TP_EXPOSURE_EXPCOMP"), ADDSET_TC_EXPCOMP, false);
    appendBehavList(mi, M("TP_EXPOSURE_COMPRHIGHLIGHTS"), ADDSET_TC_HLCOMPAMOUNT, false);
    appendBehavList(mi, M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), ADDSET_TC_HLCOMPTHRESH, false);
    appendBehavList(mi, M("TP_EXPOSURE_BLACKLEVEL"), ADDSET_TC_BLACKLEVEL, false);
    appendBehavList(mi, M("TP_EXPOSURE_COMPRSHADOWS"), ADDSET_TC_SHCOMP, false);
    appendBehavList(mi, M("TP_EXPOSURE_BRIGHTNESS"), ADDSET_TC_BRIGHTNESS, false);
    appendBehavList(mi, M("TP_EXPOSURE_CONTRAST"), ADDSET_TC_CONTRAST, false);
    appendBehavList(mi, M("TP_EXPOSURE_SATURATION"), ADDSET_TC_SATURATION, false);

    mi = behModel->append ();
    mi->set_value(behavColumns.label, M("TP_EPD_LABEL"));
    appendBehavList(mi, M("TP_EPD_STRENGTH"), ADDSET_EPD_STRENGTH, false);
    appendBehavList(mi, M("TP_EPD_GAMMA"), ADDSET_EPD_GAMMA, false);
    appendBehavList(mi, M("TP_EPD_EDGESTOPPING"), ADDSET_EPD_EDGESTOPPING, false);
    appendBehavList(mi, M("TP_EPD_SCALE"), ADDSET_EPD_SCALE, false);
    appendBehavList(mi, M("TP_EPD_REWEIGHTINGITERATES"), ADDSET_EPD_REWEIGHTINGITERATES, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_TM_FATTAL_LABEL"));
    appendBehavList (mi, M ("TP_TM_FATTAL_AMOUNT"), ADDSET_FATTAL_AMOUNT, false);
    appendBehavList (mi, M ("TP_TM_FATTAL_THRESHOLD"), ADDSET_FATTAL_THRESHOLD, false);
    appendBehavList (mi, M ("TP_TM_FATTAL_ANCHOR"), ADDSET_FATTAL_ANCHOR, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_RETINEX_LABEL"));
    appendBehavList(mi, M("TP_RETINEX_STRENGTH"), ADDSET_RETI_STR, false);
    appendBehavList(mi, M("TP_RETINEX_NEIGHBOR"), ADDSET_RETI_NEIGH, false);
    appendBehavList(mi, M("TP_RETINEX_VARIANCE"), ADDSET_RETI_VART, false);
    appendBehavList(mi, M("TP_RETINEX_GAMMA"), ADDSET_RETI_GAM, false);
    appendBehavList(mi, M("TP_RETINEX_SLOPE"), ADDSET_RETI_SLO, false);
    appendBehavList(mi, M("TP_RETINEX_GAIN"), ADDSET_RETI_GAIN, false);
    appendBehavList(mi, M("TP_RETINEX_OFFSET"), ADDSET_RETI_OFFS, false);
    appendBehavList(mi, M("TP_RETINEX_THRESHOLD"), ADDSET_RETI_LIMD, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_SHADOWSHLIGHTS_LABEL"));
    appendBehavList(mi, M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), ADDSET_SH_HIGHLIGHTS, false);
    appendBehavList(mi, M("TP_SHADOWSHLIGHTS_SHADOWS"), ADDSET_SH_SHADOWS, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_LABCURVE_LABEL"));
    appendBehavList(mi, M("TP_LABCURVE_BRIGHTNESS"), ADDSET_LC_BRIGHTNESS, false);
    appendBehavList(mi, M("TP_LABCURVE_CONTRAST"), ADDSET_LC_CONTRAST, false);
    appendBehavList(mi, M("TP_LABCURVE_CHROMATICITY"), ADDSET_LC_CHROMATICITY, false);

    mi = behModel->append();  // Used for both Resize and Post-Resize sharpening
    mi->set_value(behavColumns.label, M("TP_SHARPENING_LABEL"));
    appendBehavList (mi, M ("TP_SHARPENING_CONTRAST"), ADDSET_SHARP_CONTRAST, false);
    appendBehavList(mi, M("TP_SHARPENING_RADIUS"), ADDSET_SHARP_RADIUS, false);
    appendBehavList(mi, M("TP_SHARPENING_AMOUNT"), ADDSET_SHARP_AMOUNT, false);
    appendBehavList(mi, M("TP_SHARPENING_RLD_DAMPING"), ADDSET_SHARP_DAMPING, false);
    appendBehavList(mi, M("TP_SHARPENING_RLD_ITERATIONS"), ADDSET_SHARP_ITER, false);
    appendBehavList(mi, M("TP_SHARPENING_EDTOLERANCE"), ADDSET_SHARP_EDGETOL, false);
    appendBehavList(mi, M("TP_SHARPENING_HALOCONTROL"), ADDSET_SHARP_HALOCTRL, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_LOCALCONTRAST_LABEL"));
    appendBehavList(mi, M("TP_LOCALCONTRAST_RADIUS"), ADDSET_LOCALCONTRAST_RADIUS, false);
    appendBehavList(mi, M("TP_LOCALCONTRAST_AMOUNT"), ADDSET_LOCALCONTRAST_AMOUNT, false);
    appendBehavList(mi, M("TP_LOCALCONTRAST_DARKNESS"), ADDSET_LOCALCONTRAST_DARKNESS, false);
    appendBehavList(mi, M("TP_LOCALCONTRAST_LIGHTNESS"), ADDSET_LOCALCONTRAST_LIGHTNESS, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_SHARPENEDGE_LABEL"));
    appendBehavList(mi, M("TP_SHARPENEDGE_PASSES"), ADDSET_SHARPENEDGE_PASS, false);
    appendBehavList(mi, M("TP_SHARPENEDGE_AMOUNT"), ADDSET_SHARPENEDGE_AMOUNT, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_SHARPENMICRO_LABEL"));
    appendBehavList(mi, M("TP_SHARPENMICRO_AMOUNT"), ADDSET_SHARPENMICRO_AMOUNT, false);
    appendBehavList (mi, M ("TP_SHARPENMICRO_CONTRAST"), ADDSET_SHARPENMICRO_CONTRAST, false);
    appendBehavList(mi, M("TP_SHARPENMICRO_UNIFORMITY"), ADDSET_SHARPENMICRO_UNIFORMITY, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_DIRPYRDENOISE_LABEL"));
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_LUMINANCE_SMOOTHING"), ADDSET_DIRPYRDN_LUMA, true);
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_LUMINANCE_DETAIL"), ADDSET_DIRPYRDN_LUMDET, true);
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_CHROMINANCE_MASTER"), ADDSET_DIRPYRDN_CHROMA, true);
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_CHROMINANCE_REDGREEN"), ADDSET_DIRPYRDN_CHROMARED, true);
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_CHROMINANCE_BLUEYELLOW"), ADDSET_DIRPYRDN_CHROMABLUE, true);
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_MAIN_GAMMA"), ADDSET_DIRPYRDN_GAMMA, true);
    appendBehavList (mi, M ("TP_DIRPYRDENOISE_MEDIAN_PASSES"), ADDSET_DIRPYRDN_PASSES, true);

    mi = behModel->append();
    mi->set_value ( behavColumns.label, M ("TP_DEHAZE_LABEL") );
    appendBehavList ( mi, M ( "TP_DEHAZE_STRENGTH" ), ADDSET_DEHAZE_STRENGTH, true );

    mi = behModel->append ();
    mi->set_value(behavColumns.label, M("TP_WBALANCE_LABEL"));
    appendBehavList(mi, M("TP_WBALANCE_TEMPERATURE"), ADDSET_WB_TEMPERATURE, true);
    appendBehavList(mi, M("TP_WBALANCE_GREEN"), ADDSET_WB_GREEN, true);
    appendBehavList(mi, M("TP_WBALANCE_EQBLUERED"), ADDSET_WB_EQUAL, true);
    appendBehavList(mi, M("TP_WBALANCE_TEMPBIAS"), ADDSET_WB_TEMPBIAS, true);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_COLORAPP_LABEL"));
    appendBehavList (mi, M("TP_COLORAPP_LABEL_SCENE") + " - " + M("TP_COLORAPP_ABSOLUTELUMINANCE"), ADDSET_CAT_ADAPTSCENE, true);
    appendBehavList (mi, M("TP_COLORAPP_LABEL_VIEWING") + " - " + M("TP_COLORAPP_ABSOLUTELUMINANCE"), ADDSET_CAT_ADAPTVIEWING, true);
    appendBehavList(mi, M("TP_COLORAPP_LIGHT"), ADDSET_CAT_LIGHT, true);
    appendBehavList(mi, M("TP_COLORAPP_BRIGHT"), ADDSET_CAT_BRIGHT, true);
    appendBehavList(mi, M("TP_COLORAPP_CHROMA"), ADDSET_CAT_CHROMA, true);
    appendBehavList (mi, M ("TP_COLORAPP_CHROMA_S"), ADDSET_CAT_CHROMA_S, true);
    appendBehavList (mi, M ("TP_COLORAPP_CHROMA_M"), ADDSET_CAT_CHROMA_M, true);
    appendBehavList(mi, M("TP_COLORAPP_RSTPRO"), ADDSET_CAT_RSTPRO, true);
    appendBehavList(mi, M("TP_COLORAPP_CONTRAST"), ADDSET_CAT_CONTRAST, true);
    appendBehavList(mi, M("TP_COLORAPP_CONTRAST_Q"), ADDSET_CAT_CONTRAST_Q, true);
    appendBehavList(mi, M("TP_COLORAPP_HUE"), ADDSET_CAT_HUE, true);
    appendBehavList(mi, M("TP_COLORAPP_BADPIXSL"), ADDSET_CAT_BADPIX, true);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_VIBRANCE_LABEL"));
    appendBehavList(mi, M("TP_VIBRANCE_PASTELS"), ADDSET_VIBRANCE_PASTELS, false);
    appendBehavList(mi, M("TP_VIBRANCE_SATURATED"), ADDSET_VIBRANCE_SATURATED, false);


    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_CHMIXER_LABEL"));
    appendBehavList(mi, M("TP_CHMIXER_RED") + ", " + M("TP_CHMIXER_GREEN") + ", " + M("TP_CHMIXER_BLUE"), ADDSET_CHMIXER, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_BWMIX_LABEL"));
    appendBehavList(mi, M("TP_BWMIX_MIXC"), ADDSET_BLACKWHITE_HUES, false);
    appendBehavList(mi, M("TP_BWMIX_GAMMA"), ADDSET_BLACKWHITE_GAMMA, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_FILMSIMULATION_LABEL"));
    appendBehavList(mi, M("TP_FILMSIMULATION_STRENGTH"), ADDSET_FILMSIMULATION_STRENGTH, true);

    mi = behModel->append();
    mi->set_value ( behavColumns.label, M ("TP_SOFTLIGHT_LABEL") );
    appendBehavList ( mi, M ( "TP_SOFTLIGHT_STRENGTH" ), ADDSET_SOFTLIGHT_STRENGTH, true );

    mi = behModel->append ();
    mi->set_value(behavColumns.label, M("TP_COLORTONING_LABEL"));
    appendBehavList(mi, M("TP_COLORTONING_SPLITCOCO"), ADDSET_COLORTONING_SPLIT, true);
    appendBehavList(mi, M("TP_COLORTONING_SATURATIONTHRESHOLD"), ADDSET_COLORTONING_SATTHRESHOLD, true);
    appendBehavList(mi, M("TP_COLORTONING_SATURATEDOPACITY"), ADDSET_COLORTONING_SATOPACITY, true);
    appendBehavList(mi, M("TP_COLORTONING_BALANCE"), ADDSET_COLORTONING_BALANCE, true);
    appendBehavList(mi, M("TP_COLORTONING_STRENGTH"), ADDSET_COLORTONING_STRENGTH, true);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_ROTATE_LABEL"));
    appendBehavList(mi, M("TP_ROTATE_DEGREE"), ADDSET_ROTATE_DEGREE, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_RESIZE_LABEL"));
    appendBehavList(mi, M("TP_RESIZE_SCALE"), ADDSET_RESIZE_SCALE, true);


    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_DISTORTION_LABEL"));
    appendBehavList(mi, M("TP_DISTORTION_AMOUNT"), ADDSET_DIST_AMOUNT, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_PERSPECTIVE_LABEL"));
    appendBehavList(mi, M("TP_PERSPECTIVE_METHOD_SIMPLE") + " - " + M("TP_PERSPECTIVE_HORIZONTAL") + ", " + M("TP_PERSPECTIVE_VERTICAL"), ADDSET_PERSPECTIVE, false);
    appendBehavList(mi, M("TP_PERSPECTIVE_CAMERA_FOCAL_LENGTH") + ", " + M("TP_PERSPECTIVE_CAMERA_CROP_FACTOR"), ADDSET_PERSP_CAM_FOCAL_LENGTH, false);
    appendBehavList(mi, M("TP_PERSPECTIVE_CAMERA_FRAME") + " - " + M("TP_PERSPECTIVE_CAMERA_SHIFT_HORIZONTAL") + ", " + M("TP_PERSPECTIVE_CAMERA_SHIFT_VERTICAL"), ADDSET_PERSP_CAM_SHIFT, false);
    appendBehavList(mi, M("TP_PERSPECTIVE_CAMERA_FRAME") + " - " + M("TP_PERSPECTIVE_CAMERA_ROLL") + ", " + M("TP_PERSPECTIVE_CAMERA_YAW") + ", " + M("TP_PERSPECTIVE_CAMERA_PITCH"), ADDSET_PERSP_CAM_ANGLE, false);
    appendBehavList(mi, M("TP_PERSPECTIVE_POST_CORRECTION_ADJUSTMENT_FRAME") + " - " + M("TP_PERSPECTIVE_PROJECTION_SHIFT_HORIZONTAL") + ", " + M("TP_PERSPECTIVE_PROJECTION_SHIFT_VERTICAL"), ADDSET_PERSP_PROJ_SHIFT, false);
    appendBehavList(mi, M("TP_PERSPECTIVE_PROJECTION_ROTATE"), ADDSET_PERSP_PROJ_ROTATE, false);
    appendBehavList(mi, M("TP_PERSPECTIVE_RECOVERY_FRAME") + " - " + M("TP_PERSPECTIVE_PROJECTION_YAW") + ", " + M("TP_PERSPECTIVE_PROJECTION_PITCH"), ADDSET_PERSP_PROJ_ANGLE, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_GRADIENT_LABEL"));
    appendBehavList(mi, M("TP_GRADIENT_DEGREE"), ADDSET_GRADIENT_DEGREE, false);
    appendBehavList(mi, M("TP_GRADIENT_FEATHER"), ADDSET_GRADIENT_FEATHER, false);
    appendBehavList(mi, M("TP_GRADIENT_STRENGTH"), ADDSET_GRADIENT_STRENGTH, false);
    appendBehavList(mi, M("TP_GRADIENT_CENTER_X") + ", " + M("TP_GRADIENT_CENTER_Y"), ADDSET_GRADIENT_CENTER, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_PCVIGNETTE_LABEL"));
    appendBehavList(mi, M("TP_PCVIGNETTE_STRENGTH"), ADDSET_PCVIGNETTE_STRENGTH, false);
    appendBehavList(mi, M("TP_PCVIGNETTE_FEATHER"), ADDSET_PCVIGNETTE_FEATHER, false);
    appendBehavList(mi, M("TP_PCVIGNETTE_ROUNDNESS"), ADDSET_PCVIGNETTE_ROUNDNESS, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_CACORRECTION_LABEL"));
    appendBehavList(mi, M("TP_CACORRECTION_BLUE") + ", " + M("TP_CACORRECTION_RED"), ADDSET_CA, true);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_VIGNETTING_LABEL"));
    appendBehavList(mi, M("TP_VIGNETTING_AMOUNT"), ADDSET_VIGN_AMOUNT, false);
    appendBehavList(mi, M("TP_VIGNETTING_RADIUS"), ADDSET_VIGN_RADIUS, false);
    appendBehavList(mi, M("TP_VIGNETTING_STRENGTH"), ADDSET_VIGN_STRENGTH, false);
    appendBehavList(mi, M("TP_VIGNETTING_CENTER_X") + ", " + M("TP_VIGNETTING_CENTER_Y"), ADDSET_VIGN_CENTER, false);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_DIRPYREQUALIZER_LABEL"));
    appendBehavList(mi, M("TP_EXPOSURE_CONTRAST"), ADDSET_DIRPYREQ, true);
    appendBehavList(mi, M("TP_DIRPYREQUALIZER_THRESHOLD"), ADDSET_DIRPYREQ_THRESHOLD, true);
    appendBehavList(mi, M("TP_DIRPYREQUALIZER_SKIN"), ADDSET_DIRPYREQ_SKINPROTECT, true);

    mi = behModel->append();
    mi->set_value(behavColumns.label, M("TP_WAVELET_LABEL"));
    appendBehavList(mi, M("TP_WAVELET_LEVELS"), ADDSET_WA_THRES, true);
    appendBehavList(mi, M("TP_WAVELET_THRESHOLD"), ADDSET_WA_THRESHOLD, true);
    appendBehavList(mi, M("TP_WAVELET_THRESHOLD2"), ADDSET_WA_THRESHOLD2, true);
    appendBehavList(mi, M("TP_WAVELET_CHRO"), ADDSET_WA_CHRO, true);
    appendBehavList(mi, M("TP_WAVELET_CHR"), ADDSET_WA_CHROMA, true);
    appendBehavList(mi, M("TP_WAVELET_SKIN"), ADDSET_WA_SKINPROTECT, true);
    appendBehavList(mi, M("TP_WAVELET_EDRAD"), ADDSET_WA_EDGRAD, true);
    appendBehavList(mi, M("TP_WAVELET_EDVAL"), ADDSET_WA_EDGVAL, true);
    appendBehavList(mi, M("TP_WAVELET_RESCON"), ADDSET_WA_RESCON, true);
    appendBehavList(mi, M("TP_WAVELET_THR"), ADDSET_WA_THRR, true);
    appendBehavList(mi, M("TP_WAVELET_RESCONH"), ADDSET_WA_RESCONH, true);
    appendBehavList(mi, M("TP_WAVELET_THRH"), ADDSET_WA_THRRH, true);
    appendBehavList (mi, M ("TP_WAVELET_RADIUS"), ADDSET_WA_RADIUS, true);
    appendBehavList(mi, M("TP_WAVELET_RESCHRO"), ADDSET_WA_RESCHRO, true);
    appendBehavList(mi, M("TP_WAVELET_TMSTRENGTH"), ADDSET_WA_TMRS, true);
    appendBehavList (mi, M ("TP_WAVELET_TMEDGS"), ADDSET_WA_EDGS, true);
    appendBehavList (mi, M ("TP_WAVELET_TMSCALE"), ADDSET_WA_SCALE, true);
    appendBehavList(mi, M("TP_WAVELET_SKY"), ADDSET_WA_SKYPROTECT, true);
    appendBehavList(mi, M("TP_WAVELET_CONTRA"), ADDSET_WA_CONTRAST, true);
    appendBehavList(mi, M("TP_WAVELET_STRENGTH"), ADDSET_WA_STRENGTH, true);
    appendBehavList(mi, M("TP_WAVELET_COMPGAMMA"), ADDSET_WA_GAMMA, true);
    appendBehavList(mi, M("TP_WAVELET_EDGEDETECT"), ADDSET_WA_EDGEDETECT, true);
    appendBehavList(mi, M("TP_WAVELET_EDGEDETECTTHR"), ADDSET_WA_EDGEDETECTTHR, true);
    appendBehavList(mi, M("TP_WAVELET_EDGEDETECTTHR2"), ADDSET_WA_EDGEDETECTTHR2, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_RAW_SENSOR_BAYER_LABEL"));
    appendBehavList (mi, M ("TP_RAW_FALSECOLOR"), ADDSET_BAYER_FALSE_COLOR_SUPPRESSION, false);
    appendBehavList (mi, M ("TP_RAW_DCBITERATIONS") + ", " + M("TP_RAW_LMMSEITERATIONS"), ADDSET_BAYER_ITER, false);
    appendBehavList (mi, M ("TP_RAW_DUALDEMOSAICCONTRAST"), ADDSET_BAYER_DUALDEMOZCONTRAST, false);
    appendBehavList (mi, M ("TP_RAW_PIXELSHIFTSIGMA"), ADDSET_BAYER_PS_SIGMA, false);
    appendBehavList (mi, M ("TP_RAW_PIXELSHIFTSMOOTH"), ADDSET_BAYER_PS_SMOOTH, false);
    appendBehavList (mi, M ("TP_RAW_PIXELSHIFTEPERISO"), ADDSET_BAYER_PS_EPERISO, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_RAW_SENSOR_XTRANS_LABEL"));
    appendBehavList (mi, M ("TP_RAW_FALSECOLOR"), ADDSET_XTRANS_FALSE_COLOR_SUPPRESSION, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_PREPROCESS_LABEL"));
    appendBehavList (mi, M ("TP_PREPROCESS_GREENEQUIL"), ADDSET_PREPROCESS_GREENEQUIL, false);
    appendBehavList (mi, M ("TP_PREPROCESS_LINEDENOISE"), ADDSET_PREPROCESS_LINEDENOISE, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_EXPOS_WHITEPOINT_LABEL"));
    appendBehavList (mi, M ("TP_RAWEXPOS_LINEAR"), ADDSET_RAWEXPOS_LINEAR, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_EXPOS_BLACKPOINT_LABEL"));
    appendBehavList (mi, M ("TP_RAWEXPOS_RGB"), ADDSET_RAWEXPOS_BLACKS, false);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_FLATFIELD_LABEL"));
    appendBehavList (mi, M ("TP_FLATFIELD_CLIPCONTROL"), ADDSET_RAWFFCLIPCONTROL, true);

    mi = behModel->append ();
    mi->set_value (behavColumns.label, M("MAIN_TAB_RAW") + " - " + M("TP_RAWCACORR_LABEL"));
    appendBehavList (mi, M ("TP_RAWCACORR_CARED") + ", " + M ("TP_RAWCACORR_CABLUE"), ADDSET_RAWCACORR, true);

    behTreeView->expand_all();

    behAddAll = Gtk::manage(new Gtk::Button(M("PREFERENCES_BEHADDALL")));
    behSetAll = Gtk::manage(new Gtk::Button(M("PREFERENCES_BEHSETALL")));
    behAddAll->set_tooltip_markup(M("PREFERENCES_BEHADDALLHINT"));
    behSetAll->set_tooltip_markup(M("PREFERENCES_BEHSETALLHINT"));

    behAddAll->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::behAddAllPressed));
    behSetAll->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::behSetAllPressed));

    Gtk::HBox* buttonpanel1 = Gtk::manage(new Gtk::HBox());
    buttonpanel1->pack_end(*behSetAll, Gtk::PACK_SHRINK, 4);
    buttonpanel1->pack_end(*behAddAll, Gtk::PACK_SHRINK, 4);
    vbbeh->pack_start(*buttonpanel1, Gtk::PACK_SHRINK, 4);

    chOverwriteOutputFile = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_OVERWRITEOUTPUTFILE")));
    vbBatchProc->pack_start (*chOverwriteOutputFile, Gtk::PACK_SHRINK, 4);

    swBatchProc->add(*vbBatchProc);
    return swBatchProc;
}

void Preferences::appendBehavList(Gtk::TreeModel::iterator& parent, Glib::ustring label, int id, bool set)
{

    Gtk::TreeModel::iterator ci = behModel->append(parent->children());
    ci->set_value(behavColumns.label, label);
    ci->set_value(behavColumns.visible, true);
    ci->set_value(behavColumns.badd, !set);
    ci->set_value(behavColumns.bset, set);
    ci->set_value(behavColumns.addsetid, id);
}

void Preferences::behAddSetRadioToggled(const Glib::ustring& path, bool add)
{
    Gtk::TreeModel::iterator iter = behModel->get_iter(path);
    iter->set_value(behavColumns.badd, add);
    iter->set_value(behavColumns.bset, !add);
}

void Preferences::behAddRadioToggled(const Glib::ustring& path)
{
    behAddSetRadioToggled(path, true);
}

void Preferences::behSetRadioToggled(const Glib::ustring& path)
{
    behAddSetRadioToggled(path, false);
}


Gtk::Widget *Preferences::getDynamicProfilePanel()
{
    swDynamicProfile = Gtk::manage(new Gtk::ScrolledWindow());
    swDynamicProfile->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    dynProfilePanel = Gtk::manage(new DynamicProfilePanel());

    swDynamicProfile->add(*dynProfilePanel);
    return swDynamicProfile;
}


Gtk::Widget* Preferences::getImageProcessingPanel ()
{
    swImageProcessing = Gtk::manage(new Gtk::ScrolledWindow());
    swImageProcessing->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::VBox* vbImageProcessing = Gtk::manage (new Gtk::VBox ());

    Gtk::Frame* fpp = Gtk::manage(new Gtk::Frame(M("PREFERENCES_IMPROCPARAMS")));
    Gtk::VBox* vbpp = Gtk::manage(new Gtk::VBox());
    Gtk::Label* drlab = Gtk::manage(new Gtk::Label(M("PREFERENCES_FORRAW") + ":", Gtk::ALIGN_START));
    rprofiles = Gtk::manage(new ProfileStoreComboBox());
    const ProfileStoreEntry* dynpse = ProfileStore::getInstance()->getInternalDynamicPSE();
    rprofiles->addRow(dynpse);
    setExpandAlignProperties(rprofiles, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    rprofiles->set_size_request(50, -1);
    rpconn = rprofiles->signal_changed().connect(sigc::mem_fun(*this, &Preferences::forRAWComboChanged));
    Gtk::Label* drimg = Gtk::manage(new Gtk::Label(M("PREFERENCES_FORIMAGE") + ":", Gtk::ALIGN_START));
    iprofiles = Gtk::manage(new ProfileStoreComboBox());
    iprofiles->addRow(dynpse);
    iprofiles->set_size_request(50, -1);
    setExpandAlignProperties(iprofiles, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    ipconn = iprofiles->signal_changed().connect(sigc::mem_fun(*this, &Preferences::forImageComboChanged));
    Gtk::Table* defpt = Gtk::manage(new Gtk::Table(2, 2));
    defpt->attach(*drlab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach(*rprofiles, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    defpt->attach(*drimg, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    defpt->attach(*iprofiles, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    vbpp->pack_start(*defpt, Gtk::PACK_SHRINK, 4);
    useBundledProfiles = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_USEBUNDLEDPROFILES")));
    bpconn = useBundledProfiles->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::bundledProfilesChanged));
    vbpp->pack_start(*useBundledProfiles, Gtk::PACK_SHRINK, 4);
    fpp->add(*vbpp);
    vbImageProcessing->pack_start (*fpp, Gtk::PACK_SHRINK, 4);

    // Custom profile builder box
    Gtk::Frame* cpfrm = Gtk::manage(new Gtk::Frame(M("PREFERENCES_CUSTPROFBUILD")));
    Gtk::Label* cplab = Gtk::manage(new Gtk::Label(M("PREFERENCES_CUSTPROFBUILDPATH") + ":", Gtk::ALIGN_START));
    txtCustProfBuilderPath = Gtk::manage(new Gtk::Entry());
    txtCustProfBuilderPath->set_tooltip_markup(M("PREFERENCES_CUSTPROFBUILDHINT"));
    Gtk::Label* cpltypelab = Gtk::manage(new Gtk::Label(M("PREFERENCES_CUSTPROFBUILDKEYFORMAT") + ":", Gtk::ALIGN_START));
    custProfBuilderLabelType = Gtk::manage(new Gtk::ComboBoxText());
    custProfBuilderLabelType->append(M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_TID"));
    custProfBuilderLabelType->append(M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_NAME"));
    custProfBuilderLabelType->append(M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_TID") + "_" + M("PREFERENCES_CUSTPROFBUILDKEYFORMAT_NAME"));
    Gtk::Table* cpbt = Gtk::manage(new Gtk::Table(2, 2));
    cpbt->attach(*cplab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    cpbt->attach(*txtCustProfBuilderPath, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    cpbt->attach(*cpltypelab, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    cpbt->attach(*custProfBuilderLabelType, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    cpfrm->add(*cpbt);
    vbImageProcessing->pack_start (*cpfrm, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fdp = Gtk::manage(new Gtk::Frame(M("PREFERENCES_PROFILEHANDLING")));
    Gtk::Table* vbdp = Gtk::manage(new Gtk::Table(2, 2));
    saveParamsPreference = Gtk::manage(new Gtk::ComboBoxText());
    saveParamsPreference->append(M("PREFERENCES_PROFILESAVEINPUT"));
    saveParamsPreference->append(M("PREFERENCES_PROFILESAVECACHE"));
    saveParamsPreference->append(M("PREFERENCES_PROFILESAVEBOTH"));
    Gtk::Label *splab = Gtk::manage (new Gtk::Label (M ("PREFERENCES_PROFILESAVELOCATION") + ":", Gtk::ALIGN_START));
    vbdp->attach(*splab, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 2, 2);
    vbdp->attach(*saveParamsPreference, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    Gtk::Label* lplab = Gtk::manage (new Gtk::Label (M ("PREFERENCES_PROFILELOADPR") + ":", Gtk::ALIGN_START));
    loadParamsPreference = Gtk::manage(new Gtk::ComboBoxText());
    loadParamsPreference->append(M("PREFERENCES_PROFILEPRCACHE"));
    loadParamsPreference->append(M("PREFERENCES_PROFILEPRFILE"));
    vbdp->attach(*lplab, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, 2, 2);
    vbdp->attach(*loadParamsPreference, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL | Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    fdp->add(*vbdp);
    vbImageProcessing->pack_start (*fdp, Gtk::PACK_SHRINK, 4);

//    Gtk::Frame* fdf = Gtk::manage (new Gtk::Frame (M ("PREFERENCES_DARKFRAME")) );
//    Gtk::HBox* hb42 = Gtk::manage (new Gtk::HBox ());
//    darkFrameDir = Gtk::manage (new Gtk::FileChooserButton (M ("PREFERENCES_DIRDARKFRAMES"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    // Directories
    Gtk::Frame* cdf = Gtk::manage(new Gtk::Frame(M("PREFERENCES_DIRECTORIES")));
    Gtk::Grid* dirgrid = Gtk::manage(new Gtk::Grid());
    setExpandAlignProperties(dirgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Label *dfLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_DIRDARKFRAMES") + ":"));
    setExpandAlignProperties(dfLab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    darkFrameDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_DIRDARKFRAMES"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(darkFrameDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    dfLabel = Gtk::manage(new Gtk::Label("Found:"));
    setExpandAlignProperties(dfLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    dirgrid->attach_next_to(*dfLab, Gtk::POS_TOP, 1, 1);
    dirgrid->attach_next_to(*darkFrameDir, *dfLab, Gtk::POS_RIGHT, 1, 1);
    dirgrid->attach_next_to(*dfLabel, *darkFrameDir, Gtk::POS_RIGHT, 1, 1);

    dfconn = darkFrameDir->signal_selection_changed().connect ( sigc::mem_fun (*this, &Preferences::darkFrameChanged));

    // FLATFIELD
    Gtk::Label *ffLab = Gtk::manage(new Gtk::Label(M("PREFERENCES_FLATFIELDSDIR") + ":"));
    setExpandAlignProperties(ffLab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    flatFieldDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_FLATFIELDSDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(flatFieldDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    ffLabel = Gtk::manage(new Gtk::Label("Found:"));
    setExpandAlignProperties(ffLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    dirgrid->attach_next_to(*ffLab, *dfLab, Gtk::POS_BOTTOM, 1, 1);
    dirgrid->attach_next_to(*flatFieldDir, *ffLab, Gtk::POS_RIGHT, 1, 1);
    dirgrid->attach_next_to(*ffLabel, *flatFieldDir, Gtk::POS_RIGHT, 1, 1);

    ffconn = flatFieldDir->signal_selection_changed().connect ( sigc::mem_fun (*this, &Preferences::flatFieldChanged));

    //Cluts Dir
    Gtk::Label *clutsDirLabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_CLUTSDIR") + ":"));
    setExpandAlignProperties(clutsDirLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    clutsDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_CLUTSDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(clutsDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* clutsRestartNeeded = Gtk::manage(new Gtk::Label(Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")"));
    setExpandAlignProperties(clutsRestartNeeded, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    dirgrid->attach_next_to(*clutsDirLabel, *ffLab, Gtk::POS_BOTTOM, 1, 1);
    dirgrid->attach_next_to(*clutsDir, *clutsDirLabel, Gtk::POS_RIGHT, 1, 1);
    dirgrid->attach_next_to(*clutsRestartNeeded, *clutsDir, Gtk::POS_RIGHT, 1, 1);

    cdf->add(*dirgrid);
    vbImageProcessing->pack_start (*cdf, Gtk::PACK_SHRINK, 4 );

    // Crop
    Gtk::Frame *cropFrame = Gtk::manage(new Gtk::Frame(M("PREFERENCES_CROP")));
    Gtk::Grid *cropGrid = Gtk::manage(new Gtk::Grid());
    Gtk::Label *cropGuidesLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_CROP_GUIDES") + ": ", Gtk::ALIGN_START));
    cropGuidesCombo = Gtk::manage(new Gtk::ComboBoxText());
    cropGuidesCombo->append(M("PREFERENCES_CROP_GUIDES_NONE"));
    cropGuidesCombo->append(M("PREFERENCES_CROP_GUIDES_FRAME"));
    cropGuidesCombo->append(M("PREFERENCES_CROP_GUIDES_FULL"));
    cropAutoFitCB = Gtk::manage(new Gtk::CheckButton());
    Gtk::Label *cropAutoFitLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_CROP_AUTO_FIT"), Gtk::ALIGN_START));
    cropAutoFitLbl->set_line_wrap(true);
    cropAutoFitCB->add(*cropAutoFitLbl);
    cropGrid->attach(*cropGuidesLbl, 0, 0, 1, 1);
    cropGrid->attach(*cropGuidesCombo, 1, 0, 1, 1);
    cropGrid->attach(*cropAutoFitCB, 0, 1, 2, 1);
    cropFrame->add(*cropGrid);
    vbImageProcessing->pack_start(*cropFrame, Gtk::PACK_SHRINK, 4);

    swImageProcessing->add(*vbImageProcessing);

    return swImageProcessing;
}

Gtk::Widget* Preferences::getPerformancePanel()
{
    swPerformance = Gtk::manage(new Gtk::ScrolledWindow());
    swPerformance->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::VBox* vbPerformance = Gtk::manage ( new Gtk::VBox () );
    vbPerformance->set_spacing (4);

    Gtk::Frame* fprevdemo = Gtk::manage(new Gtk::Frame(M("PREFERENCES_PREVDEMO")));
    Gtk::HBox* hbprevdemo = Gtk::manage(new Gtk::HBox(false, 4));
    Gtk::Label* lprevdemo = Gtk::manage (new Gtk::Label (M("PREFERENCES_PREVDEMO_LABEL"), Gtk::ALIGN_START));
    cprevdemo = Gtk::manage(new Gtk::ComboBoxText());
    cprevdemo->append(M("PREFERENCES_PREVDEMO_FAST"));
    cprevdemo->append(M("PREFERENCES_PREVDEMO_SIDECAR"));
    cprevdemo->set_active(1);
    hbprevdemo->pack_start(*lprevdemo, Gtk::PACK_SHRINK);
    hbprevdemo->pack_start(*cprevdemo);
    fprevdemo->add(*hbprevdemo);
    vbPerformance->pack_start (*fprevdemo, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* ftiffserialize = Gtk::manage(new Gtk::Frame(M("PREFERENCES_SERIALIZE_TIFF_READ")));
    Gtk::HBox* htiffserialize = Gtk::manage(new Gtk::HBox(false, 4));
    ctiffserialize = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_SERIALIZE_TIFF_READ_LABEL")));
    ctiffserialize->set_tooltip_text(M("PREFERENCES_SERIALIZE_TIFF_READ_TOOLTIP"));
    htiffserialize->pack_start(*ctiffserialize);
    ftiffserialize->add(*htiffserialize);
    vbPerformance->pack_start (*ftiffserialize, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fclut = Gtk::manage(new Gtk::Frame(M("PREFERENCES_CLUTSCACHE")));
#ifdef _OPENMP
    placeSpinBox(fclut, clutCacheSizeSB, "PREFERENCES_CLUTSCACHE_LABEL", 0, 1, 5, 2, 1, 3 * omp_get_num_procs());
#else
    placeSpinBox(fclut, clutCacheSizeSB, "PREFERENCES_CLUTSCACHE_LABEL", 0, 1, 5, 2, 1, 12);
#endif
    vbPerformance->pack_start (*fclut, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* fchunksize = Gtk::manage ( new Gtk::Frame (M ("PREFERENCES_CHUNKSIZES")) );
    Gtk::VBox* chunkSizeVB = Gtk::manage ( new Gtk::VBox () );

    Gtk::HBox* measureHB = Gtk::manage ( new Gtk::HBox () );
    measureHB->set_spacing (4);
    measureCB = Gtk::manage ( new Gtk::CheckButton (M ("PREFERENCES_PERFORMANCE_MEASURE")) );
    measureCB->set_tooltip_text (M ("PREFERENCES_PERFORMANCE_MEASURE_HINT"));
    measureHB->pack_start(*measureCB, Gtk::PACK_SHRINK, 0);
    chunkSizeVB->add(*measureHB);

    placeSpinBox(chunkSizeVB, chunkSizeAMSB, "PREFERENCES_CHUNKSIZE_RAW_AMAZE", 0, 1, 5, 2, 1, 16);
    placeSpinBox(chunkSizeVB, chunkSizeCASB, "PREFERENCES_CHUNKSIZE_RAW_CA", 0, 1, 5, 2, 1, 16);
    placeSpinBox(chunkSizeVB, chunkSizeRCDSB, "PREFERENCES_CHUNKSIZE_RAW_RCD", 0, 1, 5, 2, 1, 16);
    placeSpinBox(chunkSizeVB, chunkSizeRGBSB, "PREFERENCES_CHUNKSIZE_RGB", 0, 1, 5, 2, 1, 16);
    placeSpinBox(chunkSizeVB, chunkSizeXTSB, "PREFERENCES_CHUNKSIZE_RAW_XT", 0, 1, 5, 2, 1, 16);

    fchunksize->add (*chunkSizeVB);

    vbPerformance->pack_start (*fchunksize, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* finspect = Gtk::manage ( new Gtk::Frame (M ("PREFERENCES_INSPECT_LABEL")) );
    Gtk::VBox *inspectorvb = Gtk::manage(new Gtk::VBox());
    placeSpinBox(inspectorvb, maxInspectorBuffersSB, "PREFERENCES_INSPECT_MAXBUFFERS_LABEL", 0, 1, 5, 2, 1, 12, "PREFERENCES_INSPECT_MAXBUFFERS_TOOLTIP");

    Gtk::HBox *insphb = Gtk::manage(new Gtk::HBox());
    thumbnailInspectorMode = Gtk::manage(new Gtk::ComboBoxText());
    thumbnailInspectorMode->append(M("PREFERENCES_THUMBNAIL_INSPECTOR_JPEG"));
    thumbnailInspectorMode->append(M("PREFERENCES_THUMBNAIL_INSPECTOR_RAW"));
    thumbnailInspectorMode->append(M("PREFERENCES_THUMBNAIL_INSPECTOR_RAW_IF_NO_JPEG_FULLSIZE"));
    insphb->pack_start(*Gtk::manage(new Gtk::Label(M("PREFERENCES_THUMBNAIL_INSPECTOR_MODE") + ": ")), Gtk::PACK_SHRINK, 4);
    insphb->pack_start(*thumbnailInspectorMode);
    inspectorvb->pack_start(*insphb);
    finspect->add (*inspectorvb);
    vbPerformance->pack_start (*finspect, Gtk::PACK_SHRINK, 4);

    Gtk::Frame* threadsFrame = Gtk::manage ( new Gtk::Frame (M ("PREFERENCES_PERFORMANCE_THREADS")) );
    Gtk::VBox* threadsVBox = Gtk::manage ( new Gtk::VBox (Gtk::PACK_SHRINK, 4) );

#ifdef _OPENMP
    int maxThreadNumber = omp_get_max_threads();
#else
    int maxThreadNumber = 10;
#endif


    placeSpinBox(threadsVBox, threadsSpinBtn, "PREFERENCES_PERFORMANCE_THREADS_LABEL", 0, 1, 5, 2, 0, maxThreadNumber);

    threadsFrame->add (*threadsVBox);

    vbPerformance->pack_start (*threadsFrame, Gtk::PACK_SHRINK, 4);
    swPerformance->add(*vbPerformance);

    return swPerformance;
}

Gtk::Widget* Preferences::getColorManPanel ()
{
    swColorMan = Gtk::manage(new Gtk::ScrolledWindow());
    swColorMan->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::VBox* vbColorMan = Gtk::manage (new Gtk::VBox ());
    vbColorMan->set_spacing (4);

    iccDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_ICCDIR"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(iccDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* pdlabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_ICCDIR") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(pdlabel, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    Gtk::Grid* iccdgrid = Gtk::manage(new Gtk::Grid());
    setExpandAlignProperties(iccdgrid, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    iccdgrid->set_column_spacing(4);

    Gtk::Label* monProfileRestartNeeded = Gtk::manage ( new Gtk::Label (Glib::ustring (" (") + M ("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    setExpandAlignProperties(monProfileRestartNeeded, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    iccdgrid->attach(*pdlabel, 0, 0, 1, 1);
    iccdgrid->attach(*iccDir, 1, 0, 1, 1);
    iccdgrid->attach (*monProfileRestartNeeded, 2, 0, 1, 1);

    iccDir->signal_selection_changed().connect(sigc::mem_fun(this, &Preferences::iccDirChanged));

    vbColorMan->pack_start (*iccdgrid, Gtk::PACK_SHRINK);

    //------------------------- MONITOR ----------------------

    Gtk::Frame* fmonitor = Gtk::manage(new Gtk::Frame(M("PREFERENCES_MONITOR")));
    Gtk::Grid* gmonitor = Gtk::manage(new Gtk::Grid());
    gmonitor->set_column_spacing(4);

    monProfile = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(monProfile, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* mplabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_MONPROFILE") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(mplabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    monIntent = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(monIntent, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* milabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_MONINTENT") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(milabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    monProfile->append(M("PREFERENCES_PROFILE_NONE"));
    monProfile->set_active(0);

    const std::vector<Glib::ustring> profiles = rtengine::ICCStore::getInstance()->getProfiles(rtengine::ICCStore::ProfileType::MONITOR);

    for (const auto& profile : profiles) {
        if (profile.find("file:") != 0) {
            std::string fileis_RTv4 = profile.substr(0, 4);

            if (fileis_RTv4 != "RTv4") {
            //    printf("pro=%s \n", profile.c_str());
                monProfile->append(profile);
            }
        }
    }

    // same order as the enum
    monIntent->append(M("PREFERENCES_INTENT_PERCEPTUAL"));
    monIntent->append(M("PREFERENCES_INTENT_RELATIVE"));
    monIntent->append(M("PREFERENCES_INTENT_ABSOLUTE"));
    monIntent->set_active(1);
    monIntent->set_size_request(120, -1);

    monBPC = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_CMMBPC")));
    setExpandAlignProperties(monBPC, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    monBPC->set_active(true);

    cbAutoMonProfile = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_AUTOMONPROFILE")));
    setExpandAlignProperties(cbAutoMonProfile, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    autoMonProfileConn = cbAutoMonProfile->signal_toggled().connect(sigc::mem_fun(*this, &Preferences::autoMonProfileToggled));

    int row = 0;
    gmonitor->attach(*mplabel, 0, row, 1, 1);
#if defined(__APPLE__) // monitor profile not supported on apple
    Gtk::Label *osxwarn = Gtk::manage(new Gtk::Label(M("PREFERENCES_MONPROFILE_WARNOSX"), Gtk::ALIGN_START));
    setExpandAlignProperties(osxwarn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    gmonitor->attach(*osxwarn, 1, row, 1, 1);
#else
    gmonitor->attach(*monProfile, 1, row, 1, 1);
#endif
    ++row;
    gmonitor->attach(*cbAutoMonProfile, 1, row, 1, 1);
    ++row;
    gmonitor->attach(*milabel, 0, row, 1, 1);
    gmonitor->attach(*monIntent, 1, row, 1, 1);
    ++row;
    gmonitor->attach(*monBPC, 0, row, 2, 1);

    autoMonProfileToggled();

    fmonitor->add(*gmonitor);

    vbColorMan->pack_start (*fmonitor, Gtk::PACK_SHRINK);

    //------------------------- PRINTER ----------------------

    Gtk::Frame* fprinter = Gtk::manage(new Gtk::Frame(M("PREFERENCES_PRINTER")));
    Gtk::Grid* gprinter = Gtk::manage(new Gtk::Grid());
    gprinter->set_column_spacing(4);
    prtProfile = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(prtProfile, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* pplabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_PRTPROFILE") + ":"));
    setExpandAlignProperties(pplabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    prtIntent = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(prtIntent, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* pilabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_PRTINTENT") + ":"));
    setExpandAlignProperties(pilabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    prtProfile->append(M("PREFERENCES_PROFILE_NONE"));
    prtProfile->set_active(0);

    const std::vector<Glib::ustring> prtprofiles = rtengine::ICCStore::getInstance()->getProfiles(rtengine::ICCStore::ProfileType::PRINTER);

    for (const auto& prtprofile : prtprofiles) {
        prtProfile->append(prtprofile);
    }

    // same order as the enum
    prtIntent->append(M("PREFERENCES_INTENT_PERCEPTUAL"));
    prtIntent->append(M("PREFERENCES_INTENT_RELATIVE"));
    prtIntent->append(M("PREFERENCES_INTENT_ABSOLUTE"));
    prtIntent->set_active(1);

    prtBPC = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_CMMBPC")));
    setExpandAlignProperties(prtBPC, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    prtBPC->set_active(true);

    row = 0;
    gprinter->attach(*pplabel, 0, row, 1, 1);
    gprinter->attach(*prtProfile, 1, row, 1, 1);
    ++row;
    gprinter->attach(*pilabel, 0, row, 1, 1);
    gprinter->attach(*prtIntent, 1, row, 1, 1);
    ++row;
    gprinter->attach(*prtBPC, 0, row, 2, 1);

    autoMonProfileToggled();

    fprinter->add(*gprinter);

    vbColorMan->pack_start (*fprinter, Gtk::PACK_SHRINK);
    swColorMan->add(*vbColorMan);
    return swColorMan;
}

Gtk::Widget* Preferences::getGeneralPanel()
{
    swGeneral = Gtk::manage(new Gtk::ScrolledWindow());
    swGeneral->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::Grid* vbGeneral = Gtk::manage ( new Gtk::Grid () );
    vbGeneral->set_column_spacing (4);
    vbGeneral->set_row_spacing (4);

    Gtk::Frame* fworklflow = Gtk::manage(new Gtk::Frame(M("PREFERENCES_WORKFLOW")));
    setExpandAlignProperties(fworklflow, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    Gtk::Grid* workflowGrid = Gtk::manage(new Gtk::Grid());
    workflowGrid->set_column_spacing(4);
    workflowGrid->set_row_spacing(4);
    setExpandAlignProperties(workflowGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    Gtk::Label* flayoutlab = Gtk::manage(new Gtk::Label(M("PREFERENCES_EDITORLAYOUT") + ":"));
    setExpandAlignProperties(flayoutlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    editorLayout = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(editorLayout, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    editorLayout->append(M("PREFERENCES_SINGLETAB"));
    editorLayout->append(M("PREFERENCES_SINGLETABVERTAB"));
    editorLayout->append(M("PREFERENCES_MULTITAB"));
    editorLayout->append(M("PREFERENCES_MULTITABDUALMON"));
    editorLayout->set_active(2);
    Gtk::CellRendererText* cellRenderer = dynamic_cast<Gtk::CellRendererText*>(editorLayout->get_first_cell());
    cellRenderer->property_ellipsize() = Pango::ELLIPSIZE_MIDDLE;
    cellRenderer->property_ellipsize_set() = true;
    editorLayout->signal_changed().connect(sigc::mem_fun(*this, &Preferences::layoutComboChanged));
    layoutComboChanged(); // update the tooltip
    Gtk::Label* lNextStart = Gtk::manage(new Gtk::Label(Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")"));
    setExpandAlignProperties(lNextStart, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*flayoutlab, Gtk::POS_LEFT, 1, 1);
    workflowGrid->attach_next_to(*editorLayout, *flayoutlab, Gtk::POS_RIGHT, 1, 1);
    workflowGrid->attach_next_to(*lNextStart, *editorLayout, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* curveBBoxPosL = Gtk::manage(new Gtk::Label(M("PREFERENCES_CURVEBBOXPOS") + ":"));
    setExpandAlignProperties(curveBBoxPosL, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    curveBBoxPosC = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(curveBBoxPosC, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    curveBBoxPosC->append(M("PREFERENCES_CURVEBBOXPOS_ABOVE"));
    curveBBoxPosC->append(M("PREFERENCES_CURVEBBOXPOS_RIGHT"));
    curveBBoxPosC->append(M("PREFERENCES_CURVEBBOXPOS_BELOW"));
    curveBBoxPosC->append(M("PREFERENCES_CURVEBBOXPOS_LEFT"));
    curveBBoxPosC->set_active(1);
    Gtk::Label* curveBBoxPosRestartL = Gtk::manage(new Gtk::Label(Glib::ustring("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")"));
    setExpandAlignProperties(curveBBoxPosRestartL, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*curveBBoxPosL, *flayoutlab, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*curveBBoxPosC, *editorLayout, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*curveBBoxPosRestartL, *lNextStart, Gtk::POS_BOTTOM, 1, 1);

    Gtk::Label* complexityL = Gtk::manage(new Gtk::Label(M("PREFERENCES_COMPLEXITYLOC") + ":"));
    setExpandAlignProperties(complexityL, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    complexitylocal = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(complexitylocal, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    complexitylocal->append(M("PREFERENCES_COMPLEXITY_EXP"));
    complexitylocal->append(M("PREFERENCES_COMPLEXITY_NORM"));
    complexitylocal->set_active(1);
    workflowGrid->attach_next_to(*complexityL, *curveBBoxPosL, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*complexitylocal, *curveBBoxPosC, Gtk::POS_BOTTOM, 1, 1);

    ckbHistogramPositionLeft = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_HISTOGRAMPOSITIONLEFT")));
    setExpandAlignProperties(ckbHistogramPositionLeft, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*ckbHistogramPositionLeft, *complexityL, Gtk::POS_BOTTOM, 1, 1);

    ckbFileBrowserToolbarSingleRow = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_FILEBROWSERTOOLBARSINGLEROW")));
    setExpandAlignProperties(ckbFileBrowserToolbarSingleRow, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    ckbShowFilmStripToolBar = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_SHOWFILMSTRIPTOOLBAR")));
    setExpandAlignProperties(ckbShowFilmStripToolBar, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    workflowGrid->attach_next_to(*ckbFileBrowserToolbarSingleRow, *ckbHistogramPositionLeft, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*ckbShowFilmStripToolBar, *complexitylocal, Gtk::POS_BOTTOM, 2, 1);

    Gtk::Label* hb4label = Gtk::manage(new Gtk::Label(M("PREFERENCES_TP_LABEL")));
    setExpandAlignProperties(hb4label, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    ckbHideTPVScrollbar = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_TP_VSCROLLBAR")));
    setExpandAlignProperties(ckbHideTPVScrollbar, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*hb4label, *ckbFileBrowserToolbarSingleRow, Gtk::POS_BOTTOM, 1, 1);
    workflowGrid->attach_next_to(*ckbHideTPVScrollbar, *hb4label, Gtk::POS_RIGHT, 1, 1);
    ckbAutoSaveTpOpen = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_AUTOSAVE_TP_OPEN")));
    workflowGrid->attach_next_to(*ckbAutoSaveTpOpen, *hb4label, Gtk::POS_BOTTOM, 1, 1);
    btnSaveTpOpenNow = Gtk::manage(new Gtk::Button(M("PREFERENCES_SAVE_TP_OPEN_NOW")));
    setExpandAlignProperties(btnSaveTpOpenNow, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    workflowGrid->attach_next_to(*btnSaveTpOpenNow, *ckbAutoSaveTpOpen, Gtk::POS_RIGHT, 1, 1);

    auto save_tp_open_now =
    [&]() -> void {
        parent->writeToolExpandedStatus(moptions.tpOpen);
    };
    btnSaveTpOpenNow->signal_clicked().connect(save_tp_open_now);

    ckbshowtooltiplocallab = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_SHOWTOOLTIP")));
    setExpandAlignProperties(ckbshowtooltiplocallab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    workflowGrid->attach_next_to(*ckbshowtooltiplocallab, *ckbFileBrowserToolbarSingleRow, Gtk::POS_RIGHT, 1, 1);

    fworklflow->add(*workflowGrid);

    vbGeneral->attach_next_to (*fworklflow, Gtk::POS_TOP, 2, 1);

    // ---------------------------------------------

    Gtk::Frame* flang = Gtk::manage(new Gtk::Frame(M("PREFERENCES_LANG")));
    setExpandAlignProperties(flang, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    Gtk::Grid* langGrid = Gtk::manage(new Gtk::Grid());
    langGrid->set_column_spacing(4);
    langGrid->set_row_spacing(4);
    setExpandAlignProperties(langGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);

    ckbLangAutoDetect =  Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_LANGAUTODETECT")));
    setExpandAlignProperties(ckbLangAutoDetect, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    Gtk::Label* langlab = Gtk::manage(new Gtk::Label(M("PREFERENCES_SELECTLANG") + ":"));
    setExpandAlignProperties(langlab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    languages = Gtk::manage(new Gtk::ComboBoxText());
    setExpandAlignProperties(languages, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    std::vector<Glib::ustring> langs;
    parseDir(argv0 + "/languages", langs, "");

    for (size_t i = 0; i < langs.size(); i++) {
        if ("default" != langs[i] && "README" != langs[i] && "LICENSE" != langs[i]) {
            languages->append(langs[i]);
        }
    }

    Gtk::Label* langw = Gtk::manage(new Gtk::Label(Glib::ustring(" (") + M("PREFERENCES_APPLNEXTSTARTUP") + ")"));
    setExpandAlignProperties(langw, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    langGrid->attach_next_to(*ckbLangAutoDetect, Gtk::POS_LEFT, 3, 1);
    langGrid->attach_next_to(*langlab, *ckbLangAutoDetect, Gtk::POS_BOTTOM, 1, 1);
    langGrid->attach_next_to(*languages, *langlab, Gtk::POS_RIGHT, 1, 1);
    langGrid->attach_next_to(*langw, *languages, Gtk::POS_RIGHT, 1, 1);
    flang->add(*langGrid);
    vbGeneral->attach_next_to (*flang, *fworklflow, Gtk::POS_BOTTOM, 2, 1);

    // Appearance ---------------------------------------------

    Gtk::Frame* appearanceFrame = Gtk::manage(new Gtk::Frame(M("PREFERENCES_APPEARANCE")));

    Gtk::Grid* appearanceGrid = Gtk::manage(new Gtk::Grid());
    appearanceGrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(appearanceGrid, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    Gtk::Label* themeLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_APPEARANCE_THEME") + ":"));
    setExpandAlignProperties(themeLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    Gtk::Label* themeRestartLbl = Gtk::manage ( new Gtk::Label (Glib::ustring (" (") + M ("PREFERENCES_APPLNEXTSTARTUP") + ")") );
    setExpandAlignProperties(themeRestartLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    themeCBT = Gtk::manage(new Gtk::ComboBoxText());
    themeCBT->set_active(0);
    parseThemeDir(Glib::build_filename(argv0, "themes"));
    for (size_t i = 0; i < themeFNames.size(); i++) {
        themeCBT->append(themeFNames.at(i).shortFName);
    }

    Gtk::Label* mainFontLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_APPEARANCE_MAINFONT")));
    setExpandAlignProperties(mainFontLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    mainFontFB = Gtk::manage(new Gtk::FontButton());
    mainFontFB->set_use_size(true);
    if (options.fontFamily == "default") {
        mainFontFB->set_font_name(Glib::ustring::compose("%1 %2", initialFontFamily, initialFontSize));
    } else {
        mainFontFB->set_font_name(Glib::ustring::compose("%1 %2", options.fontFamily, options.fontSize));
    }

    Gtk::Label* colorPickerFontLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_APPEARANCE_COLORPICKERFONT") + ":"));
    setExpandAlignProperties(colorPickerFontLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    colorPickerFontFB = Gtk::manage(new Gtk::FontButton());
    colorPickerFontFB->set_use_size(true);
    if (options.fontFamily == "default") {
        colorPickerFontFB->set_font_name(Glib::ustring::compose("%1 %2", initialFontFamily, initialFontSize));
    } else {
        colorPickerFontFB->set_font_name(Glib::ustring::compose("%1 %2", options.CPFontFamily, options.CPFontSize));
    }

    Gtk::Label* cropMaskColorLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_APPEARANCE_CROPMASKCOLOR") + ":"));
    setExpandAlignProperties(cropMaskColorLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    cropMaskColorCB = Gtk::manage(new Gtk::ColorButton());
    cropMaskColorCB->set_use_alpha(true);

    Gtk::Label* navGuideColorLbl = Gtk::manage(new Gtk::Label(M("PREFERENCES_APPEARANCE_NAVGUIDECOLOR") + ":"));
    setExpandAlignProperties(navGuideColorLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    navGuideColorCB = Gtk::manage(new Gtk::ColorButton());
    navGuideColorCB->set_use_alpha(true);

    pseudoHiDPI = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_APPEARANCE_PSEUDOHIDPI") + Glib::ustring (" (") + M ("PREFERENCES_APPLNEXTSTARTUP") + ")"));
    setExpandAlignProperties(pseudoHiDPI, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);

    Gtk::VSeparator *vSep = Gtk::manage(new Gtk::VSeparator());
    

    appearanceGrid->attach(*themeLbl,           0, 0, 1, 1);
    appearanceGrid->attach(*themeCBT,           1, 0, 1, 1);
    appearanceGrid->attach(*themeRestartLbl,    2, 0, 2, 1);
    appearanceGrid->attach(*vSep,               2, 1, 1, 2);
    appearanceGrid->attach(*mainFontLbl,        0, 1, 1, 1);
    appearanceGrid->attach(*mainFontFB,         1, 1, 1, 1);
    appearanceGrid->attach(*cropMaskColorLbl,   3, 1, 1, 1);
    appearanceGrid->attach(*cropMaskColorCB,    4, 1, 1, 1);
    appearanceGrid->attach(*colorPickerFontLbl, 0, 2, 1, 1);
    appearanceGrid->attach(*colorPickerFontFB,  1, 2, 1, 1);
    appearanceGrid->attach(*navGuideColorLbl,   3, 2, 1, 1);
    appearanceGrid->attach(*navGuideColorCB,    4, 2, 1, 1);
    appearanceGrid->attach(*pseudoHiDPI,        0, 3, 5, 1);

    appearanceFrame->add(*appearanceGrid);
    vbGeneral->attach_next_to(*appearanceFrame, *flang, Gtk::POS_BOTTOM, 2, 1);

    // ---------------------------------------------

    Gtk::Frame* fclip = Gtk::manage(new Gtk::Frame(M("PREFERENCES_CLIPPINGIND")));
    setExpandAlignProperties(fclip, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Grid* clipGrid = Gtk::manage(new Gtk::Grid());
    clipGrid->set_column_spacing(4);
    clipGrid->set_row_spacing(4);
    setExpandAlignProperties(clipGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    Gtk::Label* hll = Gtk::manage(new Gtk::Label(M("PREFERENCES_HLTHRESHOLD") + ": "));
    setExpandAlignProperties(hll, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    hlThresh = Gtk::manage(new Gtk::SpinButton());
    setExpandAlignProperties(hlThresh, false, false, Gtk::ALIGN_END, Gtk::ALIGN_BASELINE);
    hlThresh->set_digits(0);
    hlThresh->set_increments(1, 10);
    hlThresh->set_range(0, 255);
    clipGrid->attach_next_to(*hll, Gtk::POS_LEFT, 1, 1);
    clipGrid->attach_next_to(*hlThresh, *hll, Gtk::POS_RIGHT, 1, 1);

    Gtk::Label* shl = Gtk::manage(new Gtk::Label(M("PREFERENCES_SHTHRESHOLD") + ": "));
    setExpandAlignProperties(shl, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    shThresh = Gtk::manage(new Gtk::SpinButton());
    setExpandAlignProperties(shThresh, false, false, Gtk::ALIGN_END, Gtk::ALIGN_BASELINE);
    shThresh->show();
    shThresh->set_digits(0);
    shThresh->set_increments(1, 10);
    shThresh->set_range(0, 255);
    clipGrid->attach_next_to(*shl, *hll, Gtk::POS_BOTTOM, 1, 1);
    clipGrid->attach_next_to(*shThresh, *shl, Gtk::POS_RIGHT, 1, 1);

    fclip->add(*clipGrid);
    vbGeneral->attach_next_to (*fclip, *appearanceFrame, Gtk::POS_BOTTOM, 1, 1);

    // ---------------------------------------------

    Gtk::Frame* fnav = Gtk::manage(new Gtk::Frame(M("PREFERENCES_NAVIGATIONFRAME")));
    setExpandAlignProperties(fclip, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Grid* navigationGrid = Gtk::manage(new Gtk::Grid());
    navigationGrid->set_column_spacing(4);
    navigationGrid->set_row_spacing(4);
    setExpandAlignProperties(fclip, false, false, Gtk::ALIGN_START, Gtk::ALIGN_FILL);

    Gtk::Label* panFactorLabel = Gtk::manage(new Gtk::Label(M("PREFERENCES_PANFACTORLABEL") + ":", Gtk::ALIGN_START));
    setExpandAlignProperties(panFactorLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    panFactor = Gtk::manage(new Gtk::SpinButton());
    setExpandAlignProperties(panFactor, true, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    panFactor->set_digits(0);
    panFactor->set_increments(1, 5);
    panFactor->set_range(1, 10);
    navigationGrid->attach_next_to(*panFactorLabel, Gtk::POS_LEFT, 1, 1);
    navigationGrid->attach_next_to(*panFactor, *panFactorLabel, Gtk::POS_RIGHT, 1, 1);

    rememberZoomPanCheckbutton = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_REMEMBERZOOMPAN")));
    setExpandAlignProperties(rememberZoomPanCheckbutton, false, false, Gtk::ALIGN_START, Gtk::ALIGN_BASELINE);
    rememberZoomPanCheckbutton->set_tooltip_text(M("PREFERENCES_REMEMBERZOOMPAN_TOOLTIP"));

    navigationGrid->attach_next_to(*rememberZoomPanCheckbutton, *panFactorLabel, Gtk::POS_BOTTOM, 2, 1);

    fnav->add(*navigationGrid);
    vbGeneral->attach_next_to (*fnav, *fclip, Gtk::POS_RIGHT, 1, 1);

    // ---------------------------------------------

    Gtk::Frame* fdg = Gtk::manage(new Gtk::Frame(M("PREFERENCES_EXTERNALEDITOR")));
    setExpandAlignProperties(fdg, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    Gtk::Grid* externaleditorGrid = Gtk::manage(new Gtk::Grid());
    externaleditorGrid->set_column_spacing(4);
    externaleditorGrid->set_row_spacing(4);
    setExpandAlignProperties(externaleditorGrid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    edOther = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_EDITORCMDLINE") + ":"));
    setExpandAlignProperties(edOther, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    editorToSendTo = Gtk::manage(new Gtk::Entry());
    setExpandAlignProperties(editorToSendTo, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_BASELINE);
    Gtk::RadioButton::Group ge = edOther->get_group();

#ifdef __APPLE__
    edGimp = Gtk::manage(new Gtk::RadioButton("GIMP"));
    setExpandAlignProperties(edGimp, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    edGimp->set_group(ge);
    externaleditorGrid->attach_next_to(*edGimp, Gtk::POS_TOP, 2, 1);

    edPS = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_PSPATH") + ":"));
    setExpandAlignProperties(edPS, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    psDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(psDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edPS, *edGimp, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*psDir, *edPS, Gtk::POS_RIGHT, 1, 1);
    edPS->set_group(ge);

    externaleditorGrid->attach_next_to(*edOther, *edPS, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*editorToSendTo, *edOther, Gtk::POS_RIGHT, 1, 1);
#elif defined WIN32
    edGimp = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_GIMPPATH") + ":"));
    setExpandAlignProperties(edGimp, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    gimpDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_GIMPPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(gimpDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edGimp, Gtk::POS_TOP, 1, 1);
    externaleditorGrid->attach_next_to(*gimpDir, *edGimp, Gtk::POS_RIGHT, 1, 1);
    edGimp->set_group(ge);

    edPS = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_PSPATH") + ":"));
    setExpandAlignProperties(edPS, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    psDir = Gtk::manage(new MyFileChooserButton(M("PREFERENCES_PSPATH"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER));
    setExpandAlignProperties(psDir, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edPS, *edGimp, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*psDir, *edPS, Gtk::POS_RIGHT, 1, 1);
    edPS->set_group(ge);

    externaleditorGrid->attach_next_to(*edOther, *edPS, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*editorToSendTo, *edOther, Gtk::POS_RIGHT, 1, 1);
#else
    edGimp = Gtk::manage(new Gtk::RadioButton("GIMP"));
    setExpandAlignProperties(edGimp, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    externaleditorGrid->attach_next_to(*edGimp, Gtk::POS_TOP, 2, 1);
    edGimp->set_group(ge);

    externaleditorGrid->attach_next_to(*edOther, *edGimp, Gtk::POS_BOTTOM, 1, 1);
    externaleditorGrid->attach_next_to(*editorToSendTo, *edOther, Gtk::POS_RIGHT, 1, 1);
#endif

    fdg->add(*externaleditorGrid);
    vbGeneral->attach_next_to (*fdg, *fclip, Gtk::POS_BOTTOM, 2, 1);
    langAutoDetectConn = ckbLangAutoDetect->signal_toggled().connect(sigc::mem_fun(*this, &Preferences::langAutoDetectToggled));
    tconn = themeCBT->signal_changed().connect ( sigc::mem_fun (*this, &Preferences::themeChanged) );
    fconn = mainFontFB->signal_font_set().connect ( sigc::mem_fun (*this, &Preferences::fontChanged) );
    cpfconn = colorPickerFontFB->signal_font_set().connect ( sigc::mem_fun (*this, &Preferences::cpFontChanged) );

    swGeneral->add(*vbGeneral);
    return swGeneral;
}

Gtk::Widget* Preferences::getFileBrowserPanel()
{
    swFileBrowser = Gtk::manage(new Gtk::ScrolledWindow());
    swFileBrowser->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::VBox* vbFileBrowser = Gtk::manage ( new Gtk::VBox () );

    Gtk::Frame* fsd = Gtk::manage(new Gtk::Frame(M("PREFERENCES_STARTUPIMDIR")));

    sdcurrent  = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_DIRSOFTWARE")));
    sdlast     = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_DIRLAST")));
    sdhome     = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_DIRHOME")));
    sdother    = Gtk::manage(new Gtk::RadioButton(M("PREFERENCES_DIROTHER") + ": "));
    startupdir = Gtk::manage(new Gtk::Entry());

    Gtk::Button* sdselect = Gtk::manage(new Gtk::Button());
    sdselect->set_image (*Gtk::manage (new RTImage ("folder-open-small.png")));

    Gtk::RadioButton::Group opts = sdcurrent->get_group();
    sdlast->set_group(opts);
    sdhome->set_group(opts);
    sdother->set_group(opts);

    Gtk::VBox* vbsd = Gtk::manage(new Gtk::VBox());
    vbsd->pack_start(*sdcurrent, Gtk::PACK_SHRINK, 0);
    vbsd->pack_start(*sdlast, Gtk::PACK_SHRINK, 0);
    vbsd->pack_start(*sdhome, Gtk::PACK_SHRINK, 0);
    Gtk::HBox* otherbox = Gtk::manage(new Gtk::HBox());
    otherbox->pack_start(*sdother, Gtk::PACK_SHRINK);
    otherbox->pack_start(*startupdir);
    otherbox->pack_end(*sdselect, Gtk::PACK_SHRINK, 4);
    vbsd->pack_start(*otherbox, Gtk::PACK_SHRINK, 0);

    fsd->add(*vbsd);
    vbFileBrowser->pack_start (*fsd, Gtk::PACK_SHRINK, 4);

    sdselect->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::selectStartupDir));

//---


    Gtk::Frame* fro = Gtk::manage(new Gtk::Frame(M("PREFERENCES_FBROWSEROPTS")));
    showDateTime = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_SHOWDATETIME")));
    showBasicExif = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_SHOWBASICEXIF")));
    showExpComp = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_SHOWEXPOSURECOMPENSATION")));
    Gtk::VBox* vbro = Gtk::manage(new Gtk::VBox());
    Gtk::HBox* hbro1 = Gtk::manage(new Gtk::HBox());
    Gtk::HBox* hbro0 = Gtk::manage(new Gtk::HBox());
    overlayedFileNames = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_OVERLAY_FILENAMES")));
    filmStripOverlayedFileNames = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_OVERLAY_FILENAMES_FILMSTRIP")));
    sameThumbSize = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_FSTRIP_SAME_THUMB_HEIGHT")));
    sameThumbSize->set_tooltip_text(M("PREFERENCES_FSTRIP_SAME_THUMB_HEIGHT_HINT"));
    ckbInternalThumbIfUntouched = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_INTERNALTHUMBIFUNTOUCHED")));

    vbro->pack_start(*showDateTime, Gtk::PACK_SHRINK, 0);
    Gtk::Label* dflab = Gtk::manage(new Gtk::Label(M("PREFERENCES_DATEFORMAT") + ":", Gtk::ALIGN_START));
    dateformat = Gtk::manage(new Gtk::Entry());
    dateformat->set_tooltip_markup(M("PREFERENCES_DATEFORMATHINT"));
    dflab->set_tooltip_markup(M("PREFERENCES_DATEFORMATHINT"));
    hbro0->pack_start(*dflab, Gtk::PACK_SHRINK, 4);
    hbro0->pack_start(*dateformat, Gtk::PACK_SHRINK, 0);

    vbro->pack_start(*hbro0, Gtk::PACK_SHRINK, 0);
    hbro1->pack_start(*showBasicExif, Gtk::PACK_SHRINK, 0);
    hbro1->pack_start(*showExpComp, Gtk::PACK_SHRINK, 4);
    vbro->pack_start(*hbro1, Gtk::PACK_SHRINK, 0);
    vbro->pack_start(*overlayedFileNames, Gtk::PACK_SHRINK, 0);
    vbro->pack_start(*filmStripOverlayedFileNames, Gtk::PACK_SHRINK, 0);
    vbro->pack_start(*sameThumbSize, Gtk::PACK_SHRINK, 0);
    vbro->pack_start(*ckbInternalThumbIfUntouched, Gtk::PACK_SHRINK, 0);

    Gtk::HBox* hbrecent = Gtk::manage(new Gtk::HBox());
    Gtk::Label* labrecent = Gtk::manage (new Gtk::Label (M("PREFERENCES_MAXRECENTFOLDERS") + ":", Gtk::ALIGN_START));
    maxRecentFolders = Gtk::manage(new Gtk::SpinButton());
    hbrecent->pack_start(*labrecent, Gtk::PACK_SHRINK, 4);
    hbrecent->pack_start(*maxRecentFolders, Gtk::PACK_SHRINK, 4);
    maxRecentFolders->set_digits(0);
    maxRecentFolders->set_increments(1, 5);
    maxRecentFolders->set_range(1, 25);
    vbro->pack_start(*hbrecent, Gtk::PACK_SHRINK, 4);

    fro->add(*vbro);


    Gtk::Frame* frmnu = Gtk::manage(new Gtk::Frame(M("PREFERENCES_MENUOPTIONS")));

    Gtk::Grid* menuGrid = Gtk::manage(new Gtk::Grid());
    menuGrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(menuGrid, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    ckbmenuGroupRank = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_MENUGROUPRANK")));
    setExpandAlignProperties(ckbmenuGroupRank, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    ckbmenuGroupLabel = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_MENUGROUPLABEL")));
    ckbmenuGroupFileOperations = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_MENUGROUPFILEOPERATIONS")));
    setExpandAlignProperties(ckbmenuGroupFileOperations, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    ckbmenuGroupProfileOperations = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_MENUGROUPPROFILEOPERATIONS")));
    ckbmenuGroupExtProg = Gtk::manage(new Gtk::CheckButton(M("PREFERENCES_MENUGROUPEXTPROGS")));

    Gtk::Label* groupRestartNeeded = Gtk::manage(new Gtk::Label (Glib::ustring ("(") + M("PREFERENCES_APPLNEXTSTARTUP") + ")", Gtk::ALIGN_START));

    menuGrid->attach (*ckbmenuGroupRank, 0, 0, 1, 1);
    menuGrid->attach (*ckbmenuGroupLabel, 1, 0, 1, 1);
    menuGrid->attach (*ckbmenuGroupFileOperations, 0, 1, 1, 1);
    menuGrid->attach (*ckbmenuGroupProfileOperations, 1, 1, 1, 1);
    menuGrid->attach (*ckbmenuGroupExtProg, 0, 2, 1, 1);
    menuGrid->attach (*groupRestartNeeded, 1, 2, 1, 1);

    frmnu->add (*menuGrid);


    Gtk::Frame* fre = Gtk::manage(new Gtk::Frame(M("PREFERENCES_PARSEDEXT")));
    Gtk::VBox* vbre = Gtk::manage(new Gtk::VBox());
    Gtk::HBox* hb0 = Gtk::manage(new Gtk::HBox());
    Gtk::Label* elab = Gtk::manage (new Gtk::Label (M("PREFERENCES_PARSEDEXTADD") + ":", Gtk::ALIGN_START));
    hb0->pack_start(*elab, Gtk::PACK_SHRINK, 4);
    extension = Gtk::manage(new Gtk::Entry());
    extension->set_width_chars(5);
    extension->set_max_width_chars(5);
    hb0->pack_start(*extension);
    addExt = Gtk::manage(new Gtk::Button());
    delExt = Gtk::manage(new Gtk::Button());
    moveExtUp = Gtk::manage(new Gtk::Button());
    moveExtDown = Gtk::manage(new Gtk::Button());
    addExt->set_tooltip_text(M("PREFERENCES_PARSEDEXTADDHINT"));
    delExt->set_tooltip_text(M("PREFERENCES_PARSEDEXTDELHINT"));
    moveExtUp->set_tooltip_text(M("PREFERENCES_PARSEDEXTUPHINT"));
    moveExtDown->set_tooltip_text(M("PREFERENCES_PARSEDEXTDOWNHINT"));
    Gtk::Image* addExtImg = Gtk::manage ( new RTImage ("add-small.png") );
    Gtk::Image* delExtImg = Gtk::manage ( new RTImage ("remove-small.png") );
    Gtk::Image* moveExtUpImg = Gtk::manage(new RTImage("arrow-up-small.png"));
    Gtk::Image* moveExtDownImg = Gtk::manage(new RTImage("arrow-down-small.png"));
    addExt->add(*addExtImg);
    delExt->add(*delExtImg);
    moveExtUp->set_image(*moveExtUpImg);
    moveExtDown->set_image(*moveExtDownImg);
    hb0->pack_end(*moveExtDown, Gtk::PACK_SHRINK, 4);
    hb0->pack_end(*moveExtUp, Gtk::PACK_SHRINK, 4);
    hb0->pack_end(*delExt, Gtk::PACK_SHRINK, 4);
    hb0->pack_end(*addExt, Gtk::PACK_SHRINK, 4);
    extensions = Gtk::manage(new Gtk::TreeView());
    Gtk::ScrolledWindow* hscrollw = Gtk::manage(new Gtk::ScrolledWindow());
    hscrollw->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    hscrollw->add(*extensions);
    extensionModel = Gtk::ListStore::create(extensionColumns);
    extensions->set_model(extensionModel);
    extensions->append_column_editable("Enabled", extensionColumns.enabled);
    extensions->append_column("Extension", extensionColumns.ext);
    extensions->set_headers_visible(false);
    vbre->pack_start(*hscrollw);
    vbre->pack_start(*hb0, Gtk::PACK_SHRINK, 4);

    fre->add(*vbre);

    // Cache

    Gtk::Frame* frc = Gtk::manage (new Gtk::Frame(M("PREFERENCES_CACHEOPTS")));
    Gtk::VBox* vbc = Gtk::manage (new Gtk::VBox());
    frc->add(*vbc);

    Gtk::Grid* cacheGrid = Gtk::manage(new Gtk::Grid());
    cacheGrid->get_style_context()->add_class("grid-spacing");
    setExpandAlignProperties(cacheGrid, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    Gtk::Label* maxThumbHeightLbl = Gtk::manage (new Gtk::Label(M("PREFERENCES_CACHETHUMBHEIGHT") + ":"));
    setExpandAlignProperties(maxThumbHeightLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    maxThumbHeightSB = Gtk::manage (new Gtk::SpinButton());
    maxThumbHeightSB->set_digits (0);
    maxThumbHeightSB->set_increments (1, 10);
    maxThumbHeightSB->set_range (40, 800);

    Gtk::Label* maxCacheEntriesLbl = Gtk::manage (new Gtk::Label(M("PREFERENCES_CACHEMAXENTRIES") + ":"));
    setExpandAlignProperties(maxCacheEntriesLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    maxCacheEntriesSB = Gtk::manage (new Gtk::SpinButton());
    maxCacheEntriesSB->set_digits (0);
    maxCacheEntriesSB->set_increments (1, 10);
    maxCacheEntriesSB->set_range (10, 100000);

    // Separation is needed so that a button is not accidentally clicked when one wanted
    // to click a spinbox. Ideally, the separation wouldn't require attaching a widget, but how?
    Gtk::HSeparator *cacheSeparator = Gtk::manage (new  Gtk::HSeparator());
    cacheSeparator->get_style_context()->add_class("grid-row-separator");

    Gtk::Label* clearThumbsLbl = Gtk::manage (new Gtk::Label(M("PREFERENCES_CACHECLEAR_ALLBUTPROFILES")));
    setExpandAlignProperties(clearThumbsLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    Gtk::Button* clearThumbsBtn = Gtk::manage (new Gtk::Button(M("PREFERENCES_CACHECLEAR")));

    Gtk::Label* clearProfilesLbl = Gtk::manage (new Gtk::Label(M("PREFERENCES_CACHECLEAR_ONLYPROFILES")));
    setExpandAlignProperties(clearProfilesLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    Gtk::Button* clearProfilesBtn = Gtk::manage (new Gtk::Button(M("PREFERENCES_CACHECLEAR")));

    Gtk::Label* clearAllLbl = Gtk::manage (new Gtk::Label(M("PREFERENCES_CACHECLEAR_ALL")));
    setExpandAlignProperties(clearAllLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    Gtk::Button* clearAllBtn = Gtk::manage (new Gtk::Button(M("PREFERENCES_CACHECLEAR")));

    cacheGrid->attach (*maxThumbHeightLbl, 0, 0, 1, 1);
    cacheGrid->attach (*maxThumbHeightSB, 1, 0, 1, 1);
    cacheGrid->attach (*maxCacheEntriesLbl, 0, 1, 1, 1);
    cacheGrid->attach (*maxCacheEntriesSB, 1, 1, 1, 1);
    cacheGrid->attach (*cacheSeparator, 0, 2, 2, 1);
    cacheGrid->attach (*clearThumbsLbl, 0, 3, 1, 1);
    cacheGrid->attach (*clearThumbsBtn, 1, 3, 1, 1);
    if (moptions.saveParamsCache) {
        cacheGrid->attach (*clearProfilesLbl, 0, 4, 1, 1);
        cacheGrid->attach (*clearProfilesBtn, 1, 4, 1, 1);
        cacheGrid->attach (*clearAllLbl, 0, 5, 1, 1);
        cacheGrid->attach (*clearAllBtn, 1, 5, 1, 1);
    }

    vbc->pack_start (*cacheGrid, Gtk::PACK_SHRINK, 4);

    Gtk::Label* clearSafetyLbl = Gtk::manage (new Gtk::Label(M("PREFERENCES_CACHECLEAR_SAFETY")));
    setExpandAlignProperties(clearSafetyLbl, false, false, Gtk::ALIGN_START, Gtk::ALIGN_START);
    clearSafetyLbl->set_line_wrap(true);
    vbc->pack_start(*clearSafetyLbl, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hb6 = Gtk::manage(new Gtk::HBox());
    Gtk::VBox* vb6 = Gtk::manage(new Gtk::VBox());

    vb6->pack_start(*fro);
    vb6->pack_start(*frmnu);
    vb6->pack_end(*frc);
    hb6->pack_start(*vb6);
    hb6->pack_start(*fre);
    hb6->set_spacing(4);

    vbFileBrowser->pack_start (*hb6, Gtk::PACK_SHRINK, 4);

    addExt->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::addExtPressed));
    delExt->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::delExtPressed));
    moveExtUp->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::moveExtUpPressed));
    moveExtDown->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::moveExtDownPressed));
    extension->signal_activate().connect(sigc::mem_fun(*this, &Preferences::addExtPressed));
    clearThumbsBtn->signal_clicked().connect ( sigc::mem_fun (*this, &Preferences::clearThumbImagesPressed) );
    if (moptions.saveParamsCache) {
        clearProfilesBtn->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::clearProfilesPressed));
        clearAllBtn->signal_clicked().connect(sigc::mem_fun(*this, &Preferences::clearAllPressed));
    }

    swFileBrowser->add(*vbFileBrowser);
    return swFileBrowser;
}

Gtk::Widget* Preferences::getSoundsPanel ()
{
    swSounds = Gtk::manage(new Gtk::ScrolledWindow());
    swSounds->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);

    Gtk::VBox* vbSounds = Gtk::manage(new Gtk::VBox ());

    ckbSndEnable = Gtk::manage(new Gtk::CheckButton(M("GENERAL_ENABLE")));
    sndEnableConn = ckbSndEnable->signal_toggled().connect(sigc::mem_fun(*this, &Preferences::sndEnableToggled));

    vbSounds->pack_start (*ckbSndEnable, Gtk::PACK_SHRINK, 4);

    Gtk::HBox* hblSndHelp = Gtk::manage(new Gtk::HBox());
    Gtk::Label* lSndHelp = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_HELP"), Gtk::ALIGN_START));
    hblSndHelp->pack_start(*lSndHelp, Gtk::PACK_SHRINK, 4);
    vbSounds->pack_start (*hblSndHelp, Gtk::PACK_SHRINK, 4);

    // BatchQueueDone
    Gtk::HBox* pBatchQueueDone = Gtk::manage(new Gtk::HBox());

    Gtk::Label* lSndBatchQueueDone = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_QUEUEDONE") + Glib::ustring (":"), Gtk::ALIGN_START));
    pBatchQueueDone->pack_start (*lSndBatchQueueDone, Gtk::PACK_SHRINK, 4);

    txtSndBatchQueueDone = Gtk::manage(new Gtk::Entry());
    pBatchQueueDone->pack_end(*txtSndBatchQueueDone, Gtk::PACK_EXPAND_WIDGET, 4);

    vbSounds->pack_start (*pBatchQueueDone, Gtk::PACK_SHRINK, 4);

    // LngEditProcDone
    Gtk::HBox* pSndLngEditProcDone = Gtk::manage(new Gtk::HBox());

    Gtk::Label* lSndLngEditProcDone = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_LNGEDITPROCDONE") + Glib::ustring (":"), Gtk::ALIGN_START));
    pSndLngEditProcDone->pack_start(*lSndLngEditProcDone, Gtk::PACK_SHRINK, 4);

    txtSndLngEditProcDone = Gtk::manage(new Gtk::Entry());
    pSndLngEditProcDone->pack_start(*txtSndLngEditProcDone, Gtk::PACK_EXPAND_WIDGET, 4);

    Gtk::Label* lSndLngEditProcDoneSecs = Gtk::manage (new Gtk::Label (M("PREFERENCES_SND_THRESHOLDSECS") + Glib::ustring (":"), Gtk::ALIGN_START));
    pSndLngEditProcDone->pack_start(*lSndLngEditProcDoneSecs, Gtk::PACK_SHRINK, 12);

    spbSndLngEditProcDoneSecs = Gtk::manage(new Gtk::SpinButton());
    spbSndLngEditProcDoneSecs->set_digits(1);
    spbSndLngEditProcDoneSecs->set_increments(0.5, 1);
    spbSndLngEditProcDoneSecs->set_range(0, 10);
    pSndLngEditProcDone->pack_end(*spbSndLngEditProcDoneSecs, Gtk::PACK_SHRINK, 4);

    vbSounds->pack_start (*pSndLngEditProcDone, Gtk::PACK_SHRINK, 4);

    sndEnableToggled();

    swSounds->add(*vbSounds);
    return swSounds;
}

void Preferences::parseDir(Glib::ustring dirname, std::vector<Glib::ustring>& items, Glib::ustring ext)
{

    if (dirname.empty()) {
        return;
    }

    // process directory
    Glib::Dir* dir = nullptr;

    try {
        dir = new Glib::Dir(dirname);
    } catch (const Glib::Error& e) {
        return;
    }

    for (Glib::DirIterator i = dir->begin(); i != dir->end(); ++i) {
        Glib::ustring fname = Glib::build_filename(dirname, *i);
        Glib::ustring sname = *i;

        // ignore directories
        if (!Glib::file_test(fname, Glib::FILE_TEST_IS_DIR) && sname.size() >= ext.size() && sname.substr(sname.size() - ext.size(), ext.size()).casefold() == ext) {
            items.push_back(sname.substr(0, sname.size() - ext.size()));
        }
    }

    std::sort(items.begin(), items.end());
    delete dir;
}

void Preferences::parseThemeDir(Glib::ustring dirname)
{

    if (dirname.empty()) {
        return;
    }

    // process directory
    Glib::Dir* dir = nullptr;

    try {
        dir = new Glib::Dir(dirname);
    } catch (const Glib::Error& e) {
        return;
    }

    for (Glib::DirIterator i = dir->begin(); i != dir->end(); ++i) {
        Glib::ustring fname = Glib::build_filename(dirname, *i);
        Glib::ustring sname = *i;

        // ignore directories and filter out unsupported theme
        if (regex->match(sname, matchInfo) && !Glib::file_test(fname, Glib::FILE_TEST_IS_DIR) && sname.size() >= 4) {
            bool keepIt = false;
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

    std::sort(themeFNames.begin(), themeFNames.end(), [](const ThemeFilename & firstDir, const ThemeFilename & secondDir) {
        return firstDir.longFName < secondDir.longFName;
    });

    delete dir;
}

void Preferences::storePreferences()
{

    // With the new mechanism, we can't be sure of the availability of the DEFPROFILE_RAW & DEFPROFILE_IMG profiles,
    // because useBundledProfiles may be false. We're now using DEFPROFILE_INTERNAL instead, which is always available.

    moptions.defProfRaw = rprofiles->getFullPathFromActiveRow();

    if (moptions.defProfRaw.empty()) {
        moptions.defProfRaw = DEFPROFILE_INTERNAL;
    }

    moptions.defProfImg = iprofiles->getFullPathFromActiveRow();

    if (moptions.defProfImg.empty()) {
        moptions.defProfImg = DEFPROFILE_INTERNAL;
    }

    moptions.dateFormat = dateformat->get_text();
    moptions.panAccelFactor = (int)panFactor->get_value();
    moptions.rememberZoomAndPan = rememberZoomPanCheckbutton->get_active();
    moptions.fbShowDateTime = showDateTime->get_active();
    moptions.fbShowBasicExif = showBasicExif->get_active();
    moptions.fbShowExpComp = showExpComp->get_active();
    moptions.menuGroupRank = ckbmenuGroupRank->get_active();
    moptions.menuGroupLabel = ckbmenuGroupLabel->get_active();
    moptions.menuGroupFileOperations = ckbmenuGroupFileOperations->get_active();
    moptions.menuGroupProfileOperations = ckbmenuGroupProfileOperations->get_active();
    moptions.menuGroupExtProg = ckbmenuGroupExtProg->get_active();
    moptions.highlightThreshold = (int)hlThresh->get_value();
    moptions.shadowThreshold = (int)shThresh->get_value();
    moptions.language = languages->get_active_text();
    moptions.languageAutoDetect = ckbLangAutoDetect->get_active();
    moptions.theme = themeFNames.at (themeCBT->get_active_row_number ()).longFName;

    Gdk::RGBA cropCol = cropMaskColorCB->get_rgba();
    moptions.cutOverlayBrush[0] = cropCol.get_red();
    moptions.cutOverlayBrush[1] = cropCol.get_green();
    moptions.cutOverlayBrush[2] = cropCol.get_blue();
    moptions.cutOverlayBrush[3] = cropMaskColorCB->get_alpha() / 65535.0;

    Gdk::RGBA NavGuideCol = navGuideColorCB->get_rgba();
    moptions.navGuideBrush[0] = NavGuideCol.get_red();
    moptions.navGuideBrush[1] = NavGuideCol.get_green();
    moptions.navGuideBrush[2] = NavGuideCol.get_blue();
    moptions.navGuideBrush[3] = navGuideColorCB->get_alpha() / 65535.0;
    Pango::FontDescription fd (mainFontFB->get_font_name());


    if (newFont) {
        moptions.fontFamily = fd.get_family();
        moptions.fontSize = fd.get_size() / Pango::SCALE;
    }

    Pango::FontDescription cpfd (colorPickerFontFB->get_font_name());

    if (newCPFont) {
        moptions.CPFontFamily = cpfd.get_family();
        moptions.CPFontSize = cpfd.get_size() / Pango::SCALE;
    }

    moptions.pseudoHiDPISupport = pseudoHiDPI->get_active();

#ifdef WIN32
    moptions.gimpDir = gimpDir->get_filename();
    moptions.psDir = psDir->get_filename();
#elif defined __APPLE__
    moptions.psDir = psDir->get_filename();
#endif
    moptions.customEditorProg = editorToSendTo->get_text();

    if (edGimp->get_active()) {
        moptions.editorToSendTo = 1;
    }

#ifdef WIN32
    else if (edPS->get_active()) {
        moptions.editorToSendTo = 2;
    }

#elif defined __APPLE__
    else if (edPS->get_active()) {
        moptions.editorToSendTo = 2;
    }

#endif
    else if (edOther->get_active()) {
        moptions.editorToSendTo = 3;
    }

    moptions.CPBPath = txtCustProfBuilderPath->get_text();
    moptions.CPBKeys = CPBKeyType(custProfBuilderLabelType->get_active_row_number());

    if (!prtProfile->get_active_row_number()) {
        moptions.rtSettings.printerProfile = "";
    } else {
        moptions.rtSettings.printerProfile = prtProfile->get_active_text();
    }

    switch (prtIntent->get_active_row_number()) {
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

    moptions.rtSettings.printerBPC = prtBPC->get_active();

#if !defined(__APPLE__) // monitor profile not supported on apple

    if (!monProfile->get_active_row_number()) {
        moptions.rtSettings.monitorProfile = "";
    } else {
        moptions.rtSettings.monitorProfile = monProfile->get_active_text();
    }

    switch (monIntent->get_active_row_number()) {
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

    moptions.rtSettings.monitorBPC = monBPC->get_active();
    moptions.rtSettings.autoMonitorProfile = cbAutoMonProfile->get_active();
#endif

    moptions.rtSettings.iccDirectory = iccDir->get_filename();

    moptions.prevdemo = (prevdemo_t)cprevdemo->get_active_row_number ();
    moptions.serializeTiffRead = ctiffserialize->get_active();

    if (sdcurrent->get_active()) {
        moptions.startupDir = STARTUPDIR_CURRENT;
    } else if (sdhome->get_active()) {
        moptions.startupDir = STARTUPDIR_HOME;
    } else if (sdlast->get_active()) {
        moptions.startupDir = STARTUPDIR_LAST;
    } else if (sdother->get_active()) {
        moptions.startupDir = STARTUPDIR_CUSTOM;
        moptions.startupPath = startupdir->get_text();
    }

    moptions.parseExtensions.clear();
    moptions.parseExtensionsEnabled.clear();
    Gtk::TreeNodeChildren c = extensionModel->children();

    for (size_t i = 0; i < c.size(); i++) {
        moptions.parseExtensions.push_back(c[i][extensionColumns.ext]);
        moptions.parseExtensionsEnabled.push_back(c[i][extensionColumns.enabled]);
    }

    moptions.maxRecentFolders = (int)maxRecentFolders->get_value();
    moptions.maxThumbnailHeight = (int)maxThumbHeightSB->get_value ();
    moptions.maxCacheEntries = (int)maxCacheEntriesSB->get_value ();
    moptions.overlayedFileNames = overlayedFileNames->get_active();
    moptions.filmStripOverlayedFileNames = filmStripOverlayedFileNames->get_active();
    moptions.sameThumbSize = sameThumbSize->get_active();
    moptions.internalThumbIfUntouched = ckbInternalThumbIfUntouched->get_active();

    auto save_where = saveParamsPreference->get_active_row_number();
    moptions.saveParamsFile = save_where == 0 || save_where == 2;
    moptions.saveParamsCache = save_where == 1 || save_where == 2;
    moptions.paramsLoadLocation = (PPLoadLocation)loadParamsPreference->get_active_row_number();
    moptions.useBundledProfiles = useBundledProfiles->get_active();

    moptions.rtSettings.darkFramesPath = darkFrameDir->get_filename();
    moptions.rtSettings.flatFieldsPath = flatFieldDir->get_filename();

    moptions.clutsDir = clutsDir->get_filename();

    moptions.baBehav.resize(ADDSET_PARAM_NUM);

    for (Gtk::TreeIter sections = behModel->children().begin(); sections != behModel->children().end(); sections++)
        for (Gtk::TreeIter adjs = sections->children().begin(); adjs != sections->children().end(); adjs++) {
            moptions.baBehav[adjs->get_value(behavColumns.addsetid)] = adjs->get_value(behavColumns.badd);
        }

    int editorMode = editorLayout->get_active_row_number();
    moptions.tabbedUI = (editorMode > 1);
    moptions.multiDisplayMode = editorMode == 3 ? 1 : 0;
    moptions.mainNBVertical = editorMode == 1;

    moptions.curvebboxpos = curveBBoxPosC->get_active_row_number();
    moptions.complexity = complexitylocal->get_active_row_number();
    moptions.histogramPosition = ckbHistogramPositionLeft->get_active() ? 1 : 2;
    moptions.FileBrowserToolbarSingleRow = ckbFileBrowserToolbarSingleRow->get_active();
    moptions.showFilmStripToolBar = ckbShowFilmStripToolBar->get_active();
    moptions.hideTPVScrollbar = ckbHideTPVScrollbar->get_active();
    moptions.overwriteOutputFile = chOverwriteOutputFile->get_active();
    moptions.showtooltip = ckbshowtooltiplocallab->get_active();

    moptions.autoSaveTpOpen = ckbAutoSaveTpOpen->get_active();

    moptions.rgbDenoiseThreadLimit = threadsSpinBtn->get_value_as_int();
    moptions.clutCacheSize = clutCacheSizeSB->get_value_as_int();
    moptions.measure = measureCB->get_active();
    moptions.chunkSizeAMAZE = chunkSizeAMSB->get_value_as_int();
    moptions.chunkSizeCA = chunkSizeCASB->get_value_as_int();
    moptions.chunkSizeRCD = chunkSizeRCDSB->get_value_as_int();
    moptions.chunkSizeRGB = chunkSizeRGBSB->get_value_as_int();
    moptions.chunkSizeXT = chunkSizeXTSB->get_value_as_int();
    moptions.maxInspectorBuffers = maxInspectorBuffersSB->get_value_as_int();
    moptions.rtSettings.thumbnail_inspector_mode = static_cast<rtengine::Settings::ThumbnailInspectorMode>(thumbnailInspectorMode->get_active_row_number());

// Sounds only on Windows and Linux
#if defined(WIN32) || defined(__linux__)
    moptions.sndEnable = ckbSndEnable->get_active();
    moptions.sndBatchQueueDone = txtSndBatchQueueDone->get_text();
    moptions.sndLngEditProcDone = txtSndLngEditProcDone->get_text();
    moptions.sndLngEditProcDoneSecs = spbSndLngEditProcDoneSecs->get_value();
#endif

    moptions.cropGuides = Options::CropGuidesMode(cropGuidesCombo->get_active_row_number());
    moptions.cropAutoFit = cropAutoFitCB->get_active();
}

void Preferences::fillPreferences()
{

    tconn.block(true);
    fconn.block(true);
    cpfconn.block(true);
    sconn.block(true);
    dfconn.block(true);
    ffconn.block(true);
    rpconn.block(true);
    ipconn.block(true);
    bpconn.block(true);

    rprofiles->setActiveRowFromFullPath(moptions.defProfRaw);
    forRAWComboChanged(); // update the tooltip
    iprofiles->setActiveRowFromFullPath(moptions.defProfImg);
    forImageComboChanged(); // update the tooltip
    dateformat->set_text(moptions.dateFormat);
    panFactor->set_value(moptions.panAccelFactor);
    rememberZoomPanCheckbutton->set_active(moptions.rememberZoomAndPan);
    ctiffserialize->set_active(moptions.serializeTiffRead);

    setActiveTextOrIndex(*prtProfile, moptions.rtSettings.printerProfile, 0);

    switch (moptions.rtSettings.printerIntent) {
        default:
        case rtengine::RI_PERCEPTUAL:
            prtIntent->set_active(0);
            break;

        case rtengine::RI_RELATIVE:
            prtIntent->set_active(1);
            break;

        case rtengine::RI_ABSOLUTE:
            prtIntent->set_active(2);
            break;
    }

    prtBPC->set_active(moptions.rtSettings.printerBPC);

#if !defined(__APPLE__) // monitor profile not supported on apple
    setActiveTextOrIndex(*monProfile, moptions.rtSettings.monitorProfile, 0);

    switch (moptions.rtSettings.monitorIntent) {
        default:
        case rtengine::RI_PERCEPTUAL:
            monIntent->set_active(0);
            break;

        case rtengine::RI_RELATIVE:
            monIntent->set_active(1);
            break;

        case rtengine::RI_ABSOLUTE:
            monIntent->set_active(2);
            break;
    }

    monBPC->set_active(moptions.rtSettings.monitorBPC);
    cbAutoMonProfile->set_active(moptions.rtSettings.autoMonitorProfile);
#endif

    if (Glib::file_test(moptions.rtSettings.iccDirectory, Glib::FILE_TEST_IS_DIR)) {
        iccDir->set_current_folder(moptions.rtSettings.iccDirectory);
    }

    cprevdemo->set_active (moptions.prevdemo);
    languages->set_active_text(moptions.language);
    ckbLangAutoDetect->set_active(moptions.languageAutoDetect);
    int themeNbr = getThemeRowNumber(moptions.theme);
    themeCBT->set_active (themeNbr == -1 ? 0 : themeNbr);

    Gdk::RGBA cropCol;
    cropCol.set_rgba(moptions.cutOverlayBrush[0], moptions.cutOverlayBrush[1], moptions.cutOverlayBrush[2]);
    cropMaskColorCB->set_rgba (cropCol);
    cropMaskColorCB->set_alpha ( (unsigned short) (moptions.cutOverlayBrush[3] * 65535.0));

    Gdk::RGBA NavGuideCol;
    NavGuideCol.set_rgba(moptions.navGuideBrush[0], moptions.navGuideBrush[1], moptions.navGuideBrush[2]);
    navGuideColorCB->set_rgba (NavGuideCol);
    navGuideColorCB->set_alpha ( (unsigned short) (moptions.navGuideBrush[3] * 65535.0));

    if (options.fontFamily == "default") {
        mainFontFB->set_font_name (Glib::ustring::compose ("%1 %2", initialFontFamily, initialFontSize));
    } else {
        mainFontFB->set_font_name (Glib::ustring::compose ("%1 %2", options.fontFamily, options.fontSize));
    }

    if (options.CPFontFamily == "default") {
        colorPickerFontFB->set_font_name (Glib::ustring::compose ("%1 %2", initialFontFamily, initialFontSize));
    } else {
        colorPickerFontFB->set_font_name (Glib::ustring::compose ("%1 %2", options.CPFontFamily, options.CPFontSize));
    }

    pseudoHiDPI->set_active(options.pseudoHiDPISupport);

    showDateTime->set_active(moptions.fbShowDateTime);
    showBasicExif->set_active(moptions.fbShowBasicExif);
    showExpComp->set_active(moptions.fbShowExpComp);
    ckbmenuGroupRank->set_active(moptions.menuGroupRank);
    ckbmenuGroupLabel->set_active(moptions.menuGroupLabel);
    ckbmenuGroupFileOperations->set_active(moptions.menuGroupFileOperations);
    ckbmenuGroupProfileOperations->set_active(moptions.menuGroupProfileOperations);
    ckbmenuGroupExtProg->set_active(moptions.menuGroupExtProg);

    hlThresh->set_value(moptions.highlightThreshold);
    shThresh->set_value(moptions.shadowThreshold);

    edGimp->set_active(moptions.editorToSendTo == 1);
    edOther->set_active(moptions.editorToSendTo == 3);
#ifdef WIN32
    edPS->set_active(moptions.editorToSendTo == 2);

    if (Glib::file_test(moptions.gimpDir, Glib::FILE_TEST_IS_DIR)) {
        gimpDir->set_current_folder(moptions.gimpDir);
    } else {
        gimpDir->set_current_folder(Glib::get_home_dir());
    }

    if (Glib::file_test(moptions.psDir, Glib::FILE_TEST_IS_DIR)) {
        psDir->set_current_folder(moptions.psDir);
    } else {
        psDir->set_current_folder(Glib::get_home_dir());
    }

#elif defined __APPLE__
    edPS->set_active(moptions.editorToSendTo == 2);

    if (Glib::file_test(moptions.psDir, Glib::FILE_TEST_IS_DIR)) {
        psDir->set_current_folder(moptions.psDir);
    } else {
        psDir->set_current_folder(Glib::get_home_dir());
    }

#endif
    editorToSendTo->set_text(moptions.customEditorProg);

    txtCustProfBuilderPath->set_text(moptions.CPBPath);
    custProfBuilderLabelType->set_active(moptions.CPBKeys);


    if (moptions.startupDir == STARTUPDIR_CURRENT) {
        sdcurrent->set_active();
    } else if (moptions.startupDir == STARTUPDIR_LAST) {
        sdlast->set_active();
    } else if (moptions.startupDir == STARTUPDIR_HOME) {
        sdhome->set_active();
    } else if (moptions.startupDir == STARTUPDIR_CUSTOM) {
        sdother->set_active();
        startupdir->set_text(moptions.startupPath);
    }

    extensionModel->clear();

    for (size_t i = 0; i < moptions.parseExtensions.size(); i++) {
        Gtk::TreeRow row = * (extensionModel->append());
        row[extensionColumns.enabled] = moptions.parseExtensionsEnabled[i];
        row[extensionColumns.ext]     = moptions.parseExtensions[i];
    }

    maxRecentFolders->set_value(moptions.maxRecentFolders);
    maxThumbHeightSB->set_value (moptions.maxThumbnailHeight);
    maxCacheEntriesSB->set_value (moptions.maxCacheEntries);
    overlayedFileNames->set_active(moptions.overlayedFileNames);
    filmStripOverlayedFileNames->set_active(moptions.filmStripOverlayedFileNames);
    sameThumbSize->set_active(moptions.sameThumbSize);
    ckbInternalThumbIfUntouched->set_active(moptions.internalThumbIfUntouched);

    saveParamsPreference->set_active(moptions.saveParamsFile ? (moptions.saveParamsCache ? 2 : 0) : 1);

    loadParamsPreference->set_active(moptions.paramsLoadLocation);
    useBundledProfiles->set_active(moptions.useBundledProfiles);

    if (!moptions.tabbedUI) {
        editorLayout->set_active(moptions.mainNBVertical ? 1 : 0);
    } else {
        editorLayout->set_active(moptions.multiDisplayMode ? 3 : 2);
    }

    curveBBoxPosC->set_active(moptions.curvebboxpos);
    complexitylocal->set_active(moptions.complexity); 

    ckbHistogramPositionLeft->set_active(moptions.histogramPosition == 1);
    ckbFileBrowserToolbarSingleRow->set_active(moptions.FileBrowserToolbarSingleRow);
    ckbShowFilmStripToolBar->set_active(moptions.showFilmStripToolBar);
    ckbHideTPVScrollbar->set_active(moptions.hideTPVScrollbar);
    ckbshowtooltiplocallab->set_active(moptions.showtooltip);
    ckbAutoSaveTpOpen->set_active(moptions.autoSaveTpOpen);

    threadsSpinBtn->set_value (moptions.rgbDenoiseThreadLimit);
    clutCacheSizeSB->set_value (moptions.clutCacheSize);
    measureCB->set_active (moptions.measure);
    chunkSizeAMSB->set_value (moptions.chunkSizeAMAZE);
    chunkSizeCASB->set_value (moptions.chunkSizeCA);
    chunkSizeRGBSB->set_value (moptions.chunkSizeRGB);
    chunkSizeRCDSB->set_value (moptions.chunkSizeRCD);
    chunkSizeXTSB->set_value (moptions.chunkSizeXT);
    maxInspectorBuffersSB->set_value (moptions.maxInspectorBuffers);
    thumbnailInspectorMode->set_active(int(moptions.rtSettings.thumbnail_inspector_mode));

    darkFrameDir->set_current_folder(moptions.rtSettings.darkFramesPath);
    darkFrameChanged();

    flatFieldDir->set_current_folder(moptions.rtSettings.flatFieldsPath);
    flatFieldChanged();

    clutsDir->set_current_folder(moptions.clutsDir);

    addc.block(true);
    setc.block(true);

    moptions.baBehav.resize(ADDSET_PARAM_NUM);

    for (Gtk::TreeIter sections = behModel->children().begin(); sections != behModel->children().end(); ++sections) {
        for (Gtk::TreeIter adjs = sections->children().begin(); adjs != sections->children().end(); ++adjs) {
            const bool add = moptions.baBehav[adjs->get_value(behavColumns.addsetid)];
            adjs->set_value(behavColumns.badd, add);
            adjs->set_value(behavColumns.bset, !add);
        }
    }

    cropGuidesCombo->set_active(moptions.cropGuides);
    cropAutoFitCB->set_active(moptions.cropAutoFit);

    addc.block(false);
    setc.block(false);
    cpfconn.block(false);
    fconn.block(false);
    tconn.block(false);
    sconn.block(false);
    dfconn.block(false);
    ffconn.block(false);
    rpconn.block(true);
    ipconn.block(true);
    bpconn.block(false);

    chOverwriteOutputFile->set_active(moptions.overwriteOutputFile);

    // Sounds only on Windows and Linux
#if defined(WIN32) || defined(__linux__)
    ckbSndEnable->set_active(moptions.sndEnable);
    txtSndBatchQueueDone->set_text(moptions.sndBatchQueueDone);
    txtSndLngEditProcDone->set_text(moptions.sndLngEditProcDone);
    spbSndLngEditProcDoneSecs->set_value(moptions.sndLngEditProcDoneSecs);
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

void Preferences::autoMonProfileToggled()
{
    monProfile->set_sensitive(!cbAutoMonProfile->get_active());
}

/*
void Preferences::autocielabToggled () {
//  cbAutocielab->set_sensitive(cbAutocielab->get_active());
}
*/

void Preferences::sndEnableToggled()
{
    txtSndBatchQueueDone->set_sensitive(ckbSndEnable->get_active());
    txtSndLngEditProcDone->set_sensitive(ckbSndEnable->get_active());
    spbSndLngEditProcDoneSecs->set_sensitive(ckbSndEnable->get_active());
}

void Preferences::langAutoDetectToggled()
{
    languages->set_sensitive(!ckbLangAutoDetect->get_active());
}

void Preferences::okPressed()
{

    storePreferences();
    workflowUpdate();
    options.copyFrom(&moptions);
    options.filterOutParsedExtensions();

    try {
        Options::save();
    } catch (Options::Error &e) {
        Gtk::MessageDialog msgd(getToplevelWindow(this), e.get_msg(), true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_CLOSE, true);
        msgd.run();
    }

    dynProfilePanel->save();
    hide();
}

void Preferences::cancelPressed()
{
    // set the initial theme back
    if (themeFNames.at (themeCBT->get_active_row_number ()).longFName != options.theme) {
        RTImage::updateImages();
        switchThemeTo(options.theme);
    }

    // set the initial font back
    Pango::FontDescription fd (mainFontFB->get_font_name());

    if (fd.get_family() != options.fontFamily && (fd.get_size() / Pango::SCALE) != options.fontSize) {
        if (options.fontFamily == "default") {
            switchFontTo(initialFontFamily, initialFontSize);
        } else {
            switchFontTo(options.fontFamily, options.fontSize);
        }
    }

    // update the profileStore
    if (useBundledProfiles->get_active() != options.useBundledProfiles) {
        // we have to rescan with the old value
        bpconn.block(true);
        useBundledProfiles->set_active(false);
        bundledProfilesChanged();
        bpconn.block(false);
    }

    hide();
}

void Preferences::selectStartupDir()
{

    Gtk::FileChooserDialog dialog(getToplevelWindow(this), M("PREFERENCES_DIRSELECTDLG"), Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER);
    //dialog.set_transient_for(*this);

    //Add response buttons to the dialog:
    dialog.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(M("GENERAL_OPEN"), Gtk::RESPONSE_OK);

    int result = dialog.run();

    if (result == Gtk::RESPONSE_OK) {
        startupdir->set_text(dialog.get_filename());
    }
}

void Preferences::aboutPressed()
{

    splash = new Splash(*this);
    splash->set_transient_for(*this);
    splash->signal_delete_event().connect(sigc::mem_fun(*this, &Preferences::splashClosed));
    splash->show();
}

void Preferences::themeChanged()
{

    moptions.theme = themeFNames.at (themeCBT->get_active_row_number ()).longFName;
    RTImage::updateImages();
    switchThemeTo(moptions.theme);
}

void Preferences::forRAWComboChanged()
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

void Preferences::forImageComboChanged()
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

void Preferences::layoutComboChanged()
{
    editorLayout->set_tooltip_text(editorLayout->get_active_text());
}

void Preferences::bundledProfilesChanged()
{
    rpconn.block(true);
    ipconn.block(true);

    // parseProfiles does use options.useBundledProfiles, so we temporarily change its value
    bool currValue = options.useBundledProfiles;
    options.useBundledProfiles = useBundledProfiles->get_active();

    // rescan the file's tree
    ProfileStore::getInstance()->parseProfiles(); // This will call Preferences::updateProfileList in return

    // restoring back the old value
    options.useBundledProfiles = currValue;

    ipconn.block(false);
    rpconn.block(false);
}

void Preferences::iccDirChanged()
{
    const auto currentSelection = monProfile->get_active_text();
    const auto profiles = rtengine::ICCStore::getInstance()->getProfilesFromDir(iccDir->get_filename());

    monProfile->remove_all();

    monProfile->append(M("PREFERENCES_PROFILE_NONE"));

    for (const auto& profile : profiles) {
        monProfile->append(profile);
    }

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
    const ProfileStoreEntry* dynpse = ProfileStore::getInstance()->getInternalDynamicPSE();
    rprofiles->addRow(dynpse);
    iprofiles->addRow(dynpse);
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
        themecss->load_from_path(filename);
    } catch (Glib::Error &err) {
        printf("Error: Can't load css file \"%s\"\nMessage: %s\n", filename.c_str(), err.what().c_str());
    } catch (...) {
        printf("Error: Can't load css file \"%s\"\n", filename.c_str());
    }
}

void Preferences::fontChanged()
{
    newFont = true;
    Pango::FontDescription fd (mainFontFB->get_font_name());
    switchFontTo(fd.get_family(), fd.get_size() / Pango::SCALE);
}

void Preferences::cpFontChanged()
{

    newCPFont = true;
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
//#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
//            fontcss->load_from_data (Glib::ustring::compose ("* { font-family: %1; font-size: %2px }", newFontFamily, newFontSize * RTScalable::getScale()));
//#else
            fontcss->load_from_data (Glib::ustring::compose ("* { font-family: %1; font-size: %2pt }", newFontFamily, newFontSize * RTScalable::getScale()));
//#endif
            //GTK318
        } catch (Glib::Error &err) {
            printf("Error: \"%s\"\n", err.what().c_str());
        } catch (...) {
            printf("Error: Can't find the font named \"%s\"\n", newFontFamily.c_str());
        }
    } else {
        if (fontcss) {
            fontcss = Gtk::CssProvider::create();
            Glib::RefPtr<Gdk::Screen> screen = Gdk::Screen::get_default();
            Gtk::StyleContext::remove_provider_for_screen(screen, fontcss);
        }
    }
}

void Preferences::workflowUpdate()
{

    if (moptions.tabbedUI != options.tabbedUI) {
        parent->setEditorMode(moptions.tabbedUI);
    }

    if (moptions.hideTPVScrollbar != options.hideTPVScrollbar) {
        // Update the tool panels
        parent->updateTPVScrollbar(moptions.hideTPVScrollbar);
    }

    if (moptions.FileBrowserToolbarSingleRow != options.FileBrowserToolbarSingleRow) {
        // Update the position of the Query toolbar
        parent->updateFBQueryTB(moptions.FileBrowserToolbarSingleRow);
    }

    if (moptions.showFilmStripToolBar != options.showFilmStripToolBar) {
        // Update the visibility of FB toolbar
        parent->updateFBToolBarVisibility(moptions.showFilmStripToolBar);
    }

    if (moptions.showtooltip != options.showtooltip) {
        // Update the visibility of tooltip
        parent->updateShowtooltipVisibility(moptions.showtooltip);
    }

    if (moptions.histogramPosition != options.histogramPosition) {
        // Update the position of the Histogram
        parent->updateHistogramPosition(options.histogramPosition, moptions.histogramPosition);
    }

    if (moptions.rtSettings.printerProfile != options.rtSettings.printerProfile
            || moptions.rtSettings.printerBPC     != options.rtSettings.printerBPC
            || moptions.rtSettings.printerIntent  != options.rtSettings.printerIntent) {
        // Update the position of the Histogram
        parent->updateProfiles (moptions.rtSettings.printerProfile, rtengine::RenderingIntent(moptions.rtSettings.printerIntent), moptions.rtSettings.printerBPC);
    }

}

void Preferences::addExtPressed()
{

    Gtk::TreeNodeChildren c = extensionModel->children();

    for (size_t i = 0; i < c.size(); i++)
        if (c[i][extensionColumns.ext] == extension->get_text()) {
            return;
        }

    Gtk::TreeRow row = * (extensionModel->append());

    row[extensionColumns.enabled] = true;
    row[extensionColumns.ext]     = extension->get_text();
}

void Preferences::delExtPressed()
{

    extensionModel->erase(extensions->get_selection()->get_selected());
}

void Preferences::moveExtUpPressed()
{
    const Glib::RefPtr<Gtk::TreeSelection> selection = extensions->get_selection();

    if (!selection) {
        return;
    }

    const Gtk::TreeModel::iterator selected = selection->get_selected();

    if (!selected || selected == extensionModel->children().begin()) {
        return;
    }

    Gtk::TreeModel::iterator previous = selected;
    --previous;
    extensionModel->iter_swap(selected, previous);
}

void Preferences::moveExtDownPressed()
{
    const Glib::RefPtr<Gtk::TreeSelection> selection = extensions->get_selection();

    if (!selection) {
        return;
    }

    const Gtk::TreeModel::iterator selected = selection->get_selected();

    if (!selected) {
        return;
    }

    Gtk::TreeModel::iterator next = selected;

    if (++next) {
        extensionModel->iter_swap(selected, next);
    }
}

void Preferences::clearProfilesPressed()
{

    cacheMgr->clearProfiles();
}


void Preferences::clearThumbImagesPressed()
{

    cacheMgr->clearImages();
}

void Preferences::clearAllPressed()
{

    cacheMgr->clearAll();
}

void Preferences::darkFrameChanged()
{
    //Glib::ustring s(darkFrameDir->get_filename());
    Glib::ustring s(darkFrameDir->get_current_folder());
    //if( s.compare( rtengine::dfm.getPathname()) !=0 ){
    rtengine::dfm.init(s);
    updateDFinfos();
    //}
}

void Preferences::flatFieldChanged()
{
    //Glib::ustring s(flatFieldDir->get_filename());
    Glib::ustring s(flatFieldDir->get_current_folder());
    //if( s.compare( rtengine::ffm.getPathname()) !=0 ){
    rtengine::ffm.init(s);
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

void Preferences::behAddSetAllPressed(bool add)
{
    moptions.baBehav.assign(ADDSET_PARAM_NUM, add);

    for (Gtk::TreeIter sections = behModel->children().begin(); sections != behModel->children().end(); ++sections) {
        for (Gtk::TreeIter adjs = sections->children().begin(); adjs != sections->children().end(); ++adjs) {
            adjs->set_value(behavColumns.badd, add);
            adjs->set_value(behavColumns.bset, !add);
        }
    }
}

void Preferences::behAddAllPressed()
{
    behAddSetAllPressed(true);
}

void Preferences::behSetAllPressed()
{
    behAddSetAllPressed(false);
}
