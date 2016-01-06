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
#include "partialpastedlg.h"
#include "multilangmgr.h"
#include "paramsedited.h"

PartialPasteDlg::PartialPasteDlg (Glib::ustring title)
{

    set_modal (true);
    set_title (title);
    set_default_size(700, 600);

    everything  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EVERYTHING")));
    everything  ->set_name("partialPasteHeader");

    basic       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_BASICGROUP")));
    basic       ->set_name("partialPasteHeader");
    detail      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DETAILGROUP")));
    detail      ->set_name("partialPasteHeader");
    color       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORGROUP")));
    color       ->set_name("partialPasteHeader");
    lens        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LENSGROUP")));
    lens        ->set_name("partialPasteHeader");
    composition = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMPOSITIONGROUP")));
    composition ->set_name("partialPasteHeader");
    meta     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_METAGROUP")));
    meta     ->set_name("partialPasteHeader");
    raw         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWGROUP")));
    raw         ->set_name("partialPasteHeader");
    wav         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WAVELETGROUP")));
    wav         ->set_name("partialPasteHeader");

    // options in basic:
    wb          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WHITEBALANCE")));
    exposure    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXPOSURE")));
    sh          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHADOWSHIGHLIGHTS")));
    epd         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EPD")));
    retinex       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RETINEX")));
    pcvignette  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PCVIGNETTE")));
    gradient    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_GRADIENT")));
    labcurve    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LABCURVE")));
    colorappearance = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORAPP")));

    // options in detail:
    sharpen     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENING")));
    sharpenedge = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENEDGE")));
    sharpenmicro = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENMICRO")));
    impden      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IMPULSEDENOISE")));
    dirpyreq    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYREQUALIZER")));
    defringe    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DEFRINGE")));

    // options in wavelet:
    wavelet    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EQUALIZER"))); //TODO - rename to wavelet

    // options in color:
    icm         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ICMSETTINGS")));
    //gam         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ICMGAMMA")));
    vibrance    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_VIBRANCE")));
    chmixer     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CHANNELMIXER")));
    blackwhite   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CHANNELMIXERBW")));
    dirpyrden   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYRDENOISE")));
    hsveq       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_HSVEQUALIZER")));
    filmSimulation = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FILMSIMULATION")) );
    rgbcurves   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RGBCURVES")));
    colortoning = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORTONING")));

    // options in lens:
    distortion  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DISTORTION")));
    cacorr      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CACORRECTION")));
    vignetting  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_VIGNETTING")));
    lcp         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LENSPROFILE")));

    // options in composition:
    coarserot   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COARSETRANS")));
    finerot     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ROTATION")));
    crop        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CROP")));
    resize      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RESIZE")));
    prsharpening     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PRSHARPENING")));
    perspective = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PERSPECTIVE")));
    commonTrans = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMMONTRANSFORMPARAMS")));

    // options in meta:
    exifch      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXIFCHANGES")));
    iptc        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IPTCINFO")));

    // options in raw:
    raw_expos           = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_LINEAR")));
    raw_preser          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_PRESER")));
    raw_black           = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_BLACK")));
    raw_ca_autocorrect  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_AUTO")));
    raw_cared           = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_CARED")));
    raw_cablue          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_CABLUE")));
    raw_hotpix_filt     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_HOTPIXFILT")));
    raw_deadpix_filt    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_DEADPIXFILT")));
    raw_linenoise       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_LINEDENOISE")));
    raw_greenthresh     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_GREENEQUIL")));
    raw_method          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DMETHOD")));
    raw_ccSteps         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_FALSECOLOR")));
    raw_dcb_iterations  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DCBITERATIONS")));
    raw_dcb_enhance     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DCBENHANCE")));
    //raw_all_enhance   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_ALLENHANCE")));
    raw_lmmse_iterations = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_LMMSEITERATIONS")));

    df_file             = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEFILE")));
    df_AutoSelect       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEAUTOSELECT")));
    ff_file             = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDFILE")));
    ff_AutoSelect       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDAUTOSELECT")));
    ff_BlurRadius       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURRADIUS")));
    ff_BlurType         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURTYPE")));
    ff_ClipControl      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDCLIPCONTROL")));

    Gtk::VBox* vboxes[8];
    Gtk::HSeparator* hseps[8];

    for (int i = 0; i < 8; i++) {
        vboxes[i] = Gtk::manage (new Gtk::VBox ());
        vboxes[i]->set_border_width (6);
        hseps[i] = Gtk::manage (new Gtk::HSeparator ());
        hseps[i]->set_name("partialPasteHeaderSep");
    }

    //BASIC
    vboxes[0]->pack_start (*basic, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*hseps[0], Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*wb, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*exposure, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*sh, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*epd, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*retinex, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*pcvignette, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*gradient, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*labcurve, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*colorappearance, Gtk::PACK_SHRINK, 2);

    //DETAIL
    vboxes[1]->pack_start (*detail, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*hseps[1], Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpen, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpenedge, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpenmicro, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*impden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*dirpyrden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*defringe, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*dirpyreq, Gtk::PACK_SHRINK, 2);

    //COLOR
    vboxes[2]->pack_start (*color, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hseps[2], Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*icm, Gtk::PACK_SHRINK, 2);
    //vboxes[2]->pack_start (*gam, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*vibrance, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*chmixer, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*blackwhite, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hsveq, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*filmSimulation, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*rgbcurves, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colortoning, Gtk::PACK_SHRINK, 2);

    //LENS
    vboxes[3]->pack_start (*lens, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*hseps[3], Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*distortion, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*cacorr, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*vignetting, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*lcp, Gtk::PACK_SHRINK, 2);

    //COMPOSITION
    vboxes[4]->pack_start (*composition, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*hseps[4], Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*coarserot, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*finerot, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*crop, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*resize, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*prsharpening, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*perspective, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*commonTrans, Gtk::PACK_SHRINK, 2);

    //WAVELET
    vboxes[5]->pack_start (*wav, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*hseps[5], Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*wavelet, Gtk::PACK_SHRINK, 2);

    //RAW
    vboxes[6]->pack_start (*raw, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*hseps[6], Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_method, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_ccSteps, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_dcb_iterations, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_dcb_enhance, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_lmmse_iterations, Gtk::PACK_SHRINK, 2);
    //vboxes[6]->pack_start (*raw_all_enhance, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    vboxes[6]->pack_start (*raw_linenoise, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_greenthresh, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_hotpix_filt, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_deadpix_filt, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    vboxes[6]->pack_start (*raw_expos, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_preser, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_black, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    vboxes[6]->pack_start (*df_file, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*df_AutoSelect, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    vboxes[6]->pack_start (*ff_file, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*ff_AutoSelect, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*ff_BlurType, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*ff_BlurRadius, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*ff_ClipControl, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
    vboxes[6]->pack_start (*raw_ca_autocorrect, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_cared, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*raw_cablue, Gtk::PACK_SHRINK, 2);

    //META
    vboxes[7]->pack_start (*meta, Gtk::PACK_SHRINK, 2);
    vboxes[7]->pack_start (*hseps[7], Gtk::PACK_SHRINK, 2);
    vboxes[7]->pack_start (*exifch, Gtk::PACK_SHRINK, 2);
    vboxes[7]->pack_start (*iptc, Gtk::PACK_SHRINK, 2);

    Gtk::VBox* vbCol1 = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbCol2 = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbCol3 = Gtk::manage (new Gtk::VBox ());

    vbCol1->set_border_width (8);
    vbCol2->set_border_width (8);
    vbCol3->set_border_width (8);

    for (int i = 0; i < 3; i++) {
        vbCol1->pack_start (*vboxes[i], Gtk::PACK_SHRINK, 2);
    }

    for (int i = 3; i < 6; i++) {
        vbCol2->pack_start (*vboxes[i], Gtk::PACK_SHRINK, 2);
    }

    for (int i = 6; i < 8; i++) {
        vbCol3->pack_start (*vboxes[i], Gtk::PACK_SHRINK, 2);
    }

    Gtk::VBox* vbtop = Gtk::manage (new Gtk::VBox ());
    vbtop->pack_start (*everything, Gtk::PACK_SHRINK, 2);
    vbtop->set_border_width (8);

    Gtk::Dialog::get_content_area()->pack_start (*vbtop, Gtk::PACK_SHRINK, 2); // TODO replace with get_content_area() with GTK upgrade

    Gtk::HBox* hbmain = Gtk::manage (new Gtk::HBox ());
    hbmain->pack_start (*vbCol1);
    hbmain->pack_start (*(Gtk::manage (new Gtk::VSeparator ())));
    hbmain->pack_start (*vbCol2);
    hbmain->pack_start (*(Gtk::manage (new Gtk::VSeparator ())));
    hbmain->pack_start (*vbCol3);

    scrolledwindow = Gtk::manage ( new Gtk::ScrolledWindow() );
    scrolledwindow->set_can_focus(true);
    scrolledwindow->set_border_width(2);
    scrolledwindow->set_shadow_type(Gtk::SHADOW_NONE);
    scrolledwindow->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
    scrolledwindow->property_window_placement().set_value(Gtk::CORNER_TOP_LEFT);

    scrolledwindow->add(*hbmain);

    Gtk::Dialog::get_content_area()->pack_start (*scrolledwindow, Gtk::PACK_EXPAND_WIDGET, 2);// TODO replace with get_content_area() with GTK upgrade

    hbmain->show();
    scrolledwindow->show ();

    // This can be improved
    // there is currently no binding of subsettings to CheckButton 'everything' for its inconsistent status
    everythingConn  = everything->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::everythingToggled));
    basicConn       = basic->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::basicToggled));
    detailConn      = detail->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::detailToggled));
    colorConn       = color->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::colorToggled));
    lensConn        = lens->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::lensToggled));
    compositionConn = composition->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::compositionToggled));
    metaConn     = meta->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::metaToggled));
    rawConn         = raw->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::rawToggled));
    wavConn         = wav->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::wavToggled));

    wbConn          = wb->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    exposureConn    = exposure->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    shConn          = sh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    epdConn         = epd->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    retinexConn       = retinex->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    pcvignetteConn  = pcvignette->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    gradientConn    = gradient->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    labcurveConn    = labcurve->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    colorappearanceConn = colorappearance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));

    sharpenConn     = sharpen->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    gradsharpenConn = sharpenedge->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    microcontrastConn = sharpenmicro->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    impdenConn      = impden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    dirpyrdenConn   = dirpyrden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    dirpyreqConn    = dirpyreq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    defringeConn    = defringe->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));

    waveletConn = wavelet->signal_toggled().connect (sigc::bind (sigc::mem_fun(*wav, &Gtk::CheckButton::set_inconsistent), true));

    icmConn         = icm->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    //gamcsconn      = gam->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    vibranceConn    = vibrance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    chmixerConn     = chmixer->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    chmixerbwConn   = blackwhite->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    hsveqConn       = hsveq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    filmSimulationConn = filmSimulation->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    rgbcurvesConn   = rgbcurves->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    colortoningConn = colortoning->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));

    distortionConn  = distortion->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));
    cacorrConn      = cacorr->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));
    vignettingConn  = vignetting->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));
    lcpConn         = lcp->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));

    coarserotConn   = coarserot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    finerotConn     = finerot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    cropConn        = crop->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    resizeConn      = resize->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    prsharpeningConn      = prsharpening->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    perspectiveConn = perspective->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    commonTransConn = commonTrans->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));

    exifchConn      = exifch->signal_toggled().connect (sigc::bind (sigc::mem_fun(*meta, &Gtk::CheckButton::set_inconsistent), true));
    iptcConn        = iptc->signal_toggled().connect (sigc::bind (sigc::mem_fun(*meta, &Gtk::CheckButton::set_inconsistent), true));

    raw_methodConn         = raw_method->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_ccStepsConn         = raw_ccSteps->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_dcb_iterationsConn  = raw_dcb_iterations->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_dcb_enhanceConn     = raw_dcb_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    //raw_all_enhanceConn     = raw_all_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_lmmse_iterationsConn  = raw_lmmse_iterations->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));

    raw_exposConn           = raw_expos->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_preserConn          = raw_preser->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_blackConn           = raw_black->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_ca_autocorrectConn  = raw_ca_autocorrect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_caredConn           = raw_cared->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_cablueConn          = raw_cablue->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_hotpix_filtConn     = raw_hotpix_filt->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_deadpix_filtConn    = raw_deadpix_filt->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_linenoiseConn       = raw_linenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_greenthreshConn     = raw_greenthresh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    df_fileConn             = df_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    df_AutoSelectConn       = df_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_fileConn             = ff_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_AutoSelectConn       = ff_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_BlurRadiusConn       = ff_BlurRadius->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_BlurTypeConn         = ff_BlurType->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_ClipControlConn      = ff_ClipControl->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));

    add_button (M("GENERAL_OK"), Gtk::RESPONSE_OK);
    add_button (M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    set_response_sensitive (Gtk::RESPONSE_OK);
    set_default_response (Gtk::RESPONSE_OK);
    show_all_children ();
}

void PartialPasteDlg::everythingToggled ()
{

    basicConn.block (true);
    detailConn.block (true);
    colorConn.block (true);
    lensConn.block (true);
    compositionConn.block (true);
    metaConn.block (true);
    rawConn.block (true);
    wavConn.block (true);

    everything->set_inconsistent (false);

    //toggle group headings
    basic->set_active(everything->get_active());
    detail->set_active(everything->get_active());
    color->set_active(everything->get_active());
    lens->set_active(everything->get_active());
    composition->set_active(everything->get_active());
    meta->set_active(everything->get_active());
    raw->set_active(everything->get_active());
    wav->set_active(everything->get_active());

    //toggle group children
    PartialPasteDlg::basicToggled ();
    PartialPasteDlg::detailToggled ();
    PartialPasteDlg::colorToggled ();
    PartialPasteDlg::lensToggled ();
    PartialPasteDlg::compositionToggled ();
    PartialPasteDlg::metaToggled ();
    PartialPasteDlg::rawToggled ();
    PartialPasteDlg::wavToggled ();

    basicConn.block (false);
    detailConn.block (false);
    colorConn.block (false);
    lensConn.block (false);
    compositionConn.block (false);
    metaConn.block (false);
    rawConn.block (false);
    wavConn.block (false);
}

void PartialPasteDlg::rawToggled ()
{

    raw_methodConn.block (true);
    raw_ccStepsConn.block (true);
    raw_dcb_iterationsConn.block (true);
    raw_dcb_enhanceConn.block (true);
    //raw_all_enhanceConn.block (true);
    raw_lmmse_iterationsConn.block (true);
    raw_exposConn.block (true);
    raw_preserConn.block (true);
    raw_blackConn.block (true);
    raw_ca_autocorrectConn.block (true);
    raw_caredConn.block (true);
    raw_cablueConn.block (true);
    raw_hotpix_filtConn.block (true);
    raw_deadpix_filtConn.block (true);
    raw_linenoiseConn.block (true);
    raw_greenthreshConn.block (true);
    df_fileConn.block (true);
    df_AutoSelectConn.block (true);
    ff_fileConn.block (true);
    ff_AutoSelectConn.block (true);
    ff_BlurRadiusConn.block (true);
    ff_BlurTypeConn.block (true);
    ff_ClipControlConn.block (true);

    raw->set_inconsistent (false);

    raw_method->set_active (raw->get_active ());
    raw_ccSteps->set_active (raw->get_active ());
    raw_dcb_iterations->set_active (raw->get_active ());
    raw_dcb_enhance->set_active (raw->get_active ());
    raw_lmmse_iterations->set_active (raw->get_active ());
    //raw_all_enhance->set_active (raw->get_active ());
    raw_expos->set_active (raw->get_active ());
    raw_preser->set_active (raw->get_active ());
    raw_black->set_active (raw->get_active ());
    raw_ca_autocorrect->set_active (raw->get_active ());
    raw_cared->set_active (raw->get_active ());
    raw_cablue->set_active (raw->get_active ());
    raw_hotpix_filt->set_active (raw->get_active ());
    raw_deadpix_filt->set_active (raw->get_active ());
    raw_linenoise->set_active (raw->get_active ());
    raw_greenthresh->set_active (raw->get_active ());
    df_file->set_active (raw->get_active ());
    df_AutoSelect->set_active (raw->get_active ());
    ff_file->set_active (raw->get_active ());
    ff_AutoSelect->set_active (raw->get_active ());
    ff_BlurRadius->set_active (raw->get_active ());
    ff_BlurType->set_active (raw->get_active ());
    ff_ClipControl->set_active (raw->get_active ());

    raw_methodConn.block (false);
    raw_ccStepsConn.block (false);
    raw_dcb_iterationsConn.block (false);
    raw_dcb_enhanceConn.block (false);
    //raw_all_enhanceConn.block (false);
    raw_lmmse_iterationsConn.block (false);
    raw_exposConn.block (false);
    raw_preserConn.block (false);
    raw_blackConn.block (false);
    raw_ca_autocorrectConn.block (false);
    raw_caredConn.block (false);
    raw_cablueConn.block (false);
    raw_hotpix_filtConn.block (false);
    raw_deadpix_filtConn.block (false);
    raw_linenoiseConn.block (false);
    raw_greenthreshConn.block (false);
    df_fileConn.block (false);
    df_AutoSelectConn.block (false);
    ff_fileConn.block (false);
    ff_AutoSelectConn.block (false);
    ff_BlurRadiusConn.block (false);
    ff_BlurTypeConn.block (false);
    ff_ClipControlConn.block (false);
}

void PartialPasteDlg::basicToggled ()
{

    wbConn.block (true);
    exposureConn.block (true);
    shConn.block (true);
    epdConn.block(true);
    pcvignetteConn.block (true);
    gradientConn.block (true);
    labcurveConn.block (true);
    colorappearanceConn.block (true);

    basic->set_inconsistent (false);

    wb->set_active (basic->get_active ());
    exposure->set_active (basic->get_active ());
    sh->set_active (basic->get_active ());
    epd->set_active (basic->get_active ());
    pcvignette->set_active (basic->get_active ());
    gradient->set_active (basic->get_active ());
    labcurve->set_active (basic->get_active ());
    colorappearance->set_active (basic->get_active ());

    wbConn.block (false);
    exposureConn.block (false);
    shConn.block (false);
    epdConn.block (false);
    pcvignetteConn.block (false);
    gradientConn.block (false);
    labcurveConn.block (false);
    colorappearanceConn.block (false);
}

void PartialPasteDlg::detailToggled ()
{

    sharpenConn.block (true);
    gradsharpenConn.block(true);
    microcontrastConn.block(true);
    impdenConn.block (true);
    dirpyrdenConn.block (true);
    defringeConn.block (true);
    dirpyreqConn.block (true);

    detail->set_inconsistent (false);

    sharpen->set_active (detail->get_active ());
    sharpenedge->set_active (detail->get_active ());
    sharpenmicro->set_active (detail->get_active ());
    impden->set_active (detail->get_active ());
    dirpyrden->set_active (detail->get_active ());
    defringe->set_active (detail->get_active ());
    dirpyreq->set_active (detail->get_active ());

    sharpenConn.block (false);
    gradsharpenConn.block(false);
    microcontrastConn.block(false);
    impdenConn.block (false);
    dirpyrdenConn.block (false);
    defringeConn.block (false);
    dirpyreqConn.block (false);
}

void PartialPasteDlg::wavToggled ()
{

    waveletConn.block (true);

    wav->set_inconsistent (false);
    wavelet->set_active (wav->get_active ());

    waveletConn.block (false);
}

void PartialPasteDlg::colorToggled ()
{

    icmConn.block (true);
    //gamcsconn.block (true);
    vibranceConn.block (true);
    chmixerConn.block (true);
    chmixerbwConn.block (true);
    hsveqConn.block (true);
    filmSimulationConn.block (true);
    rgbcurvesConn.block (true);
    colortoningConn.block (true);

    color->set_inconsistent (false);
    icm->set_active (color->get_active ());
    //gam->set_active (color->get_active ());
    vibrance->set_active (color->get_active ());
    chmixer->set_active (color->get_active ());
    blackwhite->set_active (color->get_active ());
    hsveq->set_active (color->get_active ());
    filmSimulation->set_active (color->get_active ());
    rgbcurves->set_active (color->get_active ());
    colortoning->set_active(color->get_active ());

    icmConn.block (false);
    //gamcsconn.block (false);
    vibranceConn.block (false);
    chmixerbwConn.block (false);
    chmixerConn.block (false);
    hsveqConn.block (false);
    filmSimulationConn.block (false);
    rgbcurvesConn.block (false);
    colortoningConn.block (false);
}

void PartialPasteDlg::lensToggled ()
{

    distortionConn.block (true);
    cacorrConn.block (true);
    vignettingConn.block (true);
    lcpConn.block (true);

    lens->set_inconsistent (false);

    distortion->set_active (lens->get_active ());
    cacorr->set_active (lens->get_active ());
    vignetting->set_active (lens->get_active ());
    lcp->set_active (lens->get_active ());

    distortionConn.block (false);
    cacorrConn.block (false);
    vignettingConn.block (false);
    lcpConn.block (false);
}

void PartialPasteDlg::compositionToggled ()
{

    coarserotConn.block (true);
    finerotConn.block (true);
    cropConn.block (true);
    resizeConn.block (true);
    prsharpeningConn.block (true);
    perspectiveConn.block (true);
    commonTransConn.block (true);

    composition->set_inconsistent (false);

    coarserot->set_active (composition->get_active ());
    finerot->set_active (composition->get_active ());
    crop->set_active (composition->get_active ());
    resize->set_active (composition->get_active ());
    prsharpening->set_active (composition->get_active ());
    perspective->set_active (composition->get_active ());
    commonTrans->set_active (composition->get_active ());

    coarserotConn.block (false);
    finerotConn.block (false);
    cropConn.block (false);
    resizeConn.block (false);
    prsharpeningConn.block (false);
    perspectiveConn.block (false);
    commonTransConn.block (false);
}

void PartialPasteDlg::metaToggled ()
{

    exifchConn.block (true);
    iptcConn.block (true);
    meta->set_inconsistent (false);

    exifch->set_active (meta->get_active ());
    iptc->set_active (meta->get_active ());

    exifchConn.block (false);
    iptcConn.block (false);
}


/*
 * Copies the selected items from the source ProcParams+ParamsEdited(optional)
 * to the destination ProcParams.
 */
void PartialPasteDlg::applyPaste (rtengine::procparams::ProcParams* dstPP, ParamsEdited* dstPE, const rtengine::procparams::ProcParams* srcPP, const ParamsEdited* srcPE)
{

    ParamsEdited falsePE;  // falsePE is a workaround to set a group of ParamsEdited to false
    ParamsEdited filterPE(true); // Contains the initial information about the loaded values

    if (srcPE) {
        filterPE = *srcPE;
    }

    // the general section is always ignored, whichever operation we use the PartialPaste for
    filterPE.general = falsePE.general;


    // Now we filter out the filter depending on the checked items
    if (!wb->get_active ()) {
        filterPE.wb         = falsePE.wb;
    }

    if (!exposure->get_active ()) {
        filterPE.toneCurve  = falsePE.toneCurve;
    }

    if (!sh->get_active ()) {
        filterPE.sh         = falsePE.sh;
    }

    if (!epd->get_active ()) {
        filterPE.epd        = falsePE.epd;
    }

    if (!retinex->get_active ()) {
        filterPE.retinex        = falsePE.retinex;
    }

    if (!pcvignette->get_active ()) {
        filterPE.pcvignette = falsePE.pcvignette;
    }

    if (!gradient->get_active ()) {
        filterPE.gradient   = falsePE.gradient;
    }

    if (!labcurve->get_active ()) {
        filterPE.labCurve   = falsePE.labCurve;
    }

    if (!colorappearance->get_active ()) {
        filterPE.colorappearance = falsePE.colorappearance;
    }

    if (!sharpen->get_active ()) {
        filterPE.sharpening      = falsePE.sharpening;
    }

    if (!sharpenedge->get_active ()) {
        filterPE.sharpenEdge     = falsePE.sharpenEdge;
    }

    if (!sharpenmicro->get_active()) {
        filterPE.sharpenMicro    = falsePE.sharpenMicro;
    }

    if (!impden->get_active ()) {
        filterPE.impulseDenoise  = falsePE.impulseDenoise;
    }

    if (!dirpyreq->get_active ()) {
        filterPE.dirpyrequalizer = falsePE.dirpyrequalizer;
    }

    if (!defringe->get_active ()) {
        filterPE.defringe        = falsePE.defringe;
    }

    if (!dirpyrden->get_active ()) {
        filterPE.dirpyrDenoise   = falsePE.dirpyrDenoise;
    }

    if (!wavelet->get_active ()) {
        filterPE.wavelet = falsePE.wavelet;
    }

    if (!icm->get_active ()) {
        filterPE.icm          = falsePE.icm;
    }

    if (!vibrance->get_active ()) {
        filterPE.vibrance     = falsePE.vibrance;
    }

    if (!chmixer->get_active ()) {
        filterPE.chmixer      = falsePE.chmixer;
    }

    if (!blackwhite->get_active ()) {
        filterPE.blackwhite   = falsePE.blackwhite;
    }

    if (!hsveq->get_active ()) {
        filterPE.hsvequalizer = falsePE.hsvequalizer;
    }

    if (!filmSimulation->get_active ()) {
        filterPE.filmSimulation  = falsePE.filmSimulation;
    }

    if (!rgbcurves->get_active ()) {
        filterPE.rgbCurves    = falsePE.rgbCurves;
    }

    if (!colortoning->get_active ()) {
        filterPE.colorToning  = falsePE.colorToning;
    }

    if (!distortion->get_active ()) {
        filterPE.distortion   = falsePE.distortion;
    }

    if (!cacorr->get_active ()) {
        filterPE.cacorrection = falsePE.cacorrection;
    }

    if (!vignetting->get_active ()) {
        filterPE.vignetting   = falsePE.vignetting;
    }

    if (!lcp->get_active ()) {
        filterPE.lensProf     = falsePE.lensProf;
    }

    if (!coarserot->get_active ()) {
        filterPE.coarse      = falsePE.coarse;
    }

    if (!finerot->get_active ()) {
        filterPE.rotate      = falsePE.rotate;
    }

    if (!crop->get_active ()) {
        filterPE.crop        = falsePE.crop;
    }

    if (!resize->get_active ()) {
        filterPE.resize      = falsePE.resize;
    }

    if (!prsharpening->get_active ()) {
        filterPE.prsharpening      = falsePE.prsharpening;
    }

    if (!perspective->get_active ()) {
        filterPE.perspective = falsePE.perspective;
    }

    if (!commonTrans->get_active ()) {
        filterPE.commonTrans = falsePE.commonTrans;
    }

    if (!exifch->get_active ()) {
        filterPE.exif = falsePE.exif;
    }

    if (!iptc->get_active ()) {
        filterPE.iptc = falsePE.iptc;
    }

    if (!raw_method->get_active ()) {
        filterPE.raw.bayersensor.method   = falsePE.raw.bayersensor.method;
        filterPE.raw.xtranssensor.method  = falsePE.raw.xtranssensor.method;
    }

    if (!raw_ccSteps->get_active ()) {
        filterPE.raw.bayersensor.ccSteps  = falsePE.raw.bayersensor.ccSteps;
        filterPE.raw.xtranssensor.ccSteps = falsePE.raw.xtranssensor.ccSteps;
    }

    if (!raw_dcb_iterations->get_active ()) {
        filterPE.raw.bayersensor.dcbIterations   = falsePE.raw.bayersensor.dcbIterations;
    }

    if (!raw_dcb_enhance->get_active ()) {
        filterPE.raw.bayersensor.dcbEnhance      = falsePE.raw.bayersensor.dcbEnhance;
    }

    //if (!raw_all_enhance->get_active ())     filterPE.raw.bayersensor.allEnhance      = falsePE.raw.bayersensor.allEnhance;
    if (!raw_lmmse_iterations->get_active ()) {
        filterPE.raw.bayersensor.lmmseIterations = falsePE.raw.bayersensor.lmmseIterations;
    }

    if (!raw_black->get_active ()) {
        filterPE.raw.bayersensor.exBlack0        = falsePE.raw.bayersensor.exBlack0;
        filterPE.raw.bayersensor.exBlack1        = falsePE.raw.bayersensor.exBlack1;
        filterPE.raw.bayersensor.exBlack2        = falsePE.raw.bayersensor.exBlack2;
        filterPE.raw.bayersensor.exBlack3        = falsePE.raw.bayersensor.exBlack3;
        filterPE.raw.bayersensor.exTwoGreen      = falsePE.raw.bayersensor.exTwoGreen;
        filterPE.raw.xtranssensor.exBlackRed     = falsePE.raw.xtranssensor.exBlackRed;
        filterPE.raw.xtranssensor.exBlackGreen   = falsePE.raw.xtranssensor.exBlackGreen;
        filterPE.raw.xtranssensor.exBlackBlue    = falsePE.raw.xtranssensor.exBlackBlue;
    }

    if (!raw_linenoise->get_active ()) {
        filterPE.raw.bayersensor.linenoise       = falsePE.raw.bayersensor.linenoise;
    }

    if (!raw_greenthresh->get_active ()) {
        filterPE.raw.bayersensor.greenEq         = falsePE.raw.bayersensor.greenEq;
    }

    if (!raw_expos->get_active ()) {
        filterPE.raw.exPos              = falsePE.raw.exPos;
    }

    if (!raw_preser->get_active ()) {
        filterPE.raw.exPreser           = falsePE.raw.exPreser;
    }

    if (!raw_ca_autocorrect->get_active ()) {
        filterPE.raw.caCorrection       = falsePE.raw.caCorrection;
    }

    if (!raw_cared->get_active ()) {
        filterPE.raw.caRed              = falsePE.raw.caRed;
    }

    if (!raw_cablue->get_active ()) {
        filterPE.raw.caBlue             = falsePE.raw.caBlue;
    }

    if (!raw_hotpix_filt->get_active ())     {
        filterPE.raw.hotPixelFilter     = falsePE.raw.hotPixelFilter;
    }

    if (!raw_deadpix_filt->get_active ())    {
        filterPE.raw.deadPixelFilter    = falsePE.raw.deadPixelFilter;
    }

    if (!raw_deadpix_filt->get_active () && !raw_hotpix_filt->get_active ()) {
        filterPE.raw.hotDeadPixelThresh = falsePE.raw.hotDeadPixelThresh;
    }

    if (!df_file->get_active ()) {
        filterPE.raw.darkFrame          = falsePE.raw.darkFrame;
    }

    if (!df_AutoSelect->get_active ()) {
        filterPE.raw.dfAuto             = falsePE.raw.dfAuto;
    }

    if (!ff_file->get_active ()) {
        filterPE.raw.ff_file            = falsePE.raw.ff_file;
    }

    if (!ff_AutoSelect->get_active ()) {
        filterPE.raw.ff_AutoSelect      = falsePE.raw.ff_AutoSelect;
    }

    if (!ff_BlurRadius->get_active ()) {
        filterPE.raw.ff_BlurRadius      = falsePE.raw.ff_BlurRadius;
    }

    if (!ff_BlurType->get_active ()) {
        filterPE.raw.ff_BlurType        = falsePE.raw.ff_BlurType;
    }

    if (!ff_ClipControl->get_active ()) {
        filterPE.raw.ff_clipControl     = falsePE.raw.ff_clipControl;
        filterPE.raw.ff_AutoClipControl = falsePE.raw.ff_AutoClipControl;
    }

    if (dstPE) {
        *dstPE = filterPE;
    }

    // Apply the filter!
    filterPE.combine(*dstPP, *srcPP, true);
}

