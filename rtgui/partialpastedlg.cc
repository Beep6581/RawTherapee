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

PartialPasteDlg::PartialPasteDlg (Glib::ustring title) {

    set_modal (true);
    set_title (title);

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
    metaicm     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_METAICMGROUP")));
    metaicm     ->set_name("partialPasteHeader");
    raw         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWGROUP")));
    raw         ->set_name("partialPasteHeader");

    // options in basic:
    wb          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WHITEBALANCE")));
    exposure    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXPOSURE")));
    hlrec       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_HLRECONSTRUCTION")));
    sh          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHADOWSHIGHLIGHTS")));
    epd         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EPD")));
    labcurve    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LABCURVE")));

    // options in detail:
    sharpen     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENING")));
    sharpenedge = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENEDGE")));
    sharpenmicro = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENMICRO")));
    impden		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IMPULSEDENOISE")));
    dirpyreq    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYREQUALIZER")));
    defringe    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DEFRINGE")));

    // options in color:
    vibrance    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_VIBRANCE")));
    chmixer     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CHANNELMIXER")));
    dirpyrden   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYRDENOISE")));
    hsveq		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_HSVEQUALIZER")));
	rgbcurves   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RGBCURVES")));

    // options in lens:
    distortion  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DISTORTION")));
    cacorr      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CACORRECTION")));
    vignetting  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_VIGNETTING")));
    lcp         = Gtk::manage (new Gtk::CheckButton (M("TP_LENSPROFILE_LABEL")));

    // options in composition:
    coarserot   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COARSETRANS")));
    finerot     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ROTATION")));
    crop        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CROP")));
    resize      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RESIZE")));
    perspective = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PERSPECTIVE")));
    commonTrans = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMMONTRANSFORMPARAMS")));

    // options in metaicm:
    exifch      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXIFCHANGES")));
    iptc        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IPTCINFO")));
    icm         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ICMSETTINGS")));
    gam         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ICMGAMMA")));

    // options in raw:
    raw_expos			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_LINEAR")));
    raw_preser			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_PRESER")));
    raw_black			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWEXPOS_BLACK")));
    raw_ca_autocorrect	= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_AUTO")));
    raw_cared			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_CARED")));
    raw_cablue			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWCACORR_CABLUE")));
    raw_hotdeadpix_filt	= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_HOTDEADPIXFILT")));
    raw_linenoise		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_LINEDENOISE")));
    raw_greenthresh		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_PREPROCESS_GREENEQUIL")));
    raw_dmethod			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DMETHOD")));
    raw_ccSteps			= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_FALSECOLOR")));
    raw_dcb_iterations	= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DCBITERATIONS")));
    raw_dcb_enhance		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_DCBENHANCE")));
    //raw_all_enhance		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAW_ALLENHANCE")));

	df_file        		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEFILE")));
	df_AutoSelect  		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEAUTOSELECT")));
	ff_file        		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDFILE")));
	ff_AutoSelect  		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDAUTOSELECT")));
	ff_BlurRadius  		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURRADIUS")));
	ff_BlurType    		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURTYPE")));

    Gtk::VBox* vboxes[7];
    Gtk::HSeparator* hseps[7];
    for (int i=0; i<7; i++) {
        vboxes[i] = Gtk::manage (new Gtk::VBox ());
        vboxes[i]->set_border_width (6);
        hseps[i] = Gtk::manage (new Gtk::HSeparator ());
        hseps[i]->set_name("partialPasteHeaderSep");
    }
    
    vboxes[0]->pack_start (*basic, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*hseps[0], Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*wb, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*exposure, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*hlrec, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*sh, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*epd, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*labcurve, Gtk::PACK_SHRINK, 2);

    vboxes[1]->pack_start (*detail, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*hseps[1], Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpen, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpenedge, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpenmicro, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*impden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*dirpyrden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*defringe, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*dirpyreq, Gtk::PACK_SHRINK, 2);
    //vboxes[1]->pack_start (*waveq, Gtk::PACK_SHRINK, 2);

    vboxes[2]->pack_start (*color, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hseps[2], Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*vibrance, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*chmixer, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hsveq, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*rgbcurves, Gtk::PACK_SHRINK, 2);

    vboxes[3]->pack_start (*lens, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*hseps[3], Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*distortion, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*cacorr, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*vignetting, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*lcp, Gtk::PACK_SHRINK, 2);

    vboxes[4]->pack_start (*composition, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*hseps[4], Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*coarserot, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*finerot, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*crop, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*resize, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*perspective, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*commonTrans, Gtk::PACK_SHRINK, 2);

    vboxes[5]->pack_start (*metaicm, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*hseps[5], Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*exifch, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*iptc, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*icm, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*gam, Gtk::PACK_SHRINK, 2);

    vboxes[6]->pack_start (*raw, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*hseps[6], Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_dmethod, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_ccSteps, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_dcb_iterations, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_dcb_enhance, Gtk::PACK_SHRINK, 2);
	//vboxes[6]->pack_start (*raw_all_enhance, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
	vboxes[6]->pack_start (*raw_linenoise, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_greenthresh, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_hotdeadpix_filt, Gtk::PACK_SHRINK, 2);
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
	vboxes[6]->pack_start (*Gtk::manage (new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 0);
	vboxes[6]->pack_start (*raw_ca_autocorrect, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_cared, Gtk::PACK_SHRINK, 2);
	vboxes[6]->pack_start (*raw_cablue, Gtk::PACK_SHRINK, 2);

    Gtk::VBox* vbCol1 = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbCol2 = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbCol3 = Gtk::manage (new Gtk::VBox ());

    vbCol1->set_border_width (8);
    vbCol2->set_border_width (8);
    vbCol3->set_border_width (8);

    for (int i=0; i<3; i++)
        vbCol1->pack_start (*vboxes[i]);
    for (int i=3; i<6; i++)
        vbCol2->pack_start (*vboxes[i]);
    for (int i=6; i<7; i++)
        vbCol3->pack_start (*vboxes[i]);

	Gtk::VBox* vbtop = Gtk::manage (new Gtk::VBox ());
	vbtop->pack_start (*everything, Gtk::PACK_SHRINK, 2);
    vbtop->pack_start (*(Gtk::manage (new Gtk::HSeparator ())));
    vbtop->set_border_width (8);

	get_vbox()->pack_start (*vbtop);

    Gtk::HBox* hbmain = Gtk::manage (new Gtk::HBox ());
    hbmain->pack_start (*vbCol1);
    hbmain->pack_start (*(Gtk::manage (new Gtk::VSeparator ())));
    hbmain->pack_start (*vbCol2);
    hbmain->pack_start (*(Gtk::manage (new Gtk::VSeparator ())));
    hbmain->pack_start (*vbCol3);

    get_vbox()->pack_start (*hbmain);

    // This can be improved
    // there is currently no binding of subsettings to CheckButton 'everything' for its inconsistent status
    everythingConn  = everything->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::everythingToggled));
    basicConn       = basic->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::basicToggled));    
    detailConn      = detail->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::detailToggled));
    colorConn       = color->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::colorToggled));    
    lensConn        = lens->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::lensToggled));    
    compositionConn = composition->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::compositionToggled));    
    metaicmConn     = metaicm->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::metaicmToggled));    
    rawConn         = raw->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::rawToggled));

    wbConn          = wb->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    exposureConn    = exposure->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    hlrecConn       = hlrec->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    shConn          = sh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    epdConn         = epd->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));
    labcurveConn    = labcurve->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));

    sharpenConn     = sharpen->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    gradsharpenConn = sharpenedge->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    microcontrastConn = sharpenmicro->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    impdenConn		= impden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    dirpyrdenConn   = dirpyrden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    dirpyreqConn	= dirpyreq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    //waveqConn	    = waveq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));
    defringeConn    = defringe->signal_toggled().connect (sigc::bind (sigc::mem_fun(*detail, &Gtk::CheckButton::set_inconsistent), true));

    vibranceConn    = vibrance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    chmixerConn     = chmixer->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    hsveqConn		= hsveq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));
    rgbcurvesConn   = rgbcurves->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));

    distortionConn  = distortion->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    
    cacorrConn      = cacorr->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    
    vignettingConn  = vignetting->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    
    lcpConn         = lcp->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    

    coarserotConn   = coarserot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    finerotConn     = finerot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    cropConn        = crop->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    resizeConn      = resize->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    perspectiveConn = perspective->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));
    commonTransConn = commonTrans->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));

    exifchConn      = exifch->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));    
    iptcConn        = iptc->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));    
    icmConn         = icm->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));    
    gamcsconn		= gam->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));

    raw_dmethodConn         = raw_dmethod->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_ccStepsConn         = raw_ccSteps->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_dcb_iterationsConn  = raw_dcb_iterations->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_dcb_enhanceConn     = raw_dcb_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    //raw_all_enhanceConn     = raw_all_enhance->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));

    raw_exposConn           = raw_expos->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_preserConn          = raw_preser->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_blackConn           = raw_black->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_ca_autocorrectConn  = raw_ca_autocorrect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_caredConn           = raw_cared->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_cablueConn          = raw_cablue->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_hotdeadpix_filtConn = raw_hotdeadpix_filt->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_linenoiseConn       = raw_linenoise->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    raw_greenthreshConn     = raw_greenthresh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    df_fileConn             = df_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    df_AutoSelectConn       = df_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_fileConn             = ff_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_AutoSelectConn       = ff_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_BlurRadiusConn       = ff_BlurRadius->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_BlurTypeConn         = ff_BlurType->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));

    add_button (Gtk::StockID("gtk-ok"), Gtk::RESPONSE_OK);
    add_button (Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    set_response_sensitive (Gtk::RESPONSE_OK);
    set_default_response (Gtk::RESPONSE_OK);
    show_all_children ();
}

void PartialPasteDlg::everythingToggled () {

	basicConn.block (true);
	detailConn.block (true);
	colorConn.block (true);
	lensConn.block (true);
	compositionConn.block (true);
	metaicmConn.block (true);
	rawConn.block (true);

	everything->set_inconsistent (false);

	//toggle group headings
    basic->set_active(everything->get_active());
    detail->set_active(everything->get_active());
    color->set_active(everything->get_active());
    lens->set_active(everything->get_active());
    composition->set_active(everything->get_active());
    metaicm->set_active(everything->get_active());
    raw->set_active(everything->get_active());

    //toggle group children
    PartialPasteDlg::basicToggled ();
    PartialPasteDlg::detailToggled ();
    PartialPasteDlg::colorToggled ();
    PartialPasteDlg::lensToggled ();
    PartialPasteDlg::compositionToggled ();
    PartialPasteDlg::metaicmToggled ();
    PartialPasteDlg::rawToggled ();

	basicConn.block (false);
	detailConn.block (false);
	colorConn.block (false);
	lensConn.block (false);
	compositionConn.block (false);
	metaicmConn.block (false);
	rawConn.block (false);
}

void PartialPasteDlg::rawToggled () {

	raw_dmethodConn.block (true);
	raw_ccStepsConn.block (true);
	raw_dcb_iterationsConn.block (true);
	raw_dcb_enhanceConn.block (true);
	//raw_all_enhanceConn.block (true);
	raw_exposConn.block (true);
	raw_preserConn.block (true);
	raw_blackConn.block (true);
	raw_ca_autocorrectConn.block (true);
	raw_caredConn.block (true);
	raw_cablueConn.block (true);
	raw_hotdeadpix_filtConn.block (true);
	raw_linenoiseConn.block (true);
	raw_greenthreshConn.block (true);
	df_fileConn.block (true);
	df_AutoSelectConn.block (true);
	ff_fileConn.block (true);
	ff_AutoSelectConn.block (true);
	ff_BlurRadiusConn.block (true);
	ff_BlurTypeConn.block (true);

    raw->set_inconsistent (false);

    raw_dmethod->set_active (raw->get_active ());
    raw_ccSteps->set_active (raw->get_active ());
    raw_dcb_iterations->set_active (raw->get_active ());
    raw_dcb_enhance->set_active (raw->get_active ());
    //raw_all_enhance->set_active (raw->get_active ());
    raw_expos->set_active (raw->get_active ());
    raw_preser->set_active (raw->get_active ());
    raw_black->set_active (raw->get_active ());
    raw_ca_autocorrect->set_active (raw->get_active ());
    raw_cared->set_active (raw->get_active ());
    raw_cablue->set_active (raw->get_active ());
    raw_hotdeadpix_filt->set_active (raw->get_active ());
    raw_linenoise->set_active (raw->get_active ());
    raw_greenthresh->set_active (raw->get_active ());
    df_file->set_active (raw->get_active ());
    df_AutoSelect->set_active (raw->get_active ());
    ff_file->set_active (raw->get_active ());
    ff_AutoSelect->set_active (raw->get_active ());
    ff_BlurRadius->set_active (raw->get_active ());
    ff_BlurType->set_active (raw->get_active ());

    raw_dmethodConn.block (false);
    raw_ccStepsConn.block (false);
    raw_dcb_iterationsConn.block (false);
    raw_dcb_enhanceConn.block (false);
    //raw_all_enhanceConn.block (false);
    raw_exposConn.block (false);
    raw_preserConn.block (false);
    raw_blackConn.block (false);
    raw_ca_autocorrectConn.block (false);
    raw_caredConn.block (false);
    raw_cablueConn.block (false);
    raw_hotdeadpix_filtConn.block (false);
    raw_linenoiseConn.block (false);
    raw_greenthreshConn.block (false);
    df_fileConn.block (false);
    df_AutoSelectConn.block (false);
    ff_fileConn.block (false);
    ff_AutoSelectConn.block (false);
    ff_BlurRadiusConn.block (false);
    ff_BlurTypeConn.block (false);
}

void PartialPasteDlg::basicToggled () {

    wbConn.block (true);
    exposureConn.block (true);
    hlrecConn.block (true);
    shConn.block (true);
    epdConn.block(true);
    labcurveConn.block (true);

    basic->set_inconsistent (false);

    wb->set_active (basic->get_active ());
    exposure->set_active (basic->get_active ());
    hlrec->set_active (basic->get_active ());
    sh->set_active (basic->get_active ());
    epd->set_active (basic->get_active ());
    labcurve->set_active (basic->get_active ());

    wbConn.block (false);
    exposureConn.block (false);
    hlrecConn.block (false);
    shConn.block (false);
    epdConn.block (false);
    labcurveConn.block (false);
}

void PartialPasteDlg::detailToggled () {

    sharpenConn.block (true);
    gradsharpenConn.block(true);
    microcontrastConn.block(true);
    impdenConn.block (true);
    dirpyrdenConn.block (true);
    defringeConn.block (true);
    dirpyreqConn.block (true);
    //waveqConn.block (true);

    detail->set_inconsistent (false);

    sharpen->set_active (detail->get_active ());
    sharpenedge->set_active (detail->get_active ());
    sharpenmicro->set_active (detail->get_active ());
    impden->set_active (detail->get_active ());
    dirpyrden->set_active (detail->get_active ());
    defringe->set_active (detail->get_active ());
    dirpyreq->set_active (detail->get_active ());
    //waveq->set_active (detail->get_active ());

    sharpenConn.block (false);
    gradsharpenConn.block(false);
    microcontrastConn.block(false);
    impdenConn.block (false);
    dirpyrdenConn.block (false);
    defringeConn.block (false);
    dirpyreqConn.block (false);
    //waveqConn.block (false);
}

void PartialPasteDlg::colorToggled () {

    vibranceConn.block (true);
    chmixerConn.block (true);
    hsveqConn.block (true);
    rgbcurvesConn.block (true);
    gamcsconn.block (true);
    color->set_inconsistent (false);

    vibrance->set_active (color->get_active ());
    chmixer->set_active (color->get_active ());
    hsveq->set_active (color->get_active ());
    rgbcurves->set_active (color->get_active ());
    icm->set_active (color->get_active ());

    vibranceConn.block (false);
    chmixerConn.block (false);
    hsveqConn.block (false);
    rgbcurvesConn.block (false);
    gamcsconn.block (false);

}

void PartialPasteDlg::lensToggled () {

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

void PartialPasteDlg::compositionToggled () {

    coarserotConn.block (true);
    finerotConn.block (true);
    cropConn.block (true);
    resizeConn.block (true);
    perspectiveConn.block (true);
    commonTransConn.block (true);

    composition->set_inconsistent (false);

    coarserot->set_active (composition->get_active ());
    finerot->set_active (composition->get_active ());
    crop->set_active (composition->get_active ());
    resize->set_active (composition->get_active ());
    perspective->set_active (composition->get_active ());
    commonTrans->set_active (composition->get_active ());

    coarserotConn.block (false);
    finerotConn.block (false);
    cropConn.block (false);
    resizeConn.block (false);
    perspectiveConn.block (false);
    commonTransConn.block (false);
}

void PartialPasteDlg::metaicmToggled () {

    exifchConn.block (true);
    iptcConn.block (true);
    icmConn.block (true);
    gamcsconn.block (true);
    metaicm->set_inconsistent (false);

    exifch->set_active (metaicm->get_active ());
    iptc->set_active (metaicm->get_active ());
    icm->set_active (metaicm->get_active ());
    gam->set_active (metaicm->get_active ());

    exifchConn.block (false);
    iptcConn.block (false);
    icmConn.block (false);
    gamcsconn.block (false);

}


/*
 * Copies the selected items from the source ProcParams+ParamsEdited(optional)
 * to the destination ProcParams.
 */
void PartialPasteDlg::applyPaste (rtengine::procparams::ProcParams* dstPP, ParamsEdited* dstPE, const rtengine::procparams::ProcParams* srcPP, const ParamsEdited* srcPE) {

    ParamsEdited falsePE;  // falsePE is a workaround to set a group of ParamsEdited to false
    ParamsEdited filterPE; // Contains the initial information about the loaded values
    if (srcPE) {
        filterPE = *srcPE;
    }
    else {
        // By default, everything has to be copied
        filterPE.set(true);
    }


    // the general section is always ignored, whichever operation we use the PartialPaste for
    filterPE.general = falsePE.general;


    // Now we filter out the filter depending on the checked items
    if (!wb->get_active ())          filterPE.wb         = falsePE.wb;
    if (!exposure->get_active ())    filterPE.toneCurve  = falsePE.toneCurve;
    if (!hlrec->get_active ())       filterPE.hlrecovery = falsePE.hlrecovery;
    if (!sh->get_active ())          filterPE.sh         = falsePE.sh;
    if (!epd->get_active ())         filterPE.edgePreservingDecompositionUI = falsePE.edgePreservingDecompositionUI;
    if (!labcurve->get_active ())    filterPE.labCurve   = falsePE.labCurve;

    if (!sharpen->get_active ())     filterPE.sharpening      = falsePE.sharpening;
    if (!sharpenedge->get_active ()) filterPE.sharpenEdge     = falsePE.sharpenEdge;
    if (!sharpenmicro->get_active()) filterPE.sharpenMicro    = falsePE.sharpenMicro;
    if (!impden->get_active ())      filterPE.impulseDenoise  = falsePE.impulseDenoise;
    if (!dirpyreq->get_active ())    filterPE.dirpyrequalizer = falsePE.dirpyrequalizer;
    if (!defringe->get_active ())    filterPE.defringe        = falsePE.defringe;
    if (!dirpyrden->get_active ())   filterPE.dirpyrDenoise   = falsePE.dirpyrDenoise;

    if (!vibrance->get_active ())    filterPE.vibrance     = falsePE.vibrance;
    if (!chmixer->get_active ())     filterPE.chmixer      = falsePE.chmixer;
    if (!hsveq->get_active ())       filterPE.hsvequalizer = falsePE.hsvequalizer;
    if (!rgbcurves->get_active ())   filterPE.rgbCurves    = falsePE.rgbCurves;
    if (!icm->get_active ())         filterPE.icm          = falsePE.icm;

    if (!distortion->get_active ())  filterPE.distortion   = falsePE.distortion;
    if (!cacorr->get_active ())      filterPE.cacorrection = falsePE.cacorrection;
    if (!vignetting->get_active ())  filterPE.vignetting   = falsePE.vignetting;
    if (!lcp->get_active ())         filterPE.lensProf     = falsePE.lensProf;

    if (!coarserot->get_active ())   filterPE.coarse      = falsePE.coarse;
    if (!finerot->get_active ())     filterPE.rotate      = falsePE.rotate;
    if (!crop->get_active ())        filterPE.crop        = falsePE.crop;
    if (!resize->get_active ())      filterPE.resize      = falsePE.resize;
    if (!perspective->get_active ()) filterPE.perspective = falsePE.perspective;
    if (!commonTrans->get_active ()) filterPE.commonTrans = falsePE.commonTrans;

    if (!exifch->get_active ())      filterPE.exif = falsePE.exif;
    if (!iptc->get_active ())        filterPE.iptc = falsePE.iptc;
    if (!icm->get_active ())         filterPE.icm  = falsePE.icm;

    if (!raw_dmethod->get_active ())           filterPE.raw.dmethod            = falsePE.raw.dmethod;
    if (!raw_ccSteps->get_active ())           filterPE.raw.ccSteps            = falsePE.raw.ccSteps;
    if (!raw_dcb_iterations->get_active ())    filterPE.raw.dcbIterations      = falsePE.raw.dcbIterations;
    if (!raw_dcb_enhance->get_active ())       filterPE.raw.dcbEnhance         = falsePE.raw.dcbEnhance;
    //if (!raw_all_enhance->get_active ())       filterPE.raw.allEnhance         = falsePE.raw.allEnhance;

    if (!raw_expos->get_active ())             filterPE.raw.exPos              = falsePE.raw.exPos;
    if (!raw_preser->get_active ())            filterPE.raw.exPreser           = falsePE.raw.exPreser;
    if (!raw_black->get_active ()) {           filterPE.raw.exBlackzero        = falsePE.raw.exBlackzero;
                                               filterPE.raw.exBlackone         = falsePE.raw.exBlackone;
                                               filterPE.raw.exBlacktwo         = falsePE.raw.exBlacktwo;
                                               filterPE.raw.exBlackthree       = falsePE.raw.exBlackthree;
                                               filterPE.raw.exTwoGreen         = falsePE.raw.exTwoGreen; }
    if (!raw_ca_autocorrect->get_active ())    filterPE.raw.caCorrection       = falsePE.raw.caCorrection;
    if (!raw_cared->get_active ())             filterPE.raw.caRed              = falsePE.raw.caRed;
    if (!raw_cablue->get_active ())            filterPE.raw.caBlue             = falsePE.raw.caBlue;
    if (!raw_hotdeadpix_filt->get_active ()) { filterPE.raw.hotDeadPixelFilter = falsePE.raw.hotDeadPixelFilter;
                                               filterPE.raw.hotDeadPixelThresh = falsePE.raw.hotDeadPixelThresh; }
    if (!raw_linenoise->get_active ())         filterPE.raw.linenoise          = falsePE.raw.linenoise;
    if (!raw_greenthresh->get_active ())       filterPE.raw.greenEq            = falsePE.raw.greenEq;
    if (!df_file->get_active ())               filterPE.raw.darkFrame          = falsePE.raw.darkFrame;
    if (!df_AutoSelect->get_active ())         filterPE.raw.dfAuto             = falsePE.raw.dfAuto;
    if (!ff_file->get_active ())               filterPE.raw.ff_file            = falsePE.raw.ff_file;
    if (!ff_AutoSelect->get_active ())         filterPE.raw.ff_AutoSelect      = falsePE.raw.ff_AutoSelect;
    if (!ff_BlurRadius->get_active ())         filterPE.raw.ff_BlurRadius      = falsePE.raw.ff_BlurRadius;
    if (!ff_BlurType->get_active ())           filterPE.raw.ff_BlurType        = falsePE.raw.ff_BlurType;

    if (dstPE) *dstPE = filterPE;

    // Apply the filter!
    filterPE.combine(*dstPP, *srcPP, true);
}

