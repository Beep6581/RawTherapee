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
#include <partialpastedlg.h>
#include <multilangmgr.h>

PartialPasteDlg::PartialPasteDlg () {

    set_modal (true);
    set_title (M("PARTIALPASTE_DIALOGLABEL"));

		everything  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EVERYTHING")));
	
    basic       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_BASICGROUP")));
    luminance   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LUMINANCEGROUP")));
    color       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORGROUP")));
    lens        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LENSGROUP")));
    composition = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMPOSITIONGROUP")));
    metaicm     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_METAICMGROUP")));
    raw         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RAWGROUP")));

    // options in basic:
    wb          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WHITEBALANCE")));
    exposure    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXPOSURE")));
    hlrec       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_HLRECOVERY")));

    // options in luminance:
    sharpen     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENING")));
    impden		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IMPULSEDENOISE")));
    lumaden     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LUMADENOISE")));
    labcurve   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LABCURVE")));
    sh          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHADOWSHIGHLIGHTS")));
    dirpyreq    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYREQUALIZER")));
    waveq       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WAVELETEQUALIZER")));

    // options in color:
    colormixer  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORMIXER")));
    colorshift  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORSHIFT")));
    colorboost  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORBOOST")));
    colorden    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORDENOISE")));
    dirpyrden   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DIRPYRDENOISE")));
    hsveq		= Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_HSVEQUALIZER")));

    // options in lens:
    distortion  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DISTORTION")));
    cacorr      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CACORRECTION")));
    vignetting  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_VIGNETTING")));

    // options in composition:
    coarserot   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COARSETRANS")));
    finerot     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ROTATION")));
    crop        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_CROP")));
    resize      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_RESIZE")));

    // options in metaicm:
    exifch      = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXIFCHANGES")));
    iptc        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_IPTCINFO")));
    icm         = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_ICMSETTINGS")));

    // options in raw:
    df_file        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEFILE")));      
    df_AutoSelect  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_DARKFRAMEAUTOSELECT")));   
    ff_file        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDFILE")));      
    ff_AutoSelect  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDAUTOSELECT")));   
    ff_BlurRadius  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURRADIUS")));
    ff_BlurType    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_FLATFIELDBLURTYPE")));  
                                                                                                   
    Gtk::VBox* vboxes[7];                                                                          
    Gtk::HSeparator* hseps[7];                                                                     
    for (int i=0; i<7; i++) {                                                                      

        vboxes[i] = Gtk::manage (new Gtk::VBox ());
        vboxes[i]->set_border_width (16);
        hseps[i] = Gtk::manage (new Gtk::HSeparator ());
    }
    
    vboxes[0]->pack_start (*basic, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*hseps[0], Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*wb, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*exposure, Gtk::PACK_SHRINK, 2);
    vboxes[0]->pack_start (*hlrec, Gtk::PACK_SHRINK, 2);

    vboxes[1]->pack_start (*luminance, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*hseps[1], Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sharpen, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*impden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*lumaden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*labcurve, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sh, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*dirpyreq, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*waveq, Gtk::PACK_SHRINK, 2);

    vboxes[2]->pack_start (*color, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hseps[2], Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colormixer, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colorshift, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colorboost, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hsveq, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colorden, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*dirpyrden, Gtk::PACK_SHRINK, 2);


    vboxes[3]->pack_start (*lens, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*hseps[3], Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*distortion, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*cacorr, Gtk::PACK_SHRINK, 2);
    vboxes[3]->pack_start (*vignetting, Gtk::PACK_SHRINK, 2);

    vboxes[4]->pack_start (*composition, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*hseps[4], Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*coarserot, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*finerot, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*crop, Gtk::PACK_SHRINK, 2);
    vboxes[4]->pack_start (*resize, Gtk::PACK_SHRINK, 2);

    vboxes[5]->pack_start (*raw, Gtk::PACK_SHRINK, 2);          
    vboxes[5]->pack_start (*hseps[5], Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*df_file, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*df_AutoSelect, Gtk::PACK_SHRINK, 2);     
    vboxes[5]->pack_start (*ff_file, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*ff_AutoSelect, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*ff_BlurType, Gtk::PACK_SHRINK, 2);  
    vboxes[5]->pack_start (*ff_BlurRadius, Gtk::PACK_SHRINK, 2);

    vboxes[6]->pack_start (*metaicm, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*hseps[6], Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*exifch, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*iptc, Gtk::PACK_SHRINK, 2);
    vboxes[6]->pack_start (*icm, Gtk::PACK_SHRINK, 2);

    Gtk::VBox* vbCol1 = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbCol2 = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbCol3 = Gtk::manage (new Gtk::VBox ());

    vbCol1->set_border_width (16);
    vbCol2->set_border_width (16);
    vbCol3->set_border_width (16);

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
    luminanceConn   = luminance->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::luminanceToggled));    
    colorConn       = color->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::colorToggled));    
    lensConn        = lens->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::lensToggled));    
    compositionConn = composition->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::compositionToggled));    
    metaicmConn     = metaicm->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::metaicmToggled));
    rawConn         = raw->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::rawToggled));    

    wbConn          = wb->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    exposureConn    = exposure->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    hlrecConn       = hlrec->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    

    sharpenConn     = sharpen->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    impdenConn      = impden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    lumadenConn     = lumaden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    labcurveConn    = labcurve->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    shConn          = sh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    dirpyreqConn    = dirpyreq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    waveqConn       = waveq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    

    colormixerConn  = colormixer->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    colorshiftConn  = colorshift->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    colorboostConn  = colorboost->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    hsveqConn       = hsveq->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    colordenConn    = colorden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    dirpyrdenConn   = dirpyrden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    

    distortionConn  = distortion->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    
    cacorrConn      = cacorr->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    
    vignettingConn  = vignetting->signal_toggled().connect (sigc::bind (sigc::mem_fun(*lens, &Gtk::CheckButton::set_inconsistent), true));    

    coarserotConn   = coarserot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    finerotConn     = finerot->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    cropConn        = crop->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    
    resizeConn      = resize->signal_toggled().connect (sigc::bind (sigc::mem_fun(*composition, &Gtk::CheckButton::set_inconsistent), true));    

    exifchConn      = exifch->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));    
    iptcConn        = iptc->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));    
    icmConn         = icm->signal_toggled().connect (sigc::bind (sigc::mem_fun(*metaicm, &Gtk::CheckButton::set_inconsistent), true));    

    df_fileConn        = df_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    df_AutoSelectConn  = df_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_fileConn        = ff_file->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_AutoSelectConn  = ff_AutoSelect->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));   
    ff_BlurRadiusConn  = ff_BlurRadius->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));
    ff_BlurTypeConn    = ff_BlurType->signal_toggled().connect (sigc::bind (sigc::mem_fun(*raw, &Gtk::CheckButton::set_inconsistent), true));  

    add_button (Gtk::StockID("gtk-ok"), 1);
    add_button (Gtk::StockID("gtk-cancel"), 0);
    set_response_sensitive (1);
    set_default_response (1);
    show_all_children ();
}

void PartialPasteDlg::everythingToggled () {
		
		basicConn.block (true);
		luminanceConn.block (true);
		colorConn.block (true);
		lensConn.block (true);
		compositionConn.block (true);
		metaicmConn.block (true);
		rawConn.block (true);
		
		everything->set_inconsistent (false);
		
		//toggle group headings
    basic->set_active(everything->get_active());
    luminance->set_active(everything->get_active());
    color->set_active(everything->get_active());
    lens->set_active(everything->get_active());
    composition->set_active(everything->get_active());
    metaicm->set_active(everything->get_active());
    raw->set_active(everything->get_active());
                                 
    //toggle group children
    PartialPasteDlg::basicToggled ();      
    PartialPasteDlg::luminanceToggled ();  
    PartialPasteDlg::colorToggled ();      
    PartialPasteDlg::lensToggled ();       
    PartialPasteDlg::compositionToggled ();
    PartialPasteDlg::metaicmToggled ();    
    PartialPasteDlg::rawToggled ();        
                            
		basicConn.block (false);
		luminanceConn.block (false);
		colorConn.block (false);
		lensConn.block (false);
		compositionConn.block (false);
		metaicmConn.block (false);
		rawConn.block (false);    
}

void PartialPasteDlg::rawToggled () {              
                                                   
    df_fileConn.block (true);
    df_AutoSelectConn.block (true);
    ff_fileConn.block (true);
    ff_AutoSelectConn.block (true);                
    ff_BlurRadiusConn.block (true);                
    ff_BlurTypeConn.block (true);                  
                                                   
    raw->set_inconsistent (false);                 
    
    df_file->set_active (raw->get_active ());
    df_AutoSelect->set_active (raw->get_active ());                                               
    ff_file->set_active (raw->get_active ());
    ff_AutoSelect->set_active (raw->get_active ());
    ff_BlurRadius->set_active (raw->get_active ());
    ff_BlurType->set_active (raw->get_active ());  
                                                   
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

    basic->set_inconsistent (false);

    wb->set_active (basic->get_active ());
    exposure->set_active (basic->get_active ());
    hlrec->set_active (basic->get_active ());

    wbConn.block (false);
    exposureConn.block (false);
    hlrecConn.block (false);
}

void PartialPasteDlg::luminanceToggled () {

    sharpenConn.block (true);
    impdenConn.block (true);
    lumadenConn.block (true);
    labcurveConn.block (true);
    shConn.block (true);
    dirpyreqConn.block (true);
    waveqConn.block (true);

    luminance->set_inconsistent (false);

    sharpen->set_active (luminance->get_active ());
    impden->set_active (luminance->get_active ());
    lumaden->set_active (luminance->get_active ());
    labcurve->set_active (luminance->get_active ());
    sh->set_active (luminance->get_active ());
    dirpyreq->set_active (luminance->get_active ());
    waveq->set_active (luminance->get_active ());

    sharpenConn.block (false);
	impdenConn.block (false);
    lumadenConn.block (false);
    labcurveConn.block (false);
    shConn.block (false);
	dirpyreqConn.block (false);
    waveqConn.block (false);
}

void PartialPasteDlg::colorToggled () {

    colormixerConn.block (true);
    colorshiftConn.block (true);
    colorboostConn.block (true);
    hsveqConn.block (true);
    colordenConn.block (true);
    dirpyrdenConn.block (true);

    color->set_inconsistent (false);

    colormixer->set_active (color->get_active ());
    colorshift->set_active (color->get_active ());
    colorboost->set_active (color->get_active ());
    hsveq->set_active (color->get_active ());
    colorden->set_active (color->get_active ());
    dirpyrden->set_active (color->get_active ());

    colormixerConn.block (false);
    colorshiftConn.block (false);
    colorboostConn.block (false);
    hsveqConn.block (false);
    colordenConn.block (false);
    dirpyrdenConn.block (false);
}

void PartialPasteDlg::lensToggled () {

    distortionConn.block (true);
    cacorrConn.block (true);
    vignettingConn.block (true);

    lens->set_inconsistent (false);

    distortion->set_active (lens->get_active ());
    cacorr->set_active (lens->get_active ());
    vignetting->set_active (lens->get_active ());

    distortionConn.block (false);
    cacorrConn.block (false);
    vignettingConn.block (false);
}

void PartialPasteDlg::compositionToggled () {

    coarserotConn.block (true);
    finerotConn.block (true);
    cropConn.block (true);
    resizeConn.block (true);

    composition->set_inconsistent (false);

    coarserot->set_active (composition->get_active ());
    finerot->set_active (composition->get_active ());
    crop->set_active (composition->get_active ());
    resize->set_active (composition->get_active ());

    coarserotConn.block (false);
    finerotConn.block (false);
    cropConn.block (false);
    resizeConn.block (false);
}

void PartialPasteDlg::metaicmToggled () {

    exifchConn.block (true);
    iptcConn.block (true);
    icmConn.block (true);

    metaicm->set_inconsistent (false);

    exifch->set_active (metaicm->get_active ());
    iptc->set_active (metaicm->get_active ());
    icm->set_active (metaicm->get_active ());

    exifchConn.block (false);
    iptcConn.block (false);
    icmConn.block (false);
}


void PartialPasteDlg::applyPaste (rtengine::procparams::ProcParams* dst, const rtengine::procparams::ProcParams* src) {

    if (wb->get_active ())          dst->wb = src->wb;
    if (exposure->get_active ())    dst->toneCurve = src->toneCurve;
    if (hlrec->get_active ())       dst->hlrecovery = src->hlrecovery;

    if (sharpen->get_active ())     dst->sharpening = src->sharpening;
    if (impden->get_active ())		dst->impulseDenoise = src->impulseDenoise;
    if (lumaden->get_active ())     dst->lumaDenoise = src->lumaDenoise;
    if (labcurve->get_active ())   dst->labCurve = src->labCurve;
    if (sh->get_active ())          dst->sh = src->sh;
    if (dirpyreq->get_active ())	dst->dirpyrequalizer = src->dirpyrequalizer;
    if (waveq->get_active ())		dst->equalizer = src->equalizer;

    if (colormixer->get_active ())  dst->chmixer = src->chmixer;
    if (colorshift->get_active ())  dst->colorShift = src->colorShift;
    if (colorboost->get_active ())  dst->colorBoost = src->colorBoost;
    if (hsveq->get_active ())		dst->hsvequalizer = src->hsvequalizer;
    if (colorden->get_active ())    dst->colorDenoise = src->colorDenoise;
    if (dirpyrden->get_active ())   dst->dirpyrDenoise = src->dirpyrDenoise;

    if (distortion->get_active ())  dst->distortion = src->distortion;
    if (cacorr->get_active ())      dst->cacorrection = src->cacorrection;
    if (vignetting->get_active ())  dst->vignetting = src->vignetting;

    if (coarserot->get_active ())   dst->coarse = src->coarse;
    if (finerot->get_active ())     dst->rotate = src->rotate;
    if (crop->get_active ())        dst->crop = src->crop;
    if (resize->get_active ())      dst->resize = src->resize;

    if (exifch->get_active ())      dst->exif = src->exif;
    if (iptc->get_active ())        dst->iptc = src->iptc;
    if (icm->get_active ())         dst->icm = src->icm; 
                                                                                          
    if (df_file->get_active ())        dst->raw.dark_frame = src->raw.dark_frame;
    if (df_AutoSelect->get_active ())  dst->raw.df_autoselect = src->raw.df_autoselect;
    if (ff_file->get_active ())        dst->raw.ff_file = src->raw.ff_file;
    if (ff_AutoSelect->get_active ())  dst->raw.ff_AutoSelect = src->raw.ff_AutoSelect;   
		if (ff_BlurRadius->get_active ())  dst->raw.ff_BlurRadius = src->raw.ff_BlurRadius; 
		if (ff_BlurType->get_active ())    dst->raw.ff_BlurType = src->raw.ff_BlurType;
}

