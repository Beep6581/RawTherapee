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

    basic       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_BASICGROUP")));
    luminance   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LUMINANCEGROUP")));
    color       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORGROUP")));
    lens        = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LENSGROUP")));
    composition = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COMPOSITIONGROUP")));
    metaicm     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_METAICMGROUP")));

    // options in basic:
    wb          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_WHITEBALANCE")));
    exposure    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_EXPOSURE")));
    hlrec       = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_HLRECOVERY")));

    // options in luminance:
    sharpen     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHARPENING")));
    lumaden     = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LUMADENOISE")));
    lumacurve   = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_LUMACURVE")));
    sh          = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_SHADOWSHIGHLIGHTS")));

    // options in color:
    colormixer  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORMIXER")));
    colorshift  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORSHIFT")));
    colorboost  = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORBOOST")));
    colorden    = Gtk::manage (new Gtk::CheckButton (M("PARTIALPASTE_COLORDENOISE")));

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

    Gtk::VBox* vboxes[6];
    Gtk::HSeparator* hseps[6];
    for (int i=0; i<6; i++) {
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
    vboxes[1]->pack_start (*lumaden, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*lumacurve, Gtk::PACK_SHRINK, 2);
    vboxes[1]->pack_start (*sh, Gtk::PACK_SHRINK, 2);

    vboxes[2]->pack_start (*color, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*hseps[2], Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colormixer, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colorshift, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colorboost, Gtk::PACK_SHRINK, 2);
    vboxes[2]->pack_start (*colorden, Gtk::PACK_SHRINK, 2);


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

    vboxes[5]->pack_start (*metaicm, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*hseps[5], Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*exifch, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*iptc, Gtk::PACK_SHRINK, 2);
    vboxes[5]->pack_start (*icm, Gtk::PACK_SHRINK, 2);

    Gtk::VBox* vbleft = Gtk::manage (new Gtk::VBox ());
    Gtk::VBox* vbright = Gtk::manage (new Gtk::VBox ());

    vbleft->set_border_width (16);
    vbright->set_border_width (16);

    for (int i=0; i<3; i++)
        vbleft->pack_start (*vboxes[i]);
    for (int i=3; i<6; i++)
        vbright->pack_start (*vboxes[i]);

    Gtk::HBox* hbmain = Gtk::manage (new Gtk::HBox ());
    hbmain->pack_start (*vbleft);
    hbmain->pack_start (*(Gtk::manage (new Gtk::VSeparator ())));
    hbmain->pack_start (*vbright);

    get_vbox()->pack_start (*hbmain);

    basicConn       = basic->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::basicToggled));    
    luminanceConn   = luminance->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::luminanceToggled));    
    colorConn       = color->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::colorToggled));    
    lensConn        = lens->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::lensToggled));    
    compositionConn = composition->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::compositionToggled));    
    metaicmConn     = metaicm->signal_toggled().connect (sigc::mem_fun(*this, &PartialPasteDlg::metaicmToggled));    

    wbConn          = wb->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    exposureConn    = exposure->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    
    hlrecConn       = hlrec->signal_toggled().connect (sigc::bind (sigc::mem_fun(*basic, &Gtk::CheckButton::set_inconsistent), true));    

    sharpenConn     = sharpen->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    lumadenConn     = lumaden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    lumacurveConn   = lumacurve->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    
    shConn          = sh->signal_toggled().connect (sigc::bind (sigc::mem_fun(*luminance, &Gtk::CheckButton::set_inconsistent), true));    

    colormixerConn  = colormixer->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    colorshiftConn  = colorshift->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    colorboostConn  = colorboost->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    
    colordenConn    = colorden->signal_toggled().connect (sigc::bind (sigc::mem_fun(*color, &Gtk::CheckButton::set_inconsistent), true));    

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

    add_button (Gtk::StockID("gtk-ok"), 1);
    add_button (Gtk::StockID("gtk-cancel"), 0);
    set_response_sensitive (1);
    set_default_response (1);
    show_all_children ();
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
    lumadenConn.block (true);
    lumacurveConn.block (true);
    shConn.block (true);

    luminance->set_inconsistent (false);

    sharpen->set_active (luminance->get_active ());
    lumaden->set_active (luminance->get_active ());
    lumacurve->set_active (luminance->get_active ());
    sh->set_active (luminance->get_active ());

    sharpenConn.block (false);
    lumadenConn.block (false);
    lumacurveConn.block (false);
    shConn.block (false);
}

void PartialPasteDlg::colorToggled () {

    colormixerConn.block (true);
    colorshiftConn.block (true);
    colorboostConn.block (true);
    colordenConn.block (true);

    color->set_inconsistent (false);

    colormixer->set_active (color->get_active ());
    colorshift->set_active (color->get_active ());
    colorboost->set_active (color->get_active ());
    colorden->set_active (color->get_active ());

    colormixerConn.block (false);
    colorshiftConn.block (false);
    colorboostConn.block (false);
    colordenConn.block (false);
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
    if (lumaden->get_active ())     dst->lumaDenoise = src->lumaDenoise;
    if (lumacurve->get_active ())   dst->lumaCurve = src->lumaCurve;
    if (sh->get_active ())          dst->sh = src->sh;

    if (colormixer->get_active ())  dst->chmixer = src->chmixer;
    if (colorshift->get_active ())  dst->colorShift = src->colorShift;
    if (colorboost->get_active ())  dst->colorBoost = src->colorBoost;
    if (colorden->get_active ())    dst->colorDenoise = src->colorDenoise;

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
}

