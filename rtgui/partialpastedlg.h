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
#ifndef _PARTIALPASTEDLG_
#define _PARTIALPASTEDLG_

#include <gtkmm.h>
#include "../rtengine/rtengine.h"

class PartialPasteDlg : public Gtk::Dialog {

    public:

        Gtk::CheckButton* everything;

        // main groups:
        Gtk::CheckButton* basic;
        Gtk::CheckButton* detail;
        Gtk::CheckButton* color;
        Gtk::CheckButton* lens;
        Gtk::CheckButton* composition;
        Gtk::CheckButton* metaicm;
        Gtk::CheckButton* raw;

        // options in basic:
        Gtk::CheckButton* wb;
        Gtk::CheckButton* exposure;
        Gtk::CheckButton* sh;
        Gtk::CheckButton* epd;
        Gtk::CheckButton* pcvignette;
        Gtk::CheckButton* gradient;
        Gtk::CheckButton* labcurve;
        Gtk::CheckButton* colorappearance;

        // options in detail:
        Gtk::CheckButton* sharpen;
        Gtk::CheckButton* sharpenedge;
        Gtk::CheckButton* sharpenmicro;
        Gtk::CheckButton* impden;
        //Gtk::CheckButton* waveq;
        Gtk::CheckButton* dirpyrden;
        Gtk::CheckButton* defringe;
        Gtk::CheckButton* dirpyreq;

        // options in color:
        Gtk::CheckButton* vibrance;
        Gtk::CheckButton* chmixer;
        Gtk::CheckButton* blackwhite;
        Gtk::CheckButton* hsveq;
        Gtk::CheckButton* rgbcurves;
        Gtk::CheckButton* colortoning;
       // Gtk::CheckButton* icm;

        // options in lens:
        Gtk::CheckButton* distortion;
        Gtk::CheckButton* cacorr;
        Gtk::CheckButton* vignetting;
        Gtk::CheckButton* lcp;

        // options in composition:
        Gtk::CheckButton* coarserot;
        Gtk::CheckButton* finerot;
        Gtk::CheckButton* crop;
        Gtk::CheckButton* resize;
        Gtk::CheckButton* perspective;
        Gtk::CheckButton* commonTrans;

        // options in metaicm:
        Gtk::CheckButton* exifch;
        Gtk::CheckButton* iptc;
        Gtk::CheckButton* icm;
        Gtk::CheckButton* gam;

        // options in raw:
        Gtk::CheckButton* raw_expos;
        Gtk::CheckButton* raw_preser;
        Gtk::CheckButton* raw_black;
        Gtk::CheckButton* raw_ca_autocorrect;
        Gtk::CheckButton* raw_cared;
        Gtk::CheckButton* raw_cablue;
        Gtk::CheckButton* raw_hotdeadpix_filt;
        Gtk::CheckButton* raw_linenoise;
        Gtk::CheckButton* raw_greenthresh;
        Gtk::CheckButton* raw_dmethod;
        Gtk::CheckButton* raw_ccSteps;
        Gtk::CheckButton* raw_dcb_iterations;
        Gtk::CheckButton* raw_dcb_enhance;
        //Gtk::CheckButton* raw_all_enhance;
        Gtk::CheckButton* raw_lmmse_iterations;

        Gtk::CheckButton* df_file;
        Gtk::CheckButton* df_AutoSelect;
        Gtk::CheckButton* ff_file;
        Gtk::CheckButton* ff_AutoSelect;
        Gtk::CheckButton* ff_BlurRadius;
        Gtk::CheckButton* ff_BlurType;

        sigc::connection everythingConn, basicConn, detailConn, colorConn, lensConn, compositionConn, metaicmConn, rawConn;

        sigc::connection wbConn, exposureConn, shConn, pcvignetteConn, gradientConn, labcurveConn, colorappearanceConn;
        sigc::connection sharpenConn, gradsharpenConn, microcontrastConn, impdenConn, dirpyrdenConn, waveqConn, defringeConn, epdConn, dirpyreqConn;
        sigc::connection vibranceConn, chmixerConn, hsveqConn, rgbcurvesConn, chmixerbwConn, colortoningConn;
        sigc::connection distortionConn, cacorrConn, vignettingConn, lcpConn;
        sigc::connection coarserotConn, finerotConn, cropConn, resizeConn, perspectiveConn, commonTransConn;
        sigc::connection exifchConn, iptcConn, icmConn, gamcsconn;
        sigc::connection df_fileConn, df_AutoSelectConn, ff_fileConn, ff_AutoSelectConn, ff_BlurRadiusConn, ff_BlurTypeConn;
        sigc::connection raw_caredConn, raw_cablueConn, raw_ca_autocorrectConn, raw_hotdeadpix_filtConn, raw_linenoiseConn, raw_greenthreshConn, raw_ccStepsConn, raw_dmethodConn, raw_dcb_iterationsConn, raw_lmmse_iterationsConn, raw_dcb_enhanceConn, raw_exposConn, raw_preserConn, raw_blackConn; //,raw_all_enhanceConn

    public:
        PartialPasteDlg (Glib::ustring title);

        void applyPaste (rtengine::procparams::ProcParams* dstPP, ParamsEdited* dstPE, const rtengine::procparams::ProcParams* srcPP, const ParamsEdited* srcPE=NULL);

        void everythingToggled ();
        void basicToggled ();
        void detailToggled ();
        void colorToggled ();
        void lensToggled ();
        void compositionToggled ();
        void metaicmToggled ();
        void rawToggled ();
};

#endif

