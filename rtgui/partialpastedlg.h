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

class PartialPasteDlg : public Gtk::Dialog
{

public:

    Gtk::ScrolledWindow *scrolledwindow;

    Gtk::CheckButton* everything;

    // main groups:
    Gtk::CheckButton* basic;
    Gtk::CheckButton* detail;
    Gtk::CheckButton* color;
    Gtk::CheckButton* lens;
    Gtk::CheckButton* composition;
    Gtk::CheckButton* meta;
    Gtk::CheckButton* raw;
    Gtk::CheckButton* advanced;

    // options in basic:
    Gtk::CheckButton* wb;
    Gtk::CheckButton* exposure;
    Gtk::CheckButton* localcontrast;
    Gtk::CheckButton* sh;
    Gtk::CheckButton* epd;
    Gtk::CheckButton* fattal;
    Gtk::CheckButton* retinex;
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

    // options in wavelet
    Gtk::CheckButton* wavelet;

    // options in color:
    Gtk::CheckButton* icm;
    Gtk::CheckButton* vibrance;
    Gtk::CheckButton* chmixer;
    Gtk::CheckButton* blackwhite;
    Gtk::CheckButton* hsveq;
    Gtk::CheckButton* filmSimulation;
    Gtk::CheckButton* rgbcurves;
    Gtk::CheckButton* colortoning;

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
    Gtk::CheckButton* prsharpening;
    Gtk::CheckButton* perspective;
    Gtk::CheckButton* commonTrans;

    // options in meta:
    Gtk::CheckButton *metadata;
    Gtk::CheckButton* exifch;
    Gtk::CheckButton* iptc;


    // options in raw:
    Gtk::CheckButton* raw_expos;
    Gtk::CheckButton* raw_preser;
    Gtk::CheckButton* raw_black;
    Gtk::CheckButton* raw_ca_autocorrect;
    Gtk::CheckButton* raw_caredblue;
    Gtk::CheckButton* raw_hotpix_filt;
    Gtk::CheckButton* raw_deadpix_filt;
    Gtk::CheckButton* raw_linenoise;
    Gtk::CheckButton* raw_greenthresh;
    Gtk::CheckButton* raw_method;
    Gtk::CheckButton* raw_imagenum;
    Gtk::CheckButton* raw_ccSteps;
    Gtk::CheckButton* raw_dcb_iterations;
    Gtk::CheckButton* raw_dcb_enhance;
    Gtk::CheckButton* raw_lmmse_iterations;
    Gtk::CheckButton* raw_pixelshift;

    Gtk::CheckButton* df_file;
    Gtk::CheckButton* df_AutoSelect;
    Gtk::CheckButton* ff_file;
    Gtk::CheckButton* ff_AutoSelect;
    Gtk::CheckButton* ff_BlurRadius;
    Gtk::CheckButton* ff_BlurType;
    Gtk::CheckButton* ff_ClipControl;

    sigc::connection everythingConn, basicConn, detailConn, colorConn, lensConn, compositionConn, metaConn, rawConn, advancedConn;

    sigc::connection wbConn, exposureConn, localcontrastConn, shConn, pcvignetteConn, gradientConn, labcurveConn, colorappearanceConn;
    sigc::connection sharpenConn, gradsharpenConn, microcontrastConn, impdenConn, dirpyrdenConn, defringeConn, epdConn, fattalConn, dirpyreqConn, waveletConn, retinexConn;
    sigc::connection vibranceConn, chmixerConn, hsveqConn, rgbcurvesConn, chmixerbwConn, colortoningConn, filmSimulationConn;
    sigc::connection distortionConn, cacorrConn, vignettingConn, lcpConn;
    sigc::connection coarserotConn, finerotConn, cropConn, resizeConn, prsharpeningConn, perspectiveConn, commonTransConn;
    sigc::connection metadataConn, exifchConn, iptcConn, icmConn;
    sigc::connection df_fileConn, df_AutoSelectConn, ff_fileConn, ff_AutoSelectConn, ff_BlurRadiusConn, ff_BlurTypeConn, ff_ClipControlConn;
    sigc::connection raw_caredblueConn, raw_ca_autocorrectConn, raw_hotpix_filtConn, raw_deadpix_filtConn, raw_linenoiseConn, raw_greenthreshConn, raw_ccStepsConn, raw_methodConn, raw_imagenumConn, raw_dcb_iterationsConn, raw_lmmse_iterationsConn, raw_pixelshiftConn, raw_dcb_enhanceConn, raw_exposConn, raw_preserConn, raw_blackConn;

public:
    PartialPasteDlg (const Glib::ustring &title, Gtk::Window* parent);

    void applyPaste (rtengine::procparams::ProcParams* dstPP, ParamsEdited* dstPE, const rtengine::procparams::ProcParams* srcPP, const ParamsEdited* srcPE = nullptr);

    void everythingToggled ();
    void basicToggled ();
    void detailToggled ();
    void colorToggled ();
    void lensToggled ();
    void compositionToggled ();
    void metaToggled ();
    void rawToggled ();
    void advancedToggled ();
};

#endif

