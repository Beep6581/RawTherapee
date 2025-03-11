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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>

namespace rtengine
{
namespace procparams
{

class ProcParams;


}

}

struct ParamsEdited;

/* ==== PartialSpotWidgetListener ==== */
class PartialSpotWidget;
class PartialSpotWidgetListener
{
public:
    enum UpdateStatus {
        AllSelection = 1,
        NoSelection = 2,
        PartialSelection = 3
    };

public:
    PartialSpotWidgetListener() {};
    virtual ~PartialSpotWidgetListener() {};

    virtual void partialSpotUpdated(const UpdateStatus status) = 0;
};

/* ==== PartialSpotWidget ==== */
class PartialSpotWidget:
    public Gtk::Box
{
private:
    // Tree model to manage spot selection widget
    class SpotRow:
        public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<bool> keep;
        Gtk::TreeModelColumn<Glib::ustring> spotname;

        SpotRow()
        {
            add(keep);
            add(spotname);
        }
    };

    // Spot selection widgets
    Gtk::TreeView* const treeview;
    sigc::connection treeviewconn;
    SpotRow spotRow;
    Glib::RefPtr<Gtk::ListStore> treemodel;

    // Spot selection listener
    PartialSpotWidgetListener* selListener;

public:
    PartialSpotWidget();

    // Setter for spot selection listener
    void setPartialSpotWidgetListener(PartialSpotWidgetListener* pswl)
    {
        selListener = pswl;
    }

    // Spot selection widget management functions
    void updateSpotWidget(const rtengine::procparams::ProcParams* pp, const bool defValue);
    void enableAll();
    void disableAll();
    std::vector<bool> getSelectionStatus();

private:
    // GUI aspect management functions
    void render_keep(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);
    void render_spotname(Gtk::CellRenderer* cell, const Gtk::TreeModel::iterator& iter);

    // Event management function
    void keepToggled(const Glib::ustring &path);
};

/* ==== PartialPasteDlg ==== */
class PartialPasteDlg final:
    public Gtk::Dialog,
    public PartialSpotWidgetListener
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
    Gtk::CheckButton* locallab;

    // options in basic:
    Gtk::CheckButton* wb;
    Gtk::CheckButton* exposure;
    Gtk::CheckButton* localcontrast;
    Gtk::CheckButton* sh;
    Gtk::CheckButton* toneEqualizer;
    Gtk::CheckButton* epd;
    Gtk::CheckButton* fattal;
    Gtk::CheckButton* retinex;
    Gtk::CheckButton* pcvignette;
    Gtk::CheckButton* gradient;
    Gtk::CheckButton* labcurve;
    Gtk::CheckButton* colorappearance;

    // options in detail:
    Gtk::CheckButton* spot;
    Gtk::CheckButton* sharpen;
    Gtk::CheckButton* sharpenedge;
    Gtk::CheckButton* sharpenmicro;
    Gtk::CheckButton* impden;
    //Gtk::CheckButton* waveq;
    Gtk::CheckButton* dirpyrden;
    Gtk::CheckButton* defringe;
    Gtk::CheckButton* dirpyreq;
    Gtk::CheckButton* dehaze;

    // options in wavelet
    Gtk::CheckButton* wavelet;

    // options in color:
    Gtk::CheckButton* compressGamut;
    Gtk::CheckButton* icm;
    Gtk::CheckButton* vibrance;
    Gtk::CheckButton* chmixer;
    Gtk::CheckButton* blackwhite;
    Gtk::CheckButton* hsveq;
    Gtk::CheckButton* softlight;
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
    Gtk::CheckButton* framing;
    Gtk::CheckButton* perspective;
    Gtk::CheckButton* commonTrans;

    // options in meta:
    Gtk::CheckButton *metadata;
    Gtk::CheckButton* exifch;
    Gtk::CheckButton* iptc;

    // options in locallab:
    PartialSpotWidget* spots;

    // options in raw:
    Gtk::CheckButton* raw_expos;
    Gtk::CheckButton* raw_black;
    Gtk::CheckButton* raw_ca_autocorrect;
    Gtk::CheckButton* raw_caredblue;
    Gtk::CheckButton* raw_ca_avoid_colourshift;
    Gtk::CheckButton* raw_hotpix_filt;
    Gtk::CheckButton* raw_deadpix_filt;
    Gtk::CheckButton* raw_pdaf_lines_filter;
    Gtk::CheckButton* raw_linenoise;
    Gtk::CheckButton* raw_greenthresh;
    Gtk::CheckButton* raw_method;
    Gtk::CheckButton* raw_border;
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
    Gtk::CheckButton* ff_FromMetaData;
    Gtk::CheckButton* ff_BlurRadius;
    Gtk::CheckButton* ff_BlurType;
    Gtk::CheckButton* ff_ClipControl;

    Gtk::CheckButton* filmNegative;
    Gtk::CheckButton* captureSharpening;
    Gtk::CheckButton* raw_preprocwb;

    sigc::connection everythingConn, basicConn, detailConn, colorConn, lensConn, compositionConn, metaConn, rawConn, advancedConn;
    sigc::connection locallabConn;
    sigc::connection wbConn, exposureConn, localcontrastConn, shConn, pcvignetteConn, gradientConn, labcurveConn, colorappearanceConn;
    sigc::connection toneEqualizerConn;
    sigc::connection spotConn, sharpenConn, gradsharpenConn, microcontrastConn, impdenConn, dirpyrdenConn, defringeConn, epdConn, fattalConn, dirpyreqConn, waveletConn, retinexConn, dehazeConn;
    sigc::connection compressGamutConn, vibranceConn, chmixerConn, hsveqConn, rgbcurvesConn, chmixerbwConn, colortoningConn, filmSimulationConn, softlightConn;
    sigc::connection distortionConn, cacorrConn, vignettingConn, lcpConn;
    sigc::connection coarserotConn, finerotConn, cropConn, resizeConn, prsharpeningConn, framingConn, perspectiveConn, commonTransConn;
    sigc::connection metadataConn, exifchConn, iptcConn, icmConn;
    sigc::connection df_fileConn, df_AutoSelectConn, ff_fileConn, ff_AutoSelectConn, ff_FromMetaDataConn, ff_BlurRadiusConn, ff_BlurTypeConn, ff_ClipControlConn;
    sigc::connection raw_caredblueConn, raw_ca_autocorrectConn, raw_ca_avoid_colourshiftconn, raw_hotpix_filtConn, raw_deadpix_filtConn, raw_pdaf_lines_filterConn, raw_linenoiseConn, raw_greenthreshConn, raw_ccStepsConn, raw_methodConn, raw_borderConn, raw_imagenumConn, raw_dcb_iterationsConn, raw_lmmse_iterationsConn, raw_pixelshiftConn, raw_dcb_enhanceConn, raw_exposConn, raw_blackConn;
    sigc::connection filmNegativeConn;
    sigc::connection captureSharpeningConn;
    sigc::connection raw_preprocwbConn;

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
    void locallabToggled ();

    void updateSpotWidget(const rtengine::procparams::ProcParams* pp);
    void partialSpotUpdated(const UpdateStatus status);
};
