/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
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
 *  2017 Jacques Desmis <jdesmis@gmail.com>
 *  2018 Pierre Cabrera <pierre.cab@gmail.com>
 */
#pragma once

#include <array>

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "editcallbacks.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "toolpanel.h"
#include "options.h"
#include "thresholdadjuster.h"
#include "controlspotpanel.h"
#include "labgrid.h"

class Locallab :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public CurveListener,
    public ColorProvider,
    public ThresholdCurveProvider,
    public rtengine::LocallabListener,
    public ThresholdAdjusterListener

{
private:
    IdleRegister idle_register;

    // Expander widgets
    ControlSpotPanel* const expsettings;
    MyExpander* const expcolor;
    MyExpander* const expexpose;
    MyExpander* const expshadhigh;
    MyExpander* const expvibrance;
    MyExpander* const expsoft;
    MyExpander* const expblur;
    MyExpander* const exptonemap;
    MyExpander* const expreti;
    MyExpander* const expretitools;
    MyExpander* const expsharp;
    MyExpander* const expcontrastpyr;
    MyExpander* const expcontrast;
    MyExpander* const expcbdl;
    MyExpander* const expdenoi;
    MyExpander* const explog;
    MyExpander* const expmaskcol;
    MyExpander* const expmaskcol1;
    MyExpander* const expcurvcol;
    MyExpander* const expmaskexp;
    MyExpander* const expmasksh;
    MyExpander* const expmasklc;
    MyExpander* const expmaskcb;
    MyExpander* const expmaskreti;
    MyExpander* const expmasktm;
    MyExpander* const expmaskbl;
    MyExpander* const expmaskvib;
    MyExpander* const expgradexp;
    MyExpander* const exptoolexp;
    MyExpander* const expgradsh;
    MyExpander* const exptrcsh;
    MyExpander* const expgradvib;
    MyExpander* const expgradcol;

    sigc::connection enablecolorConn, enableexposeConn, enableshadhighConn, enablevibranceConn, enablesoftConn, enableblurConn, enabletonemapConn, enableretiConn, enablesharpConn, enablecontrastConn, enablecbdlConn, enabledenoiConn, enablelogConn;

    // Curve widgets
    // Color & Light
    CurveEditorGroup* const llCurveEditorG;
    CurveEditorGroup* const clCurveEditorG;
    CurveEditorGroup* const HCurveEditorG;
    CurveEditorGroup* const H2CurveEditorG;
    CurveEditorGroup* const rgbCurveEditorG;
    CurveEditorGroup* const maskCurveEditorG;
    CurveEditorGroup* const maskHCurveEditorG;
    CurveEditorGroup* const mask2CurveEditorG;
    CurveEditorGroup* const mask2CurveEditorGwav;
    DiagonalCurveEditor* const Lmaskshape;
    DiagonalCurveEditor* const llshape;
    DiagonalCurveEditor* const ccshape;
    DiagonalCurveEditor* const clshape;
    DiagonalCurveEditor* const lcshape;
    MyComboBoxText* const toneMethod;//put here to use toneMethod in rgbshape
    DiagonalCurveEditor* const rgbshape;
    FlatCurveEditor* const LHshape;
    FlatCurveEditor* const HHshape;
    FlatCurveEditor* const CCmaskshape;
    FlatCurveEditor* const LLmaskshape;
    FlatCurveEditor* const HHmaskshape;
    FlatCurveEditor* const LLmaskcolshapewav;
    FlatCurveEditor* const HHhmaskshape;
    // Exposure
    CurveEditorGroup* const curveEditorG;
    CurveEditorGroup* const maskexpCurveEditorG;
    CurveEditorGroup* const mask2expCurveEditorG;
    DiagonalCurveEditor* const Lmaskexpshape;
    DiagonalCurveEditor* const shapeexpos;
    FlatCurveEditor* const CCmaskexpshape;
    FlatCurveEditor* const LLmaskexpshape;
    FlatCurveEditor* const HHmaskexpshape;
    //Shadows Highlight
    CurveEditorGroup* const maskSHCurveEditorG;
    CurveEditorGroup* const mask2SHCurveEditorG;
    DiagonalCurveEditor* const LmaskSHshape;
    FlatCurveEditor* const CCmaskSHshape;
    FlatCurveEditor* const LLmaskSHshape;
    FlatCurveEditor* const HHmaskSHshape;
    // Vibrance
    CurveEditorGroup* const curveEditorGG;
    CurveEditorGroup* const maskvibCurveEditorG;
    DiagonalCurveEditor* const skinTonesCurve;
    CurveEditorGroup* const mask2vibCurveEditorG;
    DiagonalCurveEditor* const Lmaskvibshape;
    FlatCurveEditor* const CCmaskvibshape;
    FlatCurveEditor* const LLmaskvibshape;
    FlatCurveEditor* const HHmaskvibshape;
   
    //Blur and noise
    CurveEditorGroup* const maskblCurveEditorG;
    CurveEditorGroup* const mask2blCurveEditorG;
    CurveEditorGroup* const mask2blCurveEditorGwav;
    DiagonalCurveEditor* const Lmaskblshape;
    FlatCurveEditor* const CCmaskblshape;
    FlatCurveEditor* const LLmaskblshape;
    FlatCurveEditor* const HHmaskblshape;
    FlatCurveEditor* const LLmaskblshapewav;
    // TM
    CurveEditorGroup* const masktmCurveEditorG;
    CurveEditorGroup* const mask2tmCurveEditorG;
    DiagonalCurveEditor* const Lmasktmshape;
    FlatCurveEditor* const CCmasktmshape;
    FlatCurveEditor* const LLmasktmshape;
    FlatCurveEditor* const HHmasktmshape;
    // Retinex
    CurveEditorGroup* const LocalcurveEditortransT;
    CurveEditorGroup* const LocalcurveEditorgainT;
    CurveEditorGroup* const maskretiCurveEditorG;
    CurveEditorGroup* const mask2retiCurveEditorG;
    DiagonalCurveEditor* const Lmaskretishape;
    FlatCurveEditor* const cTtransshape;
    FlatCurveEditor* const cTgainshape;
    FlatCurveEditor* const CCmaskretishape;
    FlatCurveEditor* const LLmaskretishape;
    FlatCurveEditor* const HHmaskretishape;
    //local contrast
    CurveEditorGroup* const LocalcurveEditorwav;
    FlatCurveEditor* const wavshape;
    CurveEditorGroup* const LocalcurveEditorwavlev;
    FlatCurveEditor* const wavshapelev;
    CurveEditorGroup* const LocalcurveEditorwavcon;
    FlatCurveEditor* const wavshapecon;
    CurveEditorGroup* const LocalcurveEditorwavcomp;
    FlatCurveEditor* const wavshapecomp;
    CurveEditorGroup* const LocalcurveEditorwavcompre;
    FlatCurveEditor* const wavshapecompre;
    CurveEditorGroup* const masklcCurveEditorG;
    CurveEditorGroup* const mask2lcCurveEditorG;
    DiagonalCurveEditor* const Lmasklcshape;
    FlatCurveEditor* const CCmasklcshape;
    FlatCurveEditor* const LLmasklcshape;
    FlatCurveEditor* const HHmasklcshape;
    //Cbdl
    CurveEditorGroup* const maskcbCurveEditorG;
    CurveEditorGroup* const mask2cbCurveEditorG;
    DiagonalCurveEditor* const Lmaskcbshape;
    FlatCurveEditor* const CCmaskcbshape;
    FlatCurveEditor* const LLmaskcbshape;
    FlatCurveEditor* const HHmaskcbshape;

    //Denoise
    CurveEditorGroup* const LocalcurveEditorwavden;
    FlatCurveEditor* const wavshapeden;

    // Adjuster widgets
    // Color & Light
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Adjuster* const strengthgrid;
    Adjuster* const sensi;
    Adjuster* const structcol;
    Adjuster* const blurcolde;
    Adjuster* const strcol;
    Adjuster* const strcolab;
    Adjuster* const strcolh;
    Adjuster* const angcol;
    Adjuster* const blendmaskcol;
    Adjuster* const radmaskcol;
    Adjuster* const chromaskcol;
    Adjuster* const gammaskcol;
    Adjuster* const slomaskcol;
    Adjuster* const lapmaskcol;
    Adjuster* const shadmaskcol;
    Adjuster* const softradiuscol;
    Adjuster* const opacol;
    Adjuster* const conthrcol;
    Adjuster* const strumaskcol;
    Adjuster* const mercol;
    Adjuster* const merlucol;
    Adjuster* const blurcol;
    Adjuster* const contcol;
    // Exposure
    Adjuster* const expcomp;
    Adjuster* const hlcompr;
    Adjuster* const hlcomprthresh;
    Adjuster* const black;
    Adjuster* const shadex;
    Adjuster* const shcompr;
    Adjuster* const expchroma;
//    Adjuster* const warm;
    Adjuster* const sensiex;
    Adjuster* const structexp;
    Adjuster* const blurexpde;
    Adjuster* const strexp;
    Adjuster* const angexp;
    Adjuster* const blendmaskexp;
    Adjuster* const radmaskexp;
    Adjuster* const chromaskexp;
    Adjuster* const gammaskexp;
    Adjuster* const slomaskexp;
    Adjuster* const lapmaskexp;
    Adjuster* const softradiusexp;
    Adjuster* const laplacexp;
    Adjuster* const balanexp;
    Adjuster* const linear;
    Adjuster* const gamm;
    Adjuster* const fatamount;
    Adjuster* const fatdetail;
    Adjuster* const fatanchor;
    Adjuster* const strmaskexp;
    Adjuster* const angmaskexp;
    Adjuster* const fatlevel;
    //Shadow highlight
    const std::array<Adjuster*, 5> multipliersh;
    Adjuster* const highlights;
    Adjuster* const h_tonalwidth;
    Adjuster* const shadows;
    Adjuster* const s_tonalwidth;
    Adjuster* const sh_radius;
    Adjuster* const sensihs;
    Adjuster* const blendmaskSH;
    Adjuster* const radmaskSH;
    Adjuster* const blurSHde;
    Adjuster* const strSH;
    Adjuster* const angSH;
    Adjuster* const chromaskSH;
    Adjuster* const gammaskSH;
    Adjuster* const slomaskSH;
    Adjuster* const lapmaskSH;
    Adjuster* const detailSH;
    Adjuster* const fatamountSH;
    Adjuster* const fatanchorSH;
    Adjuster* const gamSH;
    Adjuster* const sloSH;
    // Vibrance
    Adjuster* const saturated;
    Adjuster* const pastels;
    Adjuster* const warm;
    Adjuster* const sensiv;
    Adjuster* const blendmaskvib;
    Adjuster* const radmaskvib;
    Adjuster* const chromaskvib;
    Adjuster* const gammaskvib;
    Adjuster* const slomaskvib;
    Adjuster* const lapmaskvib;
    Adjuster* const strvib;
    Adjuster* const strvibab;
    Adjuster* const strvibh;
    Adjuster* const angvib;
    // Soft Light
    Adjuster* const streng;
    Adjuster* const laplace;
    Adjuster* const sensisf;
    // Blur & Noise
    Adjuster* const radius;
    Adjuster* const strength;
    Adjuster* const itera;
    Adjuster* const guidbl;
    Adjuster* const epsbl;
    Adjuster* const sensibn;
    Adjuster* const blendmaskbl;
    Adjuster* const radmaskbl;
    Adjuster* const chromaskbl;
    Adjuster* const gammaskbl;
    Adjuster* const slomaskbl;
    Adjuster* const lapmaskbl;
    Adjuster* const shadmaskbl;
    Adjuster* const isogr;
    Adjuster* const strengr;
    Adjuster* const scalegr;
    Adjuster* const strumaskbl;
    // Tone Mapping
    Adjuster* const stren;
    Adjuster* const gamma;
    Adjuster* const estop;
    Adjuster* const scaltm;
    Adjuster* const rewei;
    Adjuster* const sensitm;
    Adjuster* const softradiustm;
    Adjuster* const amount;
    Adjuster* const satur;
    Adjuster* const blendmasktm;
    Adjuster* const radmasktm;
    Adjuster* const chromasktm;
    Adjuster* const gammasktm;
    Adjuster* const slomasktm;
    Adjuster* const lapmasktm;
    // Retinex
    Adjuster* const str;
    Adjuster* const chrrt;
    Adjuster* const neigh;
    Adjuster* const vart;
    Adjuster* const offs;
    Adjuster* const dehaz;
    Adjuster* const depth;
    Adjuster* const sensih;
    Adjuster* const softradiusret;

    Adjuster* const blendmaskreti;
    Adjuster* const radmaskreti;
    Adjuster* const chromaskreti;
    Adjuster* const gammaskreti;
    Adjuster* const slomaskreti;
    Adjuster* const lapmaskreti;
    Adjuster* const scalereti;
    Adjuster* const darkness;
    Adjuster* const lightnessreti;
    Adjuster* const limd;
    Adjuster* const cliptm;

    // Sharpening
    Adjuster* const sharcontrast;
    Adjuster* const sharradius;
    Adjuster* const sharamount;
    Adjuster* const shardamping;
    Adjuster* const shariter;
    Adjuster* const sharblur;
    Adjuster* const sensisha;
    // Local Contrast
    Adjuster* const lcradius;
    Adjuster* const lcamount;
    Adjuster* const lcdarkness;
    Adjuster* const lclightness;
    Adjuster* const levelwav;
    Adjuster* const residcont;
    Adjuster* const residblur;
    Adjuster* const levelblur;
    Adjuster* const clarilres;
    Adjuster* const clarisoft;
    Adjuster* const claricres;
    Adjuster* const sensilc;
    Adjuster* const residchro;
    Adjuster* const residcomp;
    Adjuster* const sigma;
    Adjuster* const offset;
    Adjuster* const threswav;
    Adjuster* const chromalev;
    Adjuster* const chromablu;
    Adjuster* const fatdet;
    Adjuster* const fatanch;
    Adjuster* const fatres;
    Adjuster* const blendmasklc;
    Adjuster* const radmasklc;
    Adjuster* const chromasklc;
    Adjuster* const strwav;
    Adjuster* const angwav;
    // Contrast by detail levels
    const std::array<Adjuster*, 6> multiplier;
    Adjuster* const chromacbdl;
    Adjuster* const threshold;
    Adjuster* const clarityml;
    Adjuster* const contresid;
    Adjuster* const blurcbdl;
    Adjuster* const sensicb;
    Adjuster* const softradiuscb;
    Adjuster* const blendmaskcb;
    Adjuster* const radmaskcb;
    Adjuster* const chromaskcb;
    Adjuster* const gammaskcb;
    Adjuster* const slomaskcb;
    Adjuster* const lapmaskcb;

    // Denoise
    Adjuster* const noiselumf;
    Adjuster* const noiselumf0;
    Adjuster* const noiselumf2;
    Adjuster* const noiselumc;
    Adjuster* const noiselumdetail;
    Adjuster* const noiselequal;
    Adjuster* const noisechrof;
    Adjuster* const noisechroc;
    Adjuster* const noisechrodetail;
    Adjuster* const adjblur;
    Adjuster* const bilateral;
    Adjuster* const sensiden;
    Adjuster* const detailthr;

    //Log encoding
    Adjuster* const sourceGray;
    Adjuster* const targetGray;
    Adjuster* const blackEv;
    Adjuster* const whiteEv;
    Adjuster* const detail;
    Adjuster* const sensilog;
    Adjuster* const baselog;
    Adjuster* const strlog;
    Adjuster* const anglog;

    // ButtonCheck widgets
    // Color & Light
    Gtk::CheckButton* const curvactiv;
    Gtk::CheckButton* const invers;
    Gtk::CheckButton* const special;
    Gtk::CheckButton* const toolcol;
    Gtk::CheckButton* const enaColorMask;
    Gtk::CheckButton* const fftColorMask;
    sigc::connection curvactivConn, inversConn, enaColorMaskConn, specialConn, toolcolConn, fftColorMaskConn;
    // Exposure
    Gtk::CheckButton* const enaExpMask;
    sigc::connection enaExpMaskConn;
    Gtk::CheckButton* const enaExpMaskaft;
    sigc::connection enaExpMaskaftConn;
    Gtk::CheckButton* const inversex;
    sigc::connection inversexConn;
    //Shadows highlight
    Gtk::CheckButton* const enaSHMask;
    sigc::connection enaSHMaskConn;
    Gtk::CheckButton* const inverssh;
    sigc::connection inversshConn;
    // Vibrance
    Gtk::CheckButton* const protectSkins;
    Gtk::CheckButton* const avoidColorShift;
    Gtk::CheckButton* const pastSatTog;
    Gtk::CheckButton* const enavibMask;
    sigc::connection enavibMaskConn;
    sigc::connection pskinsconn, ashiftconn, pastsattogconn;
    // Blur & Noise
    Gtk::CheckButton* const activlum;
    sigc::connection activlumConn;
    Gtk::CheckButton* const enablMask;
    sigc::connection enablMaskConn;
    Gtk::CheckButton* const fftwbl;
    sigc::connection fftwblConn;
    Gtk::CheckButton* const toolbl;
    sigc::connection toolblConn;
    //Tone mapping
    Gtk::CheckButton* const equiltm;
    sigc::connection equiltmConn;
    Gtk::CheckButton* const enatmMask;
    sigc::connection enatmMaskConn;
    Gtk::CheckButton* const enatmMaskaft;
    sigc::connection enatmMaskaftConn;
    // Retinex
    Gtk::CheckButton* const equilret;
    sigc::connection equilretConn;
    Gtk::CheckButton* const inversret;
    sigc::connection inversretConn;
    Gtk::CheckButton* const loglin;
    sigc::connection loglinConn;
    Gtk::CheckButton* const lumonly;
    sigc::connection lumonlyConn;
    Gtk::CheckButton* const enaretiMask;
    sigc::connection enaretiMaskConn;
    Gtk::CheckButton* const enaretiMasktmap;
    sigc::connection enaretiMasktmapConn;
    Gtk::CheckButton* const fftwreti;
    sigc::connection fftwretiConn;
    // Sharpening
    Gtk::CheckButton* const inverssha;
    sigc::connection inversshaConn;
    //local contrast
    Gtk::CheckButton* const fftwlc;
    sigc::connection fftwlcConn;
    Gtk::CheckButton* const blurlc;
    sigc::connection blurlcConn;
    Gtk::CheckButton* const origlc;
    sigc::connection origlcConn;
    Gtk::CheckButton* const enalcMask;
    sigc::connection enalcMaskConn;
    Gtk::CheckButton* const wavblur;
    sigc::connection wavblurConn;
    Gtk::CheckButton* const wavcont;
    sigc::connection wavcontConn;
    Gtk::CheckButton* const wavcomp;
    sigc::connection wavcompConn;
    Gtk::CheckButton* const wavcompre;
    sigc::connection wavcompreConn;
    Gtk::CheckButton* const wavgradl;
    sigc::connection wavgradlConn;

    //CBDL
    Gtk::CheckButton* const enacbMask;
    sigc::connection enacbMaskConn;
    //encoding log
    Gtk::CheckButton* const Autogray;
    sigc::connection AutograyConn;
    Gtk::CheckButton* const fullimage;
    sigc::connection fullimageConn;

    // ComboBox widgets
    // Color & Light
    MyComboBoxText* const qualitycurveMethod;
    sigc::connection qualitycurveMethodConn;
    MyComboBoxText* const gridMethod;
    sigc::connection gridMethodConn;
    MyComboBoxText* const merMethod;
    sigc::connection merMethodConn;
//    MyComboBoxText* const toneMethod;
    sigc::connection toneMethodConn;
    MyComboBoxText* const showmaskcolMethod;
    sigc::connection showmaskcolMethodConn;
    MyComboBoxText* const showmaskcolMethodinv;
    sigc::connection showmaskcolMethodConninv;
    MyComboBoxText* const mergecolMethod;
    sigc::connection mergecolMethodConn;

    ThresholdAdjuster* const csThresholdcol;
    //Exposure
    MyComboBoxText* const showmaskexpMethod;
    sigc::connection showmaskexpMethodConn;
    MyComboBoxText* const showmaskexpMethodinv;
    sigc::connection showmaskexpMethodConninv;
    MyComboBoxText* const expMethod;
    sigc::connection expMethodConn;
    MyComboBoxText* const exnoiseMethod;
    sigc::connection exnoiseMethodConn;
    //Shadows Highlight
    MyComboBoxText* const shMethod;
    sigc::connection shMethodConn;
    MyComboBoxText* const showmaskSHMethod;
    sigc::connection showmaskSHMethodConn;
    MyComboBoxText* const showmaskSHMethodinv;
    sigc::connection showmaskSHMethodConninv;
    //vibrance
    //Shadows Highlight
    MyComboBoxText* const showmaskvibMethod;
    sigc::connection showmaskvibMethodConn;
    // Blur & Noise
    MyComboBoxText* const blurMethod;
    sigc::connection blurMethodConn;
    //soft light
    MyComboBoxText* const softMethod;
    sigc::connection softMethodConn;
    MyComboBoxText* const showmasksoftMethod;
    sigc::connection showmasksoftMethodConn;
    //Blur and noise
    MyComboBoxText* const blMethod;
    sigc::connection blMethodConn;
    MyComboBoxText* const medMethod;
    sigc::connection medMethodConn;
    MyComboBoxText* const showmaskblMethod;
    sigc::connection showmaskblMethodConn;
    ThresholdAdjuster* const csThresholdblur;
    //TM
    MyComboBoxText* const showmasktmMethod;
    sigc::connection showmasktmMethodConn;
    // Retinex
    MyComboBoxText* const retinexMethod;
    sigc::connection retinexMethodConn;
    MyComboBoxText* const showmaskretiMethod;
    sigc::connection showmaskretiMethodConn;
    //Sharpening
    MyComboBoxText* const showmasksharMethod;
    sigc::connection showmasksharMethodConn;
    
    //local contrast
    MyComboBoxText* const localcontMethod;
    sigc::connection localcontMethodConn;
    ThresholdAdjuster* const csThreshold;
    MyComboBoxText* const showmasklcMethod;
    sigc::connection showmasklcMethodConn;


    //CBDL
    MyComboBoxText* const showmaskcbMethod;
    sigc::connection showmaskcbMethodConn;
    // ThresholdAdjuster widgets
    // Vibrance
    ThresholdAdjuster* const psThreshold;

    // Other widgets
    Gtk::ToggleButton* const autocompute;
    sigc::connection autoconn;
    Gtk::Label* const labqualcurv;
    Gtk::Button* const lumacontrastMinusButton;
    Gtk::Button* const lumaneutralButton;
    Gtk::Button* const lumacontrastPlusButton;
    sigc::connection lumacontrastMinusPressedConn, lumaneutralPressedConn, lumacontrastPlusPressedConn;
    Gtk::Frame* const gridFrame;
    Gtk::Frame* const struFrame;
    Gtk::Frame* const blurFrame;
    Gtk::Frame* const gridmerFrame;
    Gtk::Frame* const toolcolFrame;
    Gtk::Frame* const toolblFrame;
    Gtk::Frame* const mergecolFrame;
    Gtk::Frame* const merge1colFrame;
    Gtk::Frame* const pdeFrame;
    Gtk::Frame* const fatFrame;
    Gtk::Frame* const gradFramemask;
    Gtk::Frame* const fatSHFrame;
    Gtk::Frame* const gamFrame;
    Gtk::Frame* const dehaFrame;
    Gtk::Frame* const retiFrame;
    Gtk::Frame* const retitoolFrame;
    Gtk::Frame* const residFrame;
    Gtk::Frame* const sharFrame;
    Gtk::Frame* const clariFrame;
    Gtk::Frame* const gradwavFrame;
    Gtk::Frame* const blurlevelFrame;
    Gtk::Frame* const contFrame;
    Gtk::Frame* const compFrame;
    Gtk::Frame* const compreFrame;
    Gtk::Frame* const grainFrame;
    Gtk::Frame* const logFrame;
    Gtk::Frame* const logPFrame;
    Gtk::Frame* const gradlogFrame;
    ToolParamBlock* const retiBox;
    ToolParamBlock* const maskretiBox;
    ToolParamBlock* const mask7;
    LabGrid* const labgrid;
    LabGrid* const labgridmerg;
    Gtk::Label* const mMLabels;
    Gtk::Label* const transLabels;
    Gtk::Label* const transLabels2;
    // Others
    Gtk::HBox* const ctboxsoftmethod;
    Gtk::HBox* const ctboxexpmethod;
    double nextmin;
    double nextmax;
    double nextminiT;
    double nextmaxiT;
    double nextmeanT;
    double nextminT;
    double nextmaxT;
    double nextsigma;
//    bool lastAutogray;
//    Adjuster* sourceGray;

    int complexsoft;
    /**
     * Used to store the default ProcParams when setDefaults function is called
     * When an other spot is selected, this default ProcParams is used to update adjusters default values
     */
    const rtengine::procparams::ProcParams* defparams;
    /**
     * Used to store the default ParamsEdited when setDefaults function is called
     * When an other spot is selected, this default ParamsEdited is used to update adjusters default edited state
     */
    const ParamsEdited* defpedited;
    /**
     * Used to store the default ParamsEdited when setDefaults function is called
     * This ParamsEdited is updated when control spots are modified and is used to update adjusters edited state
     */
    ParamsEdited* pe;

    // Expander management functions
    void foldAllButMe(GdkEventButton* event, MyExpander *expander);
    void enableToggled(MyExpander *expander);

    // ButtonCheck event functions
    // Color & Light
    void curvactivChanged();
    void inversChanged();
    void specialChanged();
    void toolcolChanged();
    void enaColorMaskChanged();
    void fftColorMaskChanged();
    // Exposure
    void enaExpMaskChanged();
    void enaExpMaskaftChanged();
    void inversexChanged();
    //Shadows Highlight
    void enaSHMaskChanged();
    void inversshChanged();
    // Vibrance
    void protectskins_toggled();
    void avoidcolorshift_toggled();
    void pastsattog_toggled();
    void enavibMaskChanged();
    // Blur & Noise
    void activlumChanged();
    void enablMaskChanged();
    void fftwblChanged();
    void toolblChanged();
    //TM
    void equiltmChanged();
    void enatmMaskChanged();
    void enatmMaskaftChanged();
    // Retinex
    void equilretChanged();
    void loglinChanged();
    void lumonlyChanged();
    void inversretChanged();
    void enaretiMaskChanged();
    void enaretiMasktmapChanged();
    void fftwretiChanged();
    // Sharpening
    void inversshaChanged();
    // local contrast
    void fftwlcChanged();
    void blurlcChanged();
    void origlcChanged();
    void enalcMaskChanged();
    void wavblurChanged();
    void wavcontChanged();
    void wavcompChanged();
    void wavcompreChanged();
    void wavgradlChanged();
    //CBDL
    void enacbMaskChanged();
    // ComboBox event functions
    // Color & Light
    void qualitycurveMethodChanged();
    void gridMethodChanged();
    void merMethodChanged();
    void toneMethodChanged();
    void showmaskcolMethodChanged();
    void showmaskcolMethodChangedinv();
    void mergecolMethodChanged();
    //Exposure
    void showmaskexpMethodChanged();
    void showmaskexpMethodChangedinv();
    void expMethodChanged();
    void exnoiseMethodChanged();
    //Shadows Highlight
    void shMethodChanged();
    void showmaskSHMethodChanged();
    void showmaskSHMethodChangedinv();
    //vibrance
    void showmaskvibMethodChanged();
    // Blur & Noise
    void blMethodChanged();
    void medMethodChanged();
    // Soft light
    void softMethodChanged();
    void showmasksoftMethodChanged();
    //Blur
    void showmaskblMethodChanged();
    void blurMethodChanged();
    //TM
    void showmasktmMethodChanged();
    // Retinex
    void retinexMethodChanged();
    void showmaskretiMethodChanged();
    //sharp
    void showmasksharMethodChanged();
    
    //Local contrast
    void localcontMethodChanged();
    void showmasklcMethodChanged();
    //CBDL
    void showmaskcbMethodChanged();
    //log encoding
    void autocomputeToggled();
    void AutograyChanged();
    void fullimageChanged();

    // Other widgets event functions
    void lumacontrastMinusPressed();
    void lumaneutralPressed();
    void lumacontrastPlusPressed();

    // Locallab GUI management function
    void updateLocallabGUI(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited, int index);
    void updateSpecificGUIState();
    void setParamEditable(bool cond);


public:
    Locallab();
    ~Locallab();

    // FoldableToolPanel management functions
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited, int id);
    void setBatchMode(bool batchMode);
    void trimValues(rtengine::procparams::ProcParams* pp);
    void setListener(ToolPanelListener* tpl);
    void enableListener();
    void disableListener();
    void writeOptions(std::vector<int> &tpOpen);
    void updateToolState(std::vector<int> &tpOpen);
    void refChanged(double huer, double lumar, double chromar);
    void minmaxChanged(double cdma, double cdmin, double mini, double maxi, double Tmean, double Tsigma, double Tmin, double Tmax) override;
    void logencodChanged (float blackev, float whiteev, float sourceg, float targetg);
    void updateLabel();
    void updateTrans();

    // Mask visibility management functions
    struct llMaskVisibility {
        int colorMask;
        int colorMaskinv;
        int expMask;
        int expMaskinv;
        int vibMask;
        int SHMask;
        int SHMaskinv;
        int lcMask;
        int cbMask;
        int retiMask;
        int softMask;
        int tmMask;
        int blMask;
        int sharMask;
    };

    void resetMaskVisibility();
    llMaskVisibility* getMaskVisibility();

    // EditProvider management function
    void setEditProvider(EditDataProvider* provider);
    void subscribe();
    void unsubscribe();

    // FoldableToolPanel event function
    void enabledChanged();

    // Curve management function
    void autoOpenCurve();

    // Curve event function
    void curveChanged(CurveEditor* ce);

    // Adjuster event function
    void adjusterChanged(Adjuster* a, double newval);
    void adjusterAutoToggled(Adjuster* a, bool newval);

    // ThresholdAdjuster event functions
    virtual void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const;
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop);
    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop);
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight);
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight);
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR);
};
