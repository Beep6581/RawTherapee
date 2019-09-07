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
    MyExpander* const expsharp;
    MyExpander* const expcontrast;
    MyExpander* const expcbdl;
    MyExpander* const expdenoi;
    MyExpander* const expmaskcol;
    MyExpander* const expmaskexp;
    MyExpander* const expmasksh;
    MyExpander* const expmaskcb;
    MyExpander* const expmaskreti;
    MyExpander* const expmasktm;
    MyExpander* const expmaskbl;
    sigc::connection enablecolorConn, enableexposeConn, enableshadhighConn, enablevibranceConn, enablesoftConn, enableblurConn, enabletonemapConn, enableretiConn, enablesharpConn, enablecontrastConn, enablecbdlConn, enabledenoiConn;

    // Curve widgets
    // Color & Light
    CurveEditorGroup* const llCurveEditorG;
    CurveEditorGroup* const HCurveEditorG;
    CurveEditorGroup* const maskCurveEditorG;
    DiagonalCurveEditor* llshape;
    DiagonalCurveEditor* ccshape;
    FlatCurveEditor* LHshape;
    FlatCurveEditor* HHshape;
    FlatCurveEditor* CCmaskshape;
    FlatCurveEditor* LLmaskshape;
    FlatCurveEditor* HHmaskshape;
    // Exposure
    CurveEditorGroup* const curveEditorG;
    CurveEditorGroup* const maskexpCurveEditorG;
    DiagonalCurveEditor* shapeexpos;
    FlatCurveEditor* CCmaskexpshape;
    FlatCurveEditor* LLmaskexpshape;
    FlatCurveEditor* HHmaskexpshape;
    //Shadows Highlight
    CurveEditorGroup* const maskSHCurveEditorG;
    FlatCurveEditor* CCmaskSHshape;
    FlatCurveEditor* LLmaskSHshape;
    FlatCurveEditor* HHmaskSHshape;
    // Vibrance
    CurveEditorGroup* const curveEditorGG;
    DiagonalCurveEditor* skinTonesCurve;
    //Blur and noise
    CurveEditorGroup* const maskblCurveEditorG;
    FlatCurveEditor* CCmaskblshape;
    FlatCurveEditor* LLmaskblshape;
    FlatCurveEditor* HHmaskblshape;
    // TM
    CurveEditorGroup* const masktmCurveEditorG;
    FlatCurveEditor* CCmasktmshape;
    FlatCurveEditor* LLmasktmshape;
    FlatCurveEditor* HHmasktmshape;
    // Retinex
    CurveEditorGroup* const LocalcurveEditorgainT;
    CurveEditorGroup* const maskretiCurveEditorG;
    FlatCurveEditor* cTgainshape;
    FlatCurveEditor* CCmaskretishape;
    FlatCurveEditor* LLmaskretishape;
    FlatCurveEditor* HHmaskretishape;
    //local contrast
    CurveEditorGroup* const LocalcurveEditorwav;
    FlatCurveEditor* wavshape;

    //Cbdl
    CurveEditorGroup* const maskcbCurveEditorG;
    FlatCurveEditor* CCmaskcbshape;
    FlatCurveEditor* LLmaskcbshape;
    FlatCurveEditor* HHmaskcbshape;

    // Adjuster widgets
    // Color & Light
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Adjuster* const strengthgrid;
    Adjuster* const sensi;
    Adjuster* const structcol;
    Adjuster* const blurcolde;
    Adjuster* const blendmaskcol;
    Adjuster* const radmaskcol;
    Adjuster* const chromaskcol;
    Adjuster* const gammaskcol;
    Adjuster* const slomaskcol;
    Adjuster* const softradiuscol;
    // Exposure
    Adjuster* const expcomp;
    Adjuster* const hlcompr;
    Adjuster* const hlcomprthresh;
    Adjuster* const black;
    Adjuster* const shadex;
    Adjuster* const shcompr;
    Adjuster* const expchroma;
    Adjuster* const warm;
    Adjuster* const sensiex;
    Adjuster* const structexp;
    Adjuster* const blurexpde;
    Adjuster* const blendmaskexp;
    Adjuster* const radmaskexp;
    Adjuster* const chromaskexp;
    Adjuster* const gammaskexp;
    Adjuster* const slomaskexp;
    Adjuster* const softradiusexp;
    Adjuster* const laplacexp;
    Adjuster* const balanexp;
    Adjuster* const linear;
    Adjuster* const gamm;
    //Shadow highlight
    Adjuster* const highlights;
    Adjuster* const h_tonalwidth;
    Adjuster* const shadows;
    Adjuster* const s_tonalwidth;
    Adjuster* const sh_radius;
    Adjuster* const sensihs;
    Adjuster* const blendmaskSH;
    Adjuster* const radmaskSH;
    Adjuster* const blurSHde;
    Adjuster* const chromaskSH;
    Adjuster* const gammaskSH;
    Adjuster* const slomaskSH;
    // Vibrance
    Adjuster* const saturated;
    Adjuster* const pastels;
    Adjuster* const sensiv;
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
    // Retinex
    Adjuster* const str;
    Adjuster* const chrrt;
    Adjuster* const neigh;
    Adjuster* const vart;
    Adjuster* const dehaz;
    Adjuster* const sensih;
    Adjuster* const softradiusret;

    Adjuster* const blendmaskreti;
    Adjuster* const radmaskreti;
    Adjuster* const chromaskreti;
    Adjuster* const gammaskreti;
    Adjuster* const slomaskreti;
    Adjuster* const scalereti;
    Adjuster* const darkness;
    Adjuster* const lightnessreti;
    Adjuster* const limd;

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
    Adjuster* const clarilres;
    Adjuster* const clarisoft;
    Adjuster* const claricres;
    Adjuster* const sensilc;
    Adjuster* const residchro;
    // Contrast by detail levels
    Adjuster* multiplier[6];
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

    // ButtonCheck widgets
    // Color & Light
    Gtk::CheckButton* const curvactiv;
    Gtk::CheckButton* const invers;
    Gtk::CheckButton* const enaColorMask;
    sigc::connection curvactivConn, inversConn, enaColorMaskConn;
    // Exposure
    Gtk::CheckButton* const enaExpMask;
    sigc::connection enaExpMaskConn;
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
    sigc::connection pskinsconn, ashiftconn, pastsattogconn;
    // Blur & Noise
    Gtk::CheckButton* const activlum;
    sigc::connection activlumConn;
    Gtk::CheckButton* const enablMask;
    sigc::connection enablMaskConn;
    //Tone mapping
    Gtk::CheckButton* const equiltm;
    sigc::connection equiltmConn;
    Gtk::CheckButton* const enatmMask;
    sigc::connection enatmMaskConn;
    // Retinex
    Gtk::CheckButton* const equilret;
    sigc::connection equilretConn;
    Gtk::CheckButton* const inversret;
    sigc::connection inversretConn;
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

    //CBDL
    Gtk::CheckButton* const enacbMask;
    sigc::connection enacbMaskConn;

    // ComboBox widgets
    // Color & Light
    MyComboBoxText* const qualitycurveMethod;
    sigc::connection qualitycurveMethodConn;
    MyComboBoxText* const gridMethod;
    sigc::connection gridMethodConn;
    MyComboBoxText* const showmaskcolMethod;
    sigc::connection showmaskcolMethodConn;
    MyComboBoxText* const showmaskcolMethodinv;
    sigc::connection showmaskcolMethodConninv;
    //Exposure
    MyComboBoxText* const showmaskexpMethod;
    sigc::connection showmaskexpMethodConn;
    MyComboBoxText* const showmaskexpMethodinv;
    sigc::connection showmaskexpMethodConninv;
    MyComboBoxText* const expMethod;
    sigc::connection expMethodConn;
    //Shadows Highlight
    MyComboBoxText* const showmaskSHMethod;
    sigc::connection showmaskSHMethodConn;
    MyComboBoxText* const showmaskSHMethodinv;
    sigc::connection showmaskSHMethodConninv;
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
    //TM
    MyComboBoxText* const showmasktmMethod;
    sigc::connection showmasktmMethodConn;
    // Retinex
    MyComboBoxText* const retinexMethod;
    sigc::connection retinexMethodConn;
    MyComboBoxText* const showmaskretiMethod;
    sigc::connection showmaskretiMethodConn;
    //local contrast
    MyComboBoxText* const localcontMethod;
    sigc::connection localcontMethodConn;
    //CBDL
    MyComboBoxText* const showmaskcbMethod;
    sigc::connection showmaskcbMethodConn;
    // ThresholdAdjuster widgets
    // Vibrance
    ThresholdAdjuster* const psThreshold;

    // Other widgets
    Gtk::Label* const labqualcurv;
    Gtk::Button* const lumacontrastMinusButton;
    Gtk::Button* const lumaneutralButton;
    Gtk::Button* const lumacontrastPlusButton;
    sigc::connection lumacontrastMinusPressedConn, lumaneutralPressedConn, lumacontrastPlusPressedConn;
    Gtk::Frame* gridFrame;
    Gtk::Frame* pdeFrame;
    Gtk::Frame* residFrame;
    Gtk::Frame* clariFrame;
    LabGrid *labgrid;
    // Others
    Gtk::HBox* ctboxsoftmethod;
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
    void enaColorMaskChanged();
    // Exposure
    void enaExpMaskChanged();
    void inversexChanged();
    //Shadows Highlight
    void enaSHMaskChanged();
    void inversshChanged();
    // Vibrance
    void protectskins_toggled();
    void avoidcolorshift_toggled();
    void pastsattog_toggled();
    // Blur & Noise
    void activlumChanged();
    void enablMaskChanged();
    //TM
    void equiltmChanged();
    void enatmMaskChanged();
    // Retinex
    void equilretChanged();
    void inversretChanged();
    void enaretiMaskChanged();
    void enaretiMasktmapChanged();
    void fftwretiChanged();
    // Sharpening
    void inversshaChanged();
    // local contrast
    void fftwlcChanged();
    //CBDL
    void enacbMaskChanged();
    // ComboBox event functions
    // Color & Light
    void qualitycurveMethodChanged();
    void gridMethodChanged();
    void showmaskcolMethodChanged();
    void showmaskcolMethodChangedinv();
    //Exposure
    void showmaskexpMethodChanged();
    void showmaskexpMethodChangedinv();
    void expMethodChanged();
    //Shadows Highlight
    void showmaskSHMethodChanged();
    void showmaskSHMethodChangedinv();
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
    //Local contrast
    void localcontMethodChanged();
    //CBDL
    void showmaskcbMethodChanged();
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

    // Mask visibility management functions
    struct llMaskVisibility {
        int colorMask;
        int colorMaskinv;
        int expMask;
        int expMaskinv;
        int SHMask;
        int SHMaskinv;
        int cbMask;
        int retiMask;
        int softMask;
        int tmMask;
        int blMask;
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
