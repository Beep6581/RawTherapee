/*
 *  This file is part of RawTherapee.
 */

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "toolpanel.h"
#include "../rtengine/imagedata.h"
#include <memory>
#include "options.h"
#include <string>
#include "../rtengine/improcfun.h"
#include "thresholdadjuster.h"


class Locallab :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::localListener,
    public CurveListener,
    public EditSubscriber,
    public ColorProvider,
    public ThresholdCurveProvider,
    public ThresholdAdjusterListener

{
private:

    rtengine::ProcEvent EvLocenacolor;//548
    rtengine::ProcEvent EvLocenaexpose;//572
    rtengine::ProcEvent EvLocenavibrance;//563
    rtengine::ProcEvent EvLocenablur;//549
    rtengine::ProcEvent EvLocenatonemap;//550
    rtengine::ProcEvent EvLocenareti;//551
    rtengine::ProcEvent EvLocenasharp;//552
    rtengine::ProcEvent EvLocenacbdl;//553
    rtengine::ProcEvent EvLocenadenoi;//554

    rtengine::ProcEvent EvlocallablocX ;//= 494,
    rtengine::ProcEvent EvlocallabCenter;// = 495,
    rtengine::ProcEvent EvlocallabDegree;// = 496,
    rtengine::ProcEvent Evlocallablightness;// = 497,
    rtengine::ProcEvent Evlocallabcontrast;// = 498,
    rtengine::ProcEvent Evlocallabchroma;// = 499,
    rtengine::ProcEvent Evlocallabtransit;// = 500,
    rtengine::ProcEvent Evlocallabavoid;// = 501,
    rtengine::ProcEvent EvlocallablocYT;// = 502,
    rtengine::ProcEvent EvlocallablocXL;// = 503,
    rtengine::ProcEvent EvlocallabSmet;// = 504,
    rtengine::ProcEvent Evlocallabinvers;// = 505,
    rtengine::ProcEvent Evlocallabradius;// = 506,
    rtengine::ProcEvent Evlocallabinversrad;// = 507,
    rtengine::ProcEvent Evlocallabstrength;// = 508,
    rtengine::ProcEvent Evlocallabsensi;// = 509,
    rtengine::ProcEvent EvlocallabretinexMethod;//510
    rtengine::ProcEvent Evlocallabstr;// = 511,
    rtengine::ProcEvent Evlocallabneigh;// = 512,
    rtengine::ProcEvent Evlocallabvart;// = 513,
    rtengine::ProcEvent EvlocallabCTgainCurve;// = 514,
    rtengine::ProcEvent Evlocallabchrrt;// = 515,
    rtengine::ProcEvent Evlocallabinversret;// = 516,
    rtengine::ProcEvent Evlocallabsensih;// = 517,
    rtengine::ProcEvent Evlocallabnbspot;// = 518,
    rtengine::ProcEvent Evlocallabactivlum;// = 519,
    rtengine::ProcEvent Evlocallabanbspot;// = 520,
    rtengine::ProcEvent Evlocallabsharradius;// = 521,
    rtengine::ProcEvent Evlocallabsharamount;// = 522,
    rtengine::ProcEvent Evlocallabshardamping;// = 523,
    rtengine::ProcEvent Evlocallabshariter;// = 524,
    rtengine::ProcEvent Evlocallabsensis;// = 525,
    rtengine::ProcEvent Evlocallabinverssha;// = 526,
    rtengine::ProcEvent Evlocallabcircrad;// = 527,
    rtengine::ProcEvent Evlocallabthres;// = 528,
    rtengine::ProcEvent Evlocallabproxi;// = 529,
    rtengine::ProcEvent EvlocallabqualityMethod;// = 530,
    rtengine::ProcEvent Evlocallabnoiselumf;// = 531,
    rtengine::ProcEvent Evlocallabnoiselumc;// = 532,
    rtengine::ProcEvent Evlocallabnoisechrof;// = 533,
    rtengine::ProcEvent Evlocallabnoisechroc;// = 534,
    rtengine::ProcEvent EvlocallabThresho;// = 535,
    rtengine::ProcEvent EvlocallabEqualizer;// = 536,
    rtengine::ProcEvent Evlocallabsensicb;// = 537,
    rtengine::ProcEvent Evlocallabsensibn;// = 538,
    rtengine::ProcEvent Evlocallabstren;// = 539,
    rtengine::ProcEvent Evlocallabgamma;// = 540,
    rtengine::ProcEvent Evlocallabestop;// = 541,
    rtengine::ProcEvent Evlocallabscaltm;// = 542,
    rtengine::ProcEvent Evlocallabrewei;// = 543,
    rtengine::ProcEvent Evlocallabsensitm;// = 544,
    rtengine::ProcEvent EvlocallabCTgainCurverab;// = 545,
    rtengine::ProcEvent Evlocallabretrab;// = 546,
    rtengine::ProcEvent Evlocallabllshape;// = 547,
    rtengine::ProcEvent EvlocallabLHshape;// = 555,
    rtengine::ProcEvent Evlocallabcurvactiv;// = 556,
    rtengine::ProcEvent Evlocallabccshape;// = 557,
    rtengine::ProcEvent EvlocallabqualitycurveMethod;// = 558,
    rtengine::ProcEvent Evlocallabhueref;// = 559,
    rtengine::ProcEvent Evlocallabchromaref;// = 560,
    rtengine::ProcEvent Evlocallablumaref;// = 561,
    rtengine::ProcEvent EvlocallabHHshape;// = 562,
    rtengine::ProcEvent EvlocallabSkinTonesCurve;// = 564,
    rtengine::ProcEvent EvlocallabProtectSkins;// = 565,
    rtengine::ProcEvent EvlocallabAvoidColorShift;// = 566,
    rtengine::ProcEvent EvlocallabPastSatTog;// = 567,
    rtengine::ProcEvent EvlocallabPastels;// = 568,
    rtengine::ProcEvent EvlocallabSaturated;// = 569,
    rtengine::ProcEvent EvlocallabPastSatThreshold;// = 570,
    rtengine::ProcEvent Evlocallabsensiv;// = 571,
    rtengine::ProcEvent Evlocallabexpcomp;// = 573,
    rtengine::ProcEvent Evlocallabhlcompr;// = 574,
    rtengine::ProcEvent Evlocallabhlcomprthresh;// = 575,
    rtengine::ProcEvent Evlocallabblack;// = 576,
    rtengine::ProcEvent Evlocallabshcompr;// = 577,
    rtengine::ProcEvent Evlocallabsensiex;// = 578,
    rtengine::ProcEvent Evlocallabshapeexpos;// = 579,
    rtengine::ProcEvent EvlocallabCenterbuf;// = 580,
    rtengine::ProcEvent Evlocallabadjblur;// = 581,
    rtengine::ProcEvent Evlocallabcutpast;// = 582,
    rtengine::ProcEvent Evlocallabchromacbdl;// = 583,
    rtengine::ProcEvent EvlocallabblurMethod;//584
    rtengine::ProcEvent EvlocallabdustMethod;// = 585,
    rtengine::ProcEvent Evlocallablastdust;// = 586,
    rtengine::ProcEvent Evlocallabsobelref;// = 587,
    rtengine::ProcEvent Evlocallabexclumethod;// = 588,
    rtengine::ProcEvent Evlocallabsensiexclu;// = 589,
    rtengine::ProcEvent Evlocallabstruc;// = 590,
    rtengine::ProcEvent Evlocallabwarm;// = 591,
    rtengine::ProcEvent Evlocallabnoiselumdetail;// = 592,
    rtengine::ProcEvent Evlocallabnoisechrodetail;// = 593,
    rtengine::ProcEvent Evlocallabsensiden;// = 594,
    rtengine::ProcEvent Evlocallabhuerefblur;// = 595,
    rtengine::ProcEvent EvlocallabEnabled;// = 596,
    rtengine::ProcEvent EvlocallablocY;// = 597,
    rtengine::ProcEvent Evlocallabbilateral;// = 598,
    rtengine::ProcEvent Evlocallabnoiselequal;// = 599,
    rtengine::ProcEvent Evlocallabshapemethod;// = 600,
    rtengine::ProcEvent Evlocallabspotduplicated;

    IdleRegister idle_register;

    int lastObject;
    void foldAllButMe(GdkEventButton* event, MyExpander *expander);
    void enableToggled(MyExpander *expander);

//protected:

    MyExpander* const expcolor;
    MyExpander* const expexpose;
    MyExpander* const expvibrance;
    MyExpander* const expblur;
    MyExpander* const exptonemap;
    MyExpander* const expreti;
    MyExpander* const expsharp;
    MyExpander* const expcbdl;
    MyExpander* const expdenoi;
    MyExpander* const expsettings;

    CurveEditorGroup* const LocalcurveEditorgainT;
    CurveEditorGroup* const LocalcurveEditorgainTrab;
    CurveEditorGroup* const llCurveEditorG;


    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;

    Adjuster* nbspot;
    Adjuster* multiplier[5];

    Adjuster* const anbspot;
    Adjuster* const locX;
    Adjuster* const locXL;
    Adjuster* const degree;
    Adjuster* const locY;
    Adjuster* const locYT;
    Adjuster* const centerX;
    Adjuster* const centerY;
    Adjuster* const circrad;
    Adjuster* const sensiexclu;
    Adjuster* const struc;

    Adjuster* const thres;
    Adjuster* const proxi;
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Adjuster* const sensi;

    Adjuster* const expcomp;
    Adjuster* const hlcompr;
    Adjuster* const hlcomprthresh;
    Adjuster* const black;
    Adjuster* const shcompr;
    /*
    Adjuster* const lightnessex;
    Adjuster* const contrastex;
    Adjuster* const chromaex;
    */
    Adjuster* const sensiex;
    Adjuster* const radius;
    Adjuster* const strength;
    Adjuster* const sensibn;
    Adjuster* const transit;
    Adjuster* const stren;
    Adjuster* const gamma;
    Adjuster* const estop;
    Adjuster* const scaltm;
    Adjuster* const rewei;
    Adjuster* const sensitm;
    Adjuster* const str;
    Adjuster* const neigh;
    Adjuster* const vart;
    Adjuster* const chrrt;
    Adjuster* const sensih;
    Adjuster* const retrab;
    Adjuster* const chromacbdl;
    Adjuster* const threshold;
    Adjuster* const sensicb;
    Adjuster* const sharradius;
    Adjuster* const sharamount;
    Adjuster* const shardamping;
    Adjuster* const shariter;
    Adjuster* const sensisha;
    Adjuster* const noiselumdetail;
    Adjuster* const noisechrodetail;
    Adjuster* const bilateral;
    Adjuster* const sensiden;
    Adjuster* const hueref;
    Adjuster* const huerefblur;
    Adjuster* const chromaref;
    Adjuster* const lumaref;
    Adjuster* const sobelref;
    Adjuster* const centerXbuf;
    Adjuster* const centerYbuf;
//    Adjuster* const adjblur;

    MyComboBoxText*   const shapemethod;
    MyComboBoxText*   const Smethod;
    MyComboBoxText*   const Exclumethod;
    MyComboBoxText*   const retinexMethod;
    MyComboBoxText*   const qualityMethod;
    MyComboBoxText*   const qualitycurveMethod;
    MyComboBoxText*   const blurMethod;
    MyComboBoxText*   const dustMethod;

    Gtk::Frame* const excluFrame;

    Gtk::Frame* const artifFrame;
    Gtk::Frame* const shapeFrame;
    Gtk::Frame* const superFrame;
    Gtk::Frame* const dustFrame;
    Gtk::Frame* const wavFrame;

    Gtk::Label* const labmdh;
    Gtk::Label* const labqual;
    Gtk::Label* const labqualcurv;
    Gtk::Label* const labmS;
    Gtk::Label* const labmEx;
    Gtk::Label* const labmshape;

    Gtk::HBox* const ctboxS;
    Gtk::HBox* const ctboxshape;
    Gtk::HBox* const ctboxEx;

    Gtk::HBox* const dhbox;
    Gtk::HBox* const qualbox;
    Gtk::HBox* const qualcurvbox;

    Gtk::CheckButton* const avoid;
    Gtk::CheckButton* const activlum;
    Gtk::CheckButton* const invers;
    Gtk::CheckButton* const curvactiv;
    Gtk::CheckButton* const inversrad;
    Gtk::CheckButton* const inversret;
    Gtk::CheckButton* const inverssha;
    Gtk::CheckButton* const cutpast;
    Gtk::CheckButton* const lastdust;

    Gtk::CheckButton* spotduplicated;
    Gtk::Label* labspotdup;

    Gtk::Button* neutral;
    Gtk::HBox* neutrHBox;
    Gtk::Button* neutral1;
    Gtk::HBox* neutrHBox1;

    CurveEditorGroup* curveEditorG;
    CurveEditorGroup* curveEditorG2;
    DiagonalCurveEditor* shapeexpos;
    DiagonalCurveEditor* shape2;

    FlatCurveEditor* cTgainshape;
    FlatCurveEditor* cTgainshaperab;
    DiagonalCurveEditor* llshape;
    DiagonalCurveEditor* ccshape;
    Gtk::Image* irg;
    FlatCurveEditor* LHshape;
    FlatCurveEditor* HHshape;

    CurveEditorGroup* curveEditorGG;
    Adjuster* pastels;
    Adjuster* saturated;
    Adjuster* warm;
    Adjuster* adjblur;
    Adjuster* noiselequal;
    Adjuster* noiselumf;
    Adjuster* noiselumc;
    Adjuster* noisechrof;
    Adjuster* noisechroc;

    ThresholdAdjuster* psThreshold;
    Gtk::CheckButton* protectSkins;
    Gtk::CheckButton* avoidColorShift;
    Gtk::CheckButton* pastSatTog;
    DiagonalCurveEditor* skinTonesCurve;
    Adjuster* sensiv;

    bool lastProtectSkins;
    bool lastAvoidColorShift;
    bool lastPastSatTog;

    sigc::connection pskinsconn;
    sigc::connection ashiftconn;
    sigc::connection pastsattogconn;



    sigc::connection lumaneutralPressedConn;
    sigc::connection lumacontrastPlusPressedConn;
    sigc::connection lumacontrastMinusPressedConn;
    sigc::connection enablecolorConn, enableexposeConn, enablevibranceConn, enableblurConn, enabletonemapConn;
    sigc::connection enableretiConn, enablesharpConn, enablecbdlConn;
    sigc::connection enabledenoiConn;
    sigc::connection  editConn, avoidConn, inversConn, cutpastConn, lastdustConn, curvactivConn, activlumConn, inversradConn, inversretConn, inversshaConn,  neutralconn, neutralconn1;
    sigc::connection  Smethodconn, shapemethodconn, Exclumethodconn, spotduplicatedConn;
    sigc::connection retinexMethodConn;
    sigc::connection qualityMethodConn;
    sigc::connection qualitycurveMethodConn;
    sigc::connection blurMethodConn;
    sigc::connection dustMethodConn;

    bool lastspotduplicated;


    int nextdatasp[102];
    int nextlength;
    bool nextspotdup;
    std::string nextstr;
    std::string nextstr2;
    std::string nextll_str;
    std::string nextll_str2;
    std::string nextlh_str;
    std::string nextlh_str2;
    std::string nextcc_str;
    std::string nextcc_str2;
    std::string nexthh_str;
    std::string nexthh_str2;
    std::string nextsk_str;
    std::string nextsk_str2;
    std::string nextps_str;
    std::string nextps_str2;
    std::string nextex_str;
    std::string nextex_str2;

    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    double draggedlocYOffset;
    double draggedlocXOffset;
    double draggedlocYTOffset;
    double draggedlocXLOffset;
    rtengine::Coord draggedCenter;
    bool lastavoid, lastinvers, lastcutpast, lastlastdust, lastinversrad, lastinversret, lastactivlum, lastinverssha, lastcurvactiv;
    int lastanbspot;

    void editToggled();

public:

    Locallab();
    ~Locallab();

    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void setBatchMode(bool batchMode);

    void updateGeometry(const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth = -1, const int fullHeight = -1);
    void SmethodChanged();
    void shapemethodChanged();
    void ExclumethodChanged();
    void writeOptions(std::vector<int> &tpOpen);
    void updateToolState(std::vector<int> &tpOpen);

    void adjusterChanged(Adjuster* a, double newval);
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop);

    void enabledChanged();
    void setAdjusterBehavior(bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd, bool strengthadd);
    void trimValues(rtengine::procparams::ProcParams* pp);
    void avoidChanged();
    void activlumChanged();
    void inversChanged();
    void curvactivChanged();
    void inversradChanged();
    void inversretChanged();
    void inversshaChanged();
    void cutpastChanged();
    void lastdustChanged();
    void spotduplicatedChanged();
    bool spotdupComputed_();

    void spotdupChanged(bool spotchan);
    void curveChanged(CurveEditor* ce);
    void autoOpenCurve();
    void localChanged(int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str, std::string ex_str, int sp, int maxdat);
    void localretChanged(int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str, std::string ex_str, int sp, int maxdat);
    bool localComputed_();
    bool localretComputed_();
    void setEditProvider(EditDataProvider* provider);
    void retinexMethodChanged();
    void blurMethodChanged();
    void dustMethodChanged();
    void qualityMethodChanged();
    void qualitycurveMethodChanged();
    void lumaneutralPressed();
    void lumacontrastPlusPressed();
    void lumacontrastMinusPressed();
    void neutral_pressed();
    virtual void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);
    void protectskins_toggled();
    void avoidcolorshift_toggled();
    void pastsattog_toggled();
    std::vector<double> getCurvePoints(ThresholdSelector* tAdjuster) const;

    // EditSubscriber interface
    CursorShape getCursor(int objectID);
    bool mouseOver(int modifierKey);
    bool button1Pressed(int modifierKey);
    bool button1Released();
    bool drag1(int modifierKey);
    void switchOffEditMode();
};

