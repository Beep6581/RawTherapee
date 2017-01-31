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


class Locallab :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::localListener,
    public CurveListener,
    public EditSubscriber,
    public ColorProvider

{
private:
    int lastObject;
    void foldAllButMe (GdkEventButton* event, MyExpander *expander);
    void enableToggled (MyExpander *expander);

//protected:
    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;

    Adjuster* nbspot;
    Adjuster* multiplier[5];

    Adjuster* const degree;
    Adjuster* const locX;
    Adjuster* const locY;
    Adjuster* const locXL;
    Adjuster* const locYT;
    Adjuster* const centerX;
    Adjuster* const centerY;
    Adjuster* const circrad;
    Adjuster* const lightness;
    Adjuster* const contrast;
    Adjuster* const chroma;
    Adjuster* const sensi;
    Adjuster* const sensih;
    Adjuster* const radius;
    Adjuster* const strength;
    Adjuster* const transit;
    Adjuster* const str;
    Adjuster* const neigh;
    Adjuster* const vart;
    Adjuster* const chrrt;
    Adjuster* const anbspot;
    Adjuster* const sharradius;
    Adjuster* const sharamount;
    Adjuster* const shardamping;
    Adjuster* const shariter;
    Adjuster* const sensisha;
    Adjuster* const thres;
    Adjuster* const proxi;
    Adjuster* const noiselumf;
    Adjuster* const noiselumc;
    Adjuster* const noisechrof;
    Adjuster* const noisechroc;
    Adjuster* const threshold;
    Adjuster* const sensicb;
    Adjuster* const sensibn;
    Adjuster* const stren;
    Adjuster* const gamma;
    Adjuster* const estop;
    Adjuster* const scaltm;
    Adjuster* const rewei;
    Adjuster* const sensitm;
    Adjuster* const retrab;

    MyExpander* const expcolor;
    MyExpander* const expblur;
    MyExpander* const exptonemap;
    MyExpander* const expreti;
    MyExpander* const expsharp;
    MyExpander* const expcbdl;
    MyExpander* const expdenoi;


    sigc::connection lumaneutralPressedConn;
    sigc::connection lumacontrastPlusPressedConn;
    sigc::connection lumacontrastMinusPressedConn;
    sigc::connection enablecolorConn, enableblurConn, enabletonemapConn;
    sigc::connection enableretiConn, enablesharpConn, enablecbdlConn;
    sigc::connection enabledenoiConn;
    sigc::connection  editConn, avoidConn, inversConn, curvactivConn, activlumConn, inversradConn, inversretConn, inversshaConn,  neutralconn, neutralconn1;

    Gtk::HBox* const ctboxS;
    Gtk::HBox* const qualbox;
    Gtk::HBox* const qualcurvbox;

    Gtk::Frame* const artifFrame;
    Gtk::Frame* const shapeFrame;
    Gtk::Frame* const superFrame;

    Gtk::VBox* const artifVBox;
    Gtk::VBox* const shapeVBox;
    Gtk::VBox* const tmBox;
    Gtk::VBox* const retiBox;
    Gtk::VBox* const colorVBox;
    Gtk::VBox* const blurrVBox;
    Gtk::VBox* const sharpVBox;
    Gtk::VBox* const cbdlVBox;
    Gtk::VBox* const denoisVBox;
    Gtk::VBox* const superVBox;

    Gtk::CheckButton* const avoid;
    Gtk::CheckButton* const invers;
    Gtk::CheckButton* const curvactiv;
    Gtk::CheckButton* const inversrad;
    Gtk::CheckButton* const inversret;
    Gtk::CheckButton* const activlum;
    Gtk::CheckButton* const inverssha;

    Gtk::Button* neutral;
    Gtk::HBox* neutrHBox;
    Gtk::Button* neutral1;
    Gtk::HBox* neutrHBox1;

    MyComboBoxText*   const Smethod;
    sigc::connection  Smethodconn;
    MyComboBoxText*   const retinexMethod;
    sigc::connection retinexMethodConn;
    MyComboBoxText*   const qualityMethod;
    sigc::connection qualityMethodConn;
    MyComboBoxText*   const qualitycurveMethod;
    sigc::connection qualitycurveMethodConn;


    Gtk::Label* const labmdh;
    Gtk::Label* const labqual;
    Gtk::Label* const labqualcurv;
    Gtk::Label* const labmS;

    Gtk::HBox* const dhbox;

    CurveEditorGroup* const LocalcurveEditorgainT;
    FlatCurveEditor* cTgainshape;
    CurveEditorGroup* const LocalcurveEditorgainTrab;
    FlatCurveEditor* cTgainshaperab;
    CurveEditorGroup* const llCurveEditorG;
    DiagonalCurveEditor* llshape;
    DiagonalCurveEditor* ccshape;
    Gtk::Image* irg;
    FlatCurveEditor* LHshape;



    int nextdatasp[61];
    int nextlength;
    std::string nextstr;
    std::string nextstr2;
    std::string nextll_str;
    std::string nextll_str2;
    std::string nextlh_str;
    std::string nextlh_str2;
    std::string nextcc_str;
    std::string nextcc_str2;

    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    double draggedlocYOffset;
    double draggedlocXOffset;
    double draggedlocYTOffset;
    double draggedlocXLOffset;
    rtengine::Coord draggedCenter;
    bool lastavoid, lastinvers, lastinversrad, lastinversret, lastactivlum, lastinverssha, lastcurvactiv;
    int lastanbspot;

    void editToggled ();

public:

    Locallab ();
    ~Locallab ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);

    void setBatchMode   (bool batchMode);

    void updateGeometry (const int centerX_, const int centerY_, const int circrad_, const int locY_, const double degree_, const int locX_, const int locYT_, const int locXL_, const int fullWidth = -1, const int fullHeight = -1);
    void SmethodChanged      ();
    void writeOptions (std::vector<int> &tpOpen);
    void updateToolState (std::vector<int> &tpOpen);

    void adjusterChanged (Adjuster* a, double newval);
    void enabledChanged ();
    void setAdjusterBehavior (bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd, bool strengthadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
    void avoidChanged ();
    void activlumChanged ();
    void inversChanged ();
    void curvactivChanged ();
    void inversradChanged ();
    void inversretChanged ();
    void inversshaChanged ();
    void curveChanged (CurveEditor* ce);
    void autoOpenCurve ();
    void localChanged           (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat);
    void localretChanged           (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, int sp, int maxdat);
    bool localComputed_         ();
    bool localretComputed_         ();
    void setEditProvider (EditDataProvider* provider);
    void retinexMethodChanged();
    void qualityMethodChanged();
    void qualitycurveMethodChanged();
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();
    void neutral_pressed       ();
    virtual void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller);

    // EditSubscriber interface
    CursorShape getCursor (int objectID);
    bool mouseOver (int modifierKey);
    bool button1Pressed (int modifierKey);
    bool button1Released();
    bool drag1 (int modifierKey);
    void switchOffEditMode ();
};

