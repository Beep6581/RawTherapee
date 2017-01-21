/*
 *  This file is part of RawTherapee.
 */
#ifndef _LOCALLAB_H_
#define _LOCALLAB_H_

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


class Locallab : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public rtengine::localListener, public CurveListener, public EditSubscriber, public ColorProvider
{
private:
    int lastObject;
    void foldAllButMe (GdkEventButton* event, MyExpander *expander);
    void enableToggled (MyExpander *expander);

protected:
//   Gtk::CheckButton* enabled;
    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;
    Adjuster* degree;
    Adjuster* locX;
    Adjuster* locY;
    Adjuster* locXL;
    Adjuster* locYT;
    Adjuster* centerX;
    Adjuster* centerY;
    Adjuster* circrad;
    Adjuster* lightness;
    Adjuster* contrast;
    Adjuster* chroma;
    Adjuster* sensi;
    Adjuster* sensih;
    Adjuster* radius;
    Adjuster* strength;
    Adjuster* transit;
    Adjuster* str;
    Adjuster* neigh;
    Adjuster* vart;
    Adjuster* chrrt;
    Adjuster* nbspot;
    Adjuster* anbspot;
    Adjuster* maxn;
    Adjuster* sharradius;
    Adjuster* sharamount;
    Adjuster* shardamping;
    Adjuster* shariter;
    Adjuster* sensisha;
    Adjuster* thres;
    Adjuster* proxi;
    Adjuster* noiselumf;
    Adjuster* noiselumc;
    Adjuster* noisechrof;
    Adjuster* noisechroc;
    Adjuster* multiplier[5];
    Adjuster* threshold;
    Adjuster* sensicb;
    Adjuster* sensibn;
    Adjuster* stren;
    Adjuster* gamma;
    Adjuster* estop;
    Adjuster* scaltm;
    Adjuster* rewei;
    Adjuster* sensitm;
    Adjuster* retrab;

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

    Gtk::CheckButton* avoid;
    MyComboBoxText*   Smethod;
    sigc::connection  Smethodconn;
    Gtk::HBox* ctboxS;
    Gtk::CheckButton* invers;
    Gtk::CheckButton* curvactiv;
    Gtk::CheckButton* inversrad;
    Gtk::CheckButton* inversret;
    Gtk::CheckButton* activlum;
    Gtk::CheckButton* inverssha;

    Gtk::Button* neutral;
    Gtk::HBox* neutrHBox;

    Gtk::Button* neutral1;
    Gtk::HBox* neutrHBox1;

    MyComboBoxText*   retinexMethod;
    MyComboBoxText*   qualityMethod;
    Gtk::Label* labmdh;
    Gtk::HBox* dhbox;
    CurveEditorGroup* LocalcurveEditorgainT;
    FlatCurveEditor* cTgainshape;
    CurveEditorGroup* LocalcurveEditorgainTrab;
    FlatCurveEditor* cTgainshaperab;
    CurveEditorGroup* llCurveEditorG;
    DiagonalCurveEditor* llshape;
    Gtk::Image* irg;
    FlatCurveEditor* LHshape;

    int nextdatasp[60];
    int nextlength;
    std::string nextstr;
    std::string nextstr2;
    std::string nextll_str;
    std::string nextll_str2;
    std::string nextlh_str;
    std::string nextlh_str2;

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
    sigc::connection  editConn, avoidConn, inversConn, curvactivConn, activlumConn, inversradConn, inversretConn, inversshaConn, retinexMethodConn, qualityMethodConn, neutralconn, neutralconn1;

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
    void localChanged           (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, int sp, int maxdat);
    void localretChanged           (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, int sp, int maxdat);
    bool localComputed_         ();
    bool localretComputed_         ();
    void setEditProvider (EditDataProvider* provider);
    void retinexMethodChanged();
    void qualityMethodChanged();
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

#endif
