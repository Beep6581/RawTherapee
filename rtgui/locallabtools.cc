/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>frame
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
 *  2019 Pierre Cabrera <pierre.cab@gmail.com>
 */
#include "locallabtools.h"

#include "options.h"
#include "../rtengine/procparams.h"
#include "locallab.h"
#include "thresholdadjuster.h"
#include "rtimage.h"
#include "../rtengine/color.h"

#define MINRAD 1.5
#define MAXRAD 10000
#define CENTERRAD 100
#define MINCHRO 0.
#define MAXCHRO 150.
#define MAXCHROCC 100.

using namespace rtengine;
using namespace procparams;

extern Options options;

static double blurSlider2radius(double sval)
{
    // Slider range: 0 - 1000
    double radius;

    if (sval <= 100) {
        // Linear below center-temp
        radius = MINRAD + (sval / 100.0) * (CENTERRAD - MINRAD);
    } else {
        const double slope = (double)(CENTERRAD - MINRAD) / (MAXRAD - CENTERRAD);
        double x = (sval - 100) / 100; // x range: 0 - 1
        double y = x * slope + (1.0 - slope) * pow(x, 4.0);
        radius = CENTERRAD + y * (MAXRAD - CENTERRAD);
    }

    if (radius < MINRAD) {
        radius = MINRAD;
    }

    if (radius > MAXRAD) {
        radius = MAXRAD;
    }

    return radius;
}

static double blurRadius2Slider(double radius)
{
    double sval;

    if (radius <= CENTERRAD) {
        sval = ((radius - MINRAD) / (CENTERRAD - MINRAD)) * 100.0;
    } else {
        const double slope = (double)(CENTERRAD - MINRAD) / (MAXRAD - CENTERRAD);
        const double y = (radius - CENTERRAD) / (MAXRAD - CENTERRAD);
        double x = pow(y, 0.25); // Rough guess of x, will be a little lower
        double k = 0.1;
        bool add = true;

        // The y=f(x) function is a mess to invert, therefore we have this trial-refinement loop instead.
        // From tests, worst case is about 20 iterations, i.e. no problem
        for (;;) {
            double y1 = x * slope + (1.0 - slope) * pow(x, 4.0);

            if (100. * fabs(y1 - y) < 0.1) {
                break;
            }

            if (y1 < y) {
                if (!add) {
                    k /= 2;
                }

                x += k;
                add = true;
            } else {
                if (add) {
                    k /= 2;
                }

                x -= k;
                add = false;
            }
        }

        sval = 100.0 + x * 100.0;
    }

    if (sval < 0.) {
        sval = 0.;
    }

    if (sval > 10000.) {
        sval = 10000.;
    }

    return sval;
}

/* ==== LocallabTool ==== */
LocallabTool::LocallabTool(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11):
    ToolPanel(toolName, need11),

    // LocallabTool parameters
    isLocActivated(false),
    locToolListener(nullptr)
{
    const bool showtooltip = options.showtooltip;

    // Create expander title bar
    Gtk::HBox* const titleBox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const titleLabel = Gtk::manage(new Gtk::Label());
    titleLabel->set_markup(Glib::ustring("<b>") + escapeHtmlChars(UILabel) + Glib::ustring("</b>"));
    titleLabel->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    titleBox->pack_start(*titleLabel, Gtk::PACK_EXPAND_WIDGET, 0);

    Gtk::EventBox* const removeEvBox = Gtk::manage(new Gtk::EventBox()); // Glue to manage mouse clicking event on remove image
    removeEvBox->set_can_focus(false);
    removeEvBox->set_above_child(false); // To have priority over expander title bar when mouse clicking on remove image
    removeEvBox->signal_button_release_event().connect(sigc::mem_fun(this, &LocallabTool::on_remove_change));
    RTImage* const removeImage = Gtk::manage(new RTImage("cancel-small.png"));
    removeEvBox->add(*removeImage);
    titleBox->pack_end(*removeEvBox, Gtk::PACK_SHRINK, 4);

    if (need100Percent) {
        RTImage* titleImage = Gtk::manage(new RTImage("one-to-one-small.png"));

        if (showtooltip) {
            titleImage->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
        }

        titleBox->pack_end(*titleImage, Gtk::PACK_SHRINK, 0);
    }

    exp = Gtk::manage(new MyExpander(true, titleBox));
    exp->signal_button_release_event().connect_notify(sigc::mem_fun(this, &LocallabTool::foldThemAll));
    enaExpConn = exp->signal_enabled_toggled().connect(sigc::mem_fun(*this, &LocallabTool::enabledChanged));

    ToolParamBlock* const totalBox = Gtk::manage(new ToolParamBlock());

    // Create panel for specific widget tools
    totalBox->pack_start(*content, Gtk::PACK_SHRINK, 0);
    exp->add(*totalBox, false);
}

LocallabTool::~LocallabTool()
{
    idle_register.destroy();
}

void LocallabTool::addLocallabTool(bool raiseEvent)
{
    exp->set_visible(true);

    // Raise event if required
    if (raiseEvent) {
        if (listener) {
            listener->panelChanged(EvlocallabToolAdded,
                                   toolName + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabTool::removeLocallabTool(bool raiseEvent)
{
    exp->set_visible(false);

    // Inform LocallabToolListener to update Locallab tools list
    if (locToolListener) {
        locToolListener->toolRemoved(this);
    }

    if (exp->getEnabled()) {
        // Disable tool while removing it
        disableListener();
        exp->setEnabled(false);
        enableListener();

        // Raise event if required refreshing image
        if (raiseEvent && listener) {
            listener->panelChanged(EvlocallabToolRemovedWithRefresh,
                                   toolName + " (" + escapeHtmlChars(spotName) + ")");
        }
    } else {
        // Raise event if required without refreshing image
        if (raiseEvent && listener) {
            listener->panelChanged(EvlocallabToolRemovedWithoutRefresh,
                                   toolName + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

bool LocallabTool::isLocallabToolAdded()
{
    return exp->get_visible();
}

void LocallabTool::refChanged(const double huer, const double lumar, const double chromar)
{
    // Hue reference normalization (between 0 and 1)
    double normHuer = huer;
    float h = Color::huelab_to_huehsv2(normHuer);
    h += 1.f / 6.f;

    if (h > 1.f) {
        h -= 1.f;
    }

    normHuer = h;

    // Luma reference normalization (between 0 and 1)
    const double normLumar = lumar / 100.f;

    // Chroma reference normalization (between 0 and 1)
    const double normChromar = chromar / 137.4f;

    // Update mask curve backgrounds
    updateMaskBackground(normChromar, normLumar, normHuer);
}

void LocallabTool::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller)
{
    float R = 0.f;
    float G = 0.f;
    float B = 0.f;
    float x;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    switch (callerId) {
        case 1: // Mask CC shape (bottom bar) + Color CC/LC shape (left bar)
            Color::hsv2rgb01(float(valY), float(valX), 0.5f, R, G, B);

            break;

        case 2: // Mask HH shape (main curve and bottom bar)
            x = valX - 1.f / 6.f;

            if (x < 0.f) {
                x += 1.f;
            }

            x = log2lin(x, 3.f);
            Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);

            break;

        case 3: // Color LH/HH shape (main curve)
            Color::hsv2rgb01(float (valX), float (valY), 0.5f, R, G, B);

            break;

        case 4: // Color CC/LC shape (bottom bar)
            const float value = (1.f - 0.7f) * float (valX) + 0.7f;
            // Whole hue range
            // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
            Color::hsv2rgb01(float (valY), float (valX), value, R, G, B);
    }

    caller->ccRed = double (R);
    caller->ccGreen = double (G);
    caller->ccBlue = double (B);
}

void LocallabTool::disableListener()
{
    ToolPanel::disableListener();

    enaExpConn.block(true);
}
void LocallabTool::enableListener()
{
    ToolPanel::enableListener();

    enaExpConn.block(false);
}

bool LocallabTool::on_remove_change(GdkEventButton* event)
{
    if (event->button == GDK_BUTTON_PRIMARY) {
        // Remove Locallab tool raising an event
        removeLocallabTool(true);
    }

    return true; // No event propagation further (to avoid closing expander when mouse clicking on remove image)
}

void LocallabTool::foldThemAll(GdkEventButton* event)
{
    if (event->button == GDK_BUTTON_SECONDARY) {
        if (locToolListener) {
            (static_cast<Locallab*>(locToolListener))->foldAllButOne(this);
        }
    }
}

/* ==== LocallabColor ==== */
LocallabColor::LocallabColor():
    LocallabTool(this, M("TP_LOCALLAB_COLOR_TOOLNAME"), M("TP_LOCALLAB_COFR"), false),

    // Color & Light specific widgets
    curvactiv(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CURV")))),
    lightness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTNESS"), -100, 500, 1, 0))),
    contrast(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    gridFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABGRID")))),
    labgrid(Gtk::manage(new LabGrid(EvLocallabLabGridValue, M("TP_LOCALLAB_LABGRID_VALUES")))),
    gridMethod(Gtk::manage(new MyComboBoxText())),
    strengthgrid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRGRID"), 0, 100, 1, 30))),
    sensi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    structcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL1"), 0, 100, 1, 0))),
    blurcolde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    softradiuscol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), -10.0, 1000.0, 0.5, 0.))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    expgradcol(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPGRAD")))),
    strcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRLUM"), -4., 4., 0.05, 0.))),
    strcolab(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRCHRO"), -6., 6., 0.05, 0.))),
    strcolh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRHUE"), -6., 6., 0.05, 0.))),
    angcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    labqualcurv(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    qualitycurveMethod(Gtk::manage(new MyComboBoxText())),
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),
    llshape(static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    ccshape(static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "C(C)"))),
    clCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CH"))),
    clshape(static_cast<DiagonalCurveEditor*>(clCurveEditorG->addCurve(CT_Diagonal, "C(L)"))),
    lcshape(static_cast<DiagonalCurveEditor*>(clCurveEditorG->addCurve(CT_Diagonal, "L(C)"))),
    HCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_HLH"))),
    LHshape(static_cast<FlatCurveEditor*>(HCurveEditorG->addCurve(CT_Flat, "L(H)", nullptr, false, true))),
    H2CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_HLH"))),
    HHshape(static_cast<FlatCurveEditor*>(H2CurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true))),
    rgbCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_RGB"))),
    toneMethod(Gtk::manage(new MyComboBoxText())),
    rgbshape(static_cast<DiagonalCurveEditor*>(rgbCurveEditorG->addCurve(CT_Diagonal, "", toneMethod))),
    special(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_SPECIAL")))),
    expmaskcol1(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWC1")))),
    merMethod(Gtk::manage(new MyComboBoxText())),
    mask7(Gtk::manage(new ToolParamBlock())),
    mergecolMethod(Gtk::manage(new MyComboBoxText())),
    mercol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MERDCOL"), 0.0, 100.0, 0.5, 18.))),
    opacol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_OPACOL"), 0.0, 100.0, 0.5, 60.))),
    conthrcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTTHR"), 0.0, 100.0, 0.5, 0.))),
    gridmerFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABGRIDMERG")))),
    labgridmerg(Gtk::manage(new LabGrid(EvLocallabLabGridmergValue, M("TP_LOCALLAB_LABGRID_VALUES"), false))),
    merlucol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_MERLUCOL"), 0.0, 100.0, 0.5, 32., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png"))))),
    showmaskcolMethod(Gtk::manage(new MyComboBoxText())),
    showmaskcolMethodinv(Gtk::manage(new MyComboBoxText())),
    enaColorMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKCOL"))),
    CCmaskshape(static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskshape(static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskshape(static_cast<FlatCurveEditor *>(maskCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    strumaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUMASKCOL"), 0., 200., 0.1, 0.))),
    toolcol(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TOOLCOL")))),
    fftColorMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTCOL_MASK")))),
    contcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTCOL"), 0., 200., 0.5, 0.))),
    blurcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURCOL"), 0.2, 100., 0.5, 0.2))),
    blendmaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    lapmaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    shadmaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHAMASKCOL"), 0, 100, 1, 0))),
    maskHCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKH"))),
    HHhmaskshape(static_cast<FlatCurveEditor *>(maskHCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true))),
    mask2CurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskshape(static_cast<DiagonalCurveEditor*>(mask2CurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    mask2CurveEditorGwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVMASK"))),
    LLmaskcolshapewav(static_cast<FlatCurveEditor*>(mask2CurveEditorGwav->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    csThresholdcol(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLDBLUR"), 0, 9, 0, 0, 6, 5, 0, false)))
{
    float R, G, B;

    std::vector<GradientMilestone> six_shape;

    for (int i = 0; i < 6; i++) {
        const float x = static_cast<float>(i) * (1.f / 6.f);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        six_shape.emplace_back(x, R, G, B);
    }

    const LocallabParams::LocallabSpot defSpot;

    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Color & Light specific widgets
    curvactivConn = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::curvactivChanged));

    lightness->setAdjusterListener(this);

    if (showtooltip) {
        lightness->set_tooltip_text(M("TP_LOCALLAB_EXPCOLOR_TOOLTIP"));
    }

    contrast->setAdjusterListener(this);

    chroma->setAdjusterListener(this);

    gridFrame->set_label_align(0.025, 0.5);

    gridMethod->append(M("TP_LOCALLAB_GRIDONE"));
    gridMethod->append(M("TP_LOCALLAB_GRIDTWO"));
    gridMethod->set_active(0);
    gridMethodConn = gridMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::gridMethodChanged));

    strengthgrid->setAdjusterListener(this);

    sensi->setAdjusterListener(this);

    if (showtooltip) {
        sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    }

    structcol->setAdjusterListener(this);

    blurcolde->setAdjusterListener(this);

    softradiuscol->setLogScale(10, -10);
    softradiuscol->setAdjusterListener(this);

    inversConn = invers->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::inversChanged));

    setExpandAlignProperties(expgradcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    strcol->setAdjusterListener(this);

    if (showtooltip) {
        strcol->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
    }

    strcolab->setAdjusterListener(this);

    if (showtooltip) {
        strcolab->set_tooltip_text(M("TP_LOCALLAB_GRADSTRAB_TOOLTIP"));
    }

    strcolh->setAdjusterListener(this);

    if (showtooltip) {
        strcolh->set_tooltip_text(M("TP_LOCALLAB_GRADSTRHUE_TOOLTIP"));
    }

    angcol->setAdjusterListener(this);

    if (showtooltip) {
        angcol->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    }

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->set_active(0);

    if (showtooltip) {
        qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));
    }

    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::qualitycurveMethodChanged));

    llCurveEditorG->setCurveListener(this);

    llshape->setResetCurve(DiagonalCurveType(defSpot.llcurve.at(0)), defSpot.llcurve);

    if (showtooltip) {
        llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    llshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    llshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    ccshape->setResetCurve(DiagonalCurveType(defSpot.cccurve.at(0)), defSpot.cccurve);

    if (showtooltip) {
        ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    ccshape->setBottomBarColorProvider(this, 4);
    ccshape->setLeftBarColorProvider(this, 1);

    llCurveEditorG->curveListComplete();

    clCurveEditorG->setCurveListener(this);

    clshape->setResetCurve(DiagonalCurveType(defSpot.clcurve.at(0)), defSpot.clcurve);

    if (showtooltip) {
        clshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    clshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    clshape->setLeftBarColorProvider(this, 1);

    lcshape->setResetCurve(DiagonalCurveType(defSpot.lccurve.at(0)), defSpot.lccurve);

    if (showtooltip) {
        lcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    lcshape->setBottomBarColorProvider(this, 4);
    lcshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    clCurveEditorG->curveListComplete();

    HCurveEditorG->setCurveListener(this);

    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FlatCurveType(defSpot.LHcurve.at(0)), defSpot.LHcurve);

    if (showtooltip) {
        LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    LHshape->setCurveColorProvider(this, 3);
    LHshape->setBottomBarBgGradient(six_shape);

    HCurveEditorG->curveListComplete();

    H2CurveEditorG->setCurveListener(this);

    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FlatCurveType(defSpot.HHcurve.at(0)), defSpot.HHcurve);

    if (showtooltip) {
        HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    HHshape->setCurveColorProvider(this, 3);
    HHshape->setBottomBarBgGradient(six_shape);

    H2CurveEditorG->curveListComplete();

    rgbCurveEditorG->setCurveListener(this);

    toneMethod->append(M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneMethod->append(M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneMethod->append(M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneMethod->append(M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneMethod->set_active(0);
    toneMethodConn = toneMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::toneMethodChanged));

    rgbshape->setResetCurve(DiagonalCurveType(defSpot.rgbcurve.at(0)), defSpot.rgbcurve);

    if (showtooltip) {
        rgbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    rgbshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    rgbshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    rgbCurveEditorG->curveListComplete();

    specialConn  = special->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::specialChanged));

    setExpandAlignProperties(expmaskcol1, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    merMethod->append(M("TP_LOCALLAB_MRONE"));
    merMethod->append(M("TP_LOCALLAB_MRTWO"));
    merMethod->append(M("TP_LOCALLAB_MRTHR"));
    merMethod->append(M("TP_LOCALLAB_MRFOU"));
    merMethod->append(M("TP_LOCALLAB_MRFIV"));
    merMethod->set_active(0);
    merMethodConn = merMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::merMethodChanged));

    mergecolMethod->append(M("TP_LOCALLAB_MERONE"));
    mergecolMethod->append(M("TP_LOCALLAB_MERTWO"));
    mergecolMethod->append(M("TP_LOCALLAB_MERTHR"));
    mergecolMethod->append(M("TP_LOCALLAB_MERFOU"));
    mergecolMethod->append(M("TP_LOCALLAB_MERFIV"));
    mergecolMethod->append(M("TP_LOCALLAB_MERSIX"));
    mergecolMethod->append(M("TP_LOCALLAB_MERSEV"));
    mergecolMethod->append(M("TP_LOCALLAB_MERSEV0"));
    mergecolMethod->append(M("TP_LOCALLAB_MERSEV1"));
    mergecolMethod->append(M("TP_LOCALLAB_MERSEV2"));
    mergecolMethod->append(M("TP_LOCALLAB_MERHEI"));
    mergecolMethod->append(M("TP_LOCALLAB_MERNIN"));
    mergecolMethod->append(M("TP_LOCALLAB_MERTEN"));
    mergecolMethod->append(M("TP_LOCALLAB_MERELE"));
    mergecolMethod->append(M("TP_LOCALLAB_MERTWE"));
    mergecolMethod->append(M("TP_LOCALLAB_MERTHI"));
    mergecolMethod->append(M("TP_LOCALLAB_MERFOR"));
    mergecolMethod->append(M("TP_LOCALLAB_MERHUE"));
    mergecolMethod->append(M("TP_LOCALLAB_MERSAT"));
    mergecolMethod->append(M("TP_LOCALLAB_MERCOL"));
    mergecolMethod->append(M("TP_LOCALLAB_MERLUM"));
    mergecolMethod->set_active(0);
    mergecolMethodConn = mergecolMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::mergecolMethodChanged));

    mercol->setAdjusterListener(this);

    opacol->setAdjusterListener(this);

    if (showtooltip) {
        opacol->set_tooltip_text(M("TP_LOCALLAB_MERGEOPA_TOOLTIP"));
    }

    conthrcol->setAdjusterListener(this);

    if (showtooltip) {
        conthrcol->set_tooltip_text(M("TP_LOCALLAB_MERGEOPA_TOOLTIP"));
    }

    if (showtooltip) {
        gridmerFrame->set_tooltip_text(M("TP_LOCALLAB_GRIDFRAME_TOOLTIP"));
    }

    merlucol->setAdjusterListener(this);

    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWSTRUC"));
    showmaskcolMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskcolMethod->set_active(0);

    if (showtooltip) {
        showmaskcolMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskcolMethodConn  = showmaskcolMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::showmaskcolMethodChanged));

    showmaskcolMethodinv->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcolMethodinv->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcolMethodinv->set_active(0);

    if (showtooltip) {
        showmaskcolMethodinv->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskcolMethodConninv  = showmaskcolMethodinv->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::showmaskcolMethodChangedinv));

    enaColorMaskConn = enaColorMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::enaColorMaskChanged));

    maskCurveEditorG->setCurveListener(this);

    CCmaskshape->setIdentityValue(0.);
    CCmaskshape->setResetCurve(FlatCurveType(defSpot.CCmaskcurve.at(0)), defSpot.CCmaskcurve);

    if (showtooltip) {
        CCmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskshape->setBottomBarColorProvider(this, 1);

    LLmaskshape->setIdentityValue(0.);
    LLmaskshape->setResetCurve(FlatCurveType(defSpot.LLmaskcurve.at(0)), defSpot.LLmaskcurve);

    if (showtooltip) {
        LLmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    if (showtooltip) {
        LLmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskshape->setIdentityValue(0.);
    HHmaskshape->setResetCurve(FlatCurveType(defSpot.HHmaskcurve.at(0)), defSpot.HHmaskcurve);

    if (showtooltip) {
        HHmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskshape->setCurveColorProvider(this, 2);
    HHmaskshape->setBottomBarColorProvider(this, 2);

    maskCurveEditorG->curveListComplete();

    strumaskcol->setAdjusterListener(this);

    toolcolConn  = toolcol->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::toolcolChanged));

    fftColorMaskConn = fftColorMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::fftColorMaskChanged));

    contcol->setAdjusterListener(this);

    blurcol->setAdjusterListener(this);

    blendmaskcol->setAdjusterListener(this);

    radmaskcol->setLogScale(10, -10);
    radmaskcol->setAdjusterListener(this);

    if (showtooltip) {
        radmaskcol->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    lapmaskcol->setAdjusterListener(this);

    if (showtooltip) {
        lapmaskcol->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    chromaskcol->setAdjusterListener(this);

    gammaskcol->setAdjusterListener(this);

    slomaskcol->setAdjusterListener(this);

    shadmaskcol->setAdjusterListener(this);

    maskHCurveEditorG->setCurveListener(this);

    HHhmaskshape->setIdentityValue(0.);
    HHhmaskshape->setResetCurve(FlatCurveType(defSpot.HHhmaskcurve.at(0)), defSpot.HHhmaskcurve);

    if (showtooltip) {
        HHhmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    HHhmaskshape->setCurveColorProvider(this, 2);
    HHhmaskshape->setBottomBarColorProvider(this, 2);

    maskHCurveEditorG->curveListComplete();

    mask2CurveEditorG->setCurveListener(this);

    Lmaskshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskcurve.at(0)), defSpot.Lmaskcurve);

    if (showtooltip) {
        Lmaskshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmaskshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2CurveEditorG->curveListComplete();

    mask2CurveEditorGwav->setCurveListener(this);

    LLmaskcolshapewav->setIdentityValue(0.);
    LLmaskcolshapewav->setResetCurve(FlatCurveType(defSpot.LLmaskcolcurvewav.at(0)), defSpot.LLmaskcolcurvewav);

    if (showtooltip) {
        LLmaskcolshapewav->setTooltip(M("TP_LOCALLAB_LMASK_LEVEL_TOOLTIP"));
    }

    LLmaskcolshapewav->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2CurveEditorGwav->curveListComplete();

    csThresholdcol->setAdjusterListener(this);

    // Add Color & Light specific widgets to GUI
    Gtk::Frame* const superFrame = Gtk::manage(new Gtk::Frame());
    superFrame->set_label_align(0.025, 0.5);
    // superFrame->set_label_widget(*curvactiv);
    ToolParamBlock* const superBox = Gtk::manage(new ToolParamBlock());
    superBox->pack_start(*lightness);
    superBox->pack_start(*contrast);
    superBox->pack_start(*chroma);
    ToolParamBlock* const gridBox = Gtk::manage(new ToolParamBlock());
    gridBox->pack_start(*labgrid);
    gridBox->pack_start(*gridMethod);
    gridBox->pack_start(*strengthgrid);
    gridFrame->add(*gridBox);
    superBox->pack_start(*gridFrame);
    superFrame->add(*superBox);
    pack_start(*superFrame);
    pack_start(*sensi);
    pack_start(*structcol);

    if (complexsoft < 2) {
        pack_start(*blurcolde);
        pack_start(*softradiuscol);
    }

    pack_start(*invers);
    MyExpander* const expcurvcol = Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPCURV")));
    setExpandAlignProperties(expcurvcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    ToolParamBlock* const gradcolBox = Gtk::manage(new ToolParamBlock());
    gradcolBox->pack_start(*strcol);

    if (complexsoft < 2) {
        gradcolBox->pack_start(*strcolab);
    }

    if (complexsoft < 2) {
        gradcolBox->pack_start(*strcolh);
    }

    gradcolBox->pack_start(*angcol);
    expgradcol->add(*gradcolBox, false);
    pack_start(*expgradcol, false, false);
    ToolParamBlock* const curvBox = Gtk::manage(new ToolParamBlock());
    Gtk::HBox* const qualcurvbox = Gtk::manage(new Gtk::HBox());
    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);
    curvBox->pack_start(*qualcurvbox);
    curvBox->pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor

    if (complexsoft < 2) {
        curvBox->pack_start(*clCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        curvBox->pack_start(*HCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        curvBox->pack_start(*H2CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        curvBox->pack_start(*rgbCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        curvBox->pack_start(*special);
    }

    expcurvcol->add(*curvBox, false);
    pack_start(*expcurvcol, false, false);
    ToolParamBlock* const mask7Box = Gtk::manage(new ToolParamBlock());
    mask7Box->pack_start(*merMethod);
    Gtk::Frame* const merge1colFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_MERGE1COLFRA")));
    merge1colFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const mergecolBox = Gtk::manage(new ToolParamBlock());
    Gtk::HSeparator* const separatormer = Gtk::manage(new  Gtk::HSeparator());
    mergecolBox->pack_start(*separatormer, Gtk::PACK_SHRINK, 2);
    mergecolBox->pack_start(*mergecolMethod);
    mergecolBox->pack_start(*mercol);
    mergecolBox->pack_start(*opacol);
    mergecolBox->pack_start(*conthrcol);
    ToolParamBlock* const gridmerBox = Gtk::manage(new ToolParamBlock());
    gridmerFrame->set_label_align(0.025, 0.5);
    gridmerBox->pack_start(*labgridmerg);
    gridmerBox->pack_start(*merlucol);
    gridmerFrame->add(*gridmerBox);
    mergecolBox->pack_start(*gridmerFrame);

    if (complexsoft < 2) {
        merge1colFrame->add(*mergecolBox);
    }

    mask7->pack_start(*merge1colFrame);
    mask7Box->pack_start(*mask7);
    expmaskcol1->add(*mask7Box, false);

    if (complexsoft < 2) {
        pack_start(*expmaskcol1, false, false);
    }

    MyExpander* const expmaskcol = Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWC")));
    setExpandAlignProperties(expmaskcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expmaskcol->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    Gtk::Frame* const mergecolFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_MERGECOLFRA")));
    mergecolFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const maskcolBox = Gtk::manage(new ToolParamBlock());
    maskcolBox->pack_start(*showmaskcolMethod, Gtk::PACK_SHRINK, 4);
    maskcolBox->pack_start(*showmaskcolMethodinv, Gtk::PACK_SHRINK, 4);
    maskcolBox->pack_start(*enaColorMask, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*maskCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    Gtk::Frame* const struFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABSTRUM")));
    struFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const strumBox = Gtk::manage(new ToolParamBlock());

    if (complexsoft < 2) {
        strumBox->pack_start(*strumaskcol);
        strumBox->pack_start(*toolcol);
    }

    struFrame->add(*strumBox);

    if (complexsoft < 2) {
        maskcolBox->pack_start(*struFrame, Gtk::PACK_SHRINK, 0);
    }

    Gtk::Frame* const blurFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABBLURM")));
    blurFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const blurmBox = Gtk::manage(new ToolParamBlock());

    if (complexsoft < 2) {
        blurmBox->pack_start(*fftColorMask, Gtk::PACK_SHRINK, 0);
        blurmBox->pack_start(*contcol);
        blurmBox->pack_start(*blurcol);
    }

    blurFrame->add(*blurmBox);

    if (complexsoft < 2) {
        maskcolBox->pack_start(*blurFrame, Gtk::PACK_SHRINK, 0);
    }

    maskcolBox->pack_start(*blendmaskcol, Gtk::PACK_SHRINK, 0);
    Gtk::Frame* const toolcolFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_TOOLMASK")));
    toolcolFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const toolcolBox = Gtk::manage(new ToolParamBlock());
    toolcolBox->pack_start(*radmaskcol, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 1) {
        toolcolBox->pack_start(*lapmaskcol, Gtk::PACK_SHRINK, 0);
    }

    toolcolBox->pack_start(*chromaskcol, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        toolcolBox->pack_start(*gammaskcol, Gtk::PACK_SHRINK, 0);
        toolcolBox->pack_start(*slomaskcol, Gtk::PACK_SHRINK, 0);
        toolcolBox->pack_start(*shadmaskcol, Gtk::PACK_SHRINK, 0);
        toolcolBox->pack_start(*maskHCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    }

    toolcolBox->pack_start(*mask2CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor

    if (complexsoft < 1) {
        toolcolBox->pack_start(*mask2CurveEditorGwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        toolcolBox->pack_start(*csThresholdcol, Gtk::PACK_SHRINK, 0);
    }

    toolcolFrame->add(*toolcolBox);
    maskcolBox->pack_start(*toolcolFrame);
    mergecolFrame->add(*maskcolBox);
    expmaskcol->add(*mergecolFrame, false);
    pack_start(*expmaskcol, false, false);
}

LocallabColor::~LocallabColor()
{
    delete llCurveEditorG;
    delete clCurveEditorG;
    delete HCurveEditorG;
    delete H2CurveEditorG;
    delete rgbCurveEditorG;
    delete maskCurveEditorG;
    delete maskHCurveEditorG;
    delete mask2CurveEditorG;
    delete mask2CurveEditorGwav;
}

void LocallabColor::setListener(ToolPanelListener* tpl)
{
    LocallabTool::setListener(tpl);

    labgrid->setListener(tpl);
    labgridmerg->setListener(tpl);
}

void LocallabColor::resetMaskView()
{
    showmaskcolMethodConn.block(true);
    showmaskcolMethodConninv.block(true);

    showmaskcolMethod->set_active(0);
    showmaskcolMethodinv->set_active(0);

    showmaskcolMethodConn.block(false);
    showmaskcolMethodConninv.block(false);
}

void LocallabColor::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    colorMask = showmaskcolMethod->get_active_row_number();
    colorMaskinv = showmaskcolMethodinv->get_active_row_number();
}

void LocallabColor::disableListener()
{
    LocallabTool::disableListener();

    curvactivConn.block(true);
    gridMethodConn.block(true);
    inversConn.block(true);
    qualitycurveMethodConn.block(true);
    toneMethodConn.block(true);
    specialConn.block(true);
    merMethodConn.block(true);
    mergecolMethodConn.block(true);
    showmaskcolMethodConn.block(true);
    showmaskcolMethodConninv.block(true);
    enaColorMaskConn.block(true);
    toolcolConn.block(true);
    fftColorMaskConn.block(true);
}

void LocallabColor::enableListener()
{
    LocallabTool::enableListener();

    curvactivConn.block(false);
    gridMethodConn.block(false);
    inversConn.block(false);
    qualitycurveMethodConn.block(false);
    toneMethodConn.block(false);
    specialConn.block(false);
    merMethodConn.block(false);
    mergecolMethodConn.block(false);
    showmaskcolMethodConn.block(false);
    showmaskcolMethodConninv.block(false);
    enaColorMaskConn.block(false);
    toolcolConn.block(false);
    fftColorMaskConn.block(false);
}

void LocallabColor::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visicolor);
        exp->setEnabled(pp->locallab.spots.at(index).expcolor);

        curvactiv->set_active(pp->locallab.spots.at(index).curvactiv);
        lightness->setValue(pp->locallab.spots.at(index).lightness);
        contrast->setValue(pp->locallab.spots.at(index).contrast);
        chroma->setValue(pp->locallab.spots.at(index).chroma);
        labgrid->setParams(pp->locallab.spots.at(index).labgridALow / LocallabParams::LABGRIDL_CORR_MAX,
                           pp->locallab.spots.at(index).labgridBLow / LocallabParams::LABGRIDL_CORR_MAX,
                           pp->locallab.spots.at(index).labgridAHigh / LocallabParams::LABGRIDL_CORR_MAX,
                           pp->locallab.spots.at(index).labgridBHigh / LocallabParams::LABGRIDL_CORR_MAX,
                           false);

        if (pp->locallab.spots.at(index).gridMethod == "one") {
            gridMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).gridMethod == "two") {
            gridMethod->set_active(1);
        }

        strengthgrid->setValue(pp->locallab.spots.at(index).strengthgrid);
        sensi->setValue(pp->locallab.spots.at(index).sensi);
        structcol->setValue(pp->locallab.spots.at(index).structcol);

        if (complexsoft < 2) {
            blurcolde->setValue(pp->locallab.spots.at(index).blurcolde);
            softradiuscol->setValue(pp->locallab.spots.at(index).softradiuscol);
        } else {
            blurcolde->setValue(5.);
            softradiuscol->setValue(0.);
        }

        invers->set_active(pp->locallab.spots.at(index).invers);
        strcol->setValue(pp->locallab.spots.at(index).strcol);

        if (complexsoft < 2) {
            strcolab->setValue(pp->locallab.spots.at(index).strcolab);
            strcolh->setValue(pp->locallab.spots.at(index).strcolh);
        } else {
            strcolab->setValue(0.);
            strcolh->setValue(0.);
        }

        angcol->setValue(pp->locallab.spots.at(index).angcol);

        if (pp->locallab.spots.at(index).qualitycurveMethod == "none") {
            qualitycurveMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "std") {
            qualitycurveMethod->set_active(1);
        }

        llshape->setCurve(pp->locallab.spots.at(index).llcurve);
        ccshape->setCurve(pp->locallab.spots.at(index).cccurve);

        if (complexsoft < 2) {
            clshape->setCurve(pp->locallab.spots.at(index).clcurve);
            lcshape->setCurve(pp->locallab.spots.at(index).lccurve);
            LHshape->setCurve(pp->locallab.spots.at(index).LHcurve);
            HHshape->setCurve(pp->locallab.spots.at(index).HHcurve);
        } else {
            clshape->reset();
            lcshape->reset();
            LHshape->reset();
            HHshape->reset();
        }

        if (pp->locallab.spots.at(index).toneMethod == "one") {
            toneMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).toneMethod == "two") {
            toneMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).toneMethod == "thr") {
            toneMethod->set_active(2);
        } else if (pp->locallab.spots.at(index).toneMethod == "fou") {
            toneMethod->set_active(3);
        }

        if (complexsoft < 2) {
            rgbshape->setCurve(pp->locallab.spots.at(index).rgbcurve);
            special->set_active(pp->locallab.spots.at(index).special);

            if (pp->locallab.spots.at(index).merMethod == "mone") {
                merMethod->set_active(0);
            } else if (pp->locallab.spots.at(index).merMethod == "mtwo") {
                merMethod->set_active(1);
            } else if (pp->locallab.spots.at(index).merMethod == "mthr") {
                merMethod->set_active(2);
            } else if (pp->locallab.spots.at(index).merMethod == "mfou") {
                merMethod->set_active(3);
            } else if (pp->locallab.spots.at(index).merMethod == "mfiv") {
                merMethod->set_active(4);
            }
        } else {
            rgbshape->reset();
            special->set_active(false);
            merMethod->set_active(0);
        }

        if (pp->locallab.spots.at(index).mergecolMethod == "one") {
            mergecolMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "two") {
            mergecolMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "thr") {
            mergecolMethod->set_active(2);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "fou") {
            mergecolMethod->set_active(3);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "fiv") {
            mergecolMethod->set_active(4);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "six") {
            mergecolMethod->set_active(5);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "sev") {
            mergecolMethod->set_active(6);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "sev0") {
            mergecolMethod->set_active(7);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "sev1") {
            mergecolMethod->set_active(8);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "sev2") {
            mergecolMethod->set_active(9);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "hei") {
            mergecolMethod->set_active(10);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "nin") {
            mergecolMethod->set_active(11);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "ten") {
            mergecolMethod->set_active(12);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "ele") {
            mergecolMethod->set_active(13);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "twe") {
            mergecolMethod->set_active(14);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "thi") {
            mergecolMethod->set_active(15);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "for") {
            mergecolMethod->set_active(16);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "hue") {
            mergecolMethod->set_active(17);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "sat") {
            mergecolMethod->set_active(18);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "col") {
            mergecolMethod->set_active(19);
        } else if (pp->locallab.spots.at(index).mergecolMethod == "lum") {
            mergecolMethod->set_active(20);
        }

        mercol->setValue(pp->locallab.spots.at(index).mercol);
        opacol->setValue(pp->locallab.spots.at(index).opacol);
        conthrcol->setValue(pp->locallab.spots.at(index).conthrcol);
        labgridmerg->setParams(0, 0,
                               pp->locallab.spots.at(index).labgridAHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                               pp->locallab.spots.at(index).labgridBHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                               false);
        merlucol->setValue(pp->locallab.spots.at(index).merlucol);
        enaColorMask->set_active(pp->locallab.spots.at(index).enaColorMask);
        CCmaskshape->setCurve(pp->locallab.spots.at(index).CCmaskcurve);
        LLmaskshape->setCurve(pp->locallab.spots.at(index).LLmaskcurve);
        HHmaskshape->setCurve(pp->locallab.spots.at(index).HHmaskcurve);

        if (complexsoft < 2) {
            strumaskcol->setValue(pp->locallab.spots.at(index).strumaskcol);
            toolcol->set_active(pp->locallab.spots.at(index).toolcol);
            fftColorMask->set_active(pp->locallab.spots.at(index).fftColorMask);
            contcol->setValue(pp->locallab.spots.at(index).contcol);
            blurcol->setValue(pp->locallab.spots.at(index).blurcol);
        } else {
            strumaskcol->setValue(0.);
            toolcol->set_active(false);
            fftColorMask->set_active(false);
            contcol->setValue(0.);
            blurcol->setLimits(0.2, 100., 0.5, 0.2);
            blurcol->setValue(pp->locallab.spots.at(index).blurcol);
        }

        blendmaskcol->setValue(pp->locallab.spots.at(index).blendmaskcol);
        radmaskcol->setValue(pp->locallab.spots.at(index).radmaskcol);

        if (complexsoft < 2) {
            lapmaskcol->setValue(pp->locallab.spots.at(index).lapmaskcol);
        } else {
            lapmaskcol->setValue(0.);
        }

        chromaskcol->setValue(pp->locallab.spots.at(index).chromaskcol);

        if (complexsoft < 2) {
            gammaskcol->setValue(pp->locallab.spots.at(index).gammaskcol);
            slomaskcol->setValue(pp->locallab.spots.at(index).slomaskcol);
            shadmaskcol->setValue(pp->locallab.spots.at(index).shadmaskcol);
            HHhmaskshape->setCurve(pp->locallab.spots.at(index).HHhmaskcurve);
        } else {
            gammaskcol->setValue(1.);
            slomaskcol->setValue(0.);
            shadmaskcol->setValue(0.);
            HHhmaskshape->reset();
        }

        Lmaskshape->setCurve(pp->locallab.spots.at(index).Lmaskcurve);

        if (complexsoft == 0) {
            LLmaskcolshapewav->setCurve(pp->locallab.spots.at(index).LLmaskcolcurvewav);
        } else {
            LLmaskcolshapewav->reset();
        }

        csThresholdcol->setValue<int>(pp->locallab.spots.at(index).csthresholdcol);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to invers button state
    updateColorGUI1();

    // Update GUI according to merMethod combobox state
    updateColorGUI2();

    // Update GUI according to fftColorMash button state
    updateColorGUI3();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expcolor = exp->getEnabled();
        pp->locallab.spots.at(index).visicolor = exp->get_visible();

        pp->locallab.spots.at(index).curvactiv = curvactiv->get_active();
        pp->locallab.spots.at(index).lightness = lightness->getIntValue();
        pp->locallab.spots.at(index).contrast = contrast->getIntValue();
        pp->locallab.spots.at(index).chroma = chroma->getIntValue();
        labgrid->getParams(pp->locallab.spots.at(index).labgridALow,
                           pp->locallab.spots.at(index).labgridBLow,
                           pp->locallab.spots.at(index).labgridAHigh,
                           pp->locallab.spots.at(index).labgridBHigh);
        pp->locallab.spots.at(index).labgridALow *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).labgridAHigh *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).labgridBLow *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).labgridBHigh *= LocallabParams::LABGRIDL_CORR_MAX;

        if (gridMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).gridMethod = "one";
        } else if (gridMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).gridMethod = "two";
        }

        pp->locallab.spots.at(index).strengthgrid = strengthgrid->getIntValue();
        pp->locallab.spots.at(index).sensi = sensi->getIntValue();
        pp->locallab.spots.at(index).structcol = structcol->getIntValue();
        pp->locallab.spots.at(index).blurcolde = blurcolde->getIntValue();
        pp->locallab.spots.at(index).softradiuscol = softradiuscol->getValue();
        pp->locallab.spots.at(index).invers = invers->get_active();
        pp->locallab.spots.at(index).strcol = strcol->getValue();
        pp->locallab.spots.at(index).strcolab = strcolab->getValue();
        pp->locallab.spots.at(index).strcolh = strcolh->getValue();
        pp->locallab.spots.at(index).angcol = angcol->getValue();

        if (qualitycurveMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).qualitycurveMethod = "none";
        } else if (qualitycurveMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).qualitycurveMethod = "std";
        }

        pp->locallab.spots.at(index).llcurve = llshape->getCurve();
        pp->locallab.spots.at(index).cccurve = ccshape->getCurve();
        pp->locallab.spots.at(index).clcurve = clshape->getCurve();
        pp->locallab.spots.at(index).lccurve = lcshape->getCurve();
        pp->locallab.spots.at(index).LHcurve = LHshape->getCurve();
        pp->locallab.spots.at(index).HHcurve = HHshape->getCurve();

        if (toneMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).toneMethod = "one";
        } else if (toneMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).toneMethod = "two";
        } else if (toneMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).toneMethod = "thr";
        } else if (toneMethod->get_active_row_number() == 3) {
            pp->locallab.spots.at(index).toneMethod = "fou";
        }

        pp->locallab.spots.at(index).rgbcurve = rgbshape->getCurve();
        pp->locallab.spots.at(index).special = special->get_active();

        if (merMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).merMethod = "mone";
        } else if (merMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).merMethod = "mtwo";
        } else if (merMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).merMethod = "mthr";
        } else if (merMethod->get_active_row_number() == 3) {
            pp->locallab.spots.at(index).merMethod = "mfou";
        } else if (merMethod->get_active_row_number() == 4) {
            pp->locallab.spots.at(index).merMethod = "mfiv";
        }

        if (mergecolMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).mergecolMethod = "one";
        } else if (mergecolMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).mergecolMethod = "two";
        } else if (mergecolMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).mergecolMethod = "thr";
        } else if (mergecolMethod->get_active_row_number() == 3) {
            pp->locallab.spots.at(index).mergecolMethod = "fou";
        } else if (mergecolMethod->get_active_row_number() == 4) {
            pp->locallab.spots.at(index).mergecolMethod = "fiv";
        } else if (mergecolMethod->get_active_row_number() == 5) {
            pp->locallab.spots.at(index).mergecolMethod = "six";
        } else if (mergecolMethod->get_active_row_number() == 6) {
            pp->locallab.spots.at(index).mergecolMethod = "sev";
        } else if (mergecolMethod->get_active_row_number() == 7) {
            pp->locallab.spots.at(index).mergecolMethod = "sev0";
        } else if (mergecolMethod->get_active_row_number() == 8) {
            pp->locallab.spots.at(index).mergecolMethod = "sev1";
        } else if (mergecolMethod->get_active_row_number() == 9) {
            pp->locallab.spots.at(index).mergecolMethod = "sev2";
        } else if (mergecolMethod->get_active_row_number() == 10) {
            pp->locallab.spots.at(index).mergecolMethod = "hei";
        } else if (mergecolMethod->get_active_row_number() == 11) {
            pp->locallab.spots.at(index).mergecolMethod = "nin";
        } else if (mergecolMethod->get_active_row_number() == 12) {
            pp->locallab.spots.at(index).mergecolMethod = "ten";
        } else if (mergecolMethod->get_active_row_number() == 13) {
            pp->locallab.spots.at(index).mergecolMethod = "ele";
        } else if (mergecolMethod->get_active_row_number() == 14) {
            pp->locallab.spots.at(index).mergecolMethod = "twe";
        } else if (mergecolMethod->get_active_row_number() == 15) {
            pp->locallab.spots.at(index).mergecolMethod = "thi";
        } else if (mergecolMethod->get_active_row_number() == 16) {
            pp->locallab.spots.at(index).mergecolMethod = "for";
        } else if (mergecolMethod->get_active_row_number() == 17) {
            pp->locallab.spots.at(index).mergecolMethod = "hue";
        } else if (mergecolMethod->get_active_row_number() == 18) {
            pp->locallab.spots.at(index).mergecolMethod = "sat";
        } else if (mergecolMethod->get_active_row_number() == 19) {
            pp->locallab.spots.at(index).mergecolMethod = "col";
        } else if (mergecolMethod->get_active_row_number() == 20) {
            pp->locallab.spots.at(index).mergecolMethod = "lum";
        }

        pp->locallab.spots.at(index).mercol = mercol->getValue();
        pp->locallab.spots.at(index).opacol = opacol->getValue();
        pp->locallab.spots.at(index).conthrcol = conthrcol->getValue();
        labgridmerg->getParams(pp->locallab.spots.at(index).labgridALowmerg,
                               pp->locallab.spots.at(index).labgridBLowmerg,
                               pp->locallab.spots.at(index).labgridAHighmerg,
                               pp->locallab.spots.at(index).labgridBHighmerg);
        pp->locallab.spots.at(index).labgridALowmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).labgridAHighmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).labgridBLowmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).labgridBHighmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        pp->locallab.spots.at(index).merlucol = merlucol->getValue();
        pp->locallab.spots.at(index).enaColorMask = enaColorMask->get_active();
        pp->locallab.spots.at(index).CCmaskcurve = CCmaskshape->getCurve();
        pp->locallab.spots.at(index).LLmaskcurve = LLmaskshape->getCurve();
        pp->locallab.spots.at(index).HHmaskcurve = HHmaskshape->getCurve();
        pp->locallab.spots.at(index).strumaskcol = strumaskcol->getValue();
        pp->locallab.spots.at(index).toolcol = toolcol->get_active();
        pp->locallab.spots.at(index).fftColorMask = fftColorMask->get_active();
        pp->locallab.spots.at(index).contcol = contcol->getValue();
        pp->locallab.spots.at(index).blurcol = blurcol->getValue();
        pp->locallab.spots.at(index).blendmaskcol = blendmaskcol->getIntValue();
        pp->locallab.spots.at(index).radmaskcol = radmaskcol->getValue();
        pp->locallab.spots.at(index).lapmaskcol = lapmaskcol->getValue();
        pp->locallab.spots.at(index).chromaskcol = chromaskcol->getValue();
        pp->locallab.spots.at(index).gammaskcol = gammaskcol->getValue();
        pp->locallab.spots.at(index).slomaskcol = slomaskcol->getValue();
        pp->locallab.spots.at(index).shadmaskcol = shadmaskcol->getIntValue();
        pp->locallab.spots.at(index).HHhmaskcurve = HHhmaskshape->getCurve();
        pp->locallab.spots.at(index).Lmaskcurve = Lmaskshape->getCurve();
        pp->locallab.spots.at(index).LLmaskcolcurvewav = LLmaskcolshapewav->getCurve();
        pp->locallab.spots.at(index).csthresholdcol = csThresholdcol->getValue<int>();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster, labgrid and threshold adjuster widgets
        lightness->setDefault((double)defSpot.lightness);
        contrast->setDefault((double)defSpot.contrast);
        chroma->setDefault((double)defSpot.chroma);
        labgrid->setDefault(defSpot.labgridALow / LocallabParams::LABGRIDL_CORR_MAX,
                            defSpot.labgridBLow / LocallabParams::LABGRIDL_CORR_MAX,
                            defSpot.labgridAHigh / LocallabParams::LABGRIDL_CORR_MAX,
                            defSpot.labgridBHigh / LocallabParams::LABGRIDL_CORR_MAX);
        strengthgrid->setDefault((double) defSpot.strengthgrid);
        sensi->setDefault((double)defSpot.sensi);
        structcol->setDefault((double)defSpot.structcol);
        blurcolde->setDefault((double)defSpot.blurcolde);
        softradiuscol->setDefault(defSpot.softradiuscol);
        strcol->setDefault(defSpot.strcol);
        strcolab->setDefault(defSpot.strcolab);
        strcolh->setDefault(defSpot.strcolh);
        angcol->setDefault(defSpot.angcol);
        mercol->setDefault(defSpot.mercol);
        opacol->setDefault(defSpot.opacol);
        conthrcol->setDefault(defSpot.conthrcol);
        labgridmerg->setDefault(defSpot.labgridALowmerg / LocallabParams::LABGRIDL_CORR_MAX,
                                defSpot.labgridBLowmerg / LocallabParams::LABGRIDL_CORR_MAX,
                                defSpot.labgridAHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                                defSpot.labgridBHighmerg / LocallabParams::LABGRIDL_CORR_MAX);
        merlucol->setDefault(defSpot.merlucol);
        strumaskcol->setDefault(defSpot.strumaskcol);
        contcol->setDefault(defSpot.contcol);
        blurcol->setDefault(defSpot.blurcol);
        blendmaskcol->setDefault((double)defSpot.blendmaskcol);
        radmaskcol->setDefault(defSpot.radmaskcol);
        lapmaskcol->setDefault(defSpot.lapmaskcol);
        chromaskcol->setDefault(defSpot.chromaskcol);
        gammaskcol->setDefault(defSpot.gammaskcol);
        slomaskcol->setDefault(defSpot.slomaskcol);
        shadmaskcol->setDefault((double)defSpot.shadmaskcol);
        csThresholdcol->setDefault<int>(defSpot.csthresholdcol);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == lightness) {
            if (listener) {
                listener->panelChanged(Evlocallablightness,
                                       lightness->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == contrast) {
            if (listener) {
                listener->panelChanged(Evlocallabcontrast,
                                       contrast->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroma) {
            if (listener) {
                listener->panelChanged(Evlocallabchroma,
                                       chroma->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strengthgrid) {
            if (listener) {
                listener->panelChanged(EvLocallabLabstrengthgrid,
                                       strengthgrid->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensi) {
            if (listener) {
                listener->panelChanged(Evlocallabsensi,
                                       sensi->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == structcol) {
            if (listener) {
                listener->panelChanged(Evlocallabstructcol,
                                       structcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blurcolde) {
            if (listener) {
                listener->panelChanged(Evlocallabblurcolde,
                                       blurcolde->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiuscol) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiuscol,
                                       softradiuscol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strcol) {
            if (listener) {
                listener->panelChanged(Evlocallabstrcol,
                                       strcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strcolab) {
            if (listener) {
                listener->panelChanged(Evlocallabstrcolab,
                                       strcolab->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strcolh) {
            if (listener) {
                listener->panelChanged(Evlocallabstrcolh,
                                       strcolh->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == angcol) {
            if (listener) {
                listener->panelChanged(Evlocallabangcol,
                                       angcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == mercol) {
            if (listener) {
                listener->panelChanged(Evlocallabmercol,
                                       mercol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == opacol) {
            if (listener) {
                listener->panelChanged(Evlocallabopacol,
                                       opacol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == conthrcol) {
            if (listener) {
                listener->panelChanged(Evlocallabconthrcol,
                                       conthrcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == merlucol) {
            if (listener) {
                listener->panelChanged(Evlocallabmerlucol,
                                       merlucol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strumaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabstrumaskcol,
                                       strumaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == contcol) {
            if (listener) {
                listener->panelChanged(Evlocallabcontcol,
                                       contcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blurcol) {
            if (listener) {
                listener->panelChanged(Evlocallabblurcol,
                                       blurcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcol,
                                       blendmaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcol,
                                       radmaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskcol,
                                       lapmaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcol,
                                       chromaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcol,
                                       gammaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcol,
                                       slomaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shadmaskcol) {
            if (listener) {
                listener->panelChanged(Evlocallabshadmaskcol,
                                       shadmaskcol->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == csThresholdcol) {
            if (listener) {
                listener->panelChanged(EvlocallabcsThresholdcol,
                                       csThresholdcol->getHistoryString() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == llshape) {
            if (listener) {
                listener->panelChanged(Evlocallabllshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == ccshape) {
            if (listener) {
                listener->panelChanged(Evlocallabccshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == clshape) {
            if (listener) {
                listener->panelChanged(Evlocallabclshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == lcshape) {
            if (listener) {
                listener->panelChanged(Evlocallablcshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == rgbshape) {
            if (listener) {
                listener->panelChanged(Evlocallabrgbshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHhmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHhmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmaskshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskcolshapewav) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskcolshapewav,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenacolor,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenacolor,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskshape->updateLocallabBackground(normChromar);
        LLmaskshape->updateLocallabBackground(normLumar);
        HHmaskshape->updateLocallabBackground(normHuer);
        HHhmaskshape->updateLocallabBackground(normHuer);

        return false;
    }
    );
}

void LocallabColor::curvactivChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (curvactiv->get_active()) {
                listener->panelChanged(Evlocallabcurvactiv,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabcurvactiv,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::gridMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabgridMethod,
                                   gridMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabColor::inversChanged()
{
    updateColorGUI1(); // Update GUI according to invers button state

    // This event is called to transmit potentially resetted mask state
    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (invers->get_active()) {
                listener->panelChanged(Evlocallabinvers,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabinvers,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::qualitycurveMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabqualitycurveMethod,
                                   qualitycurveMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabColor::toneMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabtoneMethod,
                                   toneMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabColor::specialChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (special->get_active()) {
                listener->panelChanged(EvLocallabspecial,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabspecial,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::merMethodChanged()
{
    updateColorGUI2(); // Update GUI according to merMethod combobox state

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabmerMethod,
                                   merMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabColor::mergecolMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabmergecolMethod,
                                   mergecolMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabColor::showmaskcolMethodChanged()
{
    // If mask preview is activated, deactivate other Color & Light mask preview
    showmaskcolMethodConninv.block(true);
    showmaskcolMethodinv->set_active(0);
    showmaskcolMethodConninv.block(false);

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabColor::showmaskcolMethodChangedinv()
{
    // If mask preview is activated, deactivate other Color & Light mask preview
    showmaskcolMethodConn.block(true);
    showmaskcolMethod->set_active(0);
    showmaskcolMethodConn.block(false);

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabColor::enaColorMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaColorMask->get_active()) {
                listener->panelChanged(EvLocallabEnaColorMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaColorMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::toolcolChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (toolcol->get_active()) {
                listener->panelChanged(EvLocallabtoolcol,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabtoolcol,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::fftColorMaskChanged()
{
    updateColorGUI3(); // Update GUI according to fftColorMash button state

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fftColorMask->get_active()) {
                listener->panelChanged(EvLocallabfftColorMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabfftColorMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::updateColorGUI1()
{
    if (invers->get_active()) {
        gridFrame->hide();
        structcol->hide();
        softradiuscol->hide();
        expgradcol->hide();
        labqualcurv->hide();
        qualitycurveMethod->hide();
        clCurveEditorG->hide();
        HCurveEditorG->hide();
        expmaskcol1->hide();
        showmaskcolMethod->hide();
        // Reset hidden mask combobox
        showmaskcolMethodConn.block(true);
        showmaskcolMethod->set_active(0);
        showmaskcolMethodConn.block(false);
        showmaskcolMethodinv->show();
        contcol->hide();
        blurcol->hide();
    } else {
        gridFrame->show();
        structcol->show();
        softradiuscol->show();
        expgradcol->show();
        labqualcurv->show();
        qualitycurveMethod->show();
        clCurveEditorG->show();
        HCurveEditorG->show();
        expmaskcol1->show();
        showmaskcolMethod->show();
        showmaskcolMethodinv->hide();
        // Reset hidden mask combobox
        showmaskcolMethodConninv.block(true);
        showmaskcolMethodinv->set_active(0);
        showmaskcolMethodConninv.block(false);
        contcol->show();
        blurcol->show();
    }
}

void LocallabColor::updateColorGUI2()
{
    // Note: When a merMethod is selected, invers button is set insensitive to avoid the combobox to disappear
    switch (merMethod->get_active_row_number()) {
        case 0:
            invers->set_sensitive(true);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(true);
            special->set_sensitive(true);
            mask7->hide();
            conthrcol->hide();
            gridmerFrame->hide();
            break;

        case 1:
            invers->set_sensitive(false);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(true);
            special->set_sensitive(true);
            mask7->hide();
            conthrcol->hide();
            gridmerFrame->hide();
            break;

        case 2:
            sensi->set_sensitive(false);
            invers->set_sensitive(false);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(false);
            special->set_sensitive(false);
            mask7->show();
            conthrcol->show();
            gridmerFrame->hide();
            break;

        case 3:
            invers->set_sensitive(false);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(false);
            special->set_sensitive(false);
            mask7->show();
            conthrcol->show();
            gridmerFrame->hide();
            break;

        case 4:
            invers->set_sensitive(false);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(false);
            special->set_sensitive(false);
            mask7->show();
            conthrcol->hide();
            gridmerFrame->show();
    }
}

void LocallabColor::updateColorGUI3()
{
    const double temp = blurcol->getValue();

    if (fftColorMask->get_active()) {
        blurcol->setLimits(0.2, 1000., 0.5, 0.2);
    } else {
        blurcol->setLimits(0.2, 100., 0.5, 0.2);
    }

    blurcol->setValue(temp);
}

/* ==== LocallabExposure ==== */
LocallabExposure::LocallabExposure():
    LocallabTool(this, M("TP_LOCALLAB_EXP_TOOLNAME"), M("TP_LOCALLAB_EXPOSE"), false),

    // Exposure specific widgets
    expMethod(Gtk::manage(new MyComboBoxText())),
    pdeFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_PDEFRA")))),
    laplacexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPLACEXP"), 0.0, 100.0, 0.1, 0.))),
    linear(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LINEAR"), 0., 1., 0.01, 0.3))),
    balanexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BALANEXP"), 0.5, 1.5, 0.01, 1.0))),
    gamm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMM"), 0.2, 1.3, 0.01, 0.4))),
    ctboxexpmethod(Gtk::manage(new Gtk::HBox())),
    exnoiseMethod(Gtk::manage(new MyComboBoxText())),
    fatFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_FATFRA")))),
    fatamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATAMOUNT"), 1., 100., 1., 1.))),
    fatdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATDETAIL"), -100., 300., 1., 0.))),
    fatlevel(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATLEVEL"), 0.25, 2.5, 0.05, 1.))),
    fatanchor(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATANCHORA"), 0.1, 3.0, 0.05, 1.))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    structexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    blurexpde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    exptoolexp(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPTOOL")))),
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -2.0, 3.0, 0.05, 0.0))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 10, 0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0))),
    shadex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHADEX"), 0, 100, 1, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHADEXCOMP"), 0, 100, 1, 50))),
    expchroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_EXPCHROMA"), -50, 100, 1, 30))),
    curveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"))),
    shapeexpos(static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""))),
    expgradexp(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPGRAD")))),
    strexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -4., 4., 0.05, 0.))),
    angexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    softradiusexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), -10.0, 1000.0, 0.5, 0.))),
    inversex(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    expmaskexp(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWE")))),
    showmaskexpMethod(Gtk::manage(new MyComboBoxText())),
    showmaskexpMethodinv(Gtk::manage(new MyComboBoxText())),
    enaExpMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    enaExpMaskaft(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASKAFT")))),
    maskexpCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskexpshape(static_cast<FlatCurveEditor*>(maskexpCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskexpshape(static_cast<FlatCurveEditor*>(maskexpCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskexpshape(static_cast<FlatCurveEditor *>(maskexpCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    lapmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    gradFramemask(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRADFRA")))),
    strmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -2., 2., 0.05, 0.))),
    angmaskexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180., 180., 0.1, 0.))),
    mask2expCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskexpshape(static_cast<DiagonalCurveEditor*>(mask2expCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    const LocallabParams::LocallabSpot defSpot;

    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Exposure specific widgets
    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPOSURE_TOOLTIP"));
    }

    expMethod->append(M("TP_LOCALLAB_STD"));

    if (complexsoft == 1) {
        expMethod->append(M("TP_LOCALLAB_PDE"));
    }

    expMethod->set_active(0);
    expMethodConn = expMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::expMethodChanged));

    if (showtooltip) {
        expMethod->set_tooltip_text(M("TP_LOCALLAB_EXPMETHOD_TOOLTIP"));
    }

    pdeFrame->set_label_align(0.025, 0.5);

    if (showtooltip) {
        pdeFrame->set_tooltip_text(M("TP_LOCALLAB_PDEFRAME_TOOLTIP"));
    }

    laplacexp->setAdjusterListener(this);

    linear->setAdjusterListener(this);

    balanexp->setAdjusterListener(this);

    gamm->setAdjusterListener(this);

    exnoiseMethod->append(M("TP_LOCALLAB_NONENOISE"));
    exnoiseMethod->append(M("TP_LOCALLAB_MEDIAN"));
    exnoiseMethod->append(M("TP_LOCALLAB_WEDIANHI"));
    exnoiseMethod->set_active(0);
    exnoiseMethodConn  = exnoiseMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::exnoiseMethodChanged));

    if (showtooltip) {
        exnoiseMethod->set_tooltip_text(M("TP_LOCALLAB_EXPMETHOD_TOOLTIP"));
    }

    fatFrame->set_label_align(0.025, 0.5);

    if (showtooltip) {
        fatFrame->set_tooltip_text(M("TP_LOCALLAB_FATFRAME_TOOLTIP"));
    }

    fatamount->setAdjusterListener(this);

    fatdetail->setAdjusterListener(this);

    fatlevel->setAdjusterListener(this);

    fatanchor->setAdjusterListener(this);

    sensiex->setAdjusterListener(this);

    if (showtooltip) {
        sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    }

    structexp->setAdjusterListener(this);

    blurexpde->setAdjusterListener(this);

    setExpandAlignProperties(exptoolexp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    expcomp->setAdjusterListener(this);

    black->setAdjusterListener(this);

    hlcompr->setAdjusterListener(this);

    hlcomprthresh->setAdjusterListener(this);

    shadex->setAdjusterListener(this);

    shcompr->setAdjusterListener(this);

    expchroma->setAdjusterListener(this);

    curveEditorG->setCurveListener(this);

    shapeexpos->setResetCurve(DiagonalCurveType(defSpot.excurve.at(0)), defSpot.excurve);

    if (showtooltip) {
        shapeexpos->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
    }

    shapeexpos->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});
    shapeexpos->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});

    curveEditorG->curveListComplete();

    setExpandAlignProperties(expgradexp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    strexp->setAdjusterListener(this);

    if (showtooltip) {
        strexp->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
    }

    angexp->setAdjusterListener(this);

    if (showtooltip) {
        angexp->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    }

    softradiusexp->setLogScale(10, -10);
    softradiusexp->setAdjusterListener(this);

    inversexConn  = inversex->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::inversexChanged));

    setExpandAlignProperties(expmaskexp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expmaskexp->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWSTRUCEX"));
    showmaskexpMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskexpMethod->set_active(0);

    if (showtooltip) {
        showmaskexpMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskexpMethodConn  = showmaskexpMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::showmaskexpMethodChanged));

    showmaskexpMethodinv->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskexpMethodinv->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskexpMethodinv->set_active(0);

    if (showtooltip) {
        showmaskexpMethodinv->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskexpMethodConninv  = showmaskexpMethodinv->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::showmaskexpMethodChangedinv));

    enaExpMaskConn = enaExpMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::enaExpMaskChanged));

    enaExpMaskaftConn = enaExpMaskaft->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::enaExpMaskaftChanged));

    maskexpCurveEditorG->setCurveListener(this);

    CCmaskexpshape->setIdentityValue(0.);
    CCmaskexpshape->setResetCurve(FlatCurveType(defSpot.CCmaskexpcurve.at(0)), defSpot.CCmaskexpcurve);

    if (showtooltip) {
        CCmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskexpshape->setBottomBarColorProvider(this, 1);

    LLmaskexpshape->setIdentityValue(0.);
    LLmaskexpshape->setResetCurve(FlatCurveType(defSpot.LLmaskexpcurve.at(0)), defSpot.LLmaskexpcurve);

    if (showtooltip) {
        LLmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskexpshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});

    HHmaskexpshape->setIdentityValue(0.);
    HHmaskexpshape->setResetCurve(FlatCurveType(defSpot.HHmaskexpcurve.at(0)), defSpot.HHmaskexpcurve);

    if (showtooltip) {
        HHmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskexpshape->setCurveColorProvider(this, 2);
    HHmaskexpshape->setBottomBarColorProvider(this, 2);

    maskexpCurveEditorG->curveListComplete();

    blendmaskexp->setAdjusterListener(this);

    radmaskexp->setLogScale(10, -10);
    radmaskexp->setAdjusterListener(this);

    if (showtooltip) {
        radmaskexp->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    lapmaskexp->setAdjusterListener(this);

    if (showtooltip) {
        lapmaskexp->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    chromaskexp->setAdjusterListener(this);

    gammaskexp->setAdjusterListener(this);

    slomaskexp->setAdjusterListener(this);

    gradFramemask->set_label_align(0.025, 0.5);

    strmaskexp->setAdjusterListener(this);

    if (showtooltip) {
        strmaskexp->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
    }

    angmaskexp->setAdjusterListener(this);

    if (showtooltip) {
        angmaskexp->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    }

    mask2expCurveEditorG->setCurveListener(this);

    Lmaskexpshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskexpcurve.at(0)), defSpot.Lmaskexpcurve);

    if (showtooltip) {
        Lmaskexpshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmaskexpshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});
    Lmaskexpshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});

    mask2expCurveEditorG->curveListComplete();

    // Add Color & Light specific widgets to GUI
    if (complexsoft < 2) {
        pack_start(*expMethod);
    }

    ToolParamBlock* const pdeBox = Gtk::manage(new ToolParamBlock());
    pdeBox->pack_start(*laplacexp);
    pdeBox->pack_start(*linear);
    pdeBox->pack_start(*balanexp);
    pdeBox->pack_start(*gamm);
    Gtk::Label* const labelexpmethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_NOISEMETH") + ":"));
    ctboxexpmethod->pack_start(*labelexpmethod, Gtk::PACK_SHRINK, 4);
    ctboxexpmethod->pack_start(*exnoiseMethod);
    pdeBox->pack_start(*ctboxexpmethod);
    pdeFrame->add(*pdeBox);

    if (complexsoft < 1) {
        pack_start(*pdeFrame);
    }

    ToolParamBlock* const fatBox = Gtk::manage(new ToolParamBlock());
    fatBox->pack_start(*fatamount);
    fatBox->pack_start(*fatdetail);

    if (complexsoft < 2) {
        fatBox->pack_start(*fatlevel);
    }

    fatBox->pack_start(*fatanchor);
    fatFrame->add(*fatBox);
    pack_start(*fatFrame);
    pack_start(*sensiex);
    pack_start(*structexp);

    if (complexsoft < 2) {
        pack_start(*blurexpde);
    }

    ToolParamBlock* const toolBox = Gtk::manage(new ToolParamBlock());
    toolBox->pack_start(*expcomp);
    toolBox->pack_start(*black);

    if (complexsoft < 2) {
        toolBox->pack_start(*hlcompr);
        toolBox->pack_start(*hlcomprthresh);
        toolBox->pack_start(*shadex);
        toolBox->pack_start(*shcompr);
        toolBox->pack_start(*expchroma);
    }

    toolBox->pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    exptoolexp->add(*toolBox, false);
    pack_start(*exptoolexp);
    ToolParamBlock* const gradBox = Gtk::manage(new ToolParamBlock());
    gradBox->pack_start(*strexp);
    gradBox->pack_start(*angexp);
    expgradexp->add(*gradBox, false);
    pack_start(*expgradexp);

    if (complexsoft < 2) {
        pack_start(*softradiusexp);
        pack_start(*inversex);
    }

    ToolParamBlock* const maskexpBox = Gtk::manage(new ToolParamBlock());
    maskexpBox->pack_start(*showmaskexpMethod, Gtk::PACK_SHRINK, 4);
    maskexpBox->pack_start(*showmaskexpMethodinv, Gtk::PACK_SHRINK, 4);
    maskexpBox->pack_start(*enaExpMask, Gtk::PACK_SHRINK, 0);
    // maskexpBox->pack_start(*enaExpMaskaft, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*maskexpCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskexpBox->pack_start(*blendmaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*radmaskexp, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 1) {
        maskexpBox->pack_start(*lapmaskexp, Gtk::PACK_SHRINK, 0);
    }

    maskexpBox->pack_start(*chromaskexp, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        maskexpBox->pack_start(*gammaskexp, Gtk::PACK_SHRINK, 0);
        maskexpBox->pack_start(*slomaskexp, Gtk::PACK_SHRINK, 0);
    }

    ToolParamBlock* const gradmaskBox = Gtk::manage(new ToolParamBlock());
    gradmaskBox->pack_start(*strmaskexp);
    gradmaskBox->pack_start(*angmaskexp);
    gradFramemask->add(*gradmaskBox);

    if (complexsoft < 2) {
        maskexpBox->pack_start(*gradFramemask);
    }

    maskexpBox->pack_start(*mask2expCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmaskexp->add(*maskexpBox, false);
    pack_start(*expmaskexp, false, false);
}

LocallabExposure::~LocallabExposure()
{
    delete curveEditorG;
    delete maskexpCurveEditorG;
    delete mask2expCurveEditorG;
}

void LocallabExposure::resetMaskView()
{
    showmaskexpMethodConn.block(true);
    showmaskexpMethodConninv.block(true);

    showmaskexpMethod->set_active(0);
    showmaskexpMethodinv->set_active(0);

    showmaskexpMethodConn.block(false);
    showmaskexpMethodConninv.block(false);
}

void LocallabExposure::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    expMask = showmaskexpMethod->get_active_row_number();
    expMaskinv = showmaskexpMethodinv->get_active_row_number();
}

void LocallabExposure::disableListener()
{
    LocallabTool::disableListener();

    expMethodConn.block(true);
    exnoiseMethodConn.block(true);
    inversexConn.block(true);
    showmaskexpMethodConn.block(true);
    showmaskexpMethodConninv.block(true);
    enaExpMaskConn.block(true);
    enaExpMaskaftConn.block(true);
}

void LocallabExposure::enableListener()
{
    LocallabTool::enableListener();

    expMethodConn.block(false);
    exnoiseMethodConn.block(false);
    inversexConn.block(false);
    showmaskexpMethodConn.block(false);
    showmaskexpMethodConninv.block(false);
    enaExpMaskConn.block(false);
    enaExpMaskaftConn.block(false);
}

void LocallabExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visiexpose);
        exp->setEnabled(pp->locallab.spots.at(index).expexpose);

        if (complexsoft < 2) {
            if (pp->locallab.spots.at(index).expMethod == "std") {
                expMethod->set_active(0);
            } else if (pp->locallab.spots.at(index).expMethod == "pde") {
                expMethod->set_active(1);
            }
        } else {
            expMethod->set_active(0);
        }

        if (complexsoft == 0) {
            laplacexp->setValue(pp->locallab.spots.at(index).laplacexp);
        } else {
            laplacexp->setValue(0.);
        }

        linear->setValue(pp->locallab.spots.at(index).linear);
        balanexp->setValue(pp->locallab.spots.at(index).balanexp);
        gamm->setValue(pp->locallab.spots.at(index).gamm);

        if (pp->locallab.spots.at(index).exnoiseMethod == "one") {
            exnoiseMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).exnoiseMethod == "med") {
            exnoiseMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).exnoiseMethod == "medhi") {
            exnoiseMethod->set_active(2);
        }

        fatamount->setValue(pp->locallab.spots.at(index).fatamount);
        fatdetail->setValue(pp->locallab.spots.at(index).fatdetail);

        if (complexsoft == 0) {
            fatlevel->setValue(pp->locallab.spots.at(index).fatlevel);
        } else {
            fatlevel->setValue(1.);
        }

        fatanchor->setValue(pp->locallab.spots.at(index).fatanchor);
        sensiex->setValue(pp->locallab.spots.at(index).sensiex);

        if (complexsoft < 2) {
            structexp->setValue(pp->locallab.spots.at(index).structexp);
            blurexpde->setValue(pp->locallab.spots.at(index).blurexpde);
            expcomp->setValue(pp->locallab.spots.at(index).expcomp);
        } else {
            structexp->setValue(0.);
            blurexpde->setValue(5.);
            expcomp->setValue(0.);
        }

        black->setValue(pp->locallab.spots.at(index).black);

        if (complexsoft < 2) {
            hlcompr->setValue(pp->locallab.spots.at(index).hlcompr);
            hlcomprthresh->setValue(pp->locallab.spots.at(index).hlcomprthresh);
            shadex->setValue(pp->locallab.spots.at(index).shadex);
            shcompr->setValue(pp->locallab.spots.at(index).shcompr);
            expchroma->setValue(pp->locallab.spots.at(index).expchroma);
        } else {
            hlcompr->setValue(0.);
            hlcomprthresh->setValue(0.);
            shadex->setValue(0.);
            shcompr->setValue(0.);
            expchroma->setValue(0.);
        }

        shapeexpos->setCurve(pp->locallab.spots.at(index).excurve);
        strexp->setValue(pp->locallab.spots.at(index).strexp);
        angexp->setValue(pp->locallab.spots.at(index).angexp);

        if (complexsoft < 2) {
            softradiusexp->setValue(pp->locallab.spots.at(index).softradiusexp);
            inversex->set_active(pp->locallab.spots.at(index).inversex);
        } else {
            softradiusexp->setValue(0.);
            inversex->set_active(false);
        }

        enaExpMask->set_active(pp->locallab.spots.at(index).enaExpMask);
        enaExpMaskaft->set_active(pp->locallab.spots.at(index).enaExpMaskaft);
        CCmaskexpshape->setCurve(pp->locallab.spots.at(index).CCmaskexpcurve);
        LLmaskexpshape->setCurve(pp->locallab.spots.at(index).LLmaskexpcurve);
        HHmaskexpshape->setCurve(pp->locallab.spots.at(index).HHmaskexpcurve);
        blendmaskexp->setValue(pp->locallab.spots.at(index).blendmaskexp);
        radmaskexp->setValue(pp->locallab.spots.at(index).radmaskexp);

        if (complexsoft == 0) {
            lapmaskexp->setValue(pp->locallab.spots.at(index).lapmaskexp);
        } else {
            lapmaskexp->setValue(0.);
        }

        chromaskexp->setValue(pp->locallab.spots.at(index).chromaskexp);

        if (complexsoft < 2) {
            gammaskexp->setValue(pp->locallab.spots.at(index).gammaskexp);
            slomaskexp->setValue(pp->locallab.spots.at(index).slomaskexp);
            strmaskexp->setValue(pp->locallab.spots.at(index).strmaskexp);
            angmaskexp->setValue(pp->locallab.spots.at(index).angmaskexp);
        } else {
            gammaskexp->setValue(1.);
            slomaskexp->setValue(0.);
            strmaskexp->setValue(0.);
            angmaskexp->setValue(0.);
        }

        Lmaskexpshape->setCurve(pp->locallab.spots.at(index).Lmaskexpcurve);
    }

    // Enable all listeners
    enableListener();

    // Update shcompr sensitivity according to black and shadex value
    updateExposureGUI1();

    // Update exposure GUI according to expMethod value
    updateExposureGUI2();

    // Update exposure GUI according to inversex button state
    updateExposureGUI3();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabExposure::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expexpose = exp->getEnabled();
        pp->locallab.spots.at(index).visiexpose = exp->get_visible();

        if (expMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).expMethod = "std";
        } else if (expMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).expMethod = "pde";
        }

        pp->locallab.spots.at(index).laplacexp = laplacexp->getValue();
        pp->locallab.spots.at(index).linear = linear->getValue();
        pp->locallab.spots.at(index).balanexp = balanexp->getValue();
        pp->locallab.spots.at(index).gamm = gamm->getValue();

        if (exnoiseMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).exnoiseMethod = "none";
        } else if (exnoiseMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).exnoiseMethod = "med";
        } else if (exnoiseMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).exnoiseMethod = "medhi";
        }

        pp->locallab.spots.at(index).fatamount = fatamount->getValue();
        pp->locallab.spots.at(index).fatdetail = fatdetail->getValue();
        pp->locallab.spots.at(index).fatlevel = fatlevel->getValue();
        pp->locallab.spots.at(index).fatanchor = fatanchor->getValue();
        pp->locallab.spots.at(index).sensiex = sensiex->getIntValue();
        pp->locallab.spots.at(index).structexp = structexp->getIntValue();
        pp->locallab.spots.at(index).blurexpde = blurexpde->getIntValue();
        pp->locallab.spots.at(index).expcomp = expcomp->getValue();
        pp->locallab.spots.at(index).black = black->getIntValue();
        pp->locallab.spots.at(index).hlcompr = hlcompr->getIntValue();
        pp->locallab.spots.at(index).hlcomprthresh = hlcomprthresh->getIntValue();
        pp->locallab.spots.at(index).shadex = shadex->getIntValue();
        pp->locallab.spots.at(index).shcompr = shcompr->getIntValue();
        pp->locallab.spots.at(index).expchroma = expchroma->getIntValue();
        pp->locallab.spots.at(index).excurve = shapeexpos->getCurve();
        pp->locallab.spots.at(index).strexp = strexp->getValue();
        pp->locallab.spots.at(index).angexp = angexp->getValue();
        pp->locallab.spots.at(index).softradiusexp = softradiusexp->getValue();
        pp->locallab.spots.at(index).inversex = inversex->get_active();
        pp->locallab.spots.at(index).enaExpMask = enaExpMask->get_active();
        pp->locallab.spots.at(index).enaExpMaskaft = enaExpMaskaft->get_active();
        pp->locallab.spots.at(index).CCmaskexpcurve = CCmaskexpshape->getCurve();
        pp->locallab.spots.at(index).LLmaskexpcurve = LLmaskexpshape->getCurve();
        pp->locallab.spots.at(index).HHmaskexpcurve = HHmaskexpshape->getCurve();
        pp->locallab.spots.at(index).blendmaskexp = blendmaskexp->getIntValue();
        pp->locallab.spots.at(index).radmaskexp = radmaskexp->getValue();
        pp->locallab.spots.at(index).lapmaskexp = lapmaskexp->getValue();
        pp->locallab.spots.at(index).chromaskexp = chromaskexp->getValue();
        pp->locallab.spots.at(index).gammaskexp = gammaskexp->getValue();
        pp->locallab.spots.at(index).slomaskexp = slomaskexp->getValue();
        pp->locallab.spots.at(index).strmaskexp = strmaskexp->getValue();
        pp->locallab.spots.at(index).angmaskexp = angmaskexp->getValue();
        pp->locallab.spots.at(index).Lmaskexpcurve = Lmaskexpshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        laplacexp->setDefault(defSpot.laplacexp);
        linear->setDefault(defSpot.linear);
        balanexp->setDefault(defSpot.balanexp);
        gamm->setDefault(defSpot.gamm);
        fatamount->setDefault(defSpot.fatamount);
        fatdetail->setDefault(defSpot.fatdetail);
        fatlevel->setDefault(defSpot.fatlevel);
        fatanchor->setDefault(defSpot.fatanchor);
        sensiex->setDefault((double)defSpot.sensiex);
        structexp->setDefault((double)defSpot.structexp);
        blurexpde->setDefault((double)defSpot.blurexpde);
        expcomp->setDefault(defSpot.expcomp);
        black->setDefault((double)defSpot.black);
        hlcompr->setDefault((double)defSpot.hlcompr);
        hlcomprthresh->setDefault((double)defSpot.hlcomprthresh);
        shadex->setDefault((double)defSpot.shadex);
        shcompr->setDefault((double)defSpot.shcompr);
        expchroma->setDefault((double)defSpot.expchroma);
        strexp->setDefault(defSpot.strexp);
        angexp->setDefault(defSpot.angexp);
        softradiusexp->setDefault(defSpot.softradiusexp);
        blendmaskexp->setDefault((double)defSpot.blendmaskexp);
        radmaskexp->setDefault(defSpot.radmaskexp);
        lapmaskexp->setDefault(defSpot.lapmaskexp);
        chromaskexp->setDefault(defSpot.chromaskexp);
        gammaskexp->setDefault(defSpot.gammaskexp);
        slomaskexp->setDefault(defSpot.slomaskexp);
        strmaskexp->setDefault(defSpot.strmaskexp);
        angmaskexp->setDefault(defSpot.angmaskexp);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabExposure::adjusterChanged(Adjuster* a, double newval)
{
    // Update shcompr sensitivity according to black and shadex value
    if (a == black || a == shadex) {
        updateExposureGUI1();
    }

    if (isLocActivated && exp->getEnabled()) {
        if (a == laplacexp) {
            if (listener) {
                listener->panelChanged(Evlocallablaplacexp,
                                       laplacexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == linear) {
            if (listener) {
                listener->panelChanged(Evlocallablinear,
                                       linear->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == balanexp) {
            if (listener) {
                listener->panelChanged(Evlocallabbalanexp,
                                       balanexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamm) {
            if (listener) {
                listener->panelChanged(Evlocallabgamm,
                                       gamm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatamount) {
            if (listener) {
                listener->panelChanged(Evlocallabfatamount,
                                       fatamount->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatdetail) {
            if (listener) {
                listener->panelChanged(Evlocallabfatdetail,
                                       fatdetail->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatlevel) {
            if (listener) {
                listener->panelChanged(Evlocallabfatlevel,
                                       fatlevel->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatanchor) {
            if (listener) {
                listener->panelChanged(Evlocallabfatanchor,
                                       fatanchor->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensiex) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiex,
                                       sensiex->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == structexp) {
            if (listener) {
                listener->panelChanged(Evlocallabstructexp,
                                       structexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blurexpde) {
            if (listener) {
                listener->panelChanged(Evlocallabblurexpde,
                                       blurexpde->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == expcomp) {
            if (listener) {
                listener->panelChanged(Evlocallabexpcomp,
                                       expcomp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == black) {
            if (listener) {
                listener->panelChanged(Evlocallabblack,
                                       black->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == hlcompr) {
            if (listener) {
                listener->panelChanged(Evlocallabhlcompr,
                                       hlcompr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == hlcomprthresh) {
            if (listener) {
                listener->panelChanged(Evlocallabhlcomprthresh,
                                       hlcomprthresh->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shadex) {
            if (listener) {
                listener->panelChanged(Evlocallabshadex,
                                       shadex->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shcompr) {
            if (listener) {
                listener->panelChanged(Evlocallabshcompr,
                                       shcompr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == expchroma) {
            if (listener) {
                listener->panelChanged(Evlocallabexpchroma,
                                       expchroma->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strexp) {
            if (listener) {
                listener->panelChanged(Evlocallabstrexp,
                                       strexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == angexp) {
            if (listener) {
                listener->panelChanged(Evlocallabangexp,
                                       angexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == softradiusexp) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusexp,
                                       softradiusexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskexp,
                                       blendmaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskexp,
                                       radmaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskexp,
                                       lapmaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskexp,
                                       chromaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskexp,
                                       gammaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskexp,
                                       slomaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabstrmaskexp,
                                       strmaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == angmaskexp) {
            if (listener) {
                listener->panelChanged(Evlocallabangmaskexp,
                                       angmaskexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabExposure::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == shapeexpos) {
            if (listener) {
                listener->panelChanged(Evlocallabshapeexpos,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskexpshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskexpshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskexpshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmaskexpshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskexpshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabExposure::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenaexpose,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenaexpose,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}


void LocallabExposure::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskexpshape->updateLocallabBackground(normChromar);
        LLmaskexpshape->updateLocallabBackground(normLumar);
        HHmaskexpshape->updateLocallabBackground(normHuer);

        return false;
    }
    );
}

void LocallabExposure::expMethodChanged()
{
    // Update exposure GUI according to expMethod value
    updateExposureGUI2();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabexpMethod,
                                   expMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabExposure::exnoiseMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabexnoiseMethod,
                                   exnoiseMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabExposure::inversexChanged()
{
    // Update exposure GUI according to inversex button state
    updateExposureGUI3();

    // This event is called to transmit potentially resetted mask state
    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (inversex->get_active()) {
                listener->panelChanged(Evlocallabinversex,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabinversex,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabExposure::showmaskexpMethodChanged()
{
    // If mask preview is activated, deactivate other Exposure mask preview
    showmaskexpMethodConninv.block(true);
    showmaskexpMethodinv->set_active(0);
    showmaskexpMethodConninv.block(false);

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabExposure::showmaskexpMethodChangedinv()
{
    // If mask preview is activated, deactivate other Exposure mask preview
    showmaskexpMethodConn.block(true);
    showmaskexpMethod->set_active(0);
    showmaskexpMethodConn.block(false);

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabExposure::enaExpMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaExpMask->get_active()) {
                listener->panelChanged(EvLocallabEnaExpMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaExpMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabExposure::enaExpMaskaftChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaExpMaskaft->get_active()) {
                listener->panelChanged(EvLocallabEnaExpMaskaft,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaExpMaskaft,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabExposure::updateExposureGUI1()
{
    // Update shcompr sensitivity according to black and shadex value
    if (black->getIntValue() == 0 && shadex->getIntValue() == 0) {
        shcompr->set_sensitive(false);
    } else {
        shcompr->set_sensitive(true);
    }
}

void LocallabExposure::updateExposureGUI2()
{
    // Update exposure GUI according to expMethod value
    if (expMethod->get_active_row_number() == 0) {
        pdeFrame->set_sensitive(false);
        laplacexp->set_sensitive(false);
        linear->set_sensitive(false);
        balanexp->set_sensitive(false);
        gamm->set_sensitive(false);
        fatFrame->set_sensitive(false);
        fatamount->set_sensitive(false);
        fatdetail->set_sensitive(false);
        fatlevel->set_sensitive(false);
        fatanchor->set_sensitive(false);
        softradiusexp->set_sensitive(true);
    } else if (expMethod->get_active_row_number() == 1) {
        pdeFrame->set_sensitive(true);
        laplacexp->set_sensitive(true);
        linear->set_sensitive(true);
        balanexp->set_sensitive(true);
        gamm->set_sensitive(true);
        fatFrame->set_sensitive(true);
        fatamount->set_sensitive(true);
        fatdetail->set_sensitive(true);
        fatlevel->set_sensitive(true);
        fatanchor->set_sensitive(true);
        softradiusexp->set_sensitive(false);
    }
}

void LocallabExposure::updateExposureGUI3()
{
    // Update exposure GUI according to inversex button state
    if (inversex->get_active()) {
        expMethod->hide();

        // Manage specific case where expMethod is different from 0
        if (expMethod->get_active_row_number() > 0) {
            expMethodConn.block(true);
            expMethod->set_active(0);
            expMethodConn.block(false);

            // Update GUI accordingly
            updateExposureGUI2();
        }

        pdeFrame->hide();
        fatFrame->hide();
        structexp->hide();
        shadex->hide();
        softradiusexp->hide();
        expgradexp->hide();
        showmaskexpMethod->hide();
        // Reset hidden mask combobox
        showmaskexpMethodConn.block(true);
        showmaskexpMethod->set_active(0);
        showmaskexpMethodConn.block(false);
        showmaskexpMethodinv->show();
    } else {
        expMethod->show();
        pdeFrame->show();
        fatFrame->show();
        structexp->show();
        shadex->show();
        softradiusexp->show();
        expgradexp->show();
        showmaskexpMethodinv->hide();
        // Reset hidden mask combobox
        showmaskexpMethodConninv.block(true);
        showmaskexpMethodinv->set_active(0);
        showmaskexpMethodConninv.block(false);
        showmaskexpMethod->show();
    }
}

/* ==== LocallabShadow ==== */
LocallabShadow::LocallabShadow():
    LocallabTool(this, M("TP_LOCALLAB_SH_TOOLNAME"), M("TP_LOCALLAB_SHADHIGH"), false),

    // Shadow highlight specific widgets
    shMethod(Gtk::manage(new MyComboBoxText())),
    multipliersh([this]() -> std::array<Adjuster *, 5>
    {
    std::array<Adjuster*, 5> res = {};

    for (unsigned int i = 0; i < res.size(); ++i) {
        Glib::ustring ss = Glib::ustring::format(i);

        if (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_LOCALLAB_LUMADARKEST"));
        } else if (i == 4) {
            ss += Glib::ustring::compose(" (%1)", M("TP_LOCALLAB_LUMAWHITESEST"));
        }

        res[i] = Gtk::manage(new Adjuster(std::move(ss), -100, 100, 1, 0));
    }

    return res;
    }
    ()),
    detailSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DETAILSH"), -5, 5, 1, 0))),
    highlights(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0, 100, 1, 0))),
    h_tonalwidth(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HLTONALW"), 10, 100, 1, 70))),
    shadows(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHADOWS"), 0, 100, 1, 0))),
    s_tonalwidth(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHTONALW"), 10, 100, 1, 30))),
    sh_radius(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_RADIUS"), 0, 100, 1, 40))),
    sensihs(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    blurSHde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    gamFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GAMFRA")))),
    gamSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMSH"), 0.25, 15.0, 0.01, 2.4))),
    sloSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOSH"), 0.0, 150.0, 0.01, 12.92))),
    expgradsh(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPGRAD")))),
    strSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -4., 4., 0.05, 0.))),
    angSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    inverssh(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    expmasksh(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWS")))),
    showmaskSHMethod(Gtk::manage(new MyComboBoxText())),
    showmaskSHMethodinv(Gtk::manage(new MyComboBoxText())),
    enaSHMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskSHCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskSHshape(static_cast<FlatCurveEditor*>(maskSHCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskSHshape(static_cast<FlatCurveEditor*>(maskSHCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskSHshape(static_cast<FlatCurveEditor*>(maskSHCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    lapmaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2SHCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    LmaskSHshape(static_cast<DiagonalCurveEditor*>(mask2SHCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    fatSHFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_FATSHFRA")))),
    fatamountSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATAMOUNT"), 1., 100., 1., 1.))),
    fatanchorSH(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATANCHOR"), 1., 100., 1., 50., Gtk::manage(new RTImage("circle-black-small.png")), Gtk::manage(new RTImage("circle-white-small.png")))))
{
    const LocallabParams::LocallabSpot defSpot;

    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Shadow highlight specific widgets
    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_SHADOWHIGHLIGHT_TOOLTIP"));
    }

    shMethod->append(M("TP_LOCALLAB_SH1"));
    shMethod->append(M("TP_LOCALLAB_SH2"));
    shMethod->set_active(0);
    shMethodConn = shMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabShadow::shMethodChanged));

    for (unsigned int i = 0; i < multipliersh.size(); i++) {
        multipliersh[i]->setAdjusterListener(this);
    }

    detailSH->setAdjusterListener(this);

    highlights->setAdjusterListener(this);

    h_tonalwidth->setAdjusterListener(this);

    shadows->setAdjusterListener(this);

    s_tonalwidth->setAdjusterListener(this);

    sh_radius->setAdjusterListener(this);

    sensihs->setAdjusterListener(this);

    blurSHde->setAdjusterListener(this);

    gamFrame->set_label_align(0.025, 0.5);

    gamSH->setAdjusterListener(this);

    sloSH->setAdjusterListener(this);

    setExpandAlignProperties(expgradsh, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    strSH->setAdjusterListener(this);

    if (showtooltip) {
        strSH->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
    }

    angSH->setAdjusterListener(this);

    if (showtooltip) {
        angSH->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    }

    inversshConn = inverssh->signal_toggled().connect(sigc::mem_fun(*this, &LocallabShadow::inversshChanged));

    setExpandAlignProperties(expmasksh, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expmasksh->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskSHMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskSHMethod->set_active(0);

    if (showtooltip) {
        showmaskSHMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskSHMethodConn = showmaskSHMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabShadow::showmaskSHMethodChanged));

    showmaskSHMethodinv->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskSHMethodinv->append(M("TP_LOCALLAB_SHOWMASK"));

    showmaskSHMethodinv->set_active(0);

    if (showtooltip) {
        showmaskSHMethodinv->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskSHMethodConninv = showmaskSHMethodinv->signal_changed().connect(sigc::mem_fun(*this, &LocallabShadow::showmaskSHMethodChangedinv));

    enaSHMaskConn = enaSHMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabShadow::enaSHMaskChanged));

    maskSHCurveEditorG->setCurveListener(this);

    CCmaskSHshape->setIdentityValue(0.);
    CCmaskSHshape->setResetCurve(FlatCurveType(defSpot.CCmaskSHcurve.at(0)), defSpot.CCmaskSHcurve);

    if (showtooltip) {
        CCmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskSHshape->setBottomBarColorProvider(this, 1);

    LLmaskSHshape->setIdentityValue(0.);
    LLmaskSHshape->setResetCurve(FlatCurveType(defSpot.LLmaskSHcurve.at(0)), defSpot.LLmaskSHcurve);

    if (showtooltip) {
        LLmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskSHshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskSHshape->setIdentityValue(0.);
    HHmaskSHshape->setResetCurve(FlatCurveType(defSpot.HHmaskSHcurve.at(0)), defSpot.HHmaskSHcurve);

    if (showtooltip) {
        HHmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskSHshape->setCurveColorProvider(this, 2);
    HHmaskSHshape->setBottomBarColorProvider(this, 2);

    maskSHCurveEditorG->curveListComplete();

    blendmaskSH->setAdjusterListener(this);

    radmaskSH->setLogScale(10, -10);
    radmaskSH->setAdjusterListener(this);

    if (showtooltip) {
        radmaskSH->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    lapmaskSH->setAdjusterListener(this);

    if (showtooltip) {
        lapmaskSH->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    chromaskSH->setAdjusterListener(this);

    gammaskSH->setAdjusterListener(this);

    slomaskSH->setAdjusterListener(this);

    mask2SHCurveEditorG->setCurveListener(this);

    LmaskSHshape->setResetCurve(DiagonalCurveType(defSpot.LmaskSHcurve.at(0)), defSpot.LmaskSHcurve);

    if (showtooltip) {
        LmaskSHshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    LmaskSHshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    LmaskSHshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2SHCurveEditorG->curveListComplete();

    fatSHFrame->set_label_align(0.025, 0.5);

    fatamountSH->setAdjusterListener(this);

    fatanchorSH->setAdjusterListener(this);

    // Add Shadow highlight specific widgets to GUI
    if (complexsoft < 2) {
        pack_start(*shMethod);
    }

    for (int i = 0; i < 5; ++i) {
        pack_start(*multipliersh[i]);
    }

    pack_start(*detailSH);
    pack_start(*highlights);
    pack_start(*h_tonalwidth);
    pack_start(*shadows);
    pack_start(*s_tonalwidth);
    pack_start(*sh_radius);
    pack_start(*sensihs);

    if (complexsoft < 2) {
        pack_start(*blurSHde);
    }

    ToolParamBlock* const gammBox = Gtk::manage(new ToolParamBlock());
    gammBox->pack_start(*gamSH);
    gammBox->pack_start(*sloSH);
    gamFrame->add(*gammBox);
    pack_start(*gamFrame);
    ToolParamBlock* const gradSHBox = Gtk::manage(new ToolParamBlock());
    gradSHBox->pack_start(*strSH);
    gradSHBox->pack_start(*angSH);
    expgradsh->add(*gradSHBox, false);
    pack_start(*expgradsh);
    pack_start(*inverssh);
    ToolParamBlock* const maskSHBox = Gtk::manage(new ToolParamBlock());
    maskSHBox->pack_start(*showmaskSHMethod, Gtk::PACK_SHRINK, 4);
    maskSHBox->pack_start(*showmaskSHMethodinv, Gtk::PACK_SHRINK, 4);
    maskSHBox->pack_start(*enaSHMask, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*maskSHCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskSHBox->pack_start(*blendmaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*radmaskSH, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 1) {
        maskSHBox->pack_start(*lapmaskSH, Gtk::PACK_SHRINK, 0);
    }

    maskSHBox->pack_start(*chromaskSH, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        maskSHBox->pack_start(*gammaskSH, Gtk::PACK_SHRINK, 0);
        maskSHBox->pack_start(*slomaskSH, Gtk::PACK_SHRINK, 0);
    }

    maskSHBox->pack_start(*mask2SHCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    ToolParamBlock* const fatSHBox = Gtk::manage(new ToolParamBlock());
    fatSHBox->pack_start(*fatamountSH);
    fatSHBox->pack_start(*fatanchorSH);
    fatSHFrame->add(*fatSHBox);

    if (complexsoft < 1) {
        maskSHBox->pack_start(*fatSHFrame);
    }

    expmasksh->add(*maskSHBox, false);
    pack_start(*expmasksh, false, false);
}

LocallabShadow::~LocallabShadow()
{
    delete maskSHCurveEditorG;
    delete mask2SHCurveEditorG;
}

void LocallabShadow::resetMaskView()
{
    showmaskSHMethodConn.block(true);
    showmaskSHMethodConninv.block(true);

    showmaskSHMethod->set_active(0);
    showmaskSHMethodinv->set_active(0);

    showmaskSHMethodConn.block(false);
    showmaskSHMethodConninv.block(false);
}

void LocallabShadow::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    shMask = showmaskSHMethod->get_active_row_number();
    shMaskinv = showmaskSHMethodinv->get_active_row_number();
}

void LocallabShadow::disableListener()
{
    LocallabTool::disableListener();

    shMethodConn.block(true);
    inversshConn.block(true);
    showmaskSHMethodConn.block(true);
    showmaskSHMethodConninv.block(true);
    enaSHMaskConn.block(true);
}

void LocallabShadow::enableListener()
{
    LocallabTool::enableListener();

    shMethodConn.block(false);
    inversshConn.block(false);
    showmaskSHMethodConn.block(false);
    showmaskSHMethodConninv.block(false);
    enaSHMaskConn.block(false);
}

void LocallabShadow::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visishadhigh);
        exp->setEnabled(pp->locallab.spots.at(index).expshadhigh);

        if (complexsoft < 2) {
            if (pp->locallab.spots.at(index).shMethod == "std") {
                shMethod->set_active(0);
            } else if (pp->locallab.spots.at(index).shMethod == "tone") {
                shMethod->set_active(1);
            }
        } else {
            shMethod->set_active(1);
        }

        for (int i = 0; i < 5; i++) {
            multipliersh[i]->setValue((double)pp->locallab.spots.at(index).multsh[i]);
        }

        detailSH->setValue((double)pp->locallab.spots.at(index).detailSH);

        if (complexsoft < 2) {
            highlights->setValue((double)pp->locallab.spots.at(index).highlights);
        } else {
            highlights->setValue(0.);
        }

        h_tonalwidth->setValue((double)pp->locallab.spots.at(index).h_tonalwidth);

        if (complexsoft < 2) {
            shadows->setValue(pp->locallab.spots.at(index).shadows);
        } else {
            shadows->setValue(0.);
        }

        s_tonalwidth->setValue((double)pp->locallab.spots.at(index).s_tonalwidth);
        sh_radius->setValue((double)pp->locallab.spots.at(index).sh_radius);
        sensihs->setValue((double)pp->locallab.spots.at(index).sensihs);
        blurSHde->setValue((double)pp->locallab.spots.at(index).blurSHde);
        gamSH->setValue(pp->locallab.spots.at(index).gamSH);
        sloSH->setValue(pp->locallab.spots.at(index).sloSH);
        strSH->setValue(pp->locallab.spots.at(index).strSH);
        angSH->setValue(pp->locallab.spots.at(index).angSH);
        inverssh->set_active(pp->locallab.spots.at(index).inverssh);
        enaSHMask->set_active(pp->locallab.spots.at(index).enaSHMask);
        CCmaskSHshape->setCurve(pp->locallab.spots.at(index).CCmaskSHcurve);
        LLmaskSHshape->setCurve(pp->locallab.spots.at(index).LLmaskSHcurve);
        HHmaskSHshape->setCurve(pp->locallab.spots.at(index).HHmaskSHcurve);
        blendmaskSH->setValue((double)pp->locallab.spots.at(index).blendmaskSH);
        radmaskSH->setValue(pp->locallab.spots.at(index).radmaskSH);

        if (complexsoft == 0) {
            lapmaskSH->setValue(pp->locallab.spots.at(index).lapmaskSH);
        } else {
            lapmaskSH->setValue(0.);
        }

        chromaskSH->setValue(pp->locallab.spots.at(index).chromaskSH);

        if (complexsoft < 2) {
            gammaskSH->setValue(pp->locallab.spots.at(index).gammaskSH);
            slomaskSH->setValue(pp->locallab.spots.at(index).slomaskSH);
        } else {
            gammaskSH->setValue(1.);
            slomaskSH->setValue(0.);
        }

        LmaskSHshape->setCurve(pp->locallab.spots.at(index).LmaskSHcurve);
        fatamountSH->setValue(pp->locallab.spots.at(index).fatamountSH);
        fatanchorSH->setValue(pp->locallab.spots.at(index).fatanchorSH);
    }

    // Enable all listeners
    enableListener();

    // Update shadow highlight GUI according to inverssh button state
    updateShadowGUI1();

    // Update shadow highlight GUI according to shMethod combobox state
    updateShadowGUI2();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expshadhigh = exp->getEnabled();
        pp->locallab.spots.at(index).visishadhigh = exp->get_visible();

        if (shMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).shMethod = "std";
        } else if (shMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).shMethod = "tone";
        }

        for (int i = 0; i < 5; i++) {
            pp->locallab.spots.at(index).multsh[i] = multipliersh[i]->getIntValue();
        }

        pp->locallab.spots.at(index).detailSH = detailSH->getIntValue();
        pp->locallab.spots.at(index).highlights = highlights->getIntValue();
        pp->locallab.spots.at(index).h_tonalwidth = h_tonalwidth->getIntValue();
        pp->locallab.spots.at(index).shadows = shadows->getIntValue();
        pp->locallab.spots.at(index).s_tonalwidth = s_tonalwidth->getIntValue();
        pp->locallab.spots.at(index).sh_radius = sh_radius->getIntValue();
        pp->locallab.spots.at(index).sensihs = sensihs->getIntValue();
        pp->locallab.spots.at(index).blurSHde = blurSHde->getIntValue();
        pp->locallab.spots.at(index).gamSH = gamSH->getValue();
        pp->locallab.spots.at(index).sloSH = sloSH->getValue();
        pp->locallab.spots.at(index).strSH = strSH->getValue();
        pp->locallab.spots.at(index).angSH = angSH->getValue();
        pp->locallab.spots.at(index).inverssh = inverssh->get_active();
        pp->locallab.spots.at(index).enaSHMask = enaSHMask->get_active();
        pp->locallab.spots.at(index).LLmaskSHcurve = LLmaskSHshape->getCurve();
        pp->locallab.spots.at(index).CCmaskSHcurve = CCmaskSHshape->getCurve();
        pp->locallab.spots.at(index).HHmaskSHcurve = HHmaskSHshape->getCurve();
        pp->locallab.spots.at(index).blendmaskSH = blendmaskSH->getIntValue();
        pp->locallab.spots.at(index).radmaskSH = radmaskSH->getValue();
        pp->locallab.spots.at(index).lapmaskSH = lapmaskSH->getValue();
        pp->locallab.spots.at(index).chromaskSH = chromaskSH->getValue();
        pp->locallab.spots.at(index).gammaskSH = gammaskSH->getValue();
        pp->locallab.spots.at(index).slomaskSH = slomaskSH->getValue();
        pp->locallab.spots.at(index).LmaskSHcurve = LmaskSHshape->getCurve();
        pp->locallab.spots.at(index).fatamountSH = fatamountSH->getValue();
        pp->locallab.spots.at(index).fatanchorSH = fatanchorSH->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        for (int i = 0; i < 5; i++) {
            multipliersh[i]->setDefault(defSpot.multsh[i]);
        }

        detailSH->setDefault((double)defSpot.detailSH);
        highlights->setDefault((double)defSpot.highlights);
        h_tonalwidth->setDefault((double)defSpot.h_tonalwidth);
        shadows->setDefault((double)defSpot.shadows);
        s_tonalwidth->setDefault((double)defSpot.s_tonalwidth);
        sh_radius->setDefault((double)defSpot.sh_radius);
        sensihs->setDefault((double)defSpot.sensihs);
        blurSHde->setDefault((double)defSpot.blurSHde);
        gamSH->setDefault(defSpot.gamSH);
        sloSH->setDefault(defSpot.sloSH);
        strSH->setDefault(defSpot.strSH);
        angSH->setDefault(defSpot.angSH);
        blendmaskSH->setDefault((double)defSpot.blendmaskSH);
        radmaskSH->setDefault(defSpot.radmaskSH);
        lapmaskSH->setDefault(defSpot.lapmaskSH);
        chromaskSH->setDefault(defSpot.chromaskSH);
        gammaskSH->setDefault(defSpot.gammaskSH);
        slomaskSH->setDefault(defSpot.slomaskSH);
        fatamountSH->setDefault(defSpot.fatamountSH);
        fatanchorSH->setDefault(defSpot.fatanchorSH);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == multipliersh[0] || a == multipliersh[1] || a == multipliersh[2] || a == multipliersh[3] || a == multipliersh[4]) {
            if (listener) {
                listener->panelChanged(EvlocallabEqualizersh,
                                       Glib::ustring::compose("%1, %2, %3, %4, %5",
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multipliersh[0]->getIntValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multipliersh[1]->getIntValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multipliersh[2]->getIntValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multipliersh[3]->getIntValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(2), multipliersh[4]->getIntValue())) + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == detailSH) {
            if (listener) {
                listener->panelChanged(EvlocallabdetailSH,
                                       detailSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == highlights) {
            if (listener) {
                listener->panelChanged(Evlocallabhighlights,
                                       highlights->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == h_tonalwidth) {
            if (listener) {
                listener->panelChanged(Evlocallabh_tonalwidth,
                                       h_tonalwidth->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shadows) {
            if (listener) {
                listener->panelChanged(Evlocallabshadows,
                                       shadows->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == s_tonalwidth) {
            if (listener) {
                listener->panelChanged(Evlocallabs_tonalwidth,
                                       s_tonalwidth->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sh_radius) {
            if (listener) {
                listener->panelChanged(Evlocallabsh_radius,
                                       sh_radius->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensihs) {
            if (listener) {
                listener->panelChanged(Evlocallabsensihs,
                                       sensihs->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blurSHde) {
            if (listener) {
                listener->panelChanged(EvlocallabblurSHde,
                                       blurSHde->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamSH) {
            if (listener) {
                listener->panelChanged(EvlocallabgamSH,
                                       gamSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloSH) {
            if (listener) {
                listener->panelChanged(EvlocallabsloSH,
                                       sloSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strSH) {
            if (listener) {
                listener->panelChanged(EvlocallabstrSH,
                                       strSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == angSH) {
            if (listener) {
                listener->panelChanged(EvlocallabangSH,
                                       angSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabblendmaskSH,
                                       blendmaskSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabradmaskSH,
                                       radmaskSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallablapmaskSH,
                                       lapmaskSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabchromaskSH,
                                       chromaskSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabgammaskSH,
                                       gammaskSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskSH) {
            if (listener) {
                listener->panelChanged(EvlocallabslomaskSH,
                                       slomaskSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == fatamountSH) {
            if (listener) {
                listener->panelChanged(EvlocallabfatamountSH,
                                       fatamountSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }


        if (a == fatanchorSH) {
            if (listener) {
                listener->panelChanged(EvlocallabfatanchorSH,
                                       fatanchorSH->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskSHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskSHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskSHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LmaskSHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskSHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenashadhigh,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenashadhigh,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskSHshape->updateLocallabBackground(normChromar);
        LLmaskSHshape->updateLocallabBackground(normLumar);
        HHmaskSHshape->updateLocallabBackground(normHuer);

        return false;
    }
    );
}

void LocallabShadow::shMethodChanged()
{

    // Update shadow highlight GUI according to shMethod combobox state
    updateShadowGUI2();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabshMethod,
                                   shMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabShadow::inversshChanged()
{
    // Update shadow highlight GUI according to inverssh button state
    updateShadowGUI1();

    // This event is called to transmit potentially resetted mask state
    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (inverssh->get_active()) {
                listener->panelChanged(Evlocallabinverssh,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabinverssh,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::showmaskSHMethodChanged()
{
    // If mask preview is activated, deactivate other Shadow highlight mask preview
    showmaskSHMethodConninv.block(true);
    showmaskSHMethodinv->set_active(0);
    showmaskSHMethodConninv.block(false);

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabShadow::showmaskSHMethodChangedinv()
{
    // If mask preview is activated, deactivate other Shadow highlight mask preview
    showmaskSHMethodConn.block(true);
    showmaskSHMethod->set_active(0);
    showmaskSHMethodConn.block(false);

    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabShadow::enaSHMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaSHMask->get_active()) {
                listener->panelChanged(EvLocallabEnaSHMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaSHMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::updateShadowGUI1()
{
    // Update shadow highlight GUI according to inverssh button state
    if (inverssh->get_active()) {
        expgradsh->hide();
        showmaskSHMethod->hide();
        // Reset hidden mask combobox
        showmaskSHMethodConn.block(true);
        showmaskSHMethod->set_active(0);
        showmaskSHMethodConn.block(false);
        showmaskSHMethodinv->show();
    } else {
        expgradsh->show();
        showmaskSHMethod->show();
        showmaskSHMethodinv->hide();
        // Reset hidden mask combobox
        showmaskSHMethodConninv.block(true);
        showmaskSHMethodinv->set_active(0);
        showmaskSHMethodConninv.block(false);
    }
}

void LocallabShadow::updateShadowGUI2()
{
    // Update shadow highlight GUI according to shMethod combobox state
    if (shMethod->get_active_row_number() == 0) {
        for (int i = 0; i < 5; i++) {
            multipliersh[i]->hide();
        }

        detailSH->hide();
        highlights->show();
        h_tonalwidth->show();
        shadows->show();
        s_tonalwidth->show();
        sh_radius->show();
    } else if (shMethod->get_active_row_number() == 1) {
        for (int i = 0; i < 5; i++) {
            multipliersh[i]->show();
        }

        detailSH->show();
        highlights->hide();
        h_tonalwidth->hide();
        shadows->hide();
        s_tonalwidth->hide();
        sh_radius->hide();
    }
}

/* ==== LocallabVibrance ==== */
LocallabVibrance::LocallabVibrance():
    LocallabTool(this, M("TP_LOCALLAB_VIB_TOOLNAME"), M("TP_LOCALLAB_VIBRANCE"), false),

    // Vibrance specific widgets
    saturated(Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.))),
    pastels(Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.))),
    warm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    psThreshold(Gtk::manage(new ThresholdAdjuster(M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false))),
    protectSkins(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PROTECTSKINS")))),
    avoidColorShift(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_AVOIDCOLORSHIFT")))),
    pastSatTog(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PASTSATTOG")))),
    sensiv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    curveEditorGG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL"))),
    skinTonesCurve(static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")))),
    expgradvib(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPGRAD")))),
    strvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -4., 4., 0.05, 0.))),
    strvibab(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRCHRO"), -4., 4., 0.05, 0.))),
    strvibh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRHUE2"), -6., 6., 0.05, 0.))),
    angvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    expmaskvib(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWVI")))),
    showmaskvibMethod(Gtk::manage(new MyComboBoxText())),
    enavibMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskvibCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskvibshape(static_cast<FlatCurveEditor*>(maskvibCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskvibshape(static_cast<FlatCurveEditor*>(maskvibCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskvibshape(static_cast<FlatCurveEditor *>(maskvibCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    blendmaskvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    lapmaskvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.))),
    slomaskvib(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    mask2vibCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK2"))),
    Lmaskvibshape(static_cast<DiagonalCurveEditor*>(mask2vibCurveEditorG->addCurve(CT_Diagonal, "L(L)")))
{
    float R, G, B;

    const LocallabParams::LocallabSpot defSpot;

    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Vibrance specific widgets
    saturated->setAdjusterListener(this);

    if (complexsoft == 2) {
        pastels->setLabel(M("TP_LOCALLAB_PASTELS2"));
    }

    pastels->setAdjusterListener(this);

    if (showtooltip) {
        warm->set_tooltip_text(M("TP_LOCALLAB_WARM_TOOLTIP"));
    }

    warm->setAdjusterListener(this);

    if (showtooltip) {
        psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    }

    psThreshold->setAdjusterListener(this);

    pskinsConn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::protectskins_toggled));

    ashiftConn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::avoidcolorshift_toggled));

    pastsattogConn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::pastsattog_toggled));

    sensiv->setAdjusterListener(this);

    curveEditorGG->setCurveListener(this);

    if (showtooltip) {
        skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
    }

    std::vector<GradientMilestone> mskinTonesCurve;
    // -0.1 rad < Hue < 1.6 rad
    Color::hsv2rgb01(0.92f, 0.45f, 0.6f, R, G, B);
    mskinTonesCurve.emplace_back(0.0, R, G, B);
    Color::hsv2rgb01(0.14056f, 0.45f, 0.6f, R, G, B);
    mskinTonesCurve.emplace_back(0.0, R, G, B);
    skinTonesCurve->setBottomBarBgGradient(mskinTonesCurve);
    skinTonesCurve->setLeftBarBgGradient(mskinTonesCurve);
    skinTonesCurve->setRangeLabels(
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE1"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE2"),
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE3"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE4")
    );
    skinTonesCurve->setRangeDefaultMilestones(0.1, 0.4, 0.85);

    curveEditorGG->curveListComplete();

    setExpandAlignProperties(expgradvib, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        strvib->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
    }

    strvib->setAdjusterListener(this);

    if (showtooltip) {
        strvibab->set_tooltip_text(M("TP_LOCALLAB_GRADSTRAB_TOOLTIP"));
    }

    strvibab->setAdjusterListener(this);

    if (showtooltip) {
        strvibh->set_tooltip_text(M("TP_LOCALLAB_GRADSTRHUE_TOOLTIP"));
    }

    strvibh->setAdjusterListener(this);

    if (showtooltip) {
        angvib->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    }

    angvib->setAdjusterListener(this);

    setExpandAlignProperties(expmaskvib, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expmaskvib->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskvibMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskvibMethod->set_active(0);

    if (showtooltip) {
        showmaskvibMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskvibMethodConn = showmaskvibMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabVibrance::showmaskvibMethodChanged));

    enavibMaskConn = enavibMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::enavibMaskChanged));

    maskvibCurveEditorG->setCurveListener(this);

    CCmaskvibshape->setIdentityValue(0.);
    CCmaskvibshape->setResetCurve(FlatCurveType(defSpot.CCmaskvibcurve.at(0)), defSpot.CCmaskvibcurve);

    if (showtooltip) {
        CCmaskvibshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskvibshape->setBottomBarColorProvider(this, 1);

    LLmaskvibshape->setIdentityValue(0.);
    LLmaskvibshape->setResetCurve(FlatCurveType(defSpot.LLmaskvibcurve.at(0)), defSpot.LLmaskvibcurve);

    if (showtooltip) {
        LLmaskvibshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskvibshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskvibshape->setIdentityValue(0.);
    HHmaskvibshape->setResetCurve(FlatCurveType(defSpot.HHmaskvibcurve.at(0)), defSpot.HHmaskvibcurve);

    if (showtooltip) {
        HHmaskvibshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskvibshape->setCurveColorProvider(this, 2);
    HHmaskvibshape->setBottomBarColorProvider(this, 2);

    maskvibCurveEditorG->curveListComplete();

    blendmaskvib->setAdjusterListener(this);

    radmaskvib->setLogScale(10, -10);
    radmaskvib->setAdjusterListener(this);

    lapmaskvib->setAdjusterListener(this);

    chromaskvib->setAdjusterListener(this);

    gammaskvib->setAdjusterListener(this);

    slomaskvib->setAdjusterListener(this);

    mask2vibCurveEditorG->setCurveListener(this);

    Lmaskvibshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskvibcurve.at(0)), defSpot.Lmaskvibcurve);

    if (showtooltip) {
        Lmaskvibshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmaskvibshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskvibshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2vibCurveEditorG->curveListComplete();

    // Add Vibrance specific widgets to GUI
    if (complexsoft < 2) {
        pack_start(*saturated, Gtk::PACK_SHRINK, 0);
    }

    pack_start(*pastels, Gtk::PACK_SHRINK, 0);
    pack_start(*warm, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        pack_start(*psThreshold, Gtk::PACK_SHRINK, 0);
        pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);
        pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);
        pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);
    }

    pack_start(*sensiv, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        pack_start(*curveEditorGG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    }

    ToolParamBlock* const gradvibBox = Gtk::manage(new ToolParamBlock());
    gradvibBox->pack_start(*strvib);

    if (complexsoft < 2) {
        gradvibBox->pack_start(*strvibab);
        gradvibBox->pack_start(*strvibh);
    }

    gradvibBox->pack_start(*angvib);
    expgradvib->add(*gradvibBox, false);
    pack_start(*expgradvib);
    ToolParamBlock* const maskvibBox = Gtk::manage(new ToolParamBlock());
    maskvibBox->pack_start(*showmaskvibMethod, Gtk::PACK_SHRINK, 4);
    maskvibBox->pack_start(*enavibMask, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*maskvibCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskvibBox->pack_start(*blendmaskvib, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*radmaskvib, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 1) {
        maskvibBox->pack_start(*lapmaskvib, Gtk::PACK_SHRINK, 0);
    }

    maskvibBox->pack_start(*chromaskvib, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        maskvibBox->pack_start(*gammaskvib, Gtk::PACK_SHRINK, 0);
        maskvibBox->pack_start(*slomaskvib, Gtk::PACK_SHRINK, 0);
    }

    maskvibBox->pack_start(*mask2vibCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    expmaskvib->add(*maskvibBox, false);
    pack_start(*expmaskvib, false, false);
}

LocallabVibrance::~LocallabVibrance()
{
    delete curveEditorGG;
    delete maskvibCurveEditorG;
    delete mask2vibCurveEditorG;
}

void LocallabVibrance::resetMaskView()
{
    showmaskvibMethodConn.block(true);

    showmaskvibMethod->set_active(0);

    showmaskvibMethodConn.block(false);
}

void LocallabVibrance::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    vibMask = showmaskvibMethod->get_active_row_number();
}

void LocallabVibrance::disableListener()
{
    LocallabTool::disableListener();

    pskinsConn.block(true);
    ashiftConn.block(true);
    pastsattogConn.block(true);
    showmaskvibMethodConn.block(true);
    enavibMaskConn.block(true);
}

void LocallabVibrance::enableListener()
{
    LocallabTool::enableListener();

    pskinsConn.block(false);
    ashiftConn.block(false);
    pastsattogConn.block(false);
    showmaskvibMethodConn.block(false);
    enavibMaskConn.block(false);
}

void LocallabVibrance::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visivibrance);
        exp->setEnabled(pp->locallab.spots.at(index).expvibrance);

        saturated->setValue(pp->locallab.spots.at(index).saturated);
        pastels->setValue(pp->locallab.spots.at(index).pastels);
        warm->setValue(pp->locallab.spots.at(index).warm);
        psThreshold->setValue<int>(pp->locallab.spots.at(index).psthreshold);
        protectSkins->set_active(pp->locallab.spots.at(index).protectskins);
        avoidColorShift->set_active(pp->locallab.spots.at(index).avoidcolorshift);
        pastSatTog->set_active(pp->locallab.spots.at(index).pastsattog);
        sensiv->setValue(pp->locallab.spots.at(index).sensiv);

        if (complexsoft < 2) {
            skinTonesCurve->setCurve(pp->locallab.spots.at(index).skintonescurve);
        } else {
            skinTonesCurve->reset();
        }

        strvib->setValue(pp->locallab.spots.at(index).strvib);

        if (complexsoft == 0) {
            strvibab->setValue(pp->locallab.spots.at(index).strvibab);
            strvibh->setValue(pp->locallab.spots.at(index).strvibh);
        } else {
            strvibab->setValue(0.);
            strvibh->setValue(0.);
        }

        angvib->setValue(pp->locallab.spots.at(index).angvib);
        enavibMask->set_active(pp->locallab.spots.at(index).enavibMask);
        CCmaskvibshape->setCurve(pp->locallab.spots.at(index).CCmaskvibcurve);
        LLmaskvibshape->setCurve(pp->locallab.spots.at(index).LLmaskvibcurve);
        HHmaskvibshape->setCurve(pp->locallab.spots.at(index).HHmaskvibcurve);
        blendmaskvib->setValue(pp->locallab.spots.at(index).blendmaskvib);
        radmaskvib->setValue(pp->locallab.spots.at(index).radmaskvib);

        if (complexsoft == 0) {
            lapmaskvib->setValue(pp->locallab.spots.at(index).lapmaskvib);
        } else {
            lapmaskvib->setValue(0.);
        }

        chromaskvib->setValue(pp->locallab.spots.at(index).chromaskvib);

        if (complexsoft < 2) {
            gammaskvib->setValue(pp->locallab.spots.at(index).gammaskvib);
            slomaskvib->setValue(pp->locallab.spots.at(index).slomaskvib);
        } else {
            gammaskvib->setValue(1.);
            slomaskvib->setValue(0.);
        }

        Lmaskvibshape->setCurve(pp->locallab.spots.at(index).Lmaskvibcurve);
    }

    // Enable all listeners
    enableListener();

    // Update vibrance GUI according to pastsattog button state
    updateVibranceGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabVibrance::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expvibrance = exp->getEnabled();
        pp->locallab.spots.at(index).visivibrance = exp->get_visible();

        pp->locallab.spots.at(index).saturated = saturated->getIntValue();
        pp->locallab.spots.at(index).pastels = pastels->getIntValue();
        pp->locallab.spots.at(index).warm = warm->getIntValue();
        pp->locallab.spots.at(index).psthreshold = psThreshold->getValue<int>();
        pp->locallab.spots.at(index).protectskins = protectSkins->get_active();
        pp->locallab.spots.at(index).avoidcolorshift = avoidColorShift->get_active();
        pp->locallab.spots.at(index).pastsattog = pastSatTog->get_active();
        pp->locallab.spots.at(index).sensiv = sensiv->getIntValue();
        pp->locallab.spots.at(index).skintonescurve = skinTonesCurve->getCurve();
        pp->locallab.spots.at(index).strvib = strvib->getValue();
        pp->locallab.spots.at(index).strvibab = strvibab->getValue();
        pp->locallab.spots.at(index).strvibh = strvibh->getValue();
        pp->locallab.spots.at(index).angvib = angvib->getValue();
        pp->locallab.spots.at(index).enavibMask = enavibMask->get_active();
        pp->locallab.spots.at(index).CCmaskvibcurve = CCmaskvibshape->getCurve();
        pp->locallab.spots.at(index).LLmaskvibcurve = LLmaskvibshape->getCurve();
        pp->locallab.spots.at(index).HHmaskvibcurve = HHmaskvibshape->getCurve();
        pp->locallab.spots.at(index).blendmaskvib = blendmaskvib->getIntValue();
        pp->locallab.spots.at(index).radmaskvib = radmaskvib->getValue();
        pp->locallab.spots.at(index).lapmaskvib = lapmaskvib->getValue();
        pp->locallab.spots.at(index).chromaskvib = chromaskvib->getValue();
        pp->locallab.spots.at(index).gammaskvib = gammaskvib->getValue();
        pp->locallab.spots.at(index).slomaskvib = slomaskvib->getValue();
        pp->locallab.spots.at(index).Lmaskvibcurve = Lmaskvibshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabVibrance::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster and threshold adjuster widgets
        saturated->setDefault((double)defSpot.saturated);
        pastels->setDefault((double)defSpot.pastels);
        warm->setDefault((double)defSpot.warm);
        psThreshold->setDefault<int>(defSpot.psthreshold);
        sensiv->setDefault((double)defSpot.sensiv);
        strvib->setDefault(defSpot.strvib);
        strvibab->setDefault(defSpot.strvibab);
        strvibh->setDefault(defSpot.strvibh);
        angvib->setDefault(defSpot.angvib);
        blendmaskvib->setDefault((double)defSpot.blendmaskvib);
        radmaskvib->setDefault(defSpot.radmaskvib);
        lapmaskvib->setDefault(defSpot.lapmaskvib);
        chromaskvib->setDefault(defSpot.chromaskvib);
        gammaskvib->setDefault(defSpot.gammaskvib);
        slomaskvib->setDefault(defSpot.slomaskvib);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabVibrance::adjusterChanged(Adjuster* a, double newval)
{
    // Copy pastels adjuster value to saturated one according to pastSatTog button state
    if (a == pastels && pastSatTog->get_active()) {
        saturated->setValue(newval);
    }

    if (isLocActivated && exp->getEnabled()) {
        if (a == saturated && !pastSatTog->get_active()) {
            if (listener) {
                listener->panelChanged(EvlocallabSaturated,
                                       saturated->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == pastels) {
            if (listener) {
                listener->panelChanged(EvlocallabPastels,
                                       pastels->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == warm) {
            if (listener) {
                listener->panelChanged(Evlocallabwarm,
                                       warm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensiv) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiv,
                                       sensiv->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strvib) {
            if (listener) {
                listener->panelChanged(Evlocallabstrvib,
                                       strvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strvibab) {
            if (listener) {
                listener->panelChanged(Evlocallabstrvibab,
                                       strvibab->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strvibh) {
            if (listener) {
                listener->panelChanged(Evlocallabstrvibh,
                                       strvibh->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == angvib) {
            if (listener) {
                listener->panelChanged(Evlocallabangvib,
                                       angvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskvib) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskvi,
                                       blendmaskvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskvib) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskvib,
                                       radmaskvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskvib) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskvib,
                                       lapmaskvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskvib) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskvib,
                                       chromaskvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskvib) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskvib,
                                       gammaskvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskvib) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskvib,
                                       slomaskvib->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabPastSatThreshold,
                                   psThreshold->getHistoryString() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

std::vector<double> LocallabVibrance::getCurvePoints(ThresholdSelector* tAdjuster) const
{
    std::vector<double> points;
    double threshold, transitionWeighting;
    tAdjuster->getPositions<double>(transitionWeighting, threshold); // ( range -100;+100, range 0;+100 )
    transitionWeighting /= 100.; // range -1., +1.
    threshold /= 100.; // range  0., +1.

    // Initial point
    points.push_back(0.);
    points.push_back(0.);

    double p2 = 3.0 * threshold / 4.0; // same one than in ipvibrance.cc
    double s0 = threshold + (1.0 - threshold) / 4.0; // same one than in ipvibrance.cc

    // point at the beginning of the first linear transition
    points.push_back(p2);
    points.push_back(0.);

    // Y value of the chroma mean point, calculated to get a straight line between p2 and s0
    double chromaMean = (threshold / 4.0) / (s0 - p2);

    // move chromaMean up or down depending on transitionWeighting
    if (transitionWeighting > 0.0) {
        // positive values -> give more weight to Saturated
        chromaMean = (1.0 - chromaMean) * transitionWeighting + chromaMean;
    } else if (transitionWeighting < 0.0) {
        // negative values -> give more weight to Pastels
        chromaMean = chromaMean * transitionWeighting + chromaMean;
    }

    // point at the location of the Top cursor, at the end of the first linear transition and the beginning of the second one
    points.push_back(threshold);
    points.push_back(chromaMean);

    if (threshold < 1.0) {
        // point at the end of the second linear transition
        points.push_back(s0);
        points.push_back(1.0);

        // end point
        points.push_back(1.0);
        points.push_back(1.0);
    }

    return points;
}

void LocallabVibrance::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == skinTonesCurve) {
            if (listener) {
                listener->panelChanged(EvlocallabSkinTonesCurve,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCmaskvibshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskvibshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskvibshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskvibshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskvibshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskvibshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmaskvibshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskvibshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenavibrance,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenashadhigh,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskvibshape->updateLocallabBackground(normChromar);
        LLmaskvibshape->updateLocallabBackground(normLumar);
        HHmaskvibshape->updateLocallabBackground(normHuer);

        return false;
    }
    );
}

void LocallabVibrance::protectskins_toggled()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (protectSkins->get_active()) {
                listener->panelChanged(EvlocallabProtectSkins,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvlocallabProtectSkins,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::avoidcolorshift_toggled()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (avoidColorShift->get_active()) {
                listener->panelChanged(EvlocallabAvoidColorShift,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvlocallabAvoidColorShift,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::pastsattog_toggled()
{
    // Update vibrance GUI according to pastsattog button state
    updateVibranceGUI();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (pastSatTog->get_active()) {
                listener->panelChanged(EvlocallabPastSatTog,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvlocallabPastSatTog,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::showmaskvibMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabVibrance::enavibMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enavibMask->get_active()) {
                listener->panelChanged(EvLocallabEnavibMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnavibMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabVibrance::updateVibranceGUI()
{
    // Update vibrance GUI according to pastsattog button state
    if (pastSatTog->get_active()) {
        // Link both slider, so we set saturated and psThresholds unsensitive
        psThreshold->set_sensitive(false);
        saturated->set_sensitive(false);
        saturated->setValue(pastels->getValue()); // Pastels and Saturated are linked
    } else {
        // Separate sliders, so we set saturated and psThresholds sensitive again
        psThreshold->set_sensitive(true);
        saturated->set_sensitive(true);
    }
}

/* ==== LocallabSoft ==== */
LocallabSoft::LocallabSoft():
    LocallabTool(this, M("TP_LOCALLAB_SOFT_TOOLNAME"), M("TP_LOCALLAB_SOFT"), false),

    // Soft light specific widgets
    softMethod(Gtk::manage(new MyComboBoxText())),
    ctboxsoftmethod(Gtk::manage(new Gtk::HBox())),
    showmasksoftMethod(Gtk::manage(new MyComboBoxText())),
    streng(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENG"), 1, 100, 1, 1))),
    laplace(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPLACE"), 0., 100., 0.5, 25.))),
    sensisf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 1, 100, 1, 15)))
{
    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Soft light specific widgets
    softMethod->append(M("TP_LOCALLAB_SOFTM"));
    softMethod->append(M("TP_LOCALLAB_RETIM"));
    softMethod->set_active(0);

    if (showtooltip) {
        softMethod->set_tooltip_markup(M("TP_LOCALLAB_SOFTMETHOD_TOOLTIP"));
    }

    softMethodConn = softMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSoft::softMethodChanged));

    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWLAPLACE"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWFOURIER"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWPOISSON"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWNORMAL"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasksoftMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmasksoftMethod->set_active(0);

    if (showtooltip) {
        showmasksoftMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKSOFT_TOOLTIP"));
    }

    showmasksoftMethodConn = showmasksoftMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSoft::showmasksoftMethodChanged));

    streng->setAdjusterListener(this);

    laplace->setAdjusterListener(this);

    sensisf->setAdjusterListener(this);

    // Add Soft light specific widgets to GUI
    if (complexsoft < 2) {
        pack_start(*softMethod);
    }

    Gtk::Label* const labelsoftmethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_SHOWDCT") + ":"));
    ctboxsoftmethod->pack_start(*labelsoftmethod, Gtk::PACK_SHRINK, 4);
    ctboxsoftmethod->pack_start(*showmasksoftMethod);
    pack_start(*ctboxsoftmethod);
    pack_start(*streng);
    pack_start(*laplace);
    pack_start(*sensisf);
}

void LocallabSoft::resetMaskView()
{
    showmasksoftMethodConn.block(true);
    showmasksoftMethod->set_active(0);
    showmasksoftMethodConn.block(false);
}

void LocallabSoft::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    softMask = showmasksoftMethod->get_active_row_number();
}

void LocallabSoft::disableListener()
{
    LocallabTool::disableListener();

    softMethodConn.block(true);
    showmasksoftMethodConn.block(true);
}

void LocallabSoft::enableListener()
{
    LocallabTool::enableListener();

    softMethodConn.block(false);
    showmasksoftMethodConn.block(false);
}

void LocallabSoft::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visisoft);
        exp->setEnabled(pp->locallab.spots.at(index).expsoft);

        if (complexsoft < 2) {
            if (pp->locallab.spots.at(index).softMethod == "soft") {
                softMethod->set_active(0);
            } else if (pp->locallab.spots.at(index).softMethod == "reti") {
                softMethod->set_active(1);
            }

            streng->setValue((double)pp->locallab.spots.at(index).streng);
        } else {
            softMethod->set_active(0);
            streng->setValue(1.);
        }

        sensisf->setValue((double)pp->locallab.spots.at(index).sensisf);
        laplace->setValue(pp->locallab.spots.at(index).laplace);
    }

    // Enable all listeners
    enableListener();

    // Update soft light GUI according to softMethod combobox
    updateSoftGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSoft::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expsoft = exp->getEnabled();
        pp->locallab.spots.at(index).visisoft = exp->get_visible();

        if (softMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).softMethod = "soft";
        } else if (softMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).softMethod = "reti";
        }

        pp->locallab.spots.at(index).streng = streng->getIntValue();
        pp->locallab.spots.at(index).sensisf = sensisf->getIntValue();
        pp->locallab.spots.at(index).laplace = laplace->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSoft::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster widgets
        streng->setDefault((double)defSpot.streng);
        laplace->setDefault(defSpot.laplace);
        sensisf->setDefault((double)defSpot.sensisf);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSoft::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == streng) {
            if (listener) {
                listener->panelChanged(Evlocallabstreng,
                                       streng->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensisf) {
            if (listener) {
                listener->panelChanged(Evlocallabsensisf,
                                       sensisf->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == laplace) {
            if (listener) {
                listener->panelChanged(Evlocallablaplace,
                                       laplace->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabSoft::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenasoft,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenasoft,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabSoft::softMethodChanged()
{
    // Update soft light GUI according to softMethod combobox
    updateSoftGUI();

    // This event is called to transmit potentially resetted mask state
    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabsoftMethod,
                                   softMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabSoft::showmasksoftMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabSoft::updateSoftGUI()
{
    // Update soft light GUI according to softMethod combobox
    if (softMethod->get_active_row_number() == 0) {
        ctboxsoftmethod->hide();
        // Reset hidden mask combobox
        showmasksoftMethodConn.block(true);
        showmasksoftMethod->set_active(0);
        showmasksoftMethodConn.block(false);
        laplace->hide();
    } else {
        ctboxsoftmethod->show();
        laplace->show();
    }
}

/* ==== LocallabBlur ==== */
LocallabBlur::LocallabBlur():
    LocallabTool(this, M("TP_LOCALLAB_BLUR_TOOLNAME"), M("TP_LOCALLAB_BLUFR"), true),

    // Blur, Noise & Denoise specific widgets
    blMethod(Gtk::manage(new MyComboBoxText())),
    fftwbl(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_FFTWBLUR")))),
    radius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADIUS"), MINRAD, MAXRAD, 0.1, 1.5, nullptr, nullptr, &blurSlider2radius, &blurRadius2Slider))),
    strength(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0))),
    grainFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GRAINFRA")))),
    isogr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ISOGR"), 20, 6400, 1, 0))),
    strengr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGR"), 0, 100, 1, 0))),
    scalegr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALEGR"), 0, 100, 1, 100))),
    medMethod(Gtk::manage(new MyComboBoxText())),
    itera(Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_MEDIAN_PASSES"), 1, 4, 1, 1))),
    guidbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GUIDBL"), 0, 1000, 1, 0))),
    strbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRBL"), 0, 100, 1, 50))),
    epsbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_EPSBL"), -10, 10, 1, 0))),
    sensibn(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 40))),
    blurMethod(Gtk::manage(new MyComboBoxText())),
    chroMethod(Gtk::manage(new MyComboBoxText())),
    activlum(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ACTIV")))),
    LocalcurveEditorwavden(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVDEN"))),
    wavshapeden(static_cast<FlatCurveEditor*>(LocalcurveEditorwavden->addCurve(CT_Flat, "", nullptr, false, false))),
    noiselumf0(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINEZERO"), MINCHRO, MAXCHRO, 0.01, 0.))),
    noiselumf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINE"), MINCHRO, MAXCHRO, 0.01, 0.))),
    noiselumf2(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINETWO"), MINCHRO, MAXCHRO, 0.01, 0.))),
    noiselumc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMCOARSE"), MINCHRO, MAXCHROCC, 0.01, 0.))),
    noiselumdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMDETAIL"), 0., 100., 0.01, 0.))),
    noiselequal(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELEQUAL"), -2, 10, 1, 7, Gtk::manage(new RTImage("circle-white-small.png")), Gtk::manage(new RTImage("circle-black-small.png"))))),
    noisechrof(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROFINE"), MINCHRO, MAXCHRO, 0.01, 0.))),
    noisechroc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROCOARSE"), MINCHRO, MAXCHROCC, 0.01, 0.))),
    noisechrodetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHRODETAIL"), 0., 100., 0.01, 0.))),
    detailthr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_DETAILTHR"), 0, 100, 1, 0))),
    adjblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJ"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-red-small.png"))))),
    bilateral(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BILATERAL"), 0, 100, 1, 0))),
    sensiden(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIDEN"), 0, 100, 1, 20))),
    expmaskbl(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWPLUS")))),
    showmaskblMethod(Gtk::manage(new MyComboBoxText())),
    enablMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskblCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskblshape(static_cast<FlatCurveEditor*>(maskblCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskblshape(static_cast<FlatCurveEditor*>(maskblCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskblshape(static_cast<FlatCurveEditor *>(maskblCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    strumaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUMASKCOL"), 0., 200., 0.1, 0.))),
    toolbl(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TOOLCOL")))),
    blendmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    toolblFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_TOOLMASK")))),
    radmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    lapmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    shadmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_HIGHMASKCOL"), 0, 100, 1, 0))),
    mask2blCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    Lmaskblshape(static_cast<DiagonalCurveEditor*>(mask2blCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    mask2blCurveEditorGwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVMASK"))),
    LLmaskblshapewav(static_cast<FlatCurveEditor*>(mask2blCurveEditorGwav->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    csThresholdblur(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLDBLUR"), 0, 9, 0, 0, 6, 5, 0, false)))
{
    const LocallabParams::LocallabSpot defSpot;

    const bool showtooltip = options.showtooltip;
    const int complexsoft = options.complexity;

    // Parameter Blur, Noise & Denoise specific widgets


    blMethod->append(M("TP_LOCALLAB_BLUR"));
    blMethod->append(M("TP_LOCALLAB_BLMED"));
    blMethod->append(M("TP_LOCALLAB_BLGUID"));
    blMethod->set_active(0);
    blMethodConn = blMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::blMethodChanged));

    fftwblConn = fftwbl->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::fftwblChanged));

    if (showtooltip) {
        radius->set_tooltip_text(M("TP_LOCALLAB_RADIUS_TOOLTIP"));
    }

    radius->setAdjusterListener(this);

    strength->setAdjusterListener(this);

    grainFrame->set_label_align(0.025, 0.5);

    isogr->setAdjusterListener(this);

    strengr->setAdjusterListener(this);

    scalegr->setAdjusterListener(this);

    medMethod->append(M("TP_LOCALLAB_MEDNONE"));
    medMethod->append(M("TP_DIRPYRDENOISE_TYPE_3X3"));
    medMethod->append(M("TP_DIRPYRDENOISE_TYPE_5X5"));
    medMethod->append(M("TP_DIRPYRDENOISE_TYPE_7X7"));
    medMethod->append(M("TP_DIRPYRDENOISE_TYPE_9X9"));
    medMethod->set_active(0);
    medMethodConn = medMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::medMethodChanged));

    itera->setAdjusterListener(this);

    guidbl->setLogScale(100, 0);
    guidbl->setAdjusterListener(this);

    strbl->setAdjusterListener(this);

    epsbl->setAdjusterListener(this);

    if (showtooltip) {
        sensibn->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    }

    sensibn->setAdjusterListener(this);

    blurMethod->append(M("TP_LOCALLAB_BLNORM"));
    blurMethod->append(M("TP_LOCALLAB_BLINV"));
    blurMethod->set_active(0);

    if (showtooltip) {
        blurMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));
    }

    blurMethodConn = blurMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::blurMethodChanged));

    chroMethod->append(M("TP_LOCALLAB_BLLO"));
    chroMethod->append(M("TP_LOCALLAB_BLCO"));
    chroMethod->append(M("TP_LOCALLAB_BLLC"));
    chroMethod->set_active(0);
    chroMethodConn = chroMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::chroMethodChanged));

    activlumConn = activlum->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::activlumChanged));

    LocalcurveEditorwavden->setCurveListener(this);

    wavshapeden->setIdentityValue(0.);
    wavshapeden->setResetCurve(FlatCurveType(defSpot.locwavcurveden.at(0)), defSpot.locwavcurveden);

    if (showtooltip) {
        wavshapeden->setTooltip(M("TP_LOCALLAB_WASDEN_TOOLTIP"));
    }

    LocalcurveEditorwavden->curveListComplete();

    noiselumf0->setAdjusterListener(this);

    noiselumf->setAdjusterListener(this);

    noiselumf2->setAdjusterListener(this);

    if (showtooltip) {
        noiselumc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    }

    noiselumc->setAdjusterListener(this);

    noiselumdetail->setAdjusterListener(this);

    noiselequal->setAdjusterListener(this);

    noisechrof->setAdjusterListener(this);

    if (showtooltip) {
        noisechroc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    }

    noisechroc->setAdjusterListener(this);

    noisechrodetail->setAdjusterListener(this);

    detailthr->setAdjusterListener(this);

    adjblur->setAdjusterListener(this);

    bilateral->setAdjusterListener(this);

    sensiden->setAdjusterListener(this);

    setExpandAlignProperties(expmaskbl, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expmaskbl->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
    }

    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskblMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
    showmaskblMethod->set_active(0);

    if (showtooltip) {
        showmaskblMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    }

    showmaskblMethodConn = showmaskblMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::showmaskblMethodChanged));

    enablMaskConn = enablMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::enablMaskChanged));

    maskblCurveEditorG->setCurveListener(this);

    CCmaskblshape->setIdentityValue(0.);
    CCmaskblshape->setResetCurve(FlatCurveType(defSpot.CCmaskblcurve.at(0)), defSpot.CCmaskblcurve);

    if (showtooltip) {
        CCmaskblshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    CCmaskblshape->setBottomBarColorProvider(this, 1);

    LLmaskblshape->setIdentityValue(0.);
    LLmaskblshape->setResetCurve(FlatCurveType(defSpot.LLmaskblcurve.at(0)), defSpot.LLmaskblcurve);

    if (showtooltip) {
        LLmaskblshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    LLmaskblshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskblshape->setIdentityValue(0.);
    HHmaskblshape->setResetCurve(FlatCurveType(defSpot.HHmaskblcurve.at(0)), defSpot.HHmaskblcurve);

    if (showtooltip) {
        HHmaskblshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
    }

    HHmaskblshape->setCurveColorProvider(this, 2);
    HHmaskblshape->setBottomBarColorProvider(this, 2);

    maskblCurveEditorG->curveListComplete();

    strumaskbl->setAdjusterListener(this);

    toolblConn = toolbl->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::toolblChanged));

    blendmaskbl->setAdjusterListener(this);

    toolblFrame->set_label_align(0.025, 0.5);

    if (showtooltip) {
        radmaskbl->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    radmaskbl->setAdjusterListener(this);

    if (showtooltip) {
        lapmaskbl->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
    }

    lapmaskbl->setAdjusterListener(this);

    chromaskbl->setAdjusterListener(this);

    gammaskbl->setAdjusterListener(this);

    slomaskbl->setAdjusterListener(this);

    shadmaskbl->setAdjusterListener(this);

    mask2blCurveEditorG->setCurveListener(this);

    Lmaskblshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskblcurve.at(0)), defSpot.Lmaskblcurve);

    if (showtooltip) {
        Lmaskblshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
    }

    Lmaskblshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskblshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    mask2blCurveEditorG->curveListComplete();

    mask2blCurveEditorGwav->setCurveListener(this);

    LLmaskblshapewav->setIdentityValue(0.);
    LLmaskblshapewav->setResetCurve(FlatCurveType(defSpot.LLmaskblcurvewav.at(0)), defSpot.LLmaskblcurvewav);

    if (showtooltip) {
        LLmaskblshapewav->setTooltip(M("TP_LOCALLAB_LMASK_LEVEL_TOOLTIP"));
    }

    LLmaskblshapewav->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2blCurveEditorGwav->curveListComplete();

    csThresholdblur->setAdjusterListener(this);

    // Add Blur, Noise & Denoise specific widgets to GUI
    MyExpander* const expblnoise = Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_BLNOI_EXP")));
    setExpandAlignProperties(expblnoise, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expblnoise->set_tooltip_markup(M("TP_LOCALLAB_BLUMETHOD_TOOLTIP"));
    }

    expblnoise->set_expanded(false);
    ToolParamBlock* const blnoisebox = Gtk::manage(new ToolParamBlock());
    blnoisebox->pack_start(*blMethod);

    if (complexsoft < 2) {
        blnoisebox->pack_start(*fftwbl, Gtk::PACK_SHRINK, 0);
    }

    blnoisebox->pack_start(*radius);
    blnoisebox->pack_start(*strength);
    ToolParamBlock* const grainBox = Gtk::manage(new ToolParamBlock());
    grainBox->pack_start(*isogr);
    grainBox->pack_start(*strengr);
    grainBox->pack_start(*scalegr);
    grainFrame->add(*grainBox);
    blnoisebox->pack_start(*grainFrame);
    blnoisebox->pack_start(*medMethod);
    blnoisebox->pack_start(*itera);
    blnoisebox->pack_start(*guidbl);
    blnoisebox->pack_start(*strbl);
    blnoisebox->pack_start(*epsbl);
    blnoisebox->pack_start(*sensibn);
    blnoisebox->pack_start(*blurMethod);
    blnoisebox->pack_start(*chroMethod);
    // blnoisebox->pack_start(*activlum);
    expblnoise->add(*blnoisebox, false);
    pack_start(*expblnoise);
    MyExpander* const expdenoise = Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI_EXP")));
    setExpandAlignProperties(expdenoise, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    if (showtooltip) {
        expdenoise->set_tooltip_markup(M("TP_LOCALLAB_DENOI_TOOLTIP"));
    }

    expdenoise->set_expanded(false);
    ToolParamBlock* const denoisebox = Gtk::manage(new ToolParamBlock());
    Gtk::Frame* const wavFrame = Gtk::manage(new Gtk::Frame());
    ToolParamBlock* const wavBox = Gtk::manage(new ToolParamBlock());
    wavBox->pack_start(*LocalcurveEditorwavden, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    // wavBox->pack_start(*noiselumf0);
    // wavBox->pack_start(*noiselumf);
    // wavBox->pack_start(*noiselumf2);
    // wavBox->pack_start(*noiselumc);
    wavBox->pack_start(*noiselumdetail);
    wavBox->pack_start(*noiselequal);
    wavBox->pack_start(*noisechrof);
    wavBox->pack_start(*noisechroc);
    wavBox->pack_start(*noisechrodetail);
    wavBox->pack_start(*detailthr);
    wavBox->pack_start(*adjblur);
    wavFrame->add(*wavBox);
    denoisebox->pack_start(*wavFrame);
    denoisebox->pack_start(*bilateral);
    denoisebox->pack_start(*sensiden);
    expdenoise->add(*denoisebox, false);
    pack_start(*expdenoise);
    ToolParamBlock* const maskblBox = Gtk::manage(new ToolParamBlock());
    maskblBox->pack_start(*showmaskblMethod, Gtk::PACK_SHRINK, 4);
    maskblBox->pack_start(*enablMask, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        maskblBox->pack_start(*maskblCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        maskblBox->pack_start(*strumaskbl, Gtk::PACK_SHRINK, 0);
        maskblBox->pack_start(*toolbl, Gtk::PACK_SHRINK, 0);
    }

    Gtk::HSeparator* const separatorstrubl = Gtk::manage(new  Gtk::HSeparator());
    maskblBox->pack_start(*separatorstrubl, Gtk::PACK_SHRINK, 2);
    maskblBox->pack_start(*blendmaskbl, Gtk::PACK_SHRINK, 0);
    ToolParamBlock* const toolblBox = Gtk::manage(new ToolParamBlock());
    toolblBox->pack_start(*radmaskbl, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 1) {
        toolblBox->pack_start(*lapmaskbl, Gtk::PACK_SHRINK, 0);
    }

    toolblBox->pack_start(*chromaskbl, Gtk::PACK_SHRINK, 0);

    if (complexsoft < 2) {
        toolblBox->pack_start(*gammaskbl, Gtk::PACK_SHRINK, 0);
        toolblBox->pack_start(*slomaskbl, Gtk::PACK_SHRINK, 0);
        toolblBox->pack_start(*shadmaskbl, Gtk::PACK_SHRINK, 0);
    }

    toolblBox->pack_start(*mask2blCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor

    if (complexsoft < 1) {
        toolblBox->pack_start(*mask2blCurveEditorGwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        toolblBox->pack_start(*csThresholdblur, Gtk::PACK_SHRINK, 0);
    }

    toolblFrame->add(*toolblBox);
    maskblBox->pack_start(*toolblFrame);
    expmaskbl->add(*maskblBox, false);
    pack_start(*expmaskbl);
}

LocallabBlur::~LocallabBlur()
{
    delete maskblCurveEditorG;
    delete mask2blCurveEditorG;
    delete mask2blCurveEditorGwav;
    delete LocalcurveEditorwavden;
}

void LocallabBlur::resetMaskView()
{
    showmaskblMethodConn.block(true);
    showmaskblMethod->set_active(0);
    showmaskblMethodConn.block(false);
}

void LocallabBlur::getMaskView(int &colorMask, int &colorMaskinv, int &expMask, int &expMaskinv, int &shMask, int &shMaskinv, int &vibMask, int &softMask, int &blMask, int &tmMask, int &retiMask, int &sharMask, int &lcMask, int &cbMask)
{
    blMask = showmaskblMethod->get_active_row_number();
}

void LocallabBlur::disableListener()
{
    LocallabTool::disableListener();

    blMethodConn.block(true);
    fftwblConn.block(true);
    medMethodConn.block(true);
    blurMethodConn.block(true);
    chroMethodConn.block(true);
    activlumConn.block(true);
    showmaskblMethodConn.block(true);
    enablMaskConn.block(true);
    toolblConn.block(true);
}

void LocallabBlur::enableListener()
{
    LocallabTool::enableListener();

    blMethodConn.block(false);
    fftwblConn.block(false);
    medMethodConn.block(false);
    blurMethodConn.block(false);
    chroMethodConn.block(false);
    activlumConn.block(false);
    showmaskblMethodConn.block(false);
    enablMaskConn.block(false);
    toolblConn.block(false);
}

void LocallabBlur::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    const int complexsoft = options.complexity;

    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->set_visible(pp->locallab.spots.at(index).visiblur);
        exp->setEnabled(pp->locallab.spots.at(index).expblur);

        if (pp->locallab.spots.at(index).blMethod == "blur") {
            blMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).blMethod == "med") {
            blMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).blMethod == "guid") {
            blMethod->set_active(2);
        }

        if (complexsoft < 2) {
            fftwbl->set_active(pp->locallab.spots.at(index).fftwbl);
        } else {
            fftwbl->set_active(false);
        }

        radius->setValue(pp->locallab.spots.at(index).radius);
        strength->setValue(pp->locallab.spots.at(index).strength);
        isogr->setValue((double)pp->locallab.spots.at(index).isogr);
        strengr->setValue((double)pp->locallab.spots.at(index).strengr);
        scalegr->setValue((double)pp->locallab.spots.at(index).scalegr);

        if (complexsoft < 2) {
            if (pp->locallab.spots.at(index).medMethod == "none") {
                medMethod->set_active(0);
            } else if (pp->locallab.spots.at(index).medMethod == "33") {
                medMethod->set_active(1);
            } else if (pp->locallab.spots.at(index).medMethod == "55") {
                medMethod->set_active(2);
            } else if (pp->locallab.spots.at(index).medMethod == "77") {
                medMethod->set_active(3);
            } else if (pp->locallab.spots.at(index).medMethod == "99") {
                medMethod->set_active(4);
            }
        } else {
            medMethod->set_active(0);
        }

        itera->setValue((double)pp->locallab.spots.at(index).itera);
        guidbl->setValue((double)pp->locallab.spots.at(index).guidbl);
        strbl->setValue((double)pp->locallab.spots.at(index).strbl);
        epsbl->setValue((double)pp->locallab.spots.at(index).epsbl);
        sensibn->setValue((double)pp->locallab.spots.at(index).sensibn);

        if (pp->locallab.spots.at(index).blurMethod == "norm") {
            blurMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).blurMethod == "inv") {
            blurMethod->set_active(1);
        }

        if (pp->locallab.spots.at(index).chroMethod == "lum") {
            chroMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).chroMethod == "chr") {
            chroMethod->set_active(1);
        } else if (pp->locallab.spots.at(index).chroMethod == "all") {
            chroMethod->set_active(2);
        }

        activlum->set_active(pp->locallab.spots.at(index).activlum);
        wavshapeden->setCurve(pp->locallab.spots.at(index).locwavcurveden);
        noiselumf0->setValue(pp->locallab.spots.at(index).noiselumf0);
        noiselumf->setValue(pp->locallab.spots.at(index).noiselumf);
        noiselumf2->setValue(pp->locallab.spots.at(index).noiselumf2);
        noiselumc->setValue(pp->locallab.spots.at(index).noiselumc);
        noiselumdetail->setValue(pp->locallab.spots.at(index).noiselumdetail);
        noiselequal->setValue((double)pp->locallab.spots.at(index).noiselequal);
        noisechrof->setValue(pp->locallab.spots.at(index).noisechrof);
        noisechroc->setValue(pp->locallab.spots.at(index).noisechroc);
        noisechrodetail->setValue(pp->locallab.spots.at(index).noisechrodetail);
        detailthr->setValue((double)pp->locallab.spots.at(index).detailthr);
        adjblur->setValue((double)pp->locallab.spots.at(index).adjblur);
        bilateral->setValue((double)pp->locallab.spots.at(index).bilateral);
        sensiden->setValue((double)pp->locallab.spots.at(index).sensiden);
        enablMask->set_active(pp->locallab.spots.at(index).enablMask);
        CCmaskblshape->setCurve(pp->locallab.spots.at(index).CCmaskblcurve);
        LLmaskblshape->setCurve(pp->locallab.spots.at(index).LLmaskblcurve);
        HHmaskblshape->setCurve(pp->locallab.spots.at(index).HHmaskblcurve);

        if (complexsoft < 2) {
            strumaskbl->setValue(pp->locallab.spots.at(index).strumaskbl);
        } else {
            strumaskbl->setValue(0.);
        }

        toolbl->set_active(pp->locallab.spots.at(index).toolbl);
        blendmaskbl->setValue((double)pp->locallab.spots.at(index).blendmaskbl);
        radmaskbl->setValue(pp->locallab.spots.at(index).radmaskbl);

        if (complexsoft == 0) {
            lapmaskbl->setValue(pp->locallab.spots.at(index).lapmaskbl);
        } else {
            lapmaskbl->setValue(0.);
        }

        chromaskbl->setValue(pp->locallab.spots.at(index).chromaskbl);

        if (complexsoft < 2) {
            gammaskbl->setValue(pp->locallab.spots.at(index).gammaskbl);
            slomaskbl->setValue(pp->locallab.spots.at(index).slomaskbl);
            shadmaskbl->setValue((double)pp->locallab.spots.at(index).shadmaskbl);
        } else {
            gammaskbl->setValue(1.);
            slomaskbl->setValue(0.);
            shadmaskbl->setValue(0.);
        }

        Lmaskblshape->setCurve(pp->locallab.spots.at(index).Lmaskblcurve);

        if (complexsoft == 0) {
            LLmaskblshapewav->setCurve(pp->locallab.spots.at(index).LLmaskblcurvewav);
        } else {
            LLmaskblshapewav->reset();
        }

        csThresholdblur->setValue<int>(pp->locallab.spots.at(index).csthresholdblur);
    }

    // Enable all listeners
    enableListener();

    // Update Blur & Noise GUI according to blMethod combobox state
    updateBlurGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expblur = exp->getEnabled();
        pp->locallab.spots.at(index).visiblur = exp->get_visible();

        if (blMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).blMethod = "blur";
        } else if (blMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).blMethod = "med";
        } else if (blMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).blMethod = "guid";
        }

        pp->locallab.spots.at(index).fftwbl = fftwbl->get_active();
        pp->locallab.spots.at(index).radius = radius->getValue();
        pp->locallab.spots.at(index).strength = strength->getIntValue();
        pp->locallab.spots.at(index).isogr = isogr->getIntValue();
        pp->locallab.spots.at(index).strengr = strengr->getIntValue();
        pp->locallab.spots.at(index).scalegr = scalegr->getIntValue();

        if (medMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).medMethod = "none";
        } else if (medMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).medMethod = "33";
        } else if (medMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).medMethod = "55";
        } else if (medMethod->get_active_row_number() == 3) {
            pp->locallab.spots.at(index).medMethod = "77";
        } else if (medMethod->get_active_row_number() == 4) {
            pp->locallab.spots.at(index).medMethod = "99";
        }

        pp->locallab.spots.at(index).itera = itera->getIntValue();
        pp->locallab.spots.at(index).guidbl = guidbl->getIntValue();
        pp->locallab.spots.at(index).strbl = strbl->getIntValue();
        pp->locallab.spots.at(index).epsbl = epsbl->getIntValue();
        pp->locallab.spots.at(index).sensibn = sensibn->getIntValue();

        if (blurMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).blurMethod = "norm";
        } else if (blurMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).blurMethod = "inv";
        }

        if (chroMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).chroMethod = "lum";
        } else if (chroMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).chroMethod = "chr";
        } else if (chroMethod->get_active_row_number() == 2) {
            pp->locallab.spots.at(index).chroMethod = "all";
        }

        pp->locallab.spots.at(index).activlum = activlum->get_active();
        pp->locallab.spots.at(index).locwavcurveden = wavshapeden->getCurve();
        pp->locallab.spots.at(index).noiselumf0 = noiselumf0->getValue();
        pp->locallab.spots.at(index).noiselumf = noiselumf->getValue();
        pp->locallab.spots.at(index).noiselumf2 = noiselumf2->getValue();
        pp->locallab.spots.at(index).noiselumc = noiselumc->getValue();
        pp->locallab.spots.at(index).noiselumdetail = noiselumdetail->getValue();
        pp->locallab.spots.at(index).noiselequal = noiselequal->getIntValue();
        pp->locallab.spots.at(index).noisechrof = noisechrof->getValue();
        pp->locallab.spots.at(index).noisechroc = noisechroc->getValue();
        pp->locallab.spots.at(index).noisechrodetail = noisechrodetail->getValue();
        pp->locallab.spots.at(index).detailthr = detailthr->getIntValue();
        pp->locallab.spots.at(index).adjblur = adjblur->getIntValue();
        pp->locallab.spots.at(index).bilateral = bilateral->getIntValue();
        pp->locallab.spots.at(index).sensiden = sensiden->getIntValue();
        pp->locallab.spots.at(index).enablMask = enablMask->get_active();
        pp->locallab.spots.at(index).LLmaskblcurve = LLmaskblshape->getCurve();
        pp->locallab.spots.at(index).CCmaskblcurve = CCmaskblshape->getCurve();
        pp->locallab.spots.at(index).HHmaskblcurve = HHmaskblshape->getCurve();
        pp->locallab.spots.at(index).strumaskbl = strumaskbl->getValue();
        pp->locallab.spots.at(index).toolbl = toolbl->get_active();
        pp->locallab.spots.at(index).blendmaskbl = blendmaskbl->getIntValue();
        pp->locallab.spots.at(index).radmaskbl = radmaskbl->getValue();
        pp->locallab.spots.at(index).lapmaskbl = lapmaskbl->getValue();
        pp->locallab.spots.at(index).chromaskbl = chromaskbl->getValue();
        pp->locallab.spots.at(index).gammaskbl = gammaskbl->getValue();
        pp->locallab.spots.at(index).slomaskbl = slomaskbl->getValue();
        pp->locallab.spots.at(index).shadmaskbl = shadmaskbl->getIntValue();
        pp->locallab.spots.at(index).Lmaskblcurve = Lmaskblshape->getCurve();
        pp->locallab.spots.at(index).LLmaskblcurvewav = LLmaskblshapewav->getCurve();
        pp->locallab.spots.at(index).csthresholdblur = csThresholdblur->getValue<int>();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster and threshold adjuster widgets
        radius->setDefault(defSpot.radius);
        strength->setDefault((double)defSpot.strength);
        isogr->setDefault((double)defSpot.isogr);
        strengr->setDefault((double)defSpot.strengr);
        scalegr->setDefault((double)defSpot.scalegr);
        itera->setDefault((double)defSpot.itera);
        guidbl->setDefault((double)defSpot.guidbl);
        strbl->setDefault((double)defSpot.strbl);
        epsbl->setDefault((double)defSpot.epsbl);
        sensibn->setDefault((double)defSpot.sensibn);
        noiselumf0->setDefault(defSpot.noiselumf0);
        noiselumf->setDefault(defSpot.noiselumf);
        noiselumf2->setDefault(defSpot.noiselumf2);
        noiselumc->setDefault(defSpot.noiselumc);
        noiselumdetail->setDefault(defSpot.noiselumdetail);
        noiselequal->setDefault((double)defSpot.noiselequal);
        noisechrof->setDefault(defSpot.noisechrof);
        noisechroc->setDefault(defSpot.noisechroc);
        noisechrodetail->setDefault(defSpot.noisechrodetail);
        detailthr->setDefault((double)defSpot.detailthr);
        adjblur->setDefault((double)defSpot.adjblur);
        bilateral->setDefault((double)defSpot.bilateral);
        sensiden->setDefault((double)defSpot.sensiden);
        strumaskbl->setDefault(defSpot.strumaskbl);
        blendmaskbl->setDefault((double)defSpot.blendmaskbl);
        radmaskbl->setDefault(defSpot.radmaskbl);
        lapmaskbl->setDefault(defSpot.lapmaskbl);
        chromaskbl->setDefault(defSpot.chromaskbl);
        gammaskbl->setDefault(defSpot.gammaskbl);
        slomaskbl->setDefault(defSpot.slomaskbl);
        shadmaskbl->setDefault(defSpot.shadmaskbl);
        csThresholdblur->setDefault<int>(defSpot.csthresholdblur);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
        if (a == radius) {
            if (listener) {
                listener->panelChanged(Evlocallabradius,
                                       radius->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strength) {
            if (listener) {
                listener->panelChanged(Evlocallabstrength,
                                       strength->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == isogr) {
            if (listener) {
                listener->panelChanged(Evlocallabisogr,
                                       isogr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strengr) {
            if (listener) {
                listener->panelChanged(Evlocallabstrengr,
                                       strengr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == scalegr) {
            if (listener) {
                listener->panelChanged(Evlocallabscalegr,
                                       scalegr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == itera) {
            if (listener) {
                listener->panelChanged(Evlocallabitera,
                                       itera->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == guidbl) {
            if (listener) {
                listener->panelChanged(Evlocallabguidbl,
                                       guidbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strbl) {
            if (listener) {
                listener->panelChanged(Evlocallabstrbl,
                                       strbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == epsbl) {
            if (listener) {
                listener->panelChanged(Evlocallabepsbl,
                                       epsbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensibn) {
            if (listener) {
                listener->panelChanged(Evlocallabsensibn,
                                       sensibn->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumf0) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf0,
                                       noiselumf0->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumf) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf,
                                       noiselumf->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumf2) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf2,
                                       noiselumf2->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumc,
                                       noiselumc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselumdetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumdetail,
                                       noiselumdetail->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noiselequal) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselequal,
                                       noiselequal->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noisechrof) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrof,
                                       noisechrof->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noisechroc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechroc,
                                       noisechroc->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == noisechrodetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrodetail,
                                       noisechrodetail->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == detailthr) {
            if (listener) {
                listener->panelChanged(Evlocallabdetailthr,
                                       detailthr->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == adjblur) {
            if (listener) {
                listener->panelChanged(Evlocallabadjblur,
                                       adjblur->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == bilateral) {
            if (listener) {
                listener->panelChanged(Evlocallabbilateral,
                                       bilateral->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sensiden) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiden,
                                       sensiden->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == strumaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabstrumaskbl,
                                       strumaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendmaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskbl,
                                       blendmaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radmaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskbl,
                                       radmaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == lapmaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallablapmaskbl,
                                       lapmaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chromaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskbl,
                                       chromaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gammaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskbl,
                                       gammaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == slomaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskbl,
                                       slomaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == shadmaskbl) {
            if (listener) {
                listener->panelChanged(Evlocallabshadmaskbl,
                                       shadmaskbl->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR)
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabcsThresholdblur,
                                   csThresholdblur->getHistoryString() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabBlur::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == wavshapeden) {
            if (listener) {
                listener->panelChanged(EvlocallabwavCurveden,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == CCmaskblshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskblshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskblshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskblshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHmaskblshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskblshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == Lmaskblshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLmaskblshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLmaskblshapewav) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskblshapewav,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::enabledChanged()
{
    if (isLocActivated) {
        if (listener) {
            if (exp->getEnabled()) {
                listener->panelChanged(EvLocenablur,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocenablur,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::updateMaskBackground(const double normChromar, const double normLumar, const double normHuer)
{
    idle_register.add(
    [this, normHuer, normLumar, normChromar]() -> bool {
        GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

        // Update mask background
        CCmaskblshape->updateLocallabBackground(normChromar);
        LLmaskblshape->updateLocallabBackground(normLumar);
        HHmaskblshape->updateLocallabBackground(normHuer);

        return false;
    }
    );
}

void LocallabBlur::blMethodChanged()
{
    // Update Blur & Noise GUI according to blMethod combobox state
    updateBlurGUI();

    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabblMethod,
                                   blMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabBlur::fftwblChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (fftwbl->get_active()) {
                listener->panelChanged(Evlocallabfftwbl,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabfftwbl,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::medMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabmedMethod,
                                   medMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabBlur::blurMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabblurMethod,
                                   blurMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabBlur::chroMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabchroMethod,
                                   chroMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabBlur::activlumChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (activlum->get_active()) {
                listener->panelChanged(Evlocallabactivlum,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabactivlum,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::showmaskblMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabBlur::enablMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enablMask->get_active()) {
                listener->panelChanged(EvLocallabEnablMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnablMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::toolblChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (toolbl->get_active()) {
                listener->panelChanged(Evlocallabtoolbl,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(Evlocallabtoolbl,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabBlur::updateBlurGUI()
{
    if (blMethod->get_active_row_number() == 0) {
        fftwbl->show();
        radius->show();
        strength->show();
        grainFrame->show();
        medMethod->hide();
        itera->hide();
        guidbl->hide();
        strbl->hide();
        epsbl->hide();
        activlum->show();
    } else if (blMethod->get_active_row_number() == 1) {
        fftwbl->hide();
        radius->hide();
        strength->hide();
        grainFrame->hide();
        medMethod->show();
        itera->show();
        guidbl->hide();
        strbl->hide();
        epsbl->hide();
        activlum->show();
    } else if (blMethod->get_active_row_number() == 2) {
        fftwbl->hide();
        radius->hide();
        strength->hide();
        grainFrame->hide();
        medMethod->hide();
        itera->hide();
        guidbl->show();
        strbl->show();
        epsbl->show();
        activlum->hide();
    }
}
