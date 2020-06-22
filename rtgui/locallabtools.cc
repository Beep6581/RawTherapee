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
#define MINEXP -1.5
#define MAXEXP 1.5

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
LocallabTool::LocallabTool(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11, bool needMode):
    ToolPanel(toolName, need11),

    // LocallabTool parameters
    needMode(needMode),
    isLocActivated(false),
    locToolListener(nullptr),

    // LocallabTool generic widgets
    complexity(Gtk::manage(new MyComboBoxText()))
{
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

    if (needMode) {
        complexity->append(M("TP_LOCALLAB_MODE_EXPERT"));
        complexity->append(M("TP_LOCALLAB_MODE_NORMAL"));
        complexity->set_active(0);
        complexity->setPreferredWidth(100, -1);
        complexityConn = complexity->signal_changed().connect(sigc::mem_fun(*this, &LocallabTool::complexityModeChanged));
        titleBox->pack_end(*complexity, Gtk::PACK_SHRINK, 2);
    }

    Gtk::VSeparator* const separator = Gtk::manage(new Gtk::VSeparator());
    titleBox->pack_end(*separator, Gtk::PACK_SHRINK, 0);

    if (need100Percent) {
        RTImage* const titleImage = Gtk::manage(new RTImage("one-to-one-small.png"));
        titleImage->set_tooltip_text(M("TP_GENERAL_11SCALE_TOOLTIP"));
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
    if (raiseEvent) { // Note: Event is only raised when a tool is added by user
        // By default, activate newly added tool
        enaExpConn.block(true);
        exp->setEnabled(true);
        enaExpConn.block(false);

        if (needMode) {
            // Set complexity mode according to chosen default one
            complexityConn.block(true);
            complexity->set_active(options.complexity);
            complexityConn.block(false);

            // Update GUI accordingly
            if (complexity->get_active_row_number() == Normal) {
                convertParamToNormal();
                updateGUIToMode(Normal);
            } else {
                updateGUIToMode(Expert);
            }
        }

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

    if (exp->getEnabled() || isMaskViewActive()) {
        // Disable tool while removing it
        disableListener();
        exp->setEnabled(false);
        enableListener();
        // Note: Mask views are all reset when removing tool (in toolpanelcoord.cc)

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

    if (needMode) {
        complexityConn.block(true);
    }
}
void LocallabTool::enableListener()
{
    ToolPanel::enableListener();

    enaExpConn.block(false);

    if (needMode) {
        complexityConn.block(false);
    }
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

void LocallabTool::complexityModeChanged()
{
    if (complexity->get_active_row_number() == Normal) { // New selected mode is Normal one
        // Convert tool widget parameters
        convertParamToNormal();

        // Update GUI based on new mode
        updateGUIToMode(Normal);

        // Raise event with refreshing
        if (listener) {
            listener->panelChanged(EvlocallabcomplexityWithRefresh,
                                   M("TP_LOCALLAB_MODE_NORMAL") + " (" + escapeHtmlChars(spotName) + ")");
        }
    } else { // New selected mode is Expert one
        // Update GUI based on new mode
        updateGUIToMode(Expert);

        // Raise event without refreshing
        if (listener) {
            listener->panelChanged(EvlocallabcomplexityWithoutRefresh,
                                   M("TP_LOCALLAB_MODE_EXPERT") + " (" + escapeHtmlChars(spotName) + ")");
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
    softradiuscol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 1000.0, 0.5, 0.))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    expgradcol(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPGRAD")))),
    strcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRLUM"), -4., 4., 0.05, 0.))),
    strcolab(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRCHRO"), -6., 6., 0.05, 0.))),
    strcolh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTRHUE"), -6., 6., 0.05, 0.))),
    angcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    expcurvcol(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPCURV")))),
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
    expmaskcol(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWC")))),
    showmaskcolMethod(Gtk::manage(new MyComboBoxText())),
    showmaskcolMethodinv(Gtk::manage(new MyComboBoxText())),
    enaColorMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASKCOL"))),
    CCmaskshape(static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskshape(static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskshape(static_cast<FlatCurveEditor *>(maskCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    struFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABSTRUM")))),
    strumaskcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUMASKCOL"), 0., 200., 0.1, 0.))),
    toolcol(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TOOLCOL")))),
    blurFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABBLURM")))),
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

    // Parameter Color & Light specific widgets
    curvactivConn = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::curvactivChanged));

    lightness->setAdjusterListener(this);

    contrast->setAdjusterListener(this);

    chroma->setAdjusterListener(this);

    gridFrame->set_label_align(0.025, 0.5);

    gridMethod->append(M("TP_LOCALLAB_GRIDONE"));
    gridMethod->append(M("TP_LOCALLAB_GRIDTWO"));
    gridMethod->set_active(0);
    gridMethodConn = gridMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::gridMethodChanged));

    strengthgrid->setAdjusterListener(this);

    sensi->setAdjusterListener(this);

    structcol->setAdjusterListener(this);

    blurcolde->setAdjusterListener(this);

    softradiuscol->setLogScale(10, 0);
    softradiuscol->setAdjusterListener(this);

    inversConn = invers->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::inversChanged));
    invers->set_tooltip_text(M("TP_LOCALLAB_INVERS_TOOLTIP"));

    setExpandAlignProperties(expgradcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    strcol->setAdjusterListener(this);

    strcolab->setAdjusterListener(this);
    strcolab->set_tooltip_text(M("TP_LOCALLAB_GRADSTRAB_TOOLTIP"));

    strcolh->setAdjusterListener(this);
    strcolh->set_tooltip_text(M("TP_LOCALLAB_GRADSTRHUE_TOOLTIP"));

    angcol->setAdjusterListener(this);

    setExpandAlignProperties(expcurvcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->set_active(0);
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::qualitycurveMethodChanged));

    llCurveEditorG->setCurveListener(this);

    llshape->setResetCurve(DiagonalCurveType(defSpot.llcurve.at(0)), defSpot.llcurve);
    llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    llshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    llshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    ccshape->setResetCurve(DiagonalCurveType(defSpot.cccurve.at(0)), defSpot.cccurve);
    ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    ccshape->setBottomBarColorProvider(this, 4);
    ccshape->setLeftBarColorProvider(this, 1);

    llCurveEditorG->curveListComplete();

    clCurveEditorG->setCurveListener(this);

    clshape->setResetCurve(DiagonalCurveType(defSpot.clcurve.at(0)), defSpot.clcurve);
    clshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    clshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    clshape->setLeftBarColorProvider(this, 1);

    lcshape->setResetCurve(DiagonalCurveType(defSpot.lccurve.at(0)), defSpot.lccurve);
    lcshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    lcshape->setBottomBarColorProvider(this, 4);
    lcshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    clCurveEditorG->curveListComplete();

    HCurveEditorG->setCurveListener(this);

    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FlatCurveType(defSpot.LHcurve.at(0)), defSpot.LHcurve);
    LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshape->setCurveColorProvider(this, 3);
    LHshape->setBottomBarBgGradient(six_shape);

    HCurveEditorG->curveListComplete();

    H2CurveEditorG->setCurveListener(this);

    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FlatCurveType(defSpot.HHcurve.at(0)), defSpot.HHcurve);
    HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
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
    rgbshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    rgbshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    rgbshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    rgbCurveEditorG->curveListComplete();

    specialConn  = special->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::specialChanged));

    setExpandAlignProperties(expmaskcol1, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    merMethod->append(M("TP_LOCALLAB_MRONE"));
//    merMethod->append(M("TP_LOCALLAB_MRTWO"));
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

    conthrcol->setAdjusterListener(this);

    merlucol->setAdjusterListener(this);

    setExpandAlignProperties(expmaskcol, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWSTRUC"));
    showmaskcolMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskcolMethod->set_active(0);
    showmaskcolMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcolMethodConn  = showmaskcolMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::showmaskcolMethodChanged));

    showmaskcolMethodinv->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskcolMethodinv->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskcolMethodinv->set_active(0);
    showmaskcolMethodinv->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskcolMethodConninv  = showmaskcolMethodinv->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::showmaskcolMethodChangedinv));

    enaColorMaskConn = enaColorMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::enaColorMaskChanged));

    maskCurveEditorG->setCurveListener(this);

    CCmaskshape->setIdentityValue(0.);
    CCmaskshape->setResetCurve(FlatCurveType(defSpot.CCmaskcurve.at(0)), defSpot.CCmaskcurve);
    CCmaskshape->setBottomBarColorProvider(this, 1);

    LLmaskshape->setIdentityValue(0.);
    LLmaskshape->setResetCurve(FlatCurveType(defSpot.LLmaskcurve.at(0)), defSpot.LLmaskcurve);
    LLmaskshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskshape->setIdentityValue(0.);
    HHmaskshape->setResetCurve(FlatCurveType(defSpot.HHmaskcurve.at(0)), defSpot.HHmaskcurve);
    HHmaskshape->setCurveColorProvider(this, 2);
    HHmaskshape->setBottomBarColorProvider(this, 2);

    maskCurveEditorG->curveListComplete();

    struFrame->set_label_align(0.025, 0.5);

    strumaskcol->setAdjusterListener(this);

    toolcolConn  = toolcol->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::toolcolChanged));

    blurFrame->set_label_align(0.025, 0.5);

    fftColorMaskConn = fftColorMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::fftColorMaskChanged));

    contcol->setAdjusterListener(this);

    blurcol->setAdjusterListener(this);

    blendmaskcol->setAdjusterListener(this);

    radmaskcol->setLogScale(10, -10);
    radmaskcol->setAdjusterListener(this);

    lapmaskcol->setAdjusterListener(this);

    chromaskcol->setAdjusterListener(this);

    gammaskcol->setAdjusterListener(this);

    slomaskcol->setAdjusterListener(this);

    shadmaskcol->setAdjusterListener(this);

    maskHCurveEditorG->setCurveListener(this);

    HHhmaskshape->setIdentityValue(0.);
    HHhmaskshape->setResetCurve(FlatCurveType(defSpot.HHhmaskcurve.at(0)), defSpot.HHhmaskcurve);
//    HHhmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    HHhmaskshape->setCurveColorProvider(this, 2);
    HHhmaskshape->setBottomBarColorProvider(this, 2);

    maskHCurveEditorG->curveListComplete();

    mask2CurveEditorG->setCurveListener(this);

    Lmaskshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskcurve.at(0)), defSpot.Lmaskcurve);
    Lmaskshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2CurveEditorG->curveListComplete();

    mask2CurveEditorGwav->setCurveListener(this);

    LLmaskcolshapewav->setIdentityValue(0.);
    LLmaskcolshapewav->setResetCurve(FlatCurveType(defSpot.LLmaskcolcurvewav.at(0)), defSpot.LLmaskcolcurvewav);
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
//    pack_start(*sensi);
    pack_start(*structcol);
    pack_start(*blurcolde);
    pack_start(*softradiuscol);
    pack_start(*invers);
    ToolParamBlock* const gradcolBox = Gtk::manage(new ToolParamBlock());
    gradcolBox->pack_start(*strcol);
    gradcolBox->pack_start(*strcolab);
    gradcolBox->pack_start(*strcolh);
    gradcolBox->pack_start(*angcol);
    expgradcol->add(*gradcolBox, false);
    pack_start(*expgradcol, false, false);
    ToolParamBlock* const curvBox = Gtk::manage(new ToolParamBlock());
    Gtk::HBox* const qualcurvbox = Gtk::manage(new Gtk::HBox());
    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);
    curvBox->pack_start(*qualcurvbox);
    curvBox->pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    curvBox->pack_start(*clCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    curvBox->pack_start(*HCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    curvBox->pack_start(*H2CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    curvBox->pack_start(*rgbCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    curvBox->pack_start(*special);
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
    merge1colFrame->add(*mergecolBox);
    mask7->pack_start(*merge1colFrame);
    mask7Box->pack_start(*mask7);
    expmaskcol1->add(*mask7Box, false);
    pack_start(*expmaskcol1, false, false);
    Gtk::Frame* const mergecolFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_MERGECOLFRA")));
    mergecolFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const maskcolBox = Gtk::manage(new ToolParamBlock());
    maskcolBox->pack_start(*showmaskcolMethod, Gtk::PACK_SHRINK, 4);
    maskcolBox->pack_start(*showmaskcolMethodinv, Gtk::PACK_SHRINK, 4);
    maskcolBox->pack_start(*enaColorMask, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*maskCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    ToolParamBlock* const strumBox = Gtk::manage(new ToolParamBlock());
    strumBox->pack_start(*strumaskcol);
    strumBox->pack_start(*toolcol);
    struFrame->add(*strumBox);
    maskcolBox->pack_start(*struFrame, Gtk::PACK_SHRINK, 0);
    ToolParamBlock* const blurmBox = Gtk::manage(new ToolParamBlock());
    blurmBox->pack_start(*fftColorMask, Gtk::PACK_SHRINK, 0);
    blurmBox->pack_start(*contcol);
    blurmBox->pack_start(*blurcol);
    blurFrame->add(*blurmBox);
    maskcolBox->pack_start(*blurFrame, Gtk::PACK_SHRINK, 0);
    maskcolBox->pack_start(*blendmaskcol, Gtk::PACK_SHRINK, 0);
    Gtk::Frame* const toolcolFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_TOOLMASK")));
    toolcolFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const toolcolBox = Gtk::manage(new ToolParamBlock());
    toolcolBox->pack_start(*radmaskcol, Gtk::PACK_SHRINK, 0);
    toolcolBox->pack_start(*lapmaskcol, Gtk::PACK_SHRINK, 0);
    toolcolBox->pack_start(*chromaskcol, Gtk::PACK_SHRINK, 0);
    toolcolBox->pack_start(*gammaskcol, Gtk::PACK_SHRINK, 0);
    toolcolBox->pack_start(*slomaskcol, Gtk::PACK_SHRINK, 0);
    toolcolBox->pack_start(*shadmaskcol, Gtk::PACK_SHRINK, 0);
    toolcolBox->pack_start(*maskHCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolcolBox->pack_start(*mask2CurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolcolBox->pack_start(*mask2CurveEditorGwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolcolBox->pack_start(*csThresholdcol, Gtk::PACK_SHRINK, 0);
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

bool LocallabColor::isMaskViewActive()
{
    return ((showmaskcolMethod->get_active_row_number() != 0) || (showmaskcolMethodinv->get_active_row_number() != 0));
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

void LocallabColor::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPCOLOR_TOOLTIP"));
        lightness->set_tooltip_text(M("TP_LOCALLAB_LIGHTN_TOOLTIP"));
        structcol->set_tooltip_text(M("TP_LOCALLAB_STRUCT_TOOLTIP"));
        sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        strcol->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        angcol->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
        qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));
        mercol->set_tooltip_text(M("TP_LOCALLAB_MERGEMER_TOOLTIP"));
        opacol->set_tooltip_text(M("TP_LOCALLAB_MERGEOPA_TOOLTIP"));
        conthrcol->set_tooltip_text(M("TP_LOCALLAB_MERGEOPA_TOOLTIP"));
        gridmerFrame->set_tooltip_text(M("TP_LOCALLAB_GRIDFRAME_TOOLTIP"));
        expmaskcol->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        radmaskcol->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmaskcol->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        Lmaskshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        LLmaskcolshapewav->setTooltip(M("TP_LOCALLAB_LMASK_LEVEL_TOOLTIP"));
        expmaskcol1->set_tooltip_text(M("TP_LOCALLAB_EXPMERGEFILE_TOOLTIP"));
        blendmaskcol->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        struFrame->set_tooltip_text(M("TP_LOCALLAB_STRUMASK_TOOLTIP"));
        blurFrame->set_tooltip_text(M("TP_LOCALLAB_BLURMASK_TOOLTIP"));
        maskHCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_HHMASK_TOOLTIP"));
        mask2CurveEditorGwav->set_tooltip_text(M("TP_LOCALLAB_WAVMASK_TOOLTIP"));
        mask2CurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        special->set_tooltip_text(M("TP_LOCALLAB_SPECIAL_TOOLTIP"));
        } else {
        exp->set_tooltip_text("");
        lightness->set_tooltip_text("");
        structcol->set_tooltip_text("");
        sensi->set_tooltip_text("");
        angcol->set_tooltip_text(M(""));
        strcol->set_tooltip_text("");
        qualitycurveMethod->set_tooltip_text("");
        mercol->set_tooltip_text("");
        opacol->set_tooltip_text("");
        conthrcol->set_tooltip_text("");
        gridmerFrame->set_tooltip_text("");
        expmaskcol->set_tooltip_text("");
        CCmaskshape->setTooltip("");
        LLmaskshape->setTooltip("");
        HHmaskshape->setTooltip("");
        radmaskcol->set_tooltip_text("");
        lapmaskcol->set_tooltip_text("");
        Lmaskshape->setTooltip("");
        LLmaskcolshapewav->setTooltip("");
        expmaskcol1->set_tooltip_text(M(""));
        blendmaskcol->set_tooltip_text(M(""));
        struFrame->set_tooltip_text(M(""));
        blurFrame->set_tooltip_text(M(""));
        maskHCurveEditorG->set_tooltip_text(M(""));
        mask2CurveEditorGwav->set_tooltip_text(M(""));
        mask2CurveEditorG->set_tooltip_text(M(""));
        special->set_tooltip_text(M(""));
    }
}

void LocallabColor::setDefaultExpanderVisibility()
{
    expgradcol->set_expanded(false);
    expcurvcol->set_expanded(false);
    expmaskcol1->set_expanded(false);
    expmaskcol->set_expanded(false);
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
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visicolor);
        exp->setEnabled(spot.expcolor);
        complexity->set_active(spot.complexcolor);

        curvactiv->set_active(spot.curvactiv);
        lightness->setValue(spot.lightness);
        contrast->setValue(spot.contrast);
        chroma->setValue(spot.chroma);
        labgrid->setParams(spot.labgridALow / LocallabParams::LABGRIDL_CORR_MAX,
                           spot.labgridBLow / LocallabParams::LABGRIDL_CORR_MAX,
                           spot.labgridAHigh / LocallabParams::LABGRIDL_CORR_MAX,
                           spot.labgridBHigh / LocallabParams::LABGRIDL_CORR_MAX,
                           false);

        if (spot.gridMethod == "one") {
            gridMethod->set_active(0);
        } else if (spot.gridMethod == "two") {
            gridMethod->set_active(1);
        }

        strengthgrid->setValue(spot.strengthgrid);
        sensi->setValue(spot.sensi);
        structcol->setValue(spot.structcol);
        blurcolde->setValue(spot.blurcolde);
        softradiuscol->setValue(spot.softradiuscol);
        invers->set_active(spot.invers);
        strcol->setValue(spot.strcol);
        strcolab->setValue(spot.strcolab);
        strcolh->setValue(spot.strcolh);
        angcol->setValue(spot.angcol);

        if (spot.qualitycurveMethod == "none") {
            qualitycurveMethod->set_active(0);
        } else if (spot.qualitycurveMethod == "std") {
            qualitycurveMethod->set_active(1);
        }

        llshape->setCurve(spot.llcurve);
        ccshape->setCurve(spot.cccurve);
        clshape->setCurve(spot.clcurve);
        lcshape->setCurve(spot.lccurve);
        LHshape->setCurve(spot.LHcurve);
        HHshape->setCurve(spot.HHcurve);

        if (spot.toneMethod == "one") {
            toneMethod->set_active(0);
        } else if (spot.toneMethod == "two") {
            toneMethod->set_active(1);
        } else if (spot.toneMethod == "thr") {
            toneMethod->set_active(2);
        } else if (spot.toneMethod == "fou") {
            toneMethod->set_active(3);
        }

        rgbshape->setCurve(spot.rgbcurve);
        special->set_active(spot.special);

        if (spot.merMethod == "mone") {
            merMethod->set_active(0);
//        } else if (spot.merMethod == "mtwo") {
//            merMethod->set_active(1);
        } else if (spot.merMethod == "mthr") {
            merMethod->set_active(1);
        } else if (spot.merMethod == "mfou") {
            merMethod->set_active(2);
        } else if (spot.merMethod == "mfiv") {
            merMethod->set_active(3);
        }

        if (spot.mergecolMethod == "one") {
            mergecolMethod->set_active(0);
        } else if (spot.mergecolMethod == "two") {
            mergecolMethod->set_active(1);
        } else if (spot.mergecolMethod == "thr") {
            mergecolMethod->set_active(2);
        } else if (spot.mergecolMethod == "fou") {
            mergecolMethod->set_active(3);
        } else if (spot.mergecolMethod == "fiv") {
            mergecolMethod->set_active(4);
        } else if (spot.mergecolMethod == "six") {
            mergecolMethod->set_active(5);
        } else if (spot.mergecolMethod == "sev") {
            mergecolMethod->set_active(6);
        } else if (spot.mergecolMethod == "sev0") {
            mergecolMethod->set_active(7);
        } else if (spot.mergecolMethod == "sev1") {
            mergecolMethod->set_active(8);
        } else if (spot.mergecolMethod == "sev2") {
            mergecolMethod->set_active(9);
        } else if (spot.mergecolMethod == "hei") {
            mergecolMethod->set_active(10);
        } else if (spot.mergecolMethod == "nin") {
            mergecolMethod->set_active(11);
        } else if (spot.mergecolMethod == "ten") {
            mergecolMethod->set_active(12);
        } else if (spot.mergecolMethod == "ele") {
            mergecolMethod->set_active(13);
        } else if (spot.mergecolMethod == "twe") {
            mergecolMethod->set_active(14);
        } else if (spot.mergecolMethod == "thi") {
            mergecolMethod->set_active(15);
        } else if (spot.mergecolMethod == "for") {
            mergecolMethod->set_active(16);
        } else if (spot.mergecolMethod == "hue") {
            mergecolMethod->set_active(17);
        } else if (spot.mergecolMethod == "sat") {
            mergecolMethod->set_active(18);
        } else if (spot.mergecolMethod == "col") {
            mergecolMethod->set_active(19);
        } else if (spot.mergecolMethod == "lum") {
            mergecolMethod->set_active(20);
        }

        mercol->setValue(spot.mercol);
        opacol->setValue(spot.opacol);
        conthrcol->setValue(spot.conthrcol);
        labgridmerg->setParams(0, 0,
                               spot.labgridAHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                               spot.labgridBHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                               false);
        merlucol->setValue(spot.merlucol);
        enaColorMask->set_active(spot.enaColorMask);
        CCmaskshape->setCurve(spot.CCmaskcurve);
        LLmaskshape->setCurve(spot.LLmaskcurve);
        HHmaskshape->setCurve(spot.HHmaskcurve);
        strumaskcol->setValue(spot.strumaskcol);
        toolcol->set_active(spot.toolcol);
        fftColorMask->set_active(spot.fftColorMask);
        contcol->setValue(spot.contcol);
        // Update GUI according to fftColorMash button state
        // Note: Contrary to the others, shall be called before setting 'blurcol' value
        updateColorGUI3();
        blurcol->setValue(spot.blurcol);
        blendmaskcol->setValue(spot.blendmaskcol);
        radmaskcol->setValue(spot.radmaskcol);
        lapmaskcol->setValue(spot.lapmaskcol);
        chromaskcol->setValue(spot.chromaskcol);
        gammaskcol->setValue(spot.gammaskcol);
        slomaskcol->setValue(spot.slomaskcol);
        shadmaskcol->setValue(spot.shadmaskcol);
        HHhmaskshape->setCurve(spot.HHhmaskcurve);
        Lmaskshape->setCurve(spot.Lmaskcurve);
        LLmaskcolshapewav->setCurve(spot.LLmaskcolcurvewav);
        csThresholdcol->setValue<int>(spot.csthresholdcol);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update GUI according to invers button state
    updateColorGUI1();

    // Update GUI according to merMethod combobox state
    updateColorGUI2();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expcolor = exp->getEnabled();
        spot.visicolor = exp->get_visible();
        spot.complexcolor = complexity->get_active_row_number();

        spot.curvactiv = curvactiv->get_active();
        spot.lightness = lightness->getIntValue();
        spot.contrast = contrast->getIntValue();
        spot.chroma = chroma->getIntValue();
        labgrid->getParams(spot.labgridALow,
                           spot.labgridBLow,
                           spot.labgridAHigh,
                           spot.labgridBHigh);
        spot.labgridALow *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.labgridAHigh *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.labgridBLow *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.labgridBHigh *= LocallabParams::LABGRIDL_CORR_MAX;

        if (gridMethod->get_active_row_number() == 0) {
            spot.gridMethod = "one";
        } else if (gridMethod->get_active_row_number() == 1) {
            spot.gridMethod = "two";
        }

        spot.strengthgrid = strengthgrid->getIntValue();
        spot.sensi = sensi->getIntValue();
        spot.structcol = structcol->getIntValue();
        spot.blurcolde = blurcolde->getIntValue();
        spot.softradiuscol = softradiuscol->getValue();
        spot.invers = invers->get_active();
        spot.strcol = strcol->getValue();
        spot.strcolab = strcolab->getValue();
        spot.strcolh = strcolh->getValue();
        spot.angcol = angcol->getValue();

        if (qualitycurveMethod->get_active_row_number() == 0) {
            spot.qualitycurveMethod = "none";
        } else if (qualitycurveMethod->get_active_row_number() == 1) {
            spot.qualitycurveMethod = "std";
        }

        spot.llcurve = llshape->getCurve();
        spot.cccurve = ccshape->getCurve();
        spot.clcurve = clshape->getCurve();
        spot.lccurve = lcshape->getCurve();
        spot.LHcurve = LHshape->getCurve();
        spot.HHcurve = HHshape->getCurve();

        if (toneMethod->get_active_row_number() == 0) {
            spot.toneMethod = "one";
        } else if (toneMethod->get_active_row_number() == 1) {
            spot.toneMethod = "two";
        } else if (toneMethod->get_active_row_number() == 2) {
            spot.toneMethod = "thr";
        } else if (toneMethod->get_active_row_number() == 3) {
            spot.toneMethod = "fou";
        }

        spot.rgbcurve = rgbshape->getCurve();
        spot.special = special->get_active();

        if (merMethod->get_active_row_number() == 0) {
            spot.merMethod = "mone";
//        } else if (merMethod->get_active_row_number() == 1) {
//            spot.merMethod = "mtwo";
        } else if (merMethod->get_active_row_number() == 1) {
            spot.merMethod = "mthr";
        } else if (merMethod->get_active_row_number() == 2) {
            spot.merMethod = "mfou";
        } else if (merMethod->get_active_row_number() == 3) {
            spot.merMethod = "mfiv";
        }

        if (mergecolMethod->get_active_row_number() == 0) {
            spot.mergecolMethod = "one";
        } else if (mergecolMethod->get_active_row_number() == 1) {
            spot.mergecolMethod = "two";
        } else if (mergecolMethod->get_active_row_number() == 2) {
            spot.mergecolMethod = "thr";
        } else if (mergecolMethod->get_active_row_number() == 3) {
            spot.mergecolMethod = "fou";
        } else if (mergecolMethod->get_active_row_number() == 4) {
            spot.mergecolMethod = "fiv";
        } else if (mergecolMethod->get_active_row_number() == 5) {
            spot.mergecolMethod = "six";
        } else if (mergecolMethod->get_active_row_number() == 6) {
            spot.mergecolMethod = "sev";
        } else if (mergecolMethod->get_active_row_number() == 7) {
            spot.mergecolMethod = "sev0";
        } else if (mergecolMethod->get_active_row_number() == 8) {
            spot.mergecolMethod = "sev1";
        } else if (mergecolMethod->get_active_row_number() == 9) {
            spot.mergecolMethod = "sev2";
        } else if (mergecolMethod->get_active_row_number() == 10) {
            spot.mergecolMethod = "hei";
        } else if (mergecolMethod->get_active_row_number() == 11) {
            spot.mergecolMethod = "nin";
        } else if (mergecolMethod->get_active_row_number() == 12) {
            spot.mergecolMethod = "ten";
        } else if (mergecolMethod->get_active_row_number() == 13) {
            spot.mergecolMethod = "ele";
        } else if (mergecolMethod->get_active_row_number() == 14) {
            spot.mergecolMethod = "twe";
        } else if (mergecolMethod->get_active_row_number() == 15) {
            spot.mergecolMethod = "thi";
        } else if (mergecolMethod->get_active_row_number() == 16) {
            spot.mergecolMethod = "for";
        } else if (mergecolMethod->get_active_row_number() == 17) {
            spot.mergecolMethod = "hue";
        } else if (mergecolMethod->get_active_row_number() == 18) {
            spot.mergecolMethod = "sat";
        } else if (mergecolMethod->get_active_row_number() == 19) {
            spot.mergecolMethod = "col";
        } else if (mergecolMethod->get_active_row_number() == 20) {
            spot.mergecolMethod = "lum";
        }

        spot.mercol = mercol->getValue();
        spot.opacol = opacol->getValue();
        spot.conthrcol = conthrcol->getValue();
        labgridmerg->getParams(spot.labgridALowmerg,
                               spot.labgridBLowmerg,
                               spot.labgridAHighmerg,
                               spot.labgridBHighmerg);
        spot.labgridALowmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.labgridAHighmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.labgridBLowmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.labgridBHighmerg *= LocallabParams::LABGRIDL_CORR_MAX;
        spot.merlucol = merlucol->getValue();
        spot.enaColorMask = enaColorMask->get_active();
        spot.CCmaskcurve = CCmaskshape->getCurve();
        spot.LLmaskcurve = LLmaskshape->getCurve();
        spot.HHmaskcurve = HHmaskshape->getCurve();
        spot.strumaskcol = strumaskcol->getValue();
        spot.toolcol = toolcol->get_active();
        spot.fftColorMask = fftColorMask->get_active();
        spot.contcol = contcol->getValue();
        spot.blurcol = blurcol->getValue();
        spot.blendmaskcol = blendmaskcol->getIntValue();
        spot.radmaskcol = radmaskcol->getValue();
        spot.lapmaskcol = lapmaskcol->getValue();
        spot.chromaskcol = chromaskcol->getValue();
        spot.gammaskcol = gammaskcol->getValue();
        spot.slomaskcol = slomaskcol->getValue();
        spot.shadmaskcol = shadmaskcol->getIntValue();
        spot.HHhmaskcurve = HHhmaskshape->getCurve();
        spot.Lmaskcurve = Lmaskshape->getCurve();
        spot.LLmaskcolcurvewav = LLmaskcolshapewav->getCurve();
        spot.csthresholdcol = csThresholdcol->getValue<int>();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

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

void LocallabColor::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    blurcolde->setValue((double)defSpot.blurcolde);
  //  softradiuscol->setValue(defSpot.softradiuscol);
    strcolab->setValue(defSpot.strcolab);
    strcolh->setValue(defSpot.strcolh);
    clshape->setCurve(defSpot.clcurve);
    lcshape->setCurve(defSpot.lccurve);
    LHshape->setCurve(defSpot.LHcurve);
    HHshape->setCurve(defSpot.HHcurve);

    if (defSpot.toneMethod == "one") {
        toneMethod->set_active(0);
    } else if (defSpot.toneMethod == "two") {
        toneMethod->set_active(1);
    } else if (defSpot.toneMethod == "thr") {
        toneMethod->set_active(2);
    } else if (defSpot.toneMethod == "fou") {
        toneMethod->set_active(3);
    }

    rgbshape->setCurve(defSpot.rgbcurve);
    special->set_active(defSpot.special);

    if (defSpot.merMethod == "mone") {
        merMethod->set_active(0);
//    } else if (defSpot.merMethod == "mtwo") {
//        merMethod->set_active(1);
    } else if (defSpot.merMethod == "mthr") {
        merMethod->set_active(1);
    } else if (defSpot.merMethod == "mfou") {
        merMethod->set_active(2);
    } else if (defSpot.merMethod == "mfiv") {
        merMethod->set_active(3);
    }

    if (defSpot.mergecolMethod == "one") {
        mergecolMethod->set_active(0);
    } else if (defSpot.mergecolMethod == "two") {
        mergecolMethod->set_active(1);
    } else if (defSpot.mergecolMethod == "thr") {
        mergecolMethod->set_active(2);
    } else if (defSpot.mergecolMethod == "fou") {
        mergecolMethod->set_active(3);
    } else if (defSpot.mergecolMethod == "fiv") {
        mergecolMethod->set_active(4);
    } else if (defSpot.mergecolMethod == "six") {
        mergecolMethod->set_active(5);
    } else if (defSpot.mergecolMethod == "sev") {
        mergecolMethod->set_active(6);
    } else if (defSpot.mergecolMethod == "sev0") {
        mergecolMethod->set_active(7);
    } else if (defSpot.mergecolMethod == "sev1") {
        mergecolMethod->set_active(8);
    } else if (defSpot.mergecolMethod == "sev2") {
        mergecolMethod->set_active(9);
    } else if (defSpot.mergecolMethod == "hei") {
        mergecolMethod->set_active(10);
    } else if (defSpot.mergecolMethod == "nin") {
        mergecolMethod->set_active(11);
    } else if (defSpot.mergecolMethod == "ten") {
        mergecolMethod->set_active(12);
    } else if (defSpot.mergecolMethod == "ele") {
        mergecolMethod->set_active(13);
    } else if (defSpot.mergecolMethod == "twe") {
        mergecolMethod->set_active(14);
    } else if (defSpot.mergecolMethod == "thi") {
        mergecolMethod->set_active(15);
    } else if (defSpot.mergecolMethod == "for") {
        mergecolMethod->set_active(16);
    } else if (defSpot.mergecolMethod == "hue") {
        mergecolMethod->set_active(17);
    } else if (defSpot.mergecolMethod == "sat") {
        mergecolMethod->set_active(18);
    } else if (defSpot.mergecolMethod == "col") {
        mergecolMethod->set_active(19);
    } else if (defSpot.mergecolMethod == "lum") {
        mergecolMethod->set_active(20);
    }

    mercol->setValue(defSpot.mercol);
    opacol->setValue(defSpot.opacol);
    conthrcol->setValue(defSpot.conthrcol);
    labgridmerg->setParams(0, 0,
                           defSpot.labgridAHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                           defSpot.labgridBHighmerg / LocallabParams::LABGRIDL_CORR_MAX,
                           false);
    merlucol->setValue(defSpot.merlucol);
    strumaskcol->setValue(defSpot.strumaskcol);
    toolcol->set_active(defSpot.toolcol);
    fftColorMask->set_active(defSpot.fftColorMask);
    contcol->setValue(defSpot.contcol);
    lapmaskcol->setValue(defSpot.lapmaskcol);
    gammaskcol->setValue(defSpot.gammaskcol);
    slomaskcol->setValue(defSpot.slomaskcol);
    shadmaskcol->setValue((double)defSpot.shadmaskcol);
    HHhmaskshape->setCurve(defSpot.HHhmaskcurve);
    LLmaskcolshapewav->setCurve(defSpot.LLmaskcolcurvewav);
    csThresholdcol->setValue<int>(defSpot.csthresholdcol);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update GUI according to merMethod combobox state
    updateColorGUI2();
    // - Update GUI according to fftColorMash button state
    updateColorGUI3();
}

void LocallabColor::updateGUIToMode(const modeType new_type)
{
    if (new_type == Normal) {
        // Advanced widgets are hidden in Normal mode
        blurcolde->hide();
        softradiuscol->show();
        strcolab->hide();
        strcolh->hide();
        clCurveEditorG->hide();
        HCurveEditorG->hide();
        H2CurveEditorG->hide();
        rgbCurveEditorG->hide();
        special->hide();
        expmaskcol1->hide();
        struFrame->hide();
        blurFrame->hide();
        lapmaskcol->hide();
        gammaskcol->hide();
        slomaskcol->hide();
        shadmaskcol->hide();
        maskHCurveEditorG->hide();
        mask2CurveEditorGwav->hide();
        csThresholdcol->hide();
    } else {
        // Advanced widgets are shown in Expert mode
        blurcolde->show();

        if (!invers->get_active()) { // Keep widget hidden when invers is toggled
            softradiuscol->show();
        }

        strcolab->show();
        strcolh->show();

        if (!invers->get_active()) { // Keep widgets hidden when invers is toggled
            clCurveEditorG->show();
            HCurveEditorG->show();
        }

        H2CurveEditorG->show();
        rgbCurveEditorG->show();
        special->show();

        if (!invers->get_active()) { // Keep widget hidden when invers is toggled
            expmaskcol1->show();
        }

        struFrame->show();
        blurFrame->show();
        lapmaskcol->show();
        gammaskcol->show();
        slomaskcol->show();
        shadmaskcol->show();
        maskHCurveEditorG->show();
        mask2CurveEditorGwav->show();
        csThresholdcol->show();
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

    // This event is called to transmit potentially reset mask state
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
    const int mode = complexity->get_active_row_number();

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

        if (mode == Normal) { // Keep widget hidden in Normal mode
            softradiuscol->show();
        }

        expgradcol->show();
        labqualcurv->show();
        qualitycurveMethod->show();

        if (mode == Normal) { // Keep widgets hidden in Normal mode
            clCurveEditorG->show();
            HCurveEditorG->show();
            expmaskcol1->show();
        }

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
/*
        case 1:
            invers->set_sensitive(false);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(true);
            special->set_sensitive(true);
            mask7->hide();
            conthrcol->hide();
            gridmerFrame->hide();
            break;
*/
        case 1:
            invers->set_sensitive(false);
            H2CurveEditorG->set_sensitive(true);
            rgbCurveEditorG->set_sensitive(false);
            special->set_sensitive(false);
            mask7->show();
            conthrcol->show();
            gridmerFrame->hide();
            break;

        case 2:
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
    linear(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LINEAR"), 0., 1., 0.01, 0.05))),
    balanexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BALANEXP"), 0.5, 1.5, 0.01, 1.0))),
    gamm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMM"), 0.2, 1.3, 0.01, 0.4))),
    exnoiseMethod(Gtk::manage(new MyComboBoxText())),
    fatFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_FATFRA")))),
    fatamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATAMOUNT"), 1., 100., 1., 1.))),
    fatdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATDETAIL"), -100., 300., 1., 0.))),
    fatlevel(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATLEVEL"), 0.25, 2.5, 0.05, 1.))),
    fatanchor(Gtk::manage(new Adjuster(M("TP_LOCALLAB_FATANCHORA"), 0.1, 3.0, 0.05, 1.))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 60))),
    structexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    blurexpde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    exptoolexp(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPTOOL")))),
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), MINEXP, MAXEXP, 0.02, 0.))), 
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 10, 0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0))),
    shadex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHADEX"), 0, 100, 1, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHADEXCOMP"), 0, 100, 1, 50))),
    expchroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_EXPCHROMA"), -50, 100, 1, 5))),
    curveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"))),
    shapeexpos(static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""))),
    expgradexp(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_EXPGRAD")))),
    strexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADSTR"), -4., 4., 0.05, 0.))),
    angexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GRADANG"), -180, 180, 0.1, 0.))),
    softradiusexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 1000.0, 0.5, 0.))),
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

    // Parameter Exposure specific widgets
    expMethod->append(M("TP_LOCALLAB_STD"));
    expMethod->append(M("TP_LOCALLAB_PDE"));
    expMethod->set_active(0);
    expMethodConn = expMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::expMethodChanged));

    pdeFrame->set_label_align(0.025, 0.5);

    laplacexp->setAdjusterListener(this);

    linear->setAdjusterListener(this);

    balanexp->setAdjusterListener(this);

    gamm->setAdjusterListener(this);

    exnoiseMethod->append(M("TP_LOCALLAB_NONENOISE"));
    exnoiseMethod->append(M("TP_LOCALLAB_MEDIAN"));
    exnoiseMethod->append(M("TP_LOCALLAB_WEDIANHI"));
    exnoiseMethod->set_active(0);
    exnoiseMethodConn  = exnoiseMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::exnoiseMethodChanged));

    fatFrame->set_label_align(0.025, 0.5);

    fatamount->setAdjusterListener(this);

    fatdetail->setAdjusterListener(this);

    fatlevel->setAdjusterListener(this);

    fatanchor->setAdjusterListener(this);

    sensiex->setAdjusterListener(this);

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
    shapeexpos->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});
    shapeexpos->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});

    curveEditorG->curveListComplete();

    setExpandAlignProperties(expgradexp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    strexp->setAdjusterListener(this);

    angexp->setAdjusterListener(this);
    angexp->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));

    softradiusexp->setLogScale(10, 0);
    softradiusexp->setAdjusterListener(this);

    inversexConn  = inversex->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::inversexChanged));
    inversex->set_tooltip_text(M("TP_LOCALLAB_INVERS_TOOLTIP"));

    setExpandAlignProperties(expmaskexp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWSTRUCEX"));
    showmaskexpMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskexpMethod->set_active(0);
    showmaskexpMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskexpMethodConn  = showmaskexpMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::showmaskexpMethodChanged));

    showmaskexpMethodinv->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskexpMethodinv->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskexpMethodinv->set_active(0);
    showmaskexpMethodinv->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskexpMethodConninv  = showmaskexpMethodinv->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::showmaskexpMethodChangedinv));

    enaExpMaskConn = enaExpMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::enaExpMaskChanged));

    enaExpMaskaftConn = enaExpMaskaft->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::enaExpMaskaftChanged));

    maskexpCurveEditorG->setCurveListener(this);

    CCmaskexpshape->setIdentityValue(0.);
    CCmaskexpshape->setResetCurve(FlatCurveType(defSpot.CCmaskexpcurve.at(0)), defSpot.CCmaskexpcurve);
    CCmaskexpshape->setBottomBarColorProvider(this, 1);

    LLmaskexpshape->setIdentityValue(0.);
    LLmaskexpshape->setResetCurve(FlatCurveType(defSpot.LLmaskexpcurve.at(0)), defSpot.LLmaskexpcurve);
    LLmaskexpshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});

    HHmaskexpshape->setIdentityValue(0.);
    HHmaskexpshape->setResetCurve(FlatCurveType(defSpot.HHmaskexpcurve.at(0)), defSpot.HHmaskexpcurve);
    HHmaskexpshape->setCurveColorProvider(this, 2);
    HHmaskexpshape->setBottomBarColorProvider(this, 2);

    maskexpCurveEditorG->curveListComplete();

    blendmaskexp->setAdjusterListener(this);

    radmaskexp->setLogScale(10, -10);
    radmaskexp->setAdjusterListener(this);

    lapmaskexp->setAdjusterListener(this);

    chromaskexp->setAdjusterListener(this);

    gammaskexp->setAdjusterListener(this);

    slomaskexp->setAdjusterListener(this);

    gradFramemask->set_label_align(0.025, 0.5);

    strmaskexp->setAdjusterListener(this);

    angmaskexp->setAdjusterListener(this);
    angmaskexp->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));

    mask2expCurveEditorG->setCurveListener(this);

    Lmaskexpshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskexpcurve.at(0)), defSpot.Lmaskexpcurve);
    Lmaskexpshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});
    Lmaskexpshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1}});

    mask2expCurveEditorG->curveListComplete();

    // Add Color & Light specific widgets to GUI
    pack_start(*expMethod);
    ToolParamBlock* const pdeBox = Gtk::manage(new ToolParamBlock());
    pdeBox->pack_start(*laplacexp);
    pdeBox->pack_start(*linear);
    pdeBox->pack_start(*balanexp);
    pdeBox->pack_start(*gamm);
    Gtk::HBox* const ctboxexpmethod = Gtk::manage(new Gtk::HBox());
    Gtk::Label* const labelexpmethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_NOISEMETH") + ":"));
    ctboxexpmethod->pack_start(*labelexpmethod, Gtk::PACK_SHRINK, 4);
    ctboxexpmethod->pack_start(*exnoiseMethod);
    pdeBox->pack_start(*ctboxexpmethod);
    pdeFrame->add(*pdeBox);
    pack_start(*pdeFrame);
    ToolParamBlock* const fatBox = Gtk::manage(new ToolParamBlock());
    fatBox->pack_start(*fatamount);
    fatBox->pack_start(*fatdetail);
    fatBox->pack_start(*fatlevel);
    fatBox->pack_start(*fatanchor);
    fatFrame->add(*fatBox);
    pack_start(*fatFrame);
    pack_start(*expcomp);
    pack_start(*sensiex);
    pack_start(*structexp);
    pack_start(*blurexpde);
    ToolParamBlock* const toolBox = Gtk::manage(new ToolParamBlock());
//    toolBox->pack_start(*expcomp);
    toolBox->pack_start(*black);
    toolBox->pack_start(*hlcompr);
    toolBox->pack_start(*hlcomprthresh);
    toolBox->pack_start(*shadex);
    toolBox->pack_start(*shcompr);
    toolBox->pack_start(*expchroma);
    toolBox->pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    exptoolexp->add(*toolBox, false);
    pack_start(*exptoolexp);
    ToolParamBlock* const gradBox = Gtk::manage(new ToolParamBlock());
    gradBox->pack_start(*strexp);
    gradBox->pack_start(*angexp);
    expgradexp->add(*gradBox, false);
    pack_start(*expgradexp);
    pack_start(*softradiusexp);
    pack_start(*inversex);
    ToolParamBlock* const maskexpBox = Gtk::manage(new ToolParamBlock());
    maskexpBox->pack_start(*showmaskexpMethod, Gtk::PACK_SHRINK, 4);
    maskexpBox->pack_start(*showmaskexpMethodinv, Gtk::PACK_SHRINK, 4);
    maskexpBox->pack_start(*enaExpMask, Gtk::PACK_SHRINK, 0);
    // maskexpBox->pack_start(*enaExpMaskaft, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*maskexpCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskexpBox->pack_start(*blendmaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*radmaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*lapmaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*chromaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*gammaskexp, Gtk::PACK_SHRINK, 0);
    maskexpBox->pack_start(*slomaskexp, Gtk::PACK_SHRINK, 0);
    ToolParamBlock* const gradmaskBox = Gtk::manage(new ToolParamBlock());
    gradmaskBox->pack_start(*strmaskexp);
    gradmaskBox->pack_start(*angmaskexp);
    gradFramemask->add(*gradmaskBox);
    maskexpBox->pack_start(*gradFramemask);
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

bool LocallabExposure::isMaskViewActive()
{
    return ((showmaskexpMethod->get_active_row_number() != 0) || (showmaskexpMethodinv->get_active_row_number() != 0));
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

void LocallabExposure::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPOSURE_TOOLTIP"));
        expMethod->set_tooltip_text(M("TP_LOCALLAB_EXPMETHOD_TOOLTIP"));
        structexp->set_tooltip_text(M("TP_LOCALLAB_STRUCT_TOOLTIP"));
        pdeFrame->set_tooltip_text(M("TP_LOCALLAB_PDEFRAME_TOOLTIP"));
        laplacexp->set_tooltip_text(M("TP_LOCALLAB_EXPLAP_TOOLTIP"));
        balanexp->set_tooltip_text(M("TP_LOCALLAB_EXPLAPBAL_TOOLTIP"));
        gamm->set_tooltip_text(M("TP_LOCALLAB_EXPLAPGAMM_TOOLTIP"));
        linear->set_tooltip_text(M("TP_LOCALLAB_EXPLAPLIN_TOOLTIP"));
        exnoiseMethod->set_tooltip_text(M("TP_LOCALLAB_EXPNOISEMETHOD_TOOLTIP"));
        fatFrame->set_tooltip_text(M("TP_LOCALLAB_FATFRAME_TOOLTIP"));
        sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
        shapeexpos->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
        strexp->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        expmaskexp->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskexpshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        radmaskexp->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmaskexp->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        strmaskexp->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        Lmaskexpshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        blendmaskexp->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        mask2expCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        expchroma->set_tooltip_text(M("TP_LOCALLAB_EXPCHROMA_TOOLTIP"));
    } else {
        exp->set_tooltip_text("");
        expMethod->set_tooltip_text("");
        structexp->set_tooltip_text("");
        pdeFrame->set_tooltip_text("");
        exnoiseMethod->set_tooltip_text("");
        laplacexp->set_tooltip_text(M(""));
        balanexp->set_tooltip_text(M(""));
        gamm->set_tooltip_text(M(""));
        linear->set_tooltip_text(M(""));
        fatFrame->set_tooltip_text("");
        sensiex->set_tooltip_text("");
        shapeexpos->setTooltip("");
        strexp->set_tooltip_text("");
        expmaskexp->set_tooltip_text("");
        CCmaskexpshape->setTooltip("");
        LLmaskexpshape->setTooltip("");
        HHmaskexpshape->setTooltip("");
        radmaskexp->set_tooltip_text("");
        lapmaskexp->set_tooltip_text("");
        strmaskexp->set_tooltip_text("");
        Lmaskexpshape->setTooltip("");
        blendmaskexp->set_tooltip_text(M(""));
        mask2expCurveEditorG->set_tooltip_text(M(""));
        expchroma->set_tooltip_text(M(""));
    }
}

void LocallabExposure::setDefaultExpanderVisibility()
{
    exptoolexp->set_expanded(false);
    expgradexp->set_expanded(false);
    expmaskexp->set_expanded(false);
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
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visiexpose);
        exp->setEnabled(spot.expexpose);
        complexity->set_active(spot.complexexpose);

        if (spot.expMethod == "std") {
            expMethod->set_active(0);
        } else if (spot.expMethod == "pde") {
            expMethod->set_active(1);
        }

        laplacexp->setValue(spot.laplacexp);
        linear->setValue(spot.linear);
        balanexp->setValue(spot.balanexp);
        gamm->setValue(spot.gamm);

        if (spot.exnoiseMethod == "one") {
            exnoiseMethod->set_active(0);
        } else if (spot.exnoiseMethod == "med") {
            exnoiseMethod->set_active(1);
        } else if (spot.exnoiseMethod == "medhi") {
            exnoiseMethod->set_active(2);
        }

        fatamount->setValue(spot.fatamount);
        fatdetail->setValue(spot.fatdetail);
        fatlevel->setValue(spot.fatlevel);
        fatanchor->setValue(spot.fatanchor);
        sensiex->setValue(spot.sensiex);
        structexp->setValue(spot.structexp);
        blurexpde->setValue(spot.blurexpde);
        expcomp->setValue(spot.expcomp);
        black->setValue(spot.black);
        hlcompr->setValue(spot.hlcompr);
        hlcomprthresh->setValue(spot.hlcomprthresh);
        shadex->setValue(spot.shadex);
        shcompr->setValue(spot.shcompr);
        expchroma->setValue(spot.expchroma);
        shapeexpos->setCurve(spot.excurve);
        strexp->setValue(spot.strexp);
        angexp->setValue(spot.angexp);
        softradiusexp->setValue(spot.softradiusexp);
        inversex->set_active(spot.inversex);
        enaExpMask->set_active(spot.enaExpMask);
        enaExpMaskaft->set_active(spot.enaExpMaskaft);
        CCmaskexpshape->setCurve(spot.CCmaskexpcurve);
        LLmaskexpshape->setCurve(spot.LLmaskexpcurve);
        HHmaskexpshape->setCurve(spot.HHmaskexpcurve);
        blendmaskexp->setValue(spot.blendmaskexp);
        radmaskexp->setValue(spot.radmaskexp);
        lapmaskexp->setValue(spot.lapmaskexp);
        chromaskexp->setValue(spot.chromaskexp);
        gammaskexp->setValue(spot.gammaskexp);
        slomaskexp->setValue(spot.slomaskexp);
        strmaskexp->setValue(spot.strmaskexp);
        angmaskexp->setValue(spot.angmaskexp);
        Lmaskexpshape->setCurve(spot.Lmaskexpcurve);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

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
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expexpose = exp->getEnabled();
        spot.visiexpose = exp->get_visible();
        spot.complexexpose = complexity->get_active_row_number();

        if (expMethod->get_active_row_number() == 0) {
            spot.expMethod = "std";
        } else if (expMethod->get_active_row_number() == 1) {
            spot.expMethod = "pde";
        }

        spot.laplacexp = laplacexp->getValue();
        spot.linear = linear->getValue();
        spot.balanexp = balanexp->getValue();
        spot.gamm = gamm->getValue();

        if (exnoiseMethod->get_active_row_number() == 0) {
            spot.exnoiseMethod = "none";
        } else if (exnoiseMethod->get_active_row_number() == 1) {
            spot.exnoiseMethod = "med";
        } else if (exnoiseMethod->get_active_row_number() == 2) {
            spot.exnoiseMethod = "medhi";
        }

        spot.fatamount = fatamount->getValue();
        spot.fatdetail = fatdetail->getValue();
        spot.fatlevel = fatlevel->getValue();
        spot.fatanchor = fatanchor->getValue();
        spot.sensiex = sensiex->getIntValue();
        spot.structexp = structexp->getIntValue();
        spot.blurexpde = blurexpde->getIntValue();
        spot.expcomp = expcomp->getValue();
        spot.black = black->getIntValue();
        spot.hlcompr = hlcompr->getIntValue();
        spot.hlcomprthresh = hlcomprthresh->getIntValue();
        spot.shadex = shadex->getIntValue();
        spot.shcompr = shcompr->getIntValue();
        spot.expchroma = expchroma->getIntValue();
        spot.excurve = shapeexpos->getCurve();
        spot.strexp = strexp->getValue();
        spot.angexp = angexp->getValue();
        spot.softradiusexp = softradiusexp->getValue();
        spot.inversex = inversex->get_active();
        spot.enaExpMask = enaExpMask->get_active();
        spot.enaExpMaskaft = enaExpMaskaft->get_active();
        spot.CCmaskexpcurve = CCmaskexpshape->getCurve();
        spot.LLmaskexpcurve = LLmaskexpshape->getCurve();
        spot.HHmaskexpcurve = HHmaskexpshape->getCurve();
        spot.blendmaskexp = blendmaskexp->getIntValue();
        spot.radmaskexp = radmaskexp->getValue();
        spot.lapmaskexp = lapmaskexp->getValue();
        spot.chromaskexp = chromaskexp->getValue();
        spot.gammaskexp = gammaskexp->getValue();
        spot.slomaskexp = slomaskexp->getValue();
        spot.strmaskexp = strmaskexp->getValue();
        spot.angmaskexp = angmaskexp->getValue();
        spot.Lmaskexpcurve = Lmaskexpshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabExposure::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

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

void LocallabExposure::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    blurexpde->setValue((double)defSpot.blurexpde);
    lapmaskexp->setValue(defSpot.lapmaskexp);
    gammaskexp->setValue(defSpot.gammaskexp);
    slomaskexp->setValue(defSpot.slomaskexp);
    strmaskexp->setValue(defSpot.strmaskexp);
    angmaskexp->setValue(defSpot.angmaskexp);
 //   laplacexp->setValue(defSpot.laplacexp);
  //  linear->setValue(defSpot.linear);
  //  balanexp->setValue(defSpot.balanexp);
  //  gamm->setValue(defSpot.gamm);

    // Enable all listeners
    enableListener();
}

void LocallabExposure::updateGUIToMode(const modeType new_type)
{
    if (new_type == Normal) {
        // Advanced widgets are hidden in Normal mode
        lapmaskexp->hide();
        gammaskexp->hide();
        slomaskexp->hide();
        gradFramemask->hide();
        blurexpde->hide();
     //   pdeFrame->hide();
    } else {
        // Advanced widgets are shown in Expert mode
        lapmaskexp->show();
        gammaskexp->show();
        slomaskexp->show();
        gradFramemask->show();
        blurexpde->show();
      //  pdeFrame->show();
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

    // This event is called to transmit potentially reset mask state
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
        pdeFrame->hide();
        fatFrame->hide();
        softradiusexp->set_sensitive(true);
        sensiex->set_sensitive(false);
    } else if (expMethod->get_active_row_number() == 1) {
        pdeFrame->show();
        fatFrame->show();
        softradiusexp->set_sensitive(false);
        sensiex->set_sensitive(true);
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
    multipliersh([]() -> std::array<Adjuster *, 5>
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

    // Parameter Shadow highlight specific widgets
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

    gamSH->setAdjusterListener(this);

    sloSH->setAdjusterListener(this);

    setExpandAlignProperties(expgradsh, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    strSH->setAdjusterListener(this);

    angSH->setAdjusterListener(this);
    angSH->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));

    inversshConn = inverssh->signal_toggled().connect(sigc::mem_fun(*this, &LocallabShadow::inversshChanged));
    inverssh->set_tooltip_text(M("TP_LOCALLAB_INVERS_TOOLTIP"));

    setExpandAlignProperties(expmasksh, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskSHMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskSHMethod->set_active(0);
    showmaskSHMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskSHMethodConn = showmaskSHMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabShadow::showmaskSHMethodChanged));

    showmaskSHMethodinv->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskSHMethodinv->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskSHMethodinv->set_active(0);
    showmaskSHMethodinv->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskSHMethodConninv = showmaskSHMethodinv->signal_changed().connect(sigc::mem_fun(*this, &LocallabShadow::showmaskSHMethodChangedinv));

    enaSHMaskConn = enaSHMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabShadow::enaSHMaskChanged));

    maskSHCurveEditorG->setCurveListener(this);

    CCmaskSHshape->setIdentityValue(0.);
    CCmaskSHshape->setResetCurve(FlatCurveType(defSpot.CCmaskSHcurve.at(0)), defSpot.CCmaskSHcurve);
    CCmaskSHshape->setBottomBarColorProvider(this, 1);

    LLmaskSHshape->setIdentityValue(0.);
    LLmaskSHshape->setResetCurve(FlatCurveType(defSpot.LLmaskSHcurve.at(0)), defSpot.LLmaskSHcurve);
    LLmaskSHshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskSHshape->setIdentityValue(0.);
    HHmaskSHshape->setResetCurve(FlatCurveType(defSpot.HHmaskSHcurve.at(0)), defSpot.HHmaskSHcurve);
    HHmaskSHshape->setCurveColorProvider(this, 2);
    HHmaskSHshape->setBottomBarColorProvider(this, 2);

    maskSHCurveEditorG->curveListComplete();

    blendmaskSH->setAdjusterListener(this);

    radmaskSH->setLogScale(10, -10);
    radmaskSH->setAdjusterListener(this);

    lapmaskSH->setAdjusterListener(this);

    chromaskSH->setAdjusterListener(this);

    gammaskSH->setAdjusterListener(this);

    slomaskSH->setAdjusterListener(this);

    mask2SHCurveEditorG->setCurveListener(this);

    LmaskSHshape->setResetCurve(DiagonalCurveType(defSpot.LmaskSHcurve.at(0)), defSpot.LmaskSHcurve);
    LmaskSHshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    LmaskSHshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2SHCurveEditorG->curveListComplete();

    fatSHFrame->set_label_align(0.025, 0.5);

    fatamountSH->setAdjusterListener(this);

    fatanchorSH->setAdjusterListener(this);

    // Add Shadow highlight specific widgets to GUI
    pack_start(*shMethod);

    for (int i = 0; i < 5; ++i) {
        pack_start(*multipliersh[i]);
    }

    pack_start(*detailSH);
    pack_start(*highlights);
    pack_start(*h_tonalwidth);
    pack_start(*shadows);
    pack_start(*s_tonalwidth);
    pack_start(*sh_radius);
//    pack_start(*sensihs);
    pack_start(*blurSHde);
//    Gtk::Frame* const gamFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_GAMFRA")));
    gamFrame->set_label_align(0.025, 0.5);
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
    maskSHBox->pack_start(*lapmaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*chromaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*gammaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*slomaskSH, Gtk::PACK_SHRINK, 0);
    maskSHBox->pack_start(*mask2SHCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    ToolParamBlock* const fatSHBox = Gtk::manage(new ToolParamBlock());
    fatSHBox->pack_start(*fatamountSH);
    fatSHBox->pack_start(*fatanchorSH);
    fatSHFrame->add(*fatSHBox);
//    maskSHBox->pack_start(*fatSHFrame);
    expmasksh->add(*maskSHBox, false);
    pack_start(*expmasksh, false, false);
}

LocallabShadow::~LocallabShadow()
{
    delete maskSHCurveEditorG;
    delete mask2SHCurveEditorG;
}

bool LocallabShadow::isMaskViewActive()
{
    return ((showmaskSHMethod->get_active_row_number() != 0) || (showmaskSHMethodinv->get_active_row_number() != 0));
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

void LocallabShadow::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        exp->set_tooltip_text(M("TP_LOCALLAB_SHADOWHIGHLIGHT_TOOLTIP"));
        strSH->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        expmasksh->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskSHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        radmaskSH->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmaskSH->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        LmaskSHshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        for (unsigned int i = 0; i < multipliersh.size(); i++) {
            multipliersh[i]->set_tooltip_text(M("TP_LOCALLAB_MULTIPL_TOOLTIP"));
        }
        gamSH->set_tooltip_text(M("TP_LOCALLAB_SHTRC_TOOLTIP"));
        sloSH->set_tooltip_text(M("TP_LOCALLAB_SHTRC_TOOLTIP"));
        blendmaskSH->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        mask2SHCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
    } else {
        exp->set_tooltip_text("");
        strSH->set_tooltip_text("");
        expmasksh->set_tooltip_text("");
        CCmaskSHshape->setTooltip("");
        LLmaskSHshape->setTooltip("");
        HHmaskSHshape->setTooltip("");
        radmaskSH->set_tooltip_text("");
        lapmaskSH->set_tooltip_text("");
        LmaskSHshape->setTooltip("");
        for (unsigned int i = 0; i < multipliersh.size(); i++) {
            multipliersh[i]->set_tooltip_text(M(""));
        }
        gamSH->set_tooltip_text(M(""));
        sloSH->set_tooltip_text(M(""));
        blendmaskSH->set_tooltip_text(M(""));
        mask2SHCurveEditorG->set_tooltip_text(M(""));
    }
}

void LocallabShadow::setDefaultExpanderVisibility()
{
    expgradsh->set_expanded(false);
    expmasksh->set_expanded(false);
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
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visishadhigh);
        exp->setEnabled(spot.expshadhigh);
        complexity->set_active(spot.complexshadhigh);

        if (spot.shMethod == "std") {
            shMethod->set_active(0);
        } else if (spot.shMethod == "tone") {
            shMethod->set_active(1);
        }

        for (int i = 0; i < 5; i++) {
            multipliersh[i]->setValue((double)spot.multsh[i]);
        }

        detailSH->setValue((double)spot.detailSH);
        highlights->setValue((double)spot.highlights);
        h_tonalwidth->setValue((double)spot.h_tonalwidth);
        shadows->setValue(spot.shadows);
        s_tonalwidth->setValue((double)spot.s_tonalwidth);
        sh_radius->setValue((double)spot.sh_radius);
        sensihs->setValue((double)spot.sensihs);
        blurSHde->setValue((double)spot.blurSHde);
        gamSH->setValue(spot.gamSH);
        sloSH->setValue(spot.sloSH);
        strSH->setValue(spot.strSH);
        angSH->setValue(spot.angSH);
        inverssh->set_active(spot.inverssh);
        enaSHMask->set_active(spot.enaSHMask);
        CCmaskSHshape->setCurve(spot.CCmaskSHcurve);
        LLmaskSHshape->setCurve(spot.LLmaskSHcurve);
        HHmaskSHshape->setCurve(spot.HHmaskSHcurve);
        blendmaskSH->setValue((double)spot.blendmaskSH);
        radmaskSH->setValue(spot.radmaskSH);
        lapmaskSH->setValue(spot.lapmaskSH);
        chromaskSH->setValue(spot.chromaskSH);
        gammaskSH->setValue(spot.gammaskSH);
        slomaskSH->setValue(spot.slomaskSH);
        LmaskSHshape->setCurve(spot.LmaskSHcurve);
        fatamountSH->setValue(spot.fatamountSH);
        fatanchorSH->setValue(spot.fatanchorSH);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

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
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expshadhigh = exp->getEnabled();
        spot.visishadhigh = exp->get_visible();
        spot.complexshadhigh = complexity->get_active_row_number();

        if (shMethod->get_active_row_number() == 0) {
            spot.shMethod = "std";
        } else if (shMethod->get_active_row_number() == 1) {
            spot.shMethod = "tone";
        }

        for (int i = 0; i < 5; i++) {
            spot.multsh[i] = multipliersh[i]->getIntValue();
        }

        spot.detailSH = detailSH->getIntValue();
        spot.highlights = highlights->getIntValue();
        spot.h_tonalwidth = h_tonalwidth->getIntValue();
        spot.shadows = shadows->getIntValue();
        spot.s_tonalwidth = s_tonalwidth->getIntValue();
        spot.sh_radius = sh_radius->getIntValue();
        spot.sensihs = sensihs->getIntValue();
        spot.blurSHde = blurSHde->getIntValue();
        spot.gamSH = gamSH->getValue();
        spot.sloSH = sloSH->getValue();
        spot.strSH = strSH->getValue();
        spot.angSH = angSH->getValue();
        spot.inverssh = inverssh->get_active();
        spot.enaSHMask = enaSHMask->get_active();
        spot.LLmaskSHcurve = LLmaskSHshape->getCurve();
        spot.CCmaskSHcurve = CCmaskSHshape->getCurve();
        spot.HHmaskSHcurve = HHmaskSHshape->getCurve();
        spot.blendmaskSH = blendmaskSH->getIntValue();
        spot.radmaskSH = radmaskSH->getValue();
        spot.lapmaskSH = lapmaskSH->getValue();
        spot.chromaskSH = chromaskSH->getValue();
        spot.gammaskSH = gammaskSH->getValue();
        spot.slomaskSH = slomaskSH->getValue();
        spot.LmaskSHcurve = LmaskSHshape->getCurve();
        spot.fatamountSH = fatamountSH->getValue();
        spot.fatanchorSH = fatanchorSH->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

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

void LocallabShadow::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    blurSHde->setValue((double)defSpot.blurSHde);
    lapmaskSH->setValue(defSpot.lapmaskSH);
    gammaskSH->setValue(defSpot.gammaskSH);
    slomaskSH->setValue(defSpot.slomaskSH);
    fatamountSH->setValue(defSpot.fatamountSH);
    fatanchorSH->setValue(defSpot.fatanchorSH);

    // Enable all listeners
    enableListener();
}

void LocallabShadow::updateGUIToMode(const modeType new_type)
{
    if (new_type == Normal) {
        // Advanced widgets are hidden in Normal mode
        blurSHde->hide();
        lapmaskSH->hide();
        gammaskSH->hide();
        slomaskSH->hide();
        fatSHFrame->hide();
    } else {
        // Advanced widgets are shown in Expert mode
        blurSHde->show();
        lapmaskSH->show();
        gammaskSH->show();
        slomaskSH->show();
        fatSHFrame->show();
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

    // This event is called to transmit potentially reset mask state
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
        gamFrame->hide();
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
        gamFrame->show();

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

    // Parameter Vibrance specific widgets
    saturated->setAdjusterListener(this);

    pastels->setAdjusterListener(this);

    warm->setAdjusterListener(this);

    psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    psThreshold->setAdjusterListener(this);

    pskinsConn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::protectskins_toggled));

    ashiftConn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::avoidcolorshift_toggled));

    pastsattogConn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::pastsattog_toggled));

    sensiv->setAdjusterListener(this);

    curveEditorGG->setCurveListener(this);

    skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
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

    strvib->setAdjusterListener(this);

    strvibab->set_tooltip_text(M("TP_LOCALLAB_GRADSTRAB_TOOLTIP"));
    strvibab->setAdjusterListener(this);

    strvibh->set_tooltip_text(M("TP_LOCALLAB_GRADSTRHUE_TOOLTIP"));
    strvibh->setAdjusterListener(this);

    angvib->set_tooltip_text(M("TP_LOCALLAB_GRADANG_TOOLTIP"));
    angvib->setAdjusterListener(this);

    setExpandAlignProperties(expmaskvib, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWMASK"));
    showmaskvibMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskvibMethod->set_active(0);
    showmaskvibMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskvibMethodConn = showmaskvibMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabVibrance::showmaskvibMethodChanged));

    enavibMaskConn = enavibMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::enavibMaskChanged));

    maskvibCurveEditorG->setCurveListener(this);

    CCmaskvibshape->setIdentityValue(0.);
    CCmaskvibshape->setResetCurve(FlatCurveType(defSpot.CCmaskvibcurve.at(0)), defSpot.CCmaskvibcurve);
    CCmaskvibshape->setBottomBarColorProvider(this, 1);

    LLmaskvibshape->setIdentityValue(0.);
    LLmaskvibshape->setResetCurve(FlatCurveType(defSpot.LLmaskvibcurve.at(0)), defSpot.LLmaskvibcurve);
    LLmaskvibshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskvibshape->setIdentityValue(0.);
    HHmaskvibshape->setResetCurve(FlatCurveType(defSpot.HHmaskvibcurve.at(0)), defSpot.HHmaskvibcurve);
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
    Lmaskvibshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskvibshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2vibCurveEditorG->curveListComplete();

    // Add Vibrance specific widgets to GUI
    pack_start(*saturated, Gtk::PACK_SHRINK, 0);
    pack_start(*pastels, Gtk::PACK_SHRINK, 0);
    pack_start(*warm, Gtk::PACK_SHRINK, 0);
    pack_start(*psThreshold, Gtk::PACK_SHRINK, 0);
    pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);
    pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);
    pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);
//    pack_start(*sensiv, Gtk::PACK_SHRINK, 0);
    pack_start(*curveEditorGG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    ToolParamBlock* const gradvibBox = Gtk::manage(new ToolParamBlock());
    gradvibBox->pack_start(*strvib);
    gradvibBox->pack_start(*strvibab);
    gradvibBox->pack_start(*strvibh);
    gradvibBox->pack_start(*angvib);
    expgradvib->add(*gradvibBox, false);
    pack_start(*expgradvib);
    ToolParamBlock* const maskvibBox = Gtk::manage(new ToolParamBlock());
    maskvibBox->pack_start(*showmaskvibMethod, Gtk::PACK_SHRINK, 4);
    maskvibBox->pack_start(*enavibMask, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*maskvibCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskvibBox->pack_start(*blendmaskvib, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*radmaskvib, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*lapmaskvib, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*chromaskvib, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*gammaskvib, Gtk::PACK_SHRINK, 0);
    maskvibBox->pack_start(*slomaskvib, Gtk::PACK_SHRINK, 0);
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

bool LocallabVibrance::isMaskViewActive()
{
    return (showmaskvibMethod->get_active_row_number() != 0);
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

void LocallabVibrance::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        warm->set_tooltip_text(M("TP_LOCALLAB_WARM_TOOLTIP"));
        strvib->set_tooltip_text(M("TP_LOCALLAB_GRADGEN_TOOLTIP"));
        expmaskvib->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskvibshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskvibshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskvibshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        Lmaskvibshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        blendmaskvib->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        mask2vibCurveEditorG->set_tooltip_text(M("TP_LOCALLAB_CONTRASTCURVMASK_TOOLTIP"));
        } else {
        warm->set_tooltip_text("");
        strvib->set_tooltip_text("");
        expmaskvib->set_tooltip_text("");
        CCmaskvibshape->setTooltip("");
        LLmaskvibshape->setTooltip("");
        HHmaskvibshape->setTooltip("");
        Lmaskvibshape->setTooltip("");
        blendmaskvib->set_tooltip_text(M(""));
        mask2vibCurveEditorG->set_tooltip_text(M(""));
    }
}

void LocallabVibrance::setDefaultExpanderVisibility()
{
    expgradvib->set_expanded(false);
    expmaskvib->set_expanded(false);
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
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visivibrance);
        exp->setEnabled(spot.expvibrance);
        complexity->set_active(spot.complexvibrance);

        saturated->setValue(spot.saturated);
        pastels->setValue(spot.pastels);
        warm->setValue(spot.warm);
        psThreshold->setValue<int>(spot.psthreshold);
        protectSkins->set_active(spot.protectskins);
        avoidColorShift->set_active(spot.avoidcolorshift);
        pastSatTog->set_active(spot.pastsattog);
        sensiv->setValue(spot.sensiv);
        skinTonesCurve->setCurve(spot.skintonescurve);
        strvib->setValue(spot.strvib);
        strvibab->setValue(spot.strvibab);
        strvibh->setValue(spot.strvibh);
        angvib->setValue(spot.angvib);
        enavibMask->set_active(spot.enavibMask);
        CCmaskvibshape->setCurve(spot.CCmaskvibcurve);
        LLmaskvibshape->setCurve(spot.LLmaskvibcurve);
        HHmaskvibshape->setCurve(spot.HHmaskvibcurve);
        blendmaskvib->setValue(spot.blendmaskvib);
        radmaskvib->setValue(spot.radmaskvib);
        lapmaskvib->setValue(spot.lapmaskvib);
        chromaskvib->setValue(spot.chromaskvib);
        gammaskvib->setValue(spot.gammaskvib);
        slomaskvib->setValue(spot.slomaskvib);
        Lmaskvibshape->setCurve(spot.Lmaskvibcurve);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update vibrance GUI according to pastsattog button state
    updateVibranceGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabVibrance::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expvibrance = exp->getEnabled();
        spot.visivibrance = exp->get_visible();
        spot.complexvibrance = complexity->get_active_row_number();

        spot.saturated = saturated->getIntValue();
        spot.pastels = pastels->getIntValue();
        spot.warm = warm->getIntValue();
        spot.psthreshold = psThreshold->getValue<int>();
        spot.protectskins = protectSkins->get_active();
        spot.avoidcolorshift = avoidColorShift->get_active();
        spot.pastsattog = pastSatTog->get_active();
        spot.sensiv = sensiv->getIntValue();
        spot.skintonescurve = skinTonesCurve->getCurve();
        spot.strvib = strvib->getValue();
        spot.strvibab = strvibab->getValue();
        spot.strvibh = strvibh->getValue();
        spot.angvib = angvib->getValue();
        spot.enavibMask = enavibMask->get_active();
        spot.CCmaskvibcurve = CCmaskvibshape->getCurve();
        spot.LLmaskvibcurve = LLmaskvibshape->getCurve();
        spot.HHmaskvibcurve = HHmaskvibshape->getCurve();
        spot.blendmaskvib = blendmaskvib->getIntValue();
        spot.radmaskvib = radmaskvib->getValue();
        spot.lapmaskvib = lapmaskvib->getValue();
        spot.chromaskvib = chromaskvib->getValue();
        spot.gammaskvib = gammaskvib->getValue();
        spot.slomaskvib = slomaskvib->getValue();
        spot.Lmaskvibcurve = Lmaskvibshape->getCurve();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabVibrance::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

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

void LocallabVibrance::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    saturated->setValue((double)defSpot.saturated);
    psThreshold->setValue<int>(defSpot.psthreshold);
    protectSkins->set_active(defSpot.protectskins);
    avoidColorShift->set_active(defSpot.avoidcolorshift);
    pastSatTog->set_active(defSpot.pastsattog);
    skinTonesCurve->setCurve(defSpot.skintonescurve);
    strvibab->setValue(defSpot.strvibab);
    strvibh->setValue(defSpot.strvibh);
    lapmaskvib->setValue(defSpot.lapmaskvib);
    gammaskvib->setValue(defSpot.gammaskvib);
    slomaskvib->setValue(defSpot.slomaskvib);

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update vibrance GUI according to pastsattog button state
    updateVibranceGUI();
}

void LocallabVibrance::updateGUIToMode(const modeType new_type)
{
    if (new_type == Normal) {
        // Advanced widgets are hidden in Normal mode
        saturated->hide();
        pastels->setLabel(M("TP_LOCALLAB_PASTELS2"));
        psThreshold->hide();
        protectSkins->hide();
        avoidColorShift->hide();
        pastSatTog->hide();
        curveEditorGG->hide();
        strvibab->hide();
        strvibh->hide();
        lapmaskvib->hide();
        gammaskvib->hide();
        slomaskvib->hide();
    } else {
        // Advanced widgets are shown in Expert mode
        saturated->show();
        pastels->setLabel(M("TP_VIBRANCE_PASTELS"));
        psThreshold->show();
        protectSkins->show();
        avoidColorShift->show();
        pastSatTog->show();
        curveEditorGG->show();
        strvibab->show();
        strvibh->show();
        lapmaskvib->show();
        gammaskvib->show();
        slomaskvib->show();
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
    sensisf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 1, 100, 1, 30)))
{
    // Parameter Soft light specific widgets
    softMethod->append(M("TP_LOCALLAB_SOFTM"));
    softMethod->append(M("TP_LOCALLAB_RETIM"));
    softMethod->set_active(0);
    softMethodConn = softMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSoft::softMethodChanged));

    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWLAPLACE"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWFOURIER"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWPOISSON"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWNORMAL"));
    showmasksoftMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmasksoftMethod->set_active(0);
//    showmasksoftMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKSOFT_TOOLTIP"));
    showmasksoftMethodConn = showmasksoftMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSoft::showmasksoftMethodChanged));

    streng->setAdjusterListener(this);

    laplace->setAdjusterListener(this);

    sensisf->setAdjusterListener(this);

    // Add Soft light specific widgets to GUI
    pack_start(*softMethod);
    Gtk::Label* const labelsoftmethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_SHOWDCT") + ":"));
    ctboxsoftmethod->pack_start(*labelsoftmethod, Gtk::PACK_SHRINK, 4);
    ctboxsoftmethod->pack_start(*showmasksoftMethod);
    pack_start(*ctboxsoftmethod);
    pack_start(*streng);
    pack_start(*laplace);
    pack_start(*sensisf);
}

bool LocallabSoft::isMaskViewActive()
{
    return (showmasksoftMethod->get_active_row_number() != 0);
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

void LocallabSoft::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        softMethod->set_tooltip_markup(M("TP_LOCALLAB_SOFTMETHOD_TOOLTIP"));
       // ctboxsoftmethod->set_tooltip_markup(M("TP_LOCALLAB_ORRETISTEP_TOOLTIP"));
        showmasksoftMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKSOFT_TOOLTIP"));
        streng->set_tooltip_text(M("TP_LOCALLAB_ORRETISTREN_TOOLTIP"));
        laplace->set_tooltip_text(M("TP_LOCALLAB_ORRETILAP_TOOLTIP"));
    } else {
        softMethod->set_tooltip_markup(M(""));
       // ctboxsoftmethod->set_tooltip_markup(M(""));
        showmasksoftMethod->set_tooltip_markup(M(""));
        streng->set_tooltip_text(M(""));
        laplace->set_tooltip_text(M(""));
    }
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
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visisoft);
        exp->setEnabled(spot.expsoft);
        complexity->set_active(spot.complexsoft);

        if (spot.softMethod == "soft") {
            softMethod->set_active(0);
        } else if (spot.softMethod == "reti") {
            softMethod->set_active(1);
        }

        streng->setValue((double)spot.streng);
        sensisf->setValue((double)spot.sensisf);
        laplace->setValue(spot.laplace);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update soft light GUI according to softMethod combobox
    updateSoftGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSoft::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expsoft = exp->getEnabled();
        spot.visisoft = exp->get_visible();
        spot.complexsoft = complexity->get_active_row_number();

        if (softMethod->get_active_row_number() == 0) {
            spot.softMethod = "soft";
        } else if (softMethod->get_active_row_number() == 1) {
            spot.softMethod = "reti";
        }

        spot.streng = streng->getIntValue();
        spot.sensisf = sensisf->getIntValue();
        spot.laplace = laplace->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabSoft::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

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

void LocallabSoft::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    if (defSpot.softMethod == "soft") {
        softMethod->set_active(0);
    } else if (defSpot.softMethod == "reti") {
        softMethod->set_active(1);
    }

    // Enable all listeners
    enableListener();

    // Update GUI based on converted widget parameters:
    // - Update soft light GUI according to softMethod combobox
    updateSoftGUI();
}

void LocallabSoft::updateGUIToMode(const modeType new_type)
{
    if (new_type == Normal) {
        // Advanced widgets are hidden in Normal mode
        softMethod->hide();
    } else {
        // Advanced widgets are shown in Expert mode
        softMethod->show();
    }
}

void LocallabSoft::softMethodChanged()
{
    // Update soft light GUI according to softMethod combobox
    updateSoftGUI();

    // This event is called to transmit potentially reset mask state
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
    expblnoise(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_BLNOI_EXP")))),
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
    expdenoise(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_DENOI_EXP")))),
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
    adjblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJ"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-yellow-small.png")), Gtk::manage(new RTImage("circle-red-green-small.png"))))),
    bilateral(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BILATERAL"), 0, 100, 1, 0))),
    sensiden(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIDEN"), 0, 100, 1, 60))),
    expmaskbl(Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOWPLUS")))),
    showmaskblMethod(Gtk::manage(new MyComboBoxText())),
    showmaskblMethodtyp(Gtk::manage(new MyComboBoxText())),
    enablMask(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")))),
    maskblCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    CCmaskblshape(static_cast<FlatCurveEditor*>(maskblCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false))),
    LLmaskblshape(static_cast<FlatCurveEditor*>(maskblCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    HHmaskblshape(static_cast<FlatCurveEditor *>(maskblCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true))),
    strumaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUMASKCOL"), 0., 200., 0.1, 0.))),
    toolbl(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TOOLCOL")))),
    blendmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0))),
    radmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), -10.0, 1000.0, 0.1, 0.))),
    lapmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPMASKCOL"), 0.0, 100.0, 0.1, 0.))),
    chromaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.))),
    gammaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.05, 5.0, 0.01, 1.))),
    slomaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.))),
    shadmaskbl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_HIGHMASKCOL"), 0, 100, 1, 0))),
    shadmaskblsha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHAMASKCOL"), 0, 100, 1, 0))),
    mask2blCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"))),
    Lmaskblshape(static_cast<DiagonalCurveEditor*>(mask2blCurveEditorG->addCurve(CT_Diagonal, "L(L)"))),
    mask2blCurveEditorGwav(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_WAVMASK"))),
    LLmaskblshapewav(static_cast<FlatCurveEditor*>(mask2blCurveEditorGwav->addCurve(CT_Flat, "L(L)", nullptr, false, false))),
    csThresholdblur(Gtk::manage(new ThresholdAdjuster(M("TP_LOCALLAB_CSTHRESHOLDBLUR"), 0, 9, 0, 0, 6, 5, 0, false)))
{
    const LocallabParams::LocallabSpot defSpot;

    // Parameter Blur, Noise & Denoise specific widgets
    setExpandAlignProperties(expblnoise, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    blMethod->append(M("TP_LOCALLAB_BLUR"));
    blMethod->append(M("TP_LOCALLAB_BLMED"));
    blMethod->append(M("TP_LOCALLAB_BLGUID"));
    blMethod->set_active(0);
    blMethodConn = blMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::blMethodChanged));

    fftwblConn = fftwbl->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::fftwblChanged));

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

    sensibn->setAdjusterListener(this);

    blurMethod->append(M("TP_LOCALLAB_BLNORM"));
    blurMethod->append(M("TP_LOCALLAB_BLINV"));
    blurMethod->set_active(0);
    blurMethodConn = blurMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::blurMethodChanged));

    chroMethod->append(M("TP_LOCALLAB_BLLO"));
    chroMethod->append(M("TP_LOCALLAB_BLCO"));
    chroMethod->append(M("TP_LOCALLAB_BLLC"));
    chroMethod->set_active(0);
    chroMethodConn = chroMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::chroMethodChanged));

    activlumConn = activlum->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::activlumChanged));

    setExpandAlignProperties(expdenoise, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    LocalcurveEditorwavden->setCurveListener(this);

    wavshapeden->setIdentityValue(0.);
    wavshapeden->setResetCurve(FlatCurveType(defSpot.locwavcurveden.at(0)), defSpot.locwavcurveden);

    LocalcurveEditorwavden->curveListComplete();

    noiselumf0->setAdjusterListener(this);

    noiselumf->setAdjusterListener(this);

    noiselumf2->setAdjusterListener(this);

    noiselumc->setAdjusterListener(this);

    noiselumdetail->setAdjusterListener(this);

    noiselequal->setAdjusterListener(this);

    noisechrof->setAdjusterListener(this);

    noisechroc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
    noisechroc->setAdjusterListener(this);

    noisechrodetail->setAdjusterListener(this);

    detailthr->setAdjusterListener(this);

    adjblur->setAdjusterListener(this);

    bilateral->setAdjusterListener(this);

    sensiden->setAdjusterListener(this);

    setExpandAlignProperties(expmaskbl, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
    showmaskblMethod->append(M("TP_LOCALLAB_SHOWMASK"));
//    showmaskblMethod->append(M("TP_LOCALLAB_SHOWREF"));
    showmaskblMethod->set_active(0);
    showmaskblMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
    showmaskblMethodConn = showmaskblMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::showmaskblMethodChanged));

    showmaskblMethodtyp->append(M("TP_LOCALLAB_SHOWMASKTYP1"));
    showmaskblMethodtyp->append(M("TP_LOCALLAB_SHOWMASKTYP2"));
    showmaskblMethodtyp->append(M("TP_LOCALLAB_SHOWMASKTYP3"));
    showmaskblMethodtyp->set_active(0);
    showmaskblMethodtypConn = showmaskblMethodtyp->signal_changed().connect(sigc::mem_fun(*this, &LocallabBlur::showmaskblMethodtypChanged));

    enablMaskConn = enablMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::enablMaskChanged));

    maskblCurveEditorG->setCurveListener(this);

    CCmaskblshape->setIdentityValue(0.);
    CCmaskblshape->setResetCurve(FlatCurveType(defSpot.CCmaskblcurve.at(0)), defSpot.CCmaskblcurve);
    CCmaskblshape->setBottomBarColorProvider(this, 1);

    LLmaskblshape->setIdentityValue(0.);
    LLmaskblshape->setResetCurve(FlatCurveType(defSpot.LLmaskblcurve.at(0)), defSpot.LLmaskblcurve);
    LLmaskblshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    HHmaskblshape->setIdentityValue(0.);
    HHmaskblshape->setResetCurve(FlatCurveType(defSpot.HHmaskblcurve.at(0)), defSpot.HHmaskblcurve);
    HHmaskblshape->setCurveColorProvider(this, 2);
    HHmaskblshape->setBottomBarColorProvider(this, 2);

    maskblCurveEditorG->curveListComplete();

    strumaskbl->setAdjusterListener(this);

    toolblConn = toolbl->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::toolblChanged));

    blendmaskbl->setAdjusterListener(this);

    radmaskbl->setAdjusterListener(this);

    lapmaskbl->setAdjusterListener(this);

    chromaskbl->setAdjusterListener(this);

    gammaskbl->setAdjusterListener(this);

    slomaskbl->setAdjusterListener(this);

    shadmaskbl->setAdjusterListener(this);

    shadmaskblsha->setAdjusterListener(this);

    mask2blCurveEditorG->setCurveListener(this);

    Lmaskblshape->setResetCurve(DiagonalCurveType(defSpot.Lmaskblcurve.at(0)), defSpot.Lmaskblcurve);
    Lmaskblshape->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    Lmaskblshape->setLeftBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});
    mask2blCurveEditorG->curveListComplete();

    mask2blCurveEditorGwav->setCurveListener(this);

    LLmaskblshapewav->setIdentityValue(0.);
    LLmaskblshapewav->setResetCurve(FlatCurveType(defSpot.LLmaskblcurvewav.at(0)), defSpot.LLmaskblcurvewav);
    LLmaskblshapewav->setBottomBarBgGradient({{0., 0., 0., 0.}, {1., 1., 1., 1.}});

    mask2blCurveEditorGwav->curveListComplete();

    csThresholdblur->setAdjusterListener(this);

    // Add Blur, Noise & Denoise specific widgets to GUI
    ToolParamBlock* const blnoisebox = Gtk::manage(new ToolParamBlock());
    blnoisebox->pack_start(*blMethod);
    blnoisebox->pack_start(*fftwbl, Gtk::PACK_SHRINK, 0);
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
    maskblBox->pack_start(*showmaskblMethodtyp, Gtk::PACK_SHRINK, 4);
    maskblBox->pack_start(*enablMask, Gtk::PACK_SHRINK, 0);
    maskblBox->pack_start(*maskblCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    maskblBox->pack_start(*strumaskbl, Gtk::PACK_SHRINK, 0);
    maskblBox->pack_start(*toolbl, Gtk::PACK_SHRINK, 0);
    Gtk::HSeparator* const separatorstrubl = Gtk::manage(new  Gtk::HSeparator());
    maskblBox->pack_start(*separatorstrubl, Gtk::PACK_SHRINK, 2);
    maskblBox->pack_start(*blendmaskbl, Gtk::PACK_SHRINK, 0);
    Gtk::Frame* const toolblFrame = Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_TOOLMASK")));
    toolblFrame->set_label_align(0.025, 0.5);
    ToolParamBlock* const toolblBox = Gtk::manage(new ToolParamBlock());
    toolblBox->pack_start(*radmaskbl, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*lapmaskbl, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*chromaskbl, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*gammaskbl, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*slomaskbl, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*shadmaskblsha, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*shadmaskbl, Gtk::PACK_SHRINK, 0);
    toolblBox->pack_start(*mask2blCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolblBox->pack_start(*mask2blCurveEditorGwav, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    toolblBox->pack_start(*csThresholdblur, Gtk::PACK_SHRINK, 0);
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

bool LocallabBlur::isMaskViewActive()
{
    return (showmaskblMethod->get_active_row_number() != 0);
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

void LocallabBlur::updateAdviceTooltips(const bool showTooltips)
{
    if (showTooltips) {
        expblnoise->set_tooltip_markup(M("TP_LOCALLAB_BLUMETHOD_TOOLTIP"));
        radius->set_tooltip_text(M("TP_LOCALLAB_RADIUS_TOOLTIP"));
        sensibn->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
        blurMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));
        expdenoise->set_tooltip_markup(M("TP_LOCALLAB_DENOI_TOOLTIP"));
        wavshapeden->setTooltip(M("TP_LOCALLAB_WASDEN_TOOLTIP"));
        noiselumc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));
        expmaskbl->set_tooltip_markup(M("TP_LOCALLAB_MASK_TOOLTIP"));
        CCmaskblshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        LLmaskblshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        HHmaskblshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        radmaskbl->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        lapmaskbl->set_tooltip_text(M("TP_LOCALLAB_LAPRAD_TOOLTIP"));
        Lmaskblshape->setTooltip(M("TP_LOCALLAB_LMASK_LL_TOOLTIP"));
        LLmaskblshapewav->setTooltip(M("TP_LOCALLAB_LMASK_LEVEL_TOOLTIP"));
        blendmaskbl->set_tooltip_text(M("TP_LOCALLAB_BLENDMASK_TOOLTIP"));
        showmaskblMethodtyp->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKTYP_TOOLTIP"));
    } else {
        expblnoise->set_tooltip_text("");
        radius->set_tooltip_text("");
        sensibn->set_tooltip_text("");
        blurMethod->set_tooltip_text("");
        expdenoise->set_tooltip_text("");
        wavshapeden->setTooltip("");
        noiselumc->set_tooltip_text("");
        expmaskbl->set_tooltip_text("");
        CCmaskblshape->setTooltip("");
        LLmaskblshape->setTooltip("");
        HHmaskblshape->setTooltip("");
        radmaskbl->set_tooltip_text("");
        lapmaskbl->set_tooltip_text("");
        Lmaskblshape->setTooltip("");
        LLmaskblshapewav->setTooltip("");
        blendmaskbl->set_tooltip_text(M(""));
        showmaskblMethodtyp->set_tooltip_markup(M(""));
    }
}

void LocallabBlur::setDefaultExpanderVisibility()
{
    expblnoise->set_expanded(false);
    expdenoise->set_expanded(false);
    expmaskbl->set_expanded(false);
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
    showmaskblMethodtypConn.block(true);
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
    showmaskblMethodtypConn.block(false);
    enablMaskConn.block(false);
    toolblConn.block(false);
}

void LocallabBlur::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spotName = spot.name; // Update spot name according to selected spot

        exp->set_visible(spot.visiblur);
        exp->setEnabled(spot.expblur);
        complexity->set_active(spot.complexblur);

        if (spot.blMethod == "blur") {
            blMethod->set_active(0);
        } else if (spot.blMethod == "med") {
            blMethod->set_active(1);
        } else if (spot.blMethod == "guid") {
            blMethod->set_active(2);
        }

        fftwbl->set_active(spot.fftwbl);
        radius->setValue(spot.radius);
        strength->setValue(spot.strength);
        isogr->setValue((double)spot.isogr);
        strengr->setValue((double)spot.strengr);
        scalegr->setValue((double)spot.scalegr);

        if (spot.medMethod == "none") {
            medMethod->set_active(0);
        } else if (spot.medMethod == "33") {
            medMethod->set_active(1);
        } else if (spot.medMethod == "55") {
            medMethod->set_active(2);
        } else if (spot.medMethod == "77") {
            medMethod->set_active(3);
        } else if (spot.medMethod == "99") {
            medMethod->set_active(4);
        }

        itera->setValue((double)spot.itera);
        guidbl->setValue((double)spot.guidbl);
        strbl->setValue((double)spot.strbl);
        epsbl->setValue((double)spot.epsbl);
        sensibn->setValue((double)spot.sensibn);

        if (spot.blurMethod == "norm") {
            blurMethod->set_active(0);
        } else if (spot.blurMethod == "inv") {
            blurMethod->set_active(1);
        }

        if (spot.chroMethod == "lum") {
            chroMethod->set_active(0);
        } else if (spot.chroMethod == "chr") {
            chroMethod->set_active(1);
        } else if (spot.chroMethod == "all") {
            chroMethod->set_active(2);
        }


        if (spot.showmaskblMethodtyp == "blur") {
           showmaskblMethodtyp ->set_active(0);
        } else if (spot.showmaskblMethodtyp == "nois") {
            showmaskblMethodtyp->set_active(1);
        } else if (spot.showmaskblMethodtyp == "all") {
            showmaskblMethodtyp->set_active(2);
        }

        activlum->set_active(spot.activlum);
        wavshapeden->setCurve(spot.locwavcurveden);
        noiselumf0->setValue(spot.noiselumf0);
        noiselumf->setValue(spot.noiselumf);
        noiselumf2->setValue(spot.noiselumf2);
        noiselumc->setValue(spot.noiselumc);
        noiselumdetail->setValue(spot.noiselumdetail);
        noiselequal->setValue((double)spot.noiselequal);
        noisechrof->setValue(spot.noisechrof);
        noisechroc->setValue(spot.noisechroc);
        noisechrodetail->setValue(spot.noisechrodetail);
        detailthr->setValue((double)spot.detailthr);
        adjblur->setValue((double)spot.adjblur);
        bilateral->setValue((double)spot.bilateral);
        sensiden->setValue((double)spot.sensiden);
        enablMask->set_active(spot.enablMask);
        CCmaskblshape->setCurve(spot.CCmaskblcurve);
        LLmaskblshape->setCurve(spot.LLmaskblcurve);
        HHmaskblshape->setCurve(spot.HHmaskblcurve);
        strumaskbl->setValue(spot.strumaskbl);
        toolbl->set_active(spot.toolbl);
        blendmaskbl->setValue((double)spot.blendmaskbl);
        radmaskbl->setValue(spot.radmaskbl);
        lapmaskbl->setValue(spot.lapmaskbl);
        chromaskbl->setValue(spot.chromaskbl);
        gammaskbl->setValue(spot.gammaskbl);
        slomaskbl->setValue(spot.slomaskbl);
        shadmaskbl->setValue((double)spot.shadmaskbl);
        shadmaskblsha->setValue((double)spot.shadmaskblsha);
        Lmaskblshape->setCurve(spot.Lmaskblcurve);
        LLmaskblshapewav->setCurve(spot.LLmaskblcurvewav);
        csThresholdblur->setValue<int>(spot.csthresholdblur);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to complexity mode
    updateGUIToMode(static_cast<modeType>(complexity->get_active_row_number()));

    // Update Blur & Noise GUI according to blMethod combobox state
    updateBlurGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        LocallabParams::LocallabSpot& spot = pp->locallab.spots.at(index);

        spot.expblur = exp->getEnabled();
        spot.visiblur = exp->get_visible();
        spot.complexblur = complexity->get_active_row_number();

        if (blMethod->get_active_row_number() == 0) {
            spot.blMethod = "blur";
        } else if (blMethod->get_active_row_number() == 1) {
            spot.blMethod = "med";
        } else if (blMethod->get_active_row_number() == 2) {
            spot.blMethod = "guid";
        }

        spot.fftwbl = fftwbl->get_active();
        spot.radius = radius->getValue();
        spot.strength = strength->getIntValue();
        spot.isogr = isogr->getIntValue();
        spot.strengr = strengr->getIntValue();
        spot.scalegr = scalegr->getIntValue();

        if (medMethod->get_active_row_number() == 0) {
            spot.medMethod = "none";
        } else if (medMethod->get_active_row_number() == 1) {
            spot.medMethod = "33";
        } else if (medMethod->get_active_row_number() == 2) {
            spot.medMethod = "55";
        } else if (medMethod->get_active_row_number() == 3) {
            spot.medMethod = "77";
        } else if (medMethod->get_active_row_number() == 4) {
            spot.medMethod = "99";
        }

        spot.itera = itera->getIntValue();
        spot.guidbl = guidbl->getIntValue();
        spot.strbl = strbl->getIntValue();
        spot.epsbl = epsbl->getIntValue();
        spot.sensibn = sensibn->getIntValue();

        if (blurMethod->get_active_row_number() == 0) {
            spot.blurMethod = "norm";
        } else if (blurMethod->get_active_row_number() == 1) {
            spot.blurMethod = "inv";
        }

        if (chroMethod->get_active_row_number() == 0) {
            spot.chroMethod = "lum";
        } else if (chroMethod->get_active_row_number() == 1) {
            spot.chroMethod = "chr";
        } else if (chroMethod->get_active_row_number() == 2) {
            spot.chroMethod = "all";
        }


        if (showmaskblMethodtyp->get_active_row_number() == 0) {
            spot.showmaskblMethodtyp = "blur";
        } else if (showmaskblMethodtyp->get_active_row_number() == 1) {
            spot.showmaskblMethodtyp = "nois";
        } else if (showmaskblMethodtyp->get_active_row_number() == 2) {
            spot.showmaskblMethodtyp = "all";
        }

        spot.activlum = activlum->get_active();
        spot.locwavcurveden = wavshapeden->getCurve();
        spot.noiselumf0 = noiselumf0->getValue();
        spot.noiselumf = noiselumf->getValue();
        spot.noiselumf2 = noiselumf2->getValue();
        spot.noiselumc = noiselumc->getValue();
        spot.noiselumdetail = noiselumdetail->getValue();
        spot.noiselequal = noiselequal->getIntValue();
        spot.noisechrof = noisechrof->getValue();
        spot.noisechroc = noisechroc->getValue();
        spot.noisechrodetail = noisechrodetail->getValue();
        spot.detailthr = detailthr->getIntValue();
        spot.adjblur = adjblur->getIntValue();
        spot.bilateral = bilateral->getIntValue();
        spot.sensiden = sensiden->getIntValue();
        spot.enablMask = enablMask->get_active();
        spot.LLmaskblcurve = LLmaskblshape->getCurve();
        spot.CCmaskblcurve = CCmaskblshape->getCurve();
        spot.HHmaskblcurve = HHmaskblshape->getCurve();
        spot.strumaskbl = strumaskbl->getValue();
        spot.toolbl = toolbl->get_active();
        spot.blendmaskbl = blendmaskbl->getIntValue();
        spot.radmaskbl = radmaskbl->getValue();
        spot.lapmaskbl = lapmaskbl->getValue();
        spot.chromaskbl = chromaskbl->getValue();
        spot.gammaskbl = gammaskbl->getValue();
        spot.slomaskbl = slomaskbl->getValue();
        spot.shadmaskbl = shadmaskbl->getIntValue();
        spot.shadmaskblsha = shadmaskblsha->getIntValue();
        spot.Lmaskblcurve = Lmaskblshape->getCurve();
        spot.LLmaskblcurvewav = LLmaskblshapewav->getCurve();
        spot.csthresholdblur = csThresholdblur->getValue<int>();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot& defSpot = defParams->locallab.spots.at(index);

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
        shadmaskblsha->setDefault(defSpot.shadmaskblsha);
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

        if (a == shadmaskblsha) {
            if (listener) {
                listener->panelChanged(Evlocallabshadmaskblsha,
                                       shadmaskblsha->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

void LocallabBlur::convertParamToNormal()
{
    const LocallabParams::LocallabSpot defSpot;

    // Disable all listeners
    disableListener();

    // Set hidden GUI widgets in Normal mode to default spot values
    fftwbl->set_active(defSpot.fftwbl);
    strumaskbl->setValue(defSpot.strumaskbl);
    toolbl->set_active(defSpot.toolbl);
    lapmaskbl->setValue(defSpot.lapmaskbl);
//    gammaskbl->setValue(defSpot.gammaskbl);
//    slomaskbl->setValue(defSpot.slomaskbl);
    shadmaskbl->setValue((double)defSpot.shadmaskbl);
    shadmaskblsha->setValue((double)defSpot.shadmaskblsha);
    LLmaskblshapewav->setCurve(defSpot.LLmaskblcurvewav);
    csThresholdblur->setValue<int>(defSpot.csthresholdblur);

    // Enable all listeners
    enableListener();
}

void LocallabBlur::updateGUIToMode(const modeType new_type)
{
    if (new_type == Normal) {
        // Advanced widgets are hidden in Normal mode
        fftwbl->hide();
        strumaskbl->hide();
        toolbl->hide();
        lapmaskbl->hide();
        gammaskbl->show();
        slomaskbl->show();
        shadmaskbl->hide();
        shadmaskblsha->hide();
        mask2blCurveEditorGwav->hide();
        csThresholdblur->hide();
    } else {
        // Advanced widgets are shown in Expert mode
        if (blMethod->get_active_row_number() == 0) { // Keep widget hidden when blMethod is > 0
            fftwbl->show();
        }

        strumaskbl->show();
        toolbl->show();
        lapmaskbl->show();
        gammaskbl->show();
        slomaskbl->show();
        shadmaskbl->show();
        shadmaskblsha->show();
        mask2blCurveEditorGwav->show();
        csThresholdblur->show();
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

void LocallabBlur::showmaskblMethodtypChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmasktypMethod,
                               showmaskblMethodtyp->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
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
