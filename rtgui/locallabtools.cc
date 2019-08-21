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

using namespace rtengine;
using namespace procparams;

extern Options options;

/* ==== LocallabTool ==== */
LocallabTool::LocallabTool(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11, maskType usemask):
    ToolPanel(toolName, need11),

    // LocallabTool parameters
    useMask(usemask),
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

    // Create mask panel
    if (useMask) {
        maskExp = Gtk::manage(new MyExpander(false, M("TP_LOCALLAB_SHOW")));
        setExpandAlignProperties(maskExp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);

        showMaskMethod = Gtk::manage(new MyComboBoxText());
        showMaskMethod->append(M("TP_LOCALLAB_SHOWMNONE"));
        showMaskMethod->append(M("TP_LOCALLAB_SHOWMODIF"));
        showMaskMethod->append(M("TP_LOCALLAB_SHOWMODIFMASK"));
        showMaskMethod->append(M("TP_LOCALLAB_SHOWMASK"));
        showMaskMethod->append(M("TP_LOCALLAB_SHOWSTRUC"));
        showMaskMethod->append(M("TP_LOCALLAB_PREVIEWSEL"));
        showMaskMethod->set_active(0);

        if (showtooltip) {
            showMaskMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKCOL_TOOLTIP"));
        }

        showMaskMethodConn = showMaskMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabTool::showMaskMethodChanged));

        enaMask = Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ENABLE_MASK")));
        enaMaskConn = enaMask->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTool::enaMaskChanged));

        if (useMask == MaskWithTrMap) {
            enaMaskTrMap = Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_TM_MASK")));
            enaMaskTrMapConn = enaMaskTrMap->signal_toggled().connect(sigc::mem_fun(*this, &LocallabTool::enaMaskTrMapChanged));
        }

        maskCurveEditorG = new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_MASK"));
        maskCurveEditorG->setCurveListener(this);

        CCMaskShape = static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "C(C)", nullptr, false, false));
        CCMaskShape->setIdentityValue(0.);
        CCMaskShape->setResetCurve(FlatCurveType(LocallabParams::DEF_MASK_CURVE.at(0)), LocallabParams::DEF_MASK_CURVE);

        if (showtooltip) {
            CCMaskShape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        }

        CCMaskShape->setBottomBarColorProvider(this, 1);

        std::vector<GradientMilestone> bgGradient;
        bgGradient.push_back(GradientMilestone(0., 0., 0., 0.));
        bgGradient.push_back(GradientMilestone(1., 1., 1., 1.));
        LLMaskShape = static_cast<FlatCurveEditor*>(maskCurveEditorG->addCurve(CT_Flat, "L(L)", nullptr, false, false));
        LLMaskShape->setIdentityValue(0.);
        LLMaskShape->setResetCurve(FlatCurveType(LocallabParams::DEF_MASK_CURVE.at(0)), LocallabParams::DEF_MASK_CURVE);

        if (showtooltip) {
            LLMaskShape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        }

        LLMaskShape->setBottomBarBgGradient(bgGradient);

        HHMaskShape = static_cast<FlatCurveEditor *>(maskCurveEditorG->addCurve(CT_Flat, "LC(H)", nullptr, false, true));
        HHMaskShape->setIdentityValue(0.);
        HHMaskShape->setResetCurve(FlatCurveType(LocallabParams::DEF_MASK_CURVE.at(0)), LocallabParams::DEF_MASK_CURVE);

        if (showtooltip) {
            HHMaskShape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_CC_TOOLTIP"));
        }

        HHMaskShape->setCurveColorProvider(this, 2);
        HHMaskShape->setBottomBarColorProvider(this, 2);

        maskCurveEditorG->curveListComplete();

        blendMask = Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLENDMASKCOL"), -100, 100, 1, 0));
        blendMask->setAdjusterListener(this);

        radMask = Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADMASKCOL"), 0.0, 100.0, 0.1, 10.));
        radMask->setAdjusterListener(this);

        chroMask = Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMASKCOL"), -100.0, 100.0, 0.1, 0.));
        chroMask->setAdjusterListener(this);

        gamMask = Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAMMASKCOL"), 0.25, 4.0, 0.01, 1.));
        gamMask->setAdjusterListener(this);

        sloMask = Gtk::manage(new Adjuster(M("TP_LOCALLAB_SLOMASKCOL"), 0.0, 15.0, 0.1, 0.));
        sloMask->setAdjusterListener(this);

        ToolParamBlock* const maskBox = Gtk::manage(new ToolParamBlock());
        maskBox->pack_start(*showMaskMethod, Gtk::PACK_SHRINK, 4);
        maskBox->pack_start(*enaMask, Gtk::PACK_SHRINK, 0);

        if (useMask == MaskWithTrMap) {
            maskBox->pack_start(*enaMaskTrMap, Gtk::PACK_SHRINK, 0);
        }

        maskBox->pack_start(*maskCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
        maskBox->pack_start(*blendMask, Gtk::PACK_SHRINK, 0);
        maskBox->pack_start(*radMask, Gtk::PACK_SHRINK, 0);
        maskBox->pack_start(*chroMask, Gtk::PACK_SHRINK, 0);
        maskBox->pack_start(*gamMask, Gtk::PACK_SHRINK, 0);
        maskBox->pack_start(*sloMask, Gtk::PACK_SHRINK, 0);
        maskExp->add(*maskBox, false);
        totalBox->pack_start(*maskExp, Gtk::PACK_SHRINK, 0);
    }

    exp->add(*totalBox, false);
    exp->setLevel(2);
}

LocallabTool::~LocallabTool()
{
    idle_register.destroy();

    if (useMask) {
        delete maskCurveEditorG;
    }
}

void LocallabTool::addLocallabTool(bool cond)
{
    exp->set_visible(cond);
}

bool LocallabTool::isLocallabToolAdded()
{
    return exp->get_visible();
}

void LocallabTool::resetMaskView()
{
    if (useMask) {
        showMaskMethodConn.block(true);
        showMaskMethod->set_active(0);
        showMaskMethodConn.block(false);
    }
}

void LocallabTool::refChanged(const double huer, const double lumar, const double chromar)
{
    if (useMask) {
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

        // printf("nh=%f nl=%f nc=%f\n", normHuer, normLumar, normChromar);

        idle_register.add(
        [this, normHuer, normLumar, normChromar]() -> bool {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

            // Update mask background
            CCMaskShape->updateLocallabBackground(normChromar);
            LLMaskShape->updateLocallabBackground(normLumar);
            HHMaskShape->updateLocallabBackground(normHuer);

            return false;
        }
        );
    }
}

void LocallabTool::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller)
{
    if (useMask) {
        float R = 0.f;
        float G = 0.f;
        float B = 0.f;

        if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
            valY = 0.5;
        }

        switch (callerId) {
            case 1: // CCmaskshape
                Color::hsv2rgb01(float(valY), float(valX), 0.5f, R, G, B);

                break;

            case 2: // HHmaskshape
                float x = valX - 1.f / 6.f;

                if (x < 0.f) {
                    x += 1.f;
                }

                x = log2lin(x, 3.f);
                Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        }

        caller->ccRed = double (R);
        caller->ccGreen = double (G);
        caller->ccBlue = double (B);
    }
}

void LocallabTool::disableListener()
{
    ToolPanel::disableListener();

    enaExpConn.block(true);

    if (useMask) {
        enaMaskConn.block(true);

        if (useMask == MaskWithTrMap) {
            enaMaskTrMapConn.block(true);
        }

        showMaskMethodConn.block(true);
    }
}
void LocallabTool::enableListener()
{
    ToolPanel::enableListener();

    enaExpConn.block(false);

    if (useMask) {
        enaMaskConn.block(false);

        if (useMask == MaskWithTrMap) {
            enaMaskTrMapConn.block(false);
        }

        showMaskMethodConn.block(false);
    }
}

bool LocallabTool::on_remove_change(GdkEventButton* event)
{
    if (event->button == GDK_BUTTON_PRIMARY) {
        printf("Remove icon pressed\n");
        // exp->set_visible(false);
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
    LocallabTool(this, "Locallab Color&Light", M("TP_LOCALLAB_COFR"), false, MaskNormal),

    // Color & Light specific widgets
    curvactiv(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CURV")))),
    lightness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    gridFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_LABGRID")))),
    labgrid(Gtk::manage(new LabGrid(EvLocallabLabGridValue, M("TP_LOCALLAB_LABGRID_VALUES")))),
    gridMethod(Gtk::manage(new MyComboBoxText())),
    strengthgrid(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRGRID"), 0, 100, 1, 30))),
    sensi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    structcol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    blurcolde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    softradiuscol(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    labqualcurv(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    qualitycurveMethod(Gtk::manage(new MyComboBoxText())),
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),
    HCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_HLH"))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    float R, G, B;

    const bool showtooltip = options.showtooltip;

    curvactivConn = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::curvactivChanged));

    lightness->setAdjusterListener(this);

    if (showtooltip) {
        lightness->set_tooltip_text(M("TP_LOCALLAB_EXPCOLOR_TOOLTIP"));
    }

    contrast->setAdjusterListener(this);

    chroma->setAdjusterListener(this);

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

    softradiuscol->setAdjusterListener(this);

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->set_active(0);

    if (showtooltip) {
        qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));
    }

    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabColor::qualitycurveMethodChanged));

    llCurveEditorG->setCurveListener(this);
    llshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "L(L)"));
    llshape->setResetCurve(DiagonalCurveType(LocallabParams::DEF_COLOR_LCURVE.at(0)), LocallabParams::DEF_COLOR_LCURVE);

    if (showtooltip) {
        llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    std::vector<GradientMilestone> mllshape;
    mllshape.push_back(GradientMilestone(0., 0., 0., 0.));
    mllshape.push_back(GradientMilestone(1., 1., 1., 1.));
    llshape->setBottomBarBgGradient(mllshape);
    llshape->setLeftBarBgGradient(mllshape);

    ccshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "C(C)"));
    ccshape->setResetCurve(DiagonalCurveType(LocallabParams::DEF_COLOR_LCURVE.at(0)), LocallabParams::DEF_COLOR_LCURVE);

    if (showtooltip) {
        ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    std::vector<GradientMilestone> mccshape;
    mccshape.push_back(GradientMilestone(0., 0., 0., 0.));
    mccshape.push_back(GradientMilestone(1., 1., 1., 1.));
    ccshape->setBottomBarBgGradient(mccshape);
    ccshape->setLeftBarBgGradient(mccshape);

    // llCurveEditorG->newLine();
    llCurveEditorG->curveListComplete();
    HCurveEditorG->setCurveListener(this);

    LHshape = static_cast<FlatCurveEditor*>(HCurveEditorG->addCurve(CT_Flat, "L(H)", nullptr, false, true));
    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FlatCurveType(LocallabParams::DEF_COLOR_HCURVE.at(0)), LocallabParams::DEF_COLOR_HCURVE);

    if (showtooltip) {
        LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    LHshape->setCurveColorProvider(this, 3);
    std::vector<GradientMilestone> mLHshape;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        mLHshape.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    LHshape->setBottomBarBgGradient(mLHshape);

    HHshape = static_cast<FlatCurveEditor*>(HCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true));
    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FlatCurveType(LocallabParams::DEF_COLOR_HCURVE.at(0)), LocallabParams::DEF_COLOR_HCURVE);

    if (showtooltip) {
        HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    }

    HHshape->setCurveColorProvider(this, 3);
    std::vector<GradientMilestone> mHHshape;

    for (int i = 0; i < 7; i++) {
        float x = float (i) * (1.0f / 6.0);

        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        mHHshape.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    HHshape->setBottomBarBgGradient(mHHshape);

    HCurveEditorG->curveListComplete();

    inversConn = invers->signal_toggled().connect(sigc::mem_fun(*this, &LocallabColor::inversChanged));

    Gtk::Frame* const superFrame = Gtk::manage(new Gtk::Frame());
    superFrame->set_label_align(0.025, 0.5);
    superFrame->set_label_widget(*curvactiv);
    ToolParamBlock* const superBox = Gtk::manage(new ToolParamBlock());
    superBox->pack_start(*lightness);
    superBox->pack_start(*contrast);
    superBox->pack_start(*chroma);
    gridFrame->set_label_align(0.025, 0.5);
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
    pack_start(*blurcolde);
    pack_start(*softradiuscol);
    Gtk::HBox* const qualcurvbox = Gtk::manage(new Gtk::HBox());
    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);
    pack_start(*qualcurvbox);
    pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    pack_start(*HCurveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    pack_start(*invers);
}

LocallabColor::~LocallabColor()
{
    delete llCurveEditorG;
    delete HCurveEditorG;
}

void LocallabColor::setListener(ToolPanelListener* tpl)
{
    LocallabTool::setListener(tpl);

    labgrid->setListener(tpl);
}

void LocallabColor::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
{
    colorMask = showMaskMethod->get_active_row_number();
}

void LocallabColor::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller)
{
    LocallabTool::colorForValue(valX, valY, elemType, callerId, caller); // Mask curves

    float R = 0.f;
    float G = 0.f;
    float B = 0.f;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 3) { // LHshape and HHshape curves
        Color::hsv2rgb01(float (valX), float (valY), 0.5f, R, G, B);

        caller->ccRed = double (R);
        caller->ccGreen = double (G);
        caller->ccBlue = double (B);
    }
}

void LocallabColor::disableListener()
{
    LocallabTool::disableListener();

    curvactivConn.block(true);
    gridMethodConn.block(true);
    qualitycurveMethodConn.block(true);
    inversConn.block(true);
}

void LocallabColor::enableListener()
{
    LocallabTool::enableListener();

    curvactivConn.block(false);
    gridMethodConn.block(false);
    qualitycurveMethodConn.block(false);
    inversConn.block(false);
}

void LocallabColor::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

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
        blurcolde->setValue(pp->locallab.spots.at(index).blurcolde);
        softradiuscol->setValue(pp->locallab.spots.at(index).softradiuscol);

        if (pp->locallab.spots.at(index).qualitycurveMethod == "none") {
            qualitycurveMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).qualitycurveMethod == "std") {
            qualitycurveMethod->set_active(1);
        }

        llshape->setCurve(pp->locallab.spots.at(index).llcurve);
        ccshape->setCurve(pp->locallab.spots.at(index).cccurve);
        LHshape->setCurve(pp->locallab.spots.at(index).LHcurve);
        HHshape->setCurve(pp->locallab.spots.at(index).HHcurve);
        invers->set_active(pp->locallab.spots.at(index).invers);
        enaMask->set_active(pp->locallab.spots.at(index).enaColorMask);
        CCMaskShape->setCurve(pp->locallab.spots.at(index).CCmaskcurve);
        LLMaskShape->setCurve(pp->locallab.spots.at(index).LLmaskcurve);
        HHMaskShape->setCurve(pp->locallab.spots.at(index).HHmaskcurve);
        blendMask->setValue(pp->locallab.spots.at(index).blendmaskcol);
        radMask->setValue(pp->locallab.spots.at(index).radmaskcol);
        chroMask->setValue(pp->locallab.spots.at(index).chromaskcol);
        gamMask->setValue(pp->locallab.spots.at(index).gammaskcol);
        sloMask->setValue(pp->locallab.spots.at(index).slomaskcol);
    }

    // Enable all listeners
    enableListener();

    // Update GUI according to invers button state
    updateColorGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expcolor = exp->getEnabled();
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
            pp->locallab.spots.at(pp->locallab.selspot).gridMethod = "one";
        } else if (gridMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(pp->locallab.selspot).gridMethod = "two";
        }

        pp->locallab.spots.at(index).strengthgrid = strengthgrid->getIntValue();
        pp->locallab.spots.at(index).sensi = sensi->getIntValue();
        pp->locallab.spots.at(index).structcol = structcol->getIntValue();
        pp->locallab.spots.at(index).blurcolde = blurcolde->getIntValue();
        pp->locallab.spots.at(index).softradiuscol = softradiuscol->getValue();

        if (qualitycurveMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).qualitycurveMethod = "none";
        } else if (qualitycurveMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).qualitycurveMethod = "std";
        }

        pp->locallab.spots.at(index).llcurve = llshape->getCurve();
        pp->locallab.spots.at(index).cccurve = ccshape->getCurve();
        pp->locallab.spots.at(index).LHcurve = LHshape->getCurve();
        pp->locallab.spots.at(index).HHcurve = HHshape->getCurve();
        pp->locallab.spots.at(index).invers = invers->get_active();
        pp->locallab.spots.at(index).enaColorMask = enaMask->get_active();
        pp->locallab.spots.at(index).CCmaskcurve = CCMaskShape->getCurve();
        pp->locallab.spots.at(index).LLmaskcurve = LLMaskShape->getCurve();
        pp->locallab.spots.at(index).HHmaskcurve = HHMaskShape->getCurve();
        pp->locallab.spots.at(index).blendmaskcol = blendMask->getIntValue();
        pp->locallab.spots.at(index).radmaskcol = radMask->getValue();
        pp->locallab.spots.at(index).chromaskcol = chroMask->getValue();
        pp->locallab.spots.at(index).gammaskcol = gamMask->getValue();
        pp->locallab.spots.at(index).slomaskcol = sloMask->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabColor::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster and labgrid widgets
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
        blendMask->setDefault((double)defSpot.blendmaskcol);
        radMask->setDefault(defSpot.radmaskcol);
        chroMask->setDefault(defSpot.chromaskcol);
        gamMask->setDefault(defSpot.gammaskcol);
        sloMask->setDefault(defSpot.slomaskcol);
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

        if (a == blendMask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskcol,
                                       blendMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radMask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskcol,
                                       radMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroMask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskcol,
                                       chroMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamMask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskcol,
                                       gamMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloMask) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskcol,
                                       sloMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

        if (ce == CCMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskshape,
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

void LocallabColor::enaMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMask->get_active()) {
                listener->panelChanged(EvLocallabEnaColorMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaColorMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabColor::showMaskMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
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

void LocallabColor::inversChanged()
{
    updateColorGUI(); // Update GUI according to invers button state

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

void LocallabColor::gridMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvLocallabgridMethod,
                                   gridMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
        }
    }
}

void LocallabColor::updateColorGUI()
{
    if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->show();
        HCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        labqualcurv->hide();
        maskExp->hide();
        structcol->hide();
        blurcolde->show();
        gridFrame->hide();
        strengthgrid->hide();
        softradiuscol->hide();

    } else {
        sensi->show();
        llCurveEditorG->show();
        HCurveEditorG->show();
        curvactiv->hide();
        qualitycurveMethod->show();
        labqualcurv->show();
        maskExp->show();
        structcol->show();
        blurcolde->show();
        gridFrame->show();
        softradiuscol->show();
    }
}

/* ==== LocallabExposure ==== */
LocallabExposure::LocallabExposure():
    LocallabTool(this, "Locallab Exposure", M("TP_LOCALLAB_EXPOSE"), false, MaskNormal),

    // Exposure specific widgets
    expMethod(Gtk::manage(new MyComboBoxText())),
    pdeFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_PDEFRA")))),
    laplacexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPLACEXP"), 0.0, 100.0, 0.1, 20.))),
    linear(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LINEAR"), 0., 1., 0.01, 0.))),
    balanexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BALANEXP"), 0.2, 1.2, 0.01, 0.8))),
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -2.0, 4.0, 0.05, 0.0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shadex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHADEX"), 0, 100, 1, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHADEXCOMP"), 0, 100, 1, 50))),
    expchroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_EXPCHROMA"), -50, 100, 1, 30))),
    warm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., Gtk::manage(new RTImage("circle-blue-small.png")), Gtk::manage(new RTImage("circle-orange-small.png"))))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    structexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUCCOL"), 0, 100, 1, 0))),
    blurexpde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    softradiusexp(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOFTRADIUSCOL"), 0.0, 100.0, 0.1, 0.))),
    curveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"))),
    inversex(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    const bool showtooltip = options.showtooltip;

    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_EXPOSURE_TOOLTIP"));
    }

    expMethod->append(M("TP_LOCALLAB_STD"));
    expMethod->append(M("TP_LOCALLAB_PDE"));
    expMethod->set_active(0);
    expMethodConn  = expMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabExposure::expMethodChanged));

    if (showtooltip) {
        expMethod->set_tooltip_text(M("TP_LOCALLAB_EXPMETHOD_TOOLTIP"));
    }

    pdeFrame->set_label_align(0.025, 0.5);

    laplacexp->setAdjusterListener(this);

    linear->setAdjusterListener(this);

    balanexp->setAdjusterListener(this);

    expcomp->setAdjusterListener(this);

    hlcompr->setAdjusterListener(this);

    hlcomprthresh->setAdjusterListener(this);

    black->setAdjusterListener(this);

    shadex->setAdjusterListener(this);

    shcompr->setAdjusterListener(this);

    expchroma->setAdjusterListener(this);

    if (showtooltip) {
        warm->set_tooltip_text(M("TP_LOCALLAB_WARM_TOOLTIP"));
    }

    warm->setAdjusterListener(this);

    if (showtooltip) {
        sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    }

    sensiex->setAdjusterListener(this);

    structexp->setAdjusterListener(this);

    blurexpde->setAdjusterListener(this);

    softradiusexp->setAdjusterListener(this);

    curveEditorG->setCurveListener(this);

    shapeexpos = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""));
    shapeexpos->setResetCurve(DiagonalCurveType(LocallabParams::DEF_EXP_CURVE.at(0)), LocallabParams::DEF_EXP_CURVE);

    if (showtooltip) {
        shapeexpos->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
    }

    std::vector<GradientMilestone> mshapeexpos;
    mshapeexpos.push_back(GradientMilestone(0., 0., 0., 0.));
    mshapeexpos.push_back(GradientMilestone(1., 1., 1., 1.));
    shapeexpos->setBottomBarBgGradient(mshapeexpos);
    shapeexpos->setLeftBarBgGradient(mshapeexpos);

    curveEditorG->curveListComplete();

    inversexConn  = inversex->signal_toggled().connect(sigc::mem_fun(*this, &LocallabExposure::inversexChanged));

    pack_start(*expMethod);
    ToolParamBlock* const pdeBox = Gtk::manage(new ToolParamBlock());
    pdeBox->pack_start(*laplacexp);
    pdeBox->pack_start(*linear);
    pdeBox->pack_start(*balanexp);
    pdeFrame->add(*pdeBox);
    pack_start(*pdeFrame);
    pack_start(*expcomp);
    pack_start(*hlcompr);
    pack_start(*hlcomprthresh);
    pack_start(*black);
    pack_start(*shadex);
    pack_start(*shcompr);
    pack_start(*expchroma);
    pack_start(*warm);
    pack_start(*sensiex);
    pack_start(*structexp);
    pack_start(*blurexpde);
    pack_start(*softradiusexp);
    pack_start(*curveEditorG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
    pack_start(*inversex);
}

LocallabExposure::~LocallabExposure()
{
    delete curveEditorG;
}

void LocallabExposure::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
{
    expMask = showMaskMethod->get_active_row_number();
}

void LocallabExposure::disableListener()
{
    LocallabTool::disableListener();

    expMethodConn.block(true);
    inversexConn.block(true);
}

void LocallabExposure::enableListener()
{
    LocallabTool::enableListener();

    expMethodConn.block(false);
    inversexConn.block(false);
}

void LocallabExposure::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->setEnabled(pp->locallab.spots.at(index).expexpose);

        if (pp->locallab.spots.at(index).expMethod == "std") {
            expMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).expMethod == "pde") {
            expMethod->set_active(1);
        }

        laplacexp->setValue(pp->locallab.spots.at(index).laplacexp);
        linear->setValue(pp->locallab.spots.at(index).linear);
        balanexp->setValue(pp->locallab.spots.at(index).balanexp);
        expcomp->setValue(pp->locallab.spots.at(index).expcomp);
        hlcompr->setValue(pp->locallab.spots.at(index).hlcompr);
        hlcomprthresh->setValue(pp->locallab.spots.at(index).hlcomprthresh);
        black->setValue(pp->locallab.spots.at(index).black);
        shadex->setValue(pp->locallab.spots.at(index).shadex);
        shcompr->setValue(pp->locallab.spots.at(index).shcompr);
        expchroma->setValue(pp->locallab.spots.at(index).expchroma);
        warm->setValue(pp->locallab.spots.at(index).warm);
        sensiex->setValue(pp->locallab.spots.at(index).sensiex);
        structexp->setValue(pp->locallab.spots.at(index).structexp);
        blurexpde->setValue(pp->locallab.spots.at(index).blurexpde);
        shapeexpos->setCurve(pp->locallab.spots.at(index).excurve);
        softradiusexp->setValue(pp->locallab.spots.at(index).softradiusexp);
        inversex->set_active(pp->locallab.spots.at(index).inversex);
        enaMask->set_active(pp->locallab.spots.at(index).enaExpMask);
        CCMaskShape->setCurve(pp->locallab.spots.at(index).CCmaskexpcurve);
        LLMaskShape->setCurve(pp->locallab.spots.at(index).LLmaskexpcurve);
        HHMaskShape->setCurve(pp->locallab.spots.at(index).HHmaskexpcurve);
        blendMask->setValue(pp->locallab.spots.at(index).blendmaskexp);
        radMask->setValue(pp->locallab.spots.at(index).radmaskexp);
        chroMask->setValue(pp->locallab.spots.at(index).chromaskexp);
        gamMask->setValue(pp->locallab.spots.at(index).gammaskexp);
        sloMask->setValue(pp->locallab.spots.at(index).slomaskexp);
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

        if (expMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).expMethod = "std";
        } else if (expMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).expMethod = "pde";
        }

        pp->locallab.spots.at(index).laplacexp = laplacexp->getValue();
        pp->locallab.spots.at(index).linear = linear->getValue();
        pp->locallab.spots.at(index).balanexp = balanexp->getValue();
        pp->locallab.spots.at(index).expcomp = expcomp->getValue();
        pp->locallab.spots.at(index).hlcompr = hlcompr->getIntValue();
        pp->locallab.spots.at(index).hlcomprthresh = hlcomprthresh->getIntValue();
        pp->locallab.spots.at(index).black = black->getIntValue();
        pp->locallab.spots.at(index).shadex = shadex->getIntValue();
        pp->locallab.spots.at(index).shcompr = shcompr->getIntValue();
        pp->locallab.spots.at(index).expchroma = expchroma->getIntValue();
        pp->locallab.spots.at(index).warm = warm->getIntValue();
        pp->locallab.spots.at(index).sensiex = sensiex->getIntValue();
        pp->locallab.spots.at(index).structexp = structexp->getIntValue();
        pp->locallab.spots.at(index).blurexpde = blurexpde->getIntValue();
        pp->locallab.spots.at(index).softradiusexp = softradiusexp->getValue();
        pp->locallab.spots.at(index).excurve = shapeexpos->getCurve();
        pp->locallab.spots.at(index).inversex = inversex->get_active();
        pp->locallab.spots.at(index).enaExpMask = enaMask->get_active();
        pp->locallab.spots.at(index).LLmaskexpcurve = LLMaskShape->getCurve();
        pp->locallab.spots.at(index).CCmaskexpcurve = CCMaskShape->getCurve();
        pp->locallab.spots.at(index).HHmaskexpcurve = HHMaskShape->getCurve();
        pp->locallab.spots.at(index).blendmaskexp = blendMask->getIntValue();
        pp->locallab.spots.at(index).radmaskexp = radMask->getValue();
        pp->locallab.spots.at(index).chromaskexp = chroMask->getValue();
        pp->locallab.spots.at(index).gammaskexp = gamMask->getValue();
        pp->locallab.spots.at(index).slomaskexp = sloMask->getValue();
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
        expcomp->setDefault(defSpot.expcomp);
        hlcompr->setDefault((double)defSpot.hlcompr);
        hlcomprthresh->setDefault((double)defSpot.hlcomprthresh);
        black->setDefault((double)defSpot.black);
        shadex->setDefault((double)defSpot.shadex);
        shcompr->setDefault((double)defSpot.shcompr);
        expchroma->setDefault((double)defSpot.expchroma);
        warm->setDefault((double)defSpot.warm);
        sensiex->setDefault((double)defSpot.sensiex);
        structexp->setDefault((double)defSpot.structexp);
        blurexpde->setDefault((double)defSpot.blurexpde);
        softradiusexp->setDefault(defSpot.softradiusexp);
        blendMask->setDefault((double)defSpot.blendmaskexp);
        radMask->setDefault(defSpot.radmaskexp);
        chroMask->setDefault(defSpot.chromaskexp);
        gamMask->setDefault(defSpot.gammaskexp);
        sloMask->setDefault(defSpot.slomaskexp);
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

        if (a == expcomp) {
            if (listener) {
                listener->panelChanged(Evlocallabexpcomp,
                                       expcomp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

        if (a == black) {
            if (listener) {
                listener->panelChanged(Evlocallabblack,
                                       black->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

        if (a == warm) {
            if (listener) {
                listener->panelChanged(Evlocallabwarm,
                                       warm->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

        if (a == softradiusexp) {
            if (listener) {
                listener->panelChanged(Evlocallabsoftradiusexp,
                                       softradiusexp->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == blendMask) {
            if (listener) {
                listener->panelChanged(Evlocallabblendmaskexp,
                                       blendMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radMask) {
            if (listener) {
                listener->panelChanged(Evlocallabradmaskexp,
                                       radMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroMask) {
            if (listener) {
                listener->panelChanged(Evlocallabchromaskexp,
                                       chroMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamMask) {
            if (listener) {
                listener->panelChanged(Evlocallabgammaskexp,
                                       gamMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloMask) {
            if (listener) {
                listener->panelChanged(Evlocallabslomaskexp,
                                       sloMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

        if (ce == CCMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskexpshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskexpshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskexpshape,
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

void LocallabExposure::enaMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMask->get_active()) {
                listener->panelChanged(EvLocallabEnaExpMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaExpMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabExposure::showMaskMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
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

void LocallabExposure::inversexChanged()
{
    // Update exposure GUI according to inversex button state
    updateExposureGUI3();

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
        balanexp->set_sensitive(false);
        linear->set_sensitive(false);
    } else if (expMethod->get_active_row_number() == 1) {
        pdeFrame->set_sensitive(true);
        laplacexp->set_sensitive(true);
        balanexp->set_sensitive(true);
        linear->set_sensitive(true);
    }
}

void LocallabExposure::updateExposureGUI3()
{
    // Update exposure GUI according to inversex button state
    if (inversex->get_active()) {
        maskExp->hide();
        structexp->hide();
        softradiusexp->hide();
        shadex->hide();
        expMethod->hide();
        pdeFrame->hide();
    } else {
        maskExp->show();
        structexp->show();
        softradiusexp->show();
        shadex->show();
        expMethod->show();
        pdeFrame->show();
    }
}

/* ==== LocallabShadow ==== */
LocallabShadow::LocallabShadow():
    LocallabTool(this, "Locallab Shadows Highlight", M("TP_LOCALLAB_SHADHIGH"), false, MaskNormal),

    // Shadow highlight specific widgets
    highlights(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HIGHLIGHTS"), 0, 100, 1, 0))),
    h_tonalwidth(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_HLTONALW"), 10, 100, 1, 70))),
    shadows(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHADOWS"), 0, 100, 1, 0))),
    s_tonalwidth(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_SHTONALW"), 10, 100, 1, 30))),
    sh_radius(Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_RADIUS"), 0, 100, 1, 40))),
    sensihs(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    blurSHde(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BLURDE"), 2, 100, 1, 5))),
    inverssh(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS"))))
{
    const bool showtooltip = options.showtooltip;

    if (showtooltip) {
        exp->set_tooltip_text(M("TP_LOCALLAB_SHADOWHIGHLIGHT_TOOLTIP"));
    }

    highlights->setAdjusterListener(this);

    h_tonalwidth->setAdjusterListener(this);

    shadows->setAdjusterListener(this);

    s_tonalwidth->setAdjusterListener(this);

    sh_radius->setAdjusterListener(this);

    sensihs->setAdjusterListener(this);

    blurSHde->setAdjusterListener(this);

    inversshConn  = inverssh->signal_toggled().connect(sigc::mem_fun(*this, &LocallabShadow::inversshChanged));

    pack_start(*highlights);
    pack_start(*h_tonalwidth);
    pack_start(*shadows);
    pack_start(*s_tonalwidth);
    pack_start(*sh_radius);
    pack_start(*sensihs);
    pack_start(*blurSHde);
    pack_start(*inverssh);
}

void LocallabShadow::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
{
    shMask = showMaskMethod->get_active_row_number();
}

void LocallabShadow::disableListener()
{
    LocallabTool::disableListener();

    inversshConn.block(true);
}

void LocallabShadow::enableListener()
{
    LocallabTool::enableListener();

    inversshConn.block(true);
}

void LocallabShadow::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->setEnabled(pp->locallab.spots.at(index).expshadhigh);
        highlights->setValue(pp->locallab.spots.at(index).highlights);
        h_tonalwidth->setValue(pp->locallab.spots.at(index).h_tonalwidth);
        shadows->setValue(pp->locallab.spots.at(index).shadows);
        s_tonalwidth->setValue(pp->locallab.spots.at(index).s_tonalwidth);
        sh_radius->setValue(pp->locallab.spots.at(index).sh_radius);
        sensihs->setValue(pp->locallab.spots.at(index).sensihs);
        blurSHde->setValue(pp->locallab.spots.at(index).blurSHde);
        inverssh->set_active(pp->locallab.spots.at(index).inverssh);
        enaMask->set_active(pp->locallab.spots.at(index).enaSHMask);
        CCMaskShape->setCurve(pp->locallab.spots.at(index).CCmaskSHcurve);
        LLMaskShape->setCurve(pp->locallab.spots.at(index).LLmaskSHcurve);
        HHMaskShape->setCurve(pp->locallab.spots.at(index).HHmaskSHcurve);
        blendMask->setValue(pp->locallab.spots.at(index).blendmaskSH);
        radMask->setValue(pp->locallab.spots.at(index).radmaskSH);
        chroMask->setValue(pp->locallab.spots.at(index).chromaskSH);
        gamMask->setValue(pp->locallab.spots.at(index).gammaskSH);
        sloMask->setValue(pp->locallab.spots.at(index).slomaskSH);
    }

    // Enable all listeners
    enableListener();

    // Update shadow highlight GUI according to inverssh button state
    updateShadowGUI();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expshadhigh = exp->getEnabled();
        pp->locallab.spots.at(index).highlights = highlights->getIntValue();
        pp->locallab.spots.at(index).h_tonalwidth = h_tonalwidth->getIntValue();
        pp->locallab.spots.at(index).shadows = shadows->getIntValue();
        pp->locallab.spots.at(index).s_tonalwidth = s_tonalwidth->getIntValue();
        pp->locallab.spots.at(index).sh_radius = sh_radius->getIntValue();
        pp->locallab.spots.at(index).sensihs = sensihs->getIntValue();
        pp->locallab.spots.at(index).blurSHde = blurSHde->getIntValue();
        pp->locallab.spots.at(index).inverssh = inverssh->get_active();
        pp->locallab.spots.at(index).enaSHMask = enaMask->get_active();
        pp->locallab.spots.at(index).LLmaskSHcurve = LLMaskShape->getCurve();
        pp->locallab.spots.at(index).CCmaskSHcurve = CCMaskShape->getCurve();
        pp->locallab.spots.at(index).HHmaskSHcurve = HHMaskShape->getCurve();
        pp->locallab.spots.at(index).blendmaskSH = blendMask->getIntValue();
        pp->locallab.spots.at(index).radmaskSH = radMask->getValue();
        pp->locallab.spots.at(index).chromaskSH = chroMask->getValue();
        pp->locallab.spots.at(index).gammaskSH = gamMask->getValue();
        pp->locallab.spots.at(index).slomaskSH = sloMask->getValue();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default values for adjuster widgets
        highlights->setDefault((double)defSpot.highlights);
        h_tonalwidth->setDefault((double)defSpot.h_tonalwidth);
        shadows->setDefault((double)defSpot.shadows);
        s_tonalwidth->setDefault((double)defSpot.s_tonalwidth);
        sh_radius->setDefault((double)defSpot.sh_radius);
        sensihs->setDefault((double)defSpot.sensihs);
        blurSHde->setDefault((double)defSpot.blurSHde);
        blendMask->setDefault((double)defSpot.blendmaskSH);
        radMask->setDefault(defSpot.radmaskSH);
        chroMask->setDefault(defSpot.chromaskSH);
        gamMask->setDefault(defSpot.gammaskSH);
        sloMask->setDefault(defSpot.slomaskSH);
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabShadow::adjusterChanged(Adjuster* a, double newval)
{
    if (isLocActivated && exp->getEnabled()) {
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

        if (a == blendMask) {
            if (listener) {
                listener->panelChanged(EvlocallabblendmaskSH,
                                       blendMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == radMask) {
            if (listener) {
                listener->panelChanged(EvlocallabradmaskSH,
                                       radMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == chroMask) {
            if (listener) {
                listener->panelChanged(EvlocallabchromaskSH,
                                       chroMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == gamMask) {
            if (listener) {
                listener->panelChanged(EvlocallabgammaskSH,
                                       gamMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (a == sloMask) {
            if (listener) {
                listener->panelChanged(EvlocallabslomaskSH,
                                       sloMask->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::curveChanged(CurveEditor* ce)
{
    if (isLocActivated && exp->getEnabled()) {
        if (ce == CCMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabCCmaskSHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == LLMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabLLmaskSHshape,
                                       M("HISTORY_CUSTOMCURVE") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }

        if (ce == HHMaskShape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHmaskSHshape,
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

void LocallabShadow::enaMaskChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            if (enaMask->get_active()) {
                listener->panelChanged(EvLocallabEnaSHMask,
                                       M("GENERAL_ENABLED") + " (" + escapeHtmlChars(spotName) + ")");
            } else {
                listener->panelChanged(EvLocallabEnaSHMask,
                                       M("GENERAL_DISABLED") + " (" + escapeHtmlChars(spotName) + ")");
            }
        }
    }
}

void LocallabShadow::showMaskMethodChanged()
{
    // If mask preview is activated, deactivate all other tool mask preview
    if (locToolListener) {
        locToolListener->resetOtherMaskView(this);
    }

    if (listener) {
        listener->panelChanged(EvlocallabshowmaskMethod, "");
    }
}

void LocallabShadow::inversshChanged()
{
    // Update shadow highlight GUI according to inverssh button state
    updateShadowGUI();

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

void LocallabShadow::updateShadowGUI()
{
    // Update shadow highlight GUI according to inverssh button state
    if (inverssh->get_active()) {
        maskExp->hide();
    } else {
        maskExp->show();
    }
}

/* ==== LocallabVibrance ==== */
LocallabVibrance::LocallabVibrance():
    LocallabTool(this, "Locallab Vibrance", M("TP_LOCALLAB_VIBRANCE"), false, MaskNone),

    // Vibrance specific widgets
    saturated(Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.))),
    pastels(Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.))),
    psThreshold(Gtk::manage(new ThresholdAdjuster(M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false))),
    protectSkins(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PROTECTSKINS")))),
    avoidColorShift(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_AVOIDCOLORSHIFT")))),
    pastSatTog(Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PASTSATTOG")))),
    sensiv(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 15))),
    curveEditorGG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL")))
{
    float R, G, B;

    const bool showtooltip = options.showtooltip;

    saturated->setAdjusterListener(this);

    pastels->setAdjusterListener(this);

    if (showtooltip) {
        psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    }

    psThreshold->setAdjusterListener(this);

    pskinsConn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::protectskins_toggled));

    ashiftConn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::avoidcolorshift_toggled));

    pastsattogConn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &LocallabVibrance::pastsattog_toggled));

    sensiv->setAdjusterListener(this);

    curveEditorGG->setCurveListener(this);

    skinTonesCurve = static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")));

    if (showtooltip) {
        skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
    }

    std::vector<GradientMilestone> mskinTonesCurve;
    // -0.1 rad < Hue < 1.6 rad
    Color::hsv2rgb01(0.92f, 0.45f, 0.6f, R, G, B);
    mskinTonesCurve.push_back(GradientMilestone(0.0, double (R), double (G), double (B)));
    Color::hsv2rgb01(0.14056f, 0.45f, 0.6f, R, G, B);
    mskinTonesCurve.push_back(GradientMilestone(1.0, double (R), double (G), double (B)));
    skinTonesCurve->setBottomBarBgGradient(mskinTonesCurve);
    skinTonesCurve->setLeftBarBgGradient(mskinTonesCurve);
    skinTonesCurve->setRangeLabels(
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE1"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE2"),
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE3"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE4")
    );
    skinTonesCurve->setRangeDefaultMilestones(0.1, 0.4, 0.85);

    curveEditorGG->curveListComplete();

    pack_start(*saturated, Gtk::PACK_SHRINK, 0);
    pack_start(*pastels, Gtk::PACK_SHRINK, 0);
    pack_start(*psThreshold, Gtk::PACK_SHRINK, 0);
    pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);
    pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);
    pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);
    pack_start(*sensiv, Gtk::PACK_SHRINK, 0);
    pack_start(*curveEditorGG, Gtk::PACK_SHRINK, 4); // Padding is mandatory to correct behavior of curve editor
}

LocallabVibrance::~LocallabVibrance()
{
    delete curveEditorGG;
}

void LocallabVibrance::disableListener()
{
    LocallabTool::disableListener();

    pskinsConn.block(true);
    ashiftConn.block(true);
    pastsattogConn.block(true);
}

void LocallabVibrance::enableListener()
{
    LocallabTool::enableListener();

    pskinsConn.block(false);
    ashiftConn.block(false);
    pastsattogConn.block(false);
}

void LocallabVibrance::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->setEnabled(pp->locallab.spots.at(index).expvibrance);
        saturated->setValue(pp->locallab.spots.at(index).saturated);
        pastels->setValue(pp->locallab.spots.at(index).pastels);
        psThreshold->setValue<int>(pp->locallab.spots.at(index).psthreshold);
        protectSkins->set_active(pp->locallab.spots.at(index).protectskins);
        avoidColorShift->set_active(pp->locallab.spots.at(index).avoidcolorshift);
        pastSatTog->set_active(pp->locallab.spots.at(index).pastsattog);
        sensiv->setValue(pp->locallab.spots.at(index).sensiv);
        skinTonesCurve->setCurve(pp->locallab.spots.at(index).skintonescurve);
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
        pp->locallab.spots.at(index).saturated = saturated->getIntValue();
        pp->locallab.spots.at(index).pastels = pastels->getIntValue();
        pp->locallab.spots.at(index).psthreshold = psThreshold->getValue<int>();
        pp->locallab.spots.at(index).protectskins = protectSkins->get_active();
        pp->locallab.spots.at(index).avoidcolorshift = avoidColorShift->get_active();
        pp->locallab.spots.at(index).pastsattog = pastSatTog->get_active();
        pp->locallab.spots.at(index).sensiv = sensiv->getIntValue();
        pp->locallab.spots.at(index).skintonescurve = skinTonesCurve->getCurve();
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
        psThreshold->setDefault<int>(defSpot.psthreshold);
        sensiv->setDefault((double)defSpot.sensiv);
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

        if (a == sensiv) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiv,
                                       sensiv->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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
    LocallabTool(this, "Locallab Soft Light", M("TP_LOCALLAB_SOFT"), false, MaskNone),

    // Soft light specific widgets
    softMethod(Gtk::manage(new MyComboBoxText())),
    ctboxsoftmethod(Gtk::manage(new Gtk::HBox())),
    showmasksoftMethod(Gtk::manage(new MyComboBoxText())),
    streng(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENG"), 1, 100, 1, 1))),
    laplace(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LAPLACE"), 0., 100., 0.5, 25.))),
    sensisf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 1, 100, 1, 15)))
{
    const bool showtooltip = options.showtooltip;

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
    showmasksoftMethod->set_active(0);

    if (showtooltip) {
        showmasksoftMethod->set_tooltip_markup(M("TP_LOCALLAB_SHOWMASKSOFT_TOOLTIP"));
    }

    showmasksoftMethodConn  = showmasksoftMethod->signal_changed().connect(sigc::mem_fun(*this, &LocallabSoft::showmasksoftMethodChanged));

    streng->setAdjusterListener(this);

    laplace->setAdjusterListener(this);

    sensisf->setAdjusterListener(this);

    pack_start(*softMethod);
    Gtk::Label* const labelsoftmethod = Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_SHOWDCT") + ":"));
    ctboxsoftmethod->pack_start(*labelsoftmethod, Gtk::PACK_SHRINK, 4);
    ctboxsoftmethod->pack_start(*showmasksoftMethod);
    pack_start(*ctboxsoftmethod);
    pack_start(*streng);
    pack_start(*laplace);
    pack_start(*sensisf);
}

void LocallabSoft::getMaskView(int &colorMask, int &expMask, int &shMask, int &softMask, int &tmMask, int &retiMask, int &cbMask)
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
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->setEnabled(pp->locallab.spots.at(index).expsoft);

        if (pp->locallab.spots.at(index).softMethod == "soft") {
            softMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).softMethod == "reti") {
            softMethod->set_active(1);
        }

        streng->setValue(pp->locallab.spots.at(index).streng);
        sensisf->setValue(pp->locallab.spots.at(index).sensisf);
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

void LocallabSoft::resetMaskView()
{
    showmasksoftMethodConn.block(true);
    showmasksoftMethod->set_active(0);
    showmasksoftMethodConn.block(false);
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
        laplace->hide();
        ctboxsoftmethod->hide();
    } else {
        laplace->show();
        ctboxsoftmethod->show();
    }
}

/* ==== LocallabBlur ==== */
LocallabBlur::LocallabBlur():
    LocallabTool(this, "Locallab Blur & Noise", M("TP_LOCALLAB_BLUFR"), false, MaskNone),

    // Blur & Noise specific widgets
    radius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADIUS"), 1.0, 100.0, 0.1, 1.0))),
    strength(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0))),
    sensibn(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 40))),
    blurMethod(Gtk::manage(new MyComboBoxText())),
    activlum(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ACTIV"))))
{
    const bool showtooltip = options.showtooltip;

    radius->setAdjusterListener(this);

    strength->setAdjusterListener(this);

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

    activlumConn  = activlum->signal_toggled().connect(sigc::mem_fun(*this, &LocallabBlur::activlumChanged));

    pack_start(*radius);
    pack_start(*strength);
    pack_start(*sensibn);
    pack_start(*blurMethod);
    pack_start(*activlum);
}

void LocallabBlur::disableListener()
{
    LocallabTool::disableListener();

    blurMethodConn.block(true);
    activlumConn.block(true);
}

void LocallabBlur::enableListener()
{
    LocallabTool::enableListener();

    blurMethodConn.block(false);
    activlumConn.block(false);
}

void LocallabBlur::read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited)
{
    // Disable all listeners
    disableListener();

    // Update GUI to selected spot value
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        spotName = pp->locallab.spots.at(index).name; // Update spot name according to selected spot

        exp->setEnabled(pp->locallab.spots.at(index).expblur);
        radius->setValue(pp->locallab.spots.at(index).radius);
        strength->setValue(pp->locallab.spots.at(index).strength);
        sensibn->setValue(pp->locallab.spots.at(index).sensibn);

        if (pp->locallab.spots.at(index).blurMethod == "norm") {
            blurMethod->set_active(0);
        } else if (pp->locallab.spots.at(index).blurMethod == "inv") {
            blurMethod->set_active(1);
        }

        activlum->set_active(pp->locallab.spots.at(index).activlum);
    }

    // Enable all listeners
    enableListener();

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    const int index = pp->locallab.selspot;

    if (index < (int)pp->locallab.spots.size()) {
        pp->locallab.spots.at(index).expsoft = exp->getEnabled();

        pp->locallab.spots.at(index).expblur = exp->getEnabled();
        pp->locallab.spots.at(index).radius = radius->getValue();
        pp->locallab.spots.at(index).strength = strength->getIntValue();
        pp->locallab.spots.at(index).sensibn = sensibn->getIntValue();


        if (blurMethod->get_active_row_number() == 0) {
            pp->locallab.spots.at(index).blurMethod = "norm";
        } else if (blurMethod->get_active_row_number() == 1) {
            pp->locallab.spots.at(index).blurMethod = "inv";
        }

        pp->locallab.spots.at(index).activlum = activlum->get_active();
    }

    // Note: No need to manage pedited as batch mode is deactivated for Locallab
}

void LocallabBlur::setDefaults(const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited)
{
    const int index = defParams->locallab.selspot;

    if (index < (int)defParams->locallab.spots.size()) {
        const LocallabParams::LocallabSpot defSpot = defParams->locallab.spots.at(index);

        // Set default value for adjuster widgets
        radius->setDefault(defSpot.radius);
        strength->setDefault((double)defSpot.strength);
        sensibn->setDefault((double)defSpot.sensibn);
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

        if (a == sensibn) {
            if (listener) {
                listener->panelChanged(Evlocallabsensibn,
                                       sensibn->getTextValue() + " (" + escapeHtmlChars(spotName) + ")");
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

void LocallabBlur::blurMethodChanged()
{
    if (isLocActivated && exp->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabblurMethod,
                                   blurMethod->get_active_text() + " (" + escapeHtmlChars(spotName) + ")");
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
