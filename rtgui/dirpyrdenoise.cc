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
#include "dirpyrdenoise.h"
#include <iomanip>
#include <cmath>
#include "edit.h"
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;
extern Options options;

DirPyrDenoise::DirPyrDenoise () : FoldableToolPanel(this, "dirpyrdenoise", M("TP_DIRPYRDENOISE_LABEL"), true, true), lastenhance(false)
{
    std::vector<GradientMilestone> milestones;
    CurveListener::setMulti(true);
    nextnresid = 0.;
    nexthighresid = 0.;
    nextchroma = 15.;
    nextred = 0.;
    nextblue = 0.;

    std::vector<double> defaultCurve;

    Gtk::Frame* lumaFrame = Gtk::manage (new Gtk::Frame (M("TP_DIRPYRDENOISE_LUMAFR")) );
    lumaFrame->set_label_align(0.025, 0.5);

    Gtk::VBox * lumaVBox = Gtk::manage ( new Gtk::VBox());
    lumaVBox->set_spacing(2);



    ctboxL = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labmL = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_LTYPE") + ":"));
    ctboxL->pack_start (*labmL, Gtk::PACK_SHRINK, 1);

    Lmethod = Gtk::manage (new MyComboBoxText ());
    Lmethod->append (M("TP_DIRPYRDENOISE_CUR"));
    Lmethod->append (M("TP_DIRPYRDENOISE_SLI"));
    Lmethod->set_active(0);
    Lmethodconn = Lmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::LmethodChanged) );

    luma  = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_LUMA"), 0, 100, 0.01, 0));
    Ldetail  = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_LDETAIL"), 0, 100, 0.01, 50));
    NoiscurveEditorG = new CurveEditorGroup (options.lastDenoiseCurvesDir, M("TP_DIRPYRDENOISE_LCURVE"));
    //curveEditorG = new CurveEditorGroup (options.lastLabCurvesDir);
    NoiscurveEditorG->setCurveListener (this);
    rtengine::DirPyrDenoiseParams::getDefaultNoisCurve(defaultCurve);
    lshape = static_cast<FlatCurveEditor*>(NoiscurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));
    lshape->setIdentityValue(0.);
    lshape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);

    lshape->setTooltip(M("TP_DIRPYRDENOISE_CURVEEDITOR_L_TOOLTIP"));
    //lshape->setEditID(EUID_Lab_LCurve, BT_SINGLEPLANE_FLOAT);
    milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones.push_back( GradientMilestone(1., 1., 1., 1.) );
    lshape->setBottomBarBgGradient(milestones);
    //lshape->setLeftBarBgGradient(milestones);
    milestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    milestones.push_back( GradientMilestone(1., 1., 1., 1.) );
    NoiscurveEditorG->curveListComplete();
    NoiscurveEditorG->show();

    Gtk::Frame* chromaFrame = Gtk::manage (new Gtk::Frame (M("TP_DIRPYRDENOISE_CHROMAFR")) );
    chromaFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *chromaVBox = Gtk::manage ( new Gtk::VBox());
    chromaVBox->set_spacing(2);

    autochroma = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYRDENOISE_AUTO")));
    autochroma->set_active (true);
    autochroma->set_tooltip_text (M("TP_DIRPYRDENOISE_AUTO_TOOLTIP"));

    ctboxC = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labmC = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_CTYPE") + ":"));
    ctboxC->pack_start (*labmC, Gtk::PACK_SHRINK, 1);
    ctboxC->set_tooltip_markup (M("TP_DIRPYRDENOISE_CTYPE_TOOLTIP"));

    Cmethod = Gtk::manage (new MyComboBoxText ());
    Cmethod->append (M("TP_DIRPYRDENOISE_MAN"));
    Cmethod->append (M("TP_DIRPYRDENOISE_AUT"));
    Cmethod->append (M("TP_DIRPYRDENOISE_PON"));
    Cmethod->append (M("TP_DIRPYRDENOISE_PRE"));
    Cmethod->set_active(0);
    Cmethodconn = Cmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::CmethodChanged) );

    ctboxC2 = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labmC2 = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_CTYPE") + ":"));
    ctboxC2->pack_start (*labmC2, Gtk::PACK_SHRINK, 1);
    ctboxC2->set_tooltip_markup (M("TP_DIRPYRDENOISE_C2TYPE_TOOLTIP"));

    C2method = Gtk::manage (new MyComboBoxText ());
    C2method->append (M("TP_DIRPYRDENOISE_MANU"));
    C2method->append (M("TP_DIRPYRDENOISE_AUTO"));
    C2method->append (M("TP_DIRPYRDENOISE_PREV"));
    C2method->set_active(0);
    C2methodconn = C2method->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::C2methodChanged) );


    NoiseLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    NoiseLabels->set_tooltip_text(M("TP_DIRPYRDENOISE_NRESID_TOOLTIP"));

    TileLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
    PrevLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));

    chroma    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_CHROMA"), 0, 100, 0.01, 15));
    redchro    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_RED"), -100, 100, 0.1, 0));
    bluechro    = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_BLUE"), -100, 100, 0.1, 0));

    gamma   = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_GAMMA"), 1.0, 3.0, 0.01, 1.7));
    gamma->set_tooltip_text (M("TP_DIRPYRDENOISE_GAMMA_TOOLTIP"));


    Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
    hb1->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_DIRPYRDENOISE_METHOD") + ": ")), Gtk::PACK_SHRINK, 4);
    hb1->set_tooltip_markup (M("TP_DIRPYRDENOISE_METHOD_TOOLTIP"));

    dmethod = Gtk::manage (new MyComboBoxText ());
    dmethod->append (M("TP_DIRPYRDENOISE_LAB"));
    dmethod->append (M("TP_DIRPYRDENOISE_RGB"));
    dmethod->set_active(0);
    hb1->pack_end (*dmethod, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb1, Gtk::PACK_SHRINK, 4);


    dmethodconn = dmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::dmethodChanged) );

    luma->setAdjusterListener (this);
    Ldetail->setAdjusterListener (this);
    chroma->setAdjusterListener (this);
    redchro->setAdjusterListener (this);
    bluechro->setAdjusterListener (this);

    CCcurveEditorG = new CurveEditorGroup (options.lastDenoiseCurvesDir, M("TP_DIRPYRDENOISE_CCCURVE"));
    CCcurveEditorG->setCurveListener (this);
    rtengine::DirPyrDenoiseParams::getDefaultCCCurve(defaultCurve);
    ccshape = static_cast<FlatCurveEditor*>(CCcurveEditorG->addCurve(CT_Flat, "", nullptr, false, false));
    ccshape->setIdentityValue(0.);
    ccshape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);

    ccshape->setTooltip(M("TP_DIRPYRDENOISE_CURVEEDITOR_CC_TOOLTIP"));
    ccshape->setBottomBarColorProvider(this, 2);

    CCcurveEditorG->curveListComplete();


    //-----------------------------------------

    gamma->setAdjusterListener (this);

    luma->hide();
    Ldetail->show();

//  autochroma->show();
    NoiseLabels->show();
    TileLabels->show();
    PrevLabels->show();
    chroma->show();
    redchro->show();
    bluechro->show();
//  perform->show();
    gamma->show();
//  perform->set_active (true);

    enhance = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYRDENOISE_ENH")));
    enhance->set_active (false);
    enhance->set_tooltip_text (M("TP_DIRPYRDENOISE_ENH_TOOLTIP"));
    // ---- Median FIltering ----

    Gtk::Frame* medianFrame = Gtk::manage (new Gtk::Frame ());
    medianFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *medianVBox = Gtk::manage ( new Gtk::VBox());
    medianVBox->set_spacing(2);

    median = Gtk::manage (new Gtk::CheckButton (M("TP_DIRPYRDENOISE_MED") + ":"));
    median->set_active (true);
    medianFrame->set_label_widget(*median);


    Gtk::HSeparator *hsep2 = Gtk::manage (new  Gtk::HSeparator());
    hsep2->show ();

    methodmed = Gtk::manage (new MyComboBoxText ());
    methodmed->append (M("TP_DIRPYRDENOISE_LM"));
    methodmed->append (M("TP_DIRPYRDENOISE_ABM"));
    methodmed->append (M("TP_DIRPYRDENOISE_LPLABM"));
    methodmed->append (M("TP_DIRPYRDENOISE_LABM"));
    methodmed->append (M("TP_DIRPYRDENOISE_RGBM"));
    methodmed->set_active (0);
    methodmed->set_tooltip_text (M("TP_DIRPYRDENOISE_METM_TOOLTIP"));
    methodmedconn = methodmed->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::methodmedChanged) );

    rgbmethod = Gtk::manage (new MyComboBoxText ());
    rgbmethod->append (M("TP_DIRPYRDENOISE_3X3_SOFT"));
    rgbmethod->append (M("TP_DIRPYRDENOISE_3X3"));
    rgbmethod->append (M("TP_DIRPYRDENOISE_5X5_SOFT"));
    rgbmethod->set_active (0);
    rgbmethod->set_tooltip_text (M("TP_DIRPYRDENOISE_MET_TOOLTIP"));
    rgbmethodconn = rgbmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::rgbmethodChanged) );


    medmethod = Gtk::manage (new MyComboBoxText ());
    medmethod->append (M("TP_DIRPYRDENOISE_3X3_SOFT"));
    medmethod->append (M("TP_DIRPYRDENOISE_3X3"));
    medmethod->append (M("TP_DIRPYRDENOISE_5X5_SOFT"));
    medmethod->append (M("TP_DIRPYRDENOISE_5X5"));
    medmethod->append (M("TP_DIRPYRDENOISE_7X7"));
    medmethod->append (M("TP_DIRPYRDENOISE_9X9"));
    medmethod->set_active (0);
    medmethod->set_tooltip_text (M("TP_DIRPYRDENOISE_MET_TOOLTIP"));
    medmethodconn = medmethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::medmethodChanged) );

    ctboxm = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labmm = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_MEDMETHOD") + ":"));
    ctboxm->pack_start (*labmm, Gtk::PACK_SHRINK, 1);

    ctbox = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labm = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_MEDTYPE") + ":"));
    ctbox->pack_start (*labm, Gtk::PACK_SHRINK, 1);

    ctboxrgb = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* labrgb = Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_MEDTYPE") + ":"));
    ctboxrgb->pack_start (*labrgb, Gtk::PACK_SHRINK, 1);


    Gtk::HSeparator *hsep4 = Gtk::manage (new  Gtk::HSeparator());
    hsep4->show ();

    Gtk::HBox* hb11 = Gtk::manage (new Gtk::HBox ());
    hb11->pack_start (*Gtk::manage (new Gtk::Label ( M("TP_DIRPYRDENOISE_METHOD11") + ": ")), Gtk::PACK_SHRINK, 4);
    hb11->set_tooltip_markup (M("TP_DIRPYRDENOISE_METHOD11_TOOLTIP"));

    smethod = Gtk::manage (new MyComboBoxText ());
    smethod->append (M("TP_DIRPYRDENOISE_SHAL"));
//  smethod->append (M("TP_DIRPYRDENOISE_SHBI"));
    smethod->append (M("TP_DIRPYRDENOISE_SHALBI"));
//  smethod->append (M("TP_DIRPYRDENOISE_SHALAL"));
//  smethod->append (M("TP_DIRPYRDENOISE_SHBIBI"));
    smethod->set_active(1);
    hb11->pack_start (*smethod, Gtk::PACK_EXPAND_WIDGET, 4);
    pack_start( *hb11, Gtk::PACK_SHRINK, 4);
    smethodconn = smethod->signal_changed().connect ( sigc::mem_fun(*this, &DirPyrDenoise::smethodChanged) );

    passes  = Gtk::manage (new Adjuster (M("TP_DIRPYRDENOISE_PASSES"), 1.0, 3.0, 1., 1.));
    passes->set_tooltip_text (M("TP_DIRPYRDENOISE_PASSES_TOOLTIP"));
    passes->setAdjusterListener (this);
    passes->show();
    ctboxL->pack_start (*Lmethod);
    lumaVBox->pack_start (*ctboxL);
    lumaVBox->pack_start (*luma);
    lumaVBox->pack_start (*NoiscurveEditorG, Gtk::PACK_SHRINK, 4);
    lumaVBox->pack_start (*Ldetail);
    lumaFrame->add(*lumaVBox);
    pack_start (*lumaFrame);

    ctboxC->pack_start (*Cmethod);
    ctboxC2->pack_start (*C2method);

    if(options.rtSettings.leveldnautsimpl == 1) {
        chromaVBox->pack_start (*ctboxC);
    } else {
        chromaVBox->pack_start (*ctboxC2);
    }

    chromaVBox->pack_start (*NoiseLabels);
    chromaVBox->pack_start (*TileLabels);
    chromaVBox->pack_start (*PrevLabels);

    chromaVBox->pack_start (*chroma);
    chromaVBox->pack_start (*redchro);
    chromaVBox->pack_start (*bluechro);
    chromaVBox->pack_start (*CCcurveEditorG, Gtk::PACK_SHRINK, 4);
    chromaFrame->add(*chromaVBox);
    pack_start (*chromaFrame);


    pack_start (*gamma);
    //pack_start (*enhance);
    pack_start (*hsep4);

//  pack_start( *hb11, Gtk::PACK_SHRINK, 4);

//  pack_start (*hsep2);
//  pack_start (*median);

    ctboxm->pack_start (*methodmed);
    ctbox->pack_start (*medmethod);
    ctboxrgb->pack_start (*rgbmethod);
//  pack_start (*ctboxm);
//  pack_start (*ctbox);
//  pack_start (*ctboxrgb);
//  pack_start (*passes,Gtk::PACK_SHRINK, 1);

    medianVBox->pack_start (*ctboxm);
    medianVBox->pack_start (*ctbox);
    medianVBox->pack_start (*ctboxrgb);
    medianVBox->pack_start (*passes);
    medianFrame->add(*medianVBox);

    pack_start (*medianFrame);


//  pack_start (*perform);
    enhanConn = enhance->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::enhanceChanged) );
    autochromaConn = autochroma->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::autochromaChanged) );
    medianConn = median->signal_toggled().connect( sigc::mem_fun(*this, &DirPyrDenoise::medianChanged) );
    ctboxrgb->hide();

}

DirPyrDenoise::~DirPyrDenoise ()
{
    delete NoiscurveEditorG;
    delete CCcurveEditorG;

}
int chromaChangedUI (void* data)
{
    (static_cast<DirPyrDenoise*>(data))->chromaComputed_ ();
    return 0;
}
void DirPyrDenoise::chromaChanged (double autchroma, double autred, double autblue)
{
    nextchroma = autchroma;
//  printf("CHROM=%f\n",nextchroma);
    nextred = autred;
    nextblue = autblue;
    add_idle(chromaChangedUI, this);
}

bool DirPyrDenoise::chromaComputed_ ()
{

    disableListener ();
    chroma->setValue (nextchroma);
    redchro->setValue (nextred);
    bluechro->setValue (nextblue);
    enableListener ();
    updateNoiseLabel ();
    return false;
}
int TilePrevChangedUI (void* data)
{
    (static_cast<DirPyrDenoise*>(data))->TilePrevComputed_ ();
    return 0;
}


void DirPyrDenoise::noiseTilePrev (int tileX, int tileY, int prevX, int prevY, int sizeT, int sizeP)
{
    nexttileX = tileX;
    nexttileY = tileY;
    nextprevX = prevX;
    nextprevY = prevY;
    nextsizeT = sizeT;
    nextsizeP = sizeP;

    add_idle(TilePrevChangedUI, this);


}
bool DirPyrDenoise::TilePrevComputed_ ()
{

    disableListener ();
    enableListener ();
    updateTileLabel ();
    updatePrevLabel ();
    return false;
}
void DirPyrDenoise::updateTileLabel ()
{
    if (!batchMode) {
        float sT;
        float nX, nY;
        sT = nextsizeT;
        nX = nexttileX;
        nY = nexttileY;
        {
            TileLabels->set_text(
                Glib::ustring::compose(M("TP_DIRPYRDENOISE_TILELABEL"),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), sT),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), nX),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), nY))
            );
        }
    }
}
void DirPyrDenoise::updatePrevLabel ()
{
    if (!batchMode) {
        float sP;
        float pX, pY;
        sP = nextsizeP;
        pX = nextprevX;
        pY = nextprevY;
        {
            PrevLabels->set_text(
                Glib::ustring::compose(M("TP_DIRPYRDENOISE_PREVLABEL"),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), sP),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), pX),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), pY))
            );
        }
    }
}


int noiseChangedUI (void* data)
{
    (static_cast<DirPyrDenoise*>(data))->noiseComputed_ ();
    return 0;
}


void DirPyrDenoise::noiseChanged (double nresid, double highresid)
{
    nextnresid = nresid;
    nexthighresid = highresid;
    add_idle(noiseChangedUI, this);
}

bool DirPyrDenoise::noiseComputed_ ()
{

    disableListener ();
    enableListener ();
    updateNoiseLabel ();
    return false;
}

void DirPyrDenoise::updateNoiseLabel ()
{
    if (!batchMode) {
        float nois, high;
        nois = nextnresid;
        high = nexthighresid;

        if(nois == 0.f && high == 0.f) {
            NoiseLabels->set_text(M("TP_DIRPYRDENOISE_NOISELABELEMPTY"));
        } else {
            NoiseLabels->set_text(
                Glib::ustring::compose(M("TP_DIRPYRDENOISE_NOISELABEL"),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), nois),
                                       Glib::ustring::format(std::fixed, std::setprecision(0), high))
            );
        }
    }
}



void DirPyrDenoise::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
    dmethodconn.block(true);
    Lmethodconn.block(true);
    Cmethodconn.block(true);
    C2methodconn.block(true);
    smethodconn.block(true);
    autochromaConn.block(true);
    medmethodconn.block(true);
    rgbmethodconn.block(true);
    methodmedconn.block(true);


    autochromaChanged ();
    dmethod->set_active (0);

    if (pp->dirpyrDenoise.dmethod == "Lab") {
        dmethod->set_active (0);
    } else if (pp->dirpyrDenoise.dmethod == "RGB") {
        dmethod->set_active (1);
    }

    dmethodChanged ();

    Lmethod->set_active (0);

    if (pp->dirpyrDenoise.Lmethod == "CUR") {
        Lmethod->set_active (0);
    } else if (pp->dirpyrDenoise.Lmethod == "SLI") {
        Lmethod->set_active (1);
    }

    LmethodChanged();

    if(options.rtSettings.leveldnautsimpl == 1) {
        Cmethod->set_active (0);

        if (pp->dirpyrDenoise.Cmethod == "MAN") {
            Cmethod->set_active (0);
        } else if (pp->dirpyrDenoise.Cmethod == "AUT") {
            Cmethod->set_active (1);
        } else if (pp->dirpyrDenoise.Cmethod == "PON") {
            Cmethod->set_active (2);
        } else if (pp->dirpyrDenoise.Cmethod == "PRE") {
            Cmethod->set_active (3);
        }

        CmethodChanged();
    } else {
        C2method->set_active (0);

        if (pp->dirpyrDenoise.C2method == "MANU") {
            C2method->set_active (0);
        } else if (pp->dirpyrDenoise.C2method == "AUTO") {
            C2method->set_active (1);
        } else if (pp->dirpyrDenoise.C2method == "PREV") {
            C2method->set_active (2);
        }

        C2methodChanged();
    }

    smethod->set_active (0);

    if (pp->dirpyrDenoise.smethod == "shal") {
        smethod->set_active (0);
    }
//    else if (pp->dirpyrDenoise.smethod=="shbi")
//        smethod->set_active (1);
    else if (pp->dirpyrDenoise.smethod == "shalbi") {
        smethod->set_active (1);
    }

//   else if (pp->dirpyrDenoise.smethod=="shalal")
//       smethod->set_active (3);
//   else if (pp->dirpyrDenoise.smethod=="shbibi")
//       smethod->set_active (4);

    methodmed->set_active (0);

//   if (pp->dirpyrDenoise.methodmed=="none")
//       methodmed->set_active (0);
    if (pp->dirpyrDenoise.methodmed == "Lonly") {
        methodmed->set_active (0);
    } else if (pp->dirpyrDenoise.methodmed == "ab") {
        methodmed->set_active (1);
    } else if (pp->dirpyrDenoise.methodmed == "Lpab") {
        methodmed->set_active (2);
    } else if (pp->dirpyrDenoise.methodmed == "Lab") {
        methodmed->set_active (3);
    } else if (pp->dirpyrDenoise.methodmed == "RGB") {
        methodmed->set_active (4);
    }

    methodmedChanged();

    medmethod->set_active (0);

//    if (pp->dirpyrDenoise.medmethod=="none")
//        medmethod->set_active (0);
    if (pp->dirpyrDenoise.medmethod == "soft") {
        medmethod->set_active (0);
    } else if (pp->dirpyrDenoise.medmethod == "33") {
        medmethod->set_active (1);
    } else if (pp->dirpyrDenoise.medmethod == "55soft") {
        medmethod->set_active (2);
    } else if (pp->dirpyrDenoise.medmethod == "55") {
        medmethod->set_active (3);
    } else if (pp->dirpyrDenoise.medmethod == "77") {
        medmethod->set_active (4);
    } else if (pp->dirpyrDenoise.medmethod == "99") {
        medmethod->set_active (5);
    }

    medmethodChanged();

    rgbmethod->set_active (0);

//    if (pp->dirpyrDenoise.medmethod=="none")
//        medmethod->set_active (0);
    if (pp->dirpyrDenoise.rgbmethod == "soft") {
        rgbmethod->set_active (0);
    } else if (pp->dirpyrDenoise.rgbmethod == "33") {
        rgbmethod->set_active (1);
    } else if (pp->dirpyrDenoise.rgbmethod == "55soft") {
        rgbmethod->set_active (2);
    }

    rgbmethodChanged();


    if (pedited) {
        if (!pedited->dirpyrDenoise.dmethod) {
            dmethod->set_active (2);
        }

        if (!pedited->dirpyrDenoise.smethod) {
            smethod->set_active (2);
        }

        if (!pedited->dirpyrDenoise.rgbmethod) {
            rgbmethod->set_active (2);
        }

        if (!pedited->dirpyrDenoise.medmethod) {
            medmethod->set_active (6);
        }

        if (!pedited->dirpyrDenoise.methodmed) {
            methodmed->set_active (5);
        }

        if (!pedited->dirpyrDenoise.Cmethod) {
            Cmethod->set_active (4);
        }

        if (!pedited->dirpyrDenoise.C2method) {
            C2method->set_active (3);
        }

        if (!pedited->dirpyrDenoise.Lmethod) {
            Lmethod->set_active (2);
        }

        luma->setEditedState       (pedited->dirpyrDenoise.luma ? Edited : UnEdited);
        Ldetail->setEditedState    (pedited->dirpyrDenoise.Ldetail ? Edited : UnEdited);
        chroma->setEditedState     (pedited->dirpyrDenoise.chroma ? Edited : UnEdited);
        redchro->setEditedState    (pedited->dirpyrDenoise.redchro ? Edited : UnEdited);
        bluechro->setEditedState   (pedited->dirpyrDenoise.bluechro ? Edited : UnEdited);

        gamma->setEditedState      (pedited->dirpyrDenoise.gamma ? Edited : UnEdited);
        passes->setEditedState     (pedited->dirpyrDenoise.passes ? Edited : UnEdited);
        set_inconsistent           (multiImage && !pedited->dirpyrDenoise.enabled);
        enhance->set_inconsistent  (!pedited->dirpyrDenoise.enhance);
        median->set_inconsistent   (!pedited->dirpyrDenoise.median);
        ccshape->setUnChanged      (!pedited->dirpyrDenoise.cccurve);

        //      perform->set_inconsistent (!pedited->dirpyrDenoise.perform);
    }

//  perfconn.block (true);
    setEnabled(pp->dirpyrDenoise.enabled);
    enhance->set_active (pp->dirpyrDenoise.enhance);
//   perform->set_active (pp->dirpyrDenoise.perform);
    median->set_active (pp->dirpyrDenoise.median);
    autochroma->set_active (pp->dirpyrDenoise.autochroma);

//   perfconn.block (false);
    lastmedian = pp->dirpyrDenoise.median;
    lastautochroma = pp->dirpyrDenoise.autochroma;
    lastenhance = pp->dirpyrDenoise.enhance;
//  lastperform = pp->dirpyrDenoise.perform;
    luma->setValue    (pp->dirpyrDenoise.luma);
    Ldetail->setValue (pp->dirpyrDenoise.Ldetail);
    chroma->setValue  (pp->dirpyrDenoise.chroma);
    redchro->setValue  (pp->dirpyrDenoise.redchro);
    bluechro->setValue  (pp->dirpyrDenoise.bluechro);

    gamma->setValue   (pp->dirpyrDenoise.gamma);
    passes->setValue   (pp->dirpyrDenoise.passes);
    lshape->setCurve   (pp->dirpyrDenoise.lcurve);
    ccshape->setCurve   (pp->dirpyrDenoise.cccurve);

    autochromaConn.block(false);

    dmethodconn.block(false);
    Lmethodconn.block(false);
    Cmethodconn.block(false);
    C2methodconn.block(false);
    smethodconn.block(false);
    medmethodconn.block(false);
    rgbmethodconn.block(false);
    methodmedconn.block(false);
    updateNoiseLabel ();

    enableListener ();

}
void DirPyrDenoise::setEditProvider  (EditDataProvider *provider)
{
    lshape->setEditProvider(provider);
    ccshape->setEditProvider(provider);

}
void DirPyrDenoise::autoOpenCurve ()
{
    lshape->openIfNonlinear();
    ccshape->openIfNonlinear();
}
void DirPyrDenoise::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->dirpyrDenoise.luma      = luma->getValue ();
    pp->dirpyrDenoise.Ldetail   = Ldetail->getValue ();
    pp->dirpyrDenoise.chroma    = chroma->getValue ();
    pp->dirpyrDenoise.redchro   = redchro->getValue ();
    pp->dirpyrDenoise.bluechro  = bluechro->getValue ();
    pp->dirpyrDenoise.gamma     = gamma->getValue ();
    pp->dirpyrDenoise.passes    = passes->getValue ();
    pp->dirpyrDenoise.enabled   = getEnabled();
    pp->dirpyrDenoise.enhance   = enhance->get_active();
//  pp->dirpyrDenoise.perform   = perform->get_active();
    pp->dirpyrDenoise.median   = median->get_active();
    pp->dirpyrDenoise.autochroma   = autochroma->get_active();
    pp->dirpyrDenoise.lcurve  = lshape->getCurve ();
    pp->dirpyrDenoise.cccurve  = ccshape->getCurve ();

    if (pedited) {
        pedited->dirpyrDenoise.dmethod  = dmethod->get_active_row_number() != 2;
        pedited->dirpyrDenoise.Lmethod  = Lmethod->get_active_row_number() != 2;
        pedited->dirpyrDenoise.Cmethod  = Cmethod->get_active_row_number() != 4;
        pedited->dirpyrDenoise.C2method  = C2method->get_active_row_number() != 3;
        pedited->dirpyrDenoise.smethod  = smethod->get_active_row_number() != 2;
        pedited->dirpyrDenoise.medmethod  = medmethod->get_active_row_number() != 6;
        pedited->dirpyrDenoise.rgbmethod  = rgbmethod->get_active_row_number() != 2;
        pedited->dirpyrDenoise.methodmed  = methodmed->get_active_row_number() != 5;
        pedited->dirpyrDenoise.luma     = luma->getEditedState ();
        pedited->dirpyrDenoise.Ldetail  = Ldetail->getEditedState ();
        pedited->dirpyrDenoise.chroma   = chroma->getEditedState ();
        pedited->dirpyrDenoise.redchro  = redchro->getEditedState ();
        pedited->dirpyrDenoise.bluechro = bluechro->getEditedState ();
        pedited->dirpyrDenoise.gamma    = gamma->getEditedState ();
        pedited->dirpyrDenoise.passes    = passes->getEditedState ();
        pedited->dirpyrDenoise.enabled  = !get_inconsistent();
        pedited->dirpyrDenoise.enhance  = !enhance->get_inconsistent();
        pedited->dirpyrDenoise.median  = !median->get_inconsistent();
        pedited->dirpyrDenoise.autochroma  = !autochroma->get_inconsistent();
        pedited->dirpyrDenoise.lcurve    = !lshape->isUnChanged ();
        pedited->dirpyrDenoise.cccurve    = !ccshape->isUnChanged ();

        //    pedited->dirpyrDenoise.perform  = !perform->get_inconsistent();
    }

    if (dmethod->get_active_row_number() == 0) {
        pp->dirpyrDenoise.dmethod = "Lab";
    } else if (dmethod->get_active_row_number() == 1) {
        pp->dirpyrDenoise.dmethod = "RGB";
    }

    if (Lmethod->get_active_row_number() == 0) {
        pp->dirpyrDenoise.Lmethod = "CUR";
    } else if (Lmethod->get_active_row_number() == 1) {
        pp->dirpyrDenoise.Lmethod = "SLI";
    }

    if(options.rtSettings.leveldnautsimpl == 1) {
        if (Cmethod->get_active_row_number() == 0) {
            pp->dirpyrDenoise.Cmethod = "MAN";
        } else if (Cmethod->get_active_row_number() == 1) {
            pp->dirpyrDenoise.Cmethod = "AUT";
        } else if (Cmethod->get_active_row_number() == 2) {
            pp->dirpyrDenoise.Cmethod = "PON";
        } else if (Cmethod->get_active_row_number() == 3) {
            pp->dirpyrDenoise.Cmethod = "PRE";
        }
    } else {
        if (C2method->get_active_row_number() == 0) {
            pp->dirpyrDenoise.C2method = "MANU";
        } else if (C2method->get_active_row_number() == 1) {
            pp->dirpyrDenoise.C2method = "AUTO";
        } else if (C2method->get_active_row_number() == 2) {
            pp->dirpyrDenoise.C2method = "PREV";
        }
    }

    if (smethod->get_active_row_number() == 0) {
        pp->dirpyrDenoise.smethod = "shal";
    }
//   else if (smethod->get_active_row_number()==1)
//        pp->dirpyrDenoise.smethod = "shbi";
    else if (smethod->get_active_row_number() == 1) {
        pp->dirpyrDenoise.smethod = "shalbi";
    }

//    else if (smethod->get_active_row_number()==3)
//        pp->dirpyrDenoise.smethod = "shalal";
//    else if (smethod->get_active_row_number()==4)
//        pp->dirpyrDenoise.smethod = "shbibi";

    //  if (methodmed->get_active_row_number()==0)
    //   pp->dirpyrDenoise.methodmed = "none";
    if (methodmed->get_active_row_number() == 0) {
        pp->dirpyrDenoise.methodmed = "Lonly";
    } else if (methodmed->get_active_row_number() == 1) {
        pp->dirpyrDenoise.methodmed = "ab";
    } else if (methodmed->get_active_row_number() == 2) {
        pp->dirpyrDenoise.methodmed = "Lpab";
    } else if (methodmed->get_active_row_number() == 3) {
        pp->dirpyrDenoise.methodmed = "Lab";
    } else if (methodmed->get_active_row_number() == 4) {
        pp->dirpyrDenoise.methodmed = "RGB";
    }



//   if (medmethod->get_active_row_number()==0)
//       pp->dirpyrDenoise.medmethod = "none";
    if (medmethod->get_active_row_number() == 0) {
        pp->dirpyrDenoise.medmethod = "soft";
    } else if (medmethod->get_active_row_number() == 1) {
        pp->dirpyrDenoise.medmethod = "33";
    } else if (medmethod->get_active_row_number() == 2) {
        pp->dirpyrDenoise.medmethod = "55soft";
    } else if (medmethod->get_active_row_number() == 3) {
        pp->dirpyrDenoise.medmethod = "55";
    } else if (medmethod->get_active_row_number() == 4) {
        pp->dirpyrDenoise.medmethod = "77";
    } else if (medmethod->get_active_row_number() == 5) {
        pp->dirpyrDenoise.medmethod = "99";
    }

    if (rgbmethod->get_active_row_number() == 0) {
        pp->dirpyrDenoise.rgbmethod = "soft";
    } else if (rgbmethod->get_active_row_number() == 1) {
        pp->dirpyrDenoise.rgbmethod = "33";
    } else if (rgbmethod->get_active_row_number() == 2) {
        pp->dirpyrDenoise.rgbmethod = "55soft";
    }

}

void DirPyrDenoise::curveChanged (CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == lshape) {
            listener->panelChanged (EvDPDNLCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == ccshape) {
            listener->panelChanged (EvDPDNCCCurve, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void DirPyrDenoise::dmethodChanged ()
{

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDPDNmet, dmethod->get_active_text ());
    }
}
void DirPyrDenoise::LmethodChanged ()
{
    if (!batchMode) {
        if(Lmethod->get_active_row_number() == 0) {             // CUR
            luma->hide();
            NoiscurveEditorG->show();
        } else if(Lmethod->get_active_row_number() == 1) {          // SLI
            luma->show();
            NoiscurveEditorG->hide();
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDPDNLmet, Lmethod->get_active_text ());
    }
}

void DirPyrDenoise::CmethodChanged ()
{
    if (!batchMode) {
        if(Cmethod->get_active_row_number() == 0 ) { //MAN
            chroma->show();
            redchro->show();
            bluechro->show();
            chroma->set_sensitive(true);
            redchro->set_sensitive(true);
            bluechro->set_sensitive(true);
            NoiseLabels->show();
            TileLabels->hide();
            PrevLabels->hide();

        } else if(Cmethod->get_active_row_number() == 1) {          // AUT
            chroma->show();
            redchro->show();
            bluechro->show();
            NoiseLabels->show();
            TileLabels->hide();
            PrevLabels->hide();

            chroma->set_sensitive(false);
            redchro->set_sensitive(false);
            bluechro->set_sensitive(false);
        } else if(Cmethod->get_active_row_number() == 2) {          // PON
            chroma->hide();
            redchro->hide();
            bluechro->hide();
            NoiseLabels->hide();
            TileLabels->hide();
            PrevLabels->hide();

            chroma->set_sensitive(false);
            redchro->set_sensitive(false);
            bluechro->set_sensitive(false);
        } else if(Cmethod->get_active_row_number() == 3) {          // PRE
            chroma->show();
            redchro->show();
            bluechro->show();
            NoiseLabels->show();
            TileLabels->show();
            PrevLabels->show();

            chroma->set_sensitive(false);
            redchro->set_sensitive(false);
            bluechro->set_sensitive(false);
        }

    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDPDNCmet, Cmethod->get_active_text ());
    }
}

void DirPyrDenoise::C2methodChanged ()
{
    if (!batchMode) {
        if(C2method->get_active_row_number() == 0 ) { //MAN
            chroma->show();
            redchro->show();
            bluechro->show();
            chroma->set_sensitive(true);
            redchro->set_sensitive(true);
            bluechro->set_sensitive(true);
            NoiseLabels->show();
            TileLabels->hide();
            PrevLabels->hide();

        } else if(C2method->get_active_row_number() == 1) {         // AUTO
            chroma->show();
            redchro->show();
            bluechro->show();
            NoiseLabels->show();
            TileLabels->hide();
            PrevLabels->hide();

            chroma->set_sensitive(false);
            redchro->set_sensitive(false);
            bluechro->set_sensitive(false);
        } else if(C2method->get_active_row_number() == 2) {         // PREV
            chroma->show();
            redchro->show();
            bluechro->show();
            NoiseLabels->show();
            TileLabels->hide();
            PrevLabels->hide();

            chroma->set_sensitive(false);
            redchro->set_sensitive(false);
            bluechro->set_sensitive(false);
        }


    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDPDNC2met, C2method->get_active_text ());
    }
}


void DirPyrDenoise::smethodChanged ()
{

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDPDNsmet, smethod->get_active_text ());
    }
}

void DirPyrDenoise::medmethodChanged ()
{


    if (listener && (multiImage || getEnabled())  && median->get_active() ) {
        listener->panelChanged (EvDPDNmedmet, medmethod->get_active_text ());
    }
}

void DirPyrDenoise::rgbmethodChanged ()
{
    ctboxrgb->hide();

    if(methodmed->get_active_row_number() == 4) {
        ctboxrgb->show();
    }

    if (listener && (multiImage || getEnabled())  && median->get_active()) {
        listener->panelChanged (EvDPDNrgbmet, rgbmethod->get_active_text ());
    }
}



void DirPyrDenoise::methodmedChanged ()
{
    if(methodmed->get_active_row_number() == 4) {
        ctboxrgb->show();
        ctbox->hide();
    } else {
        ctboxrgb->hide();
        ctbox->show();
    }

    if (listener && (multiImage || getEnabled())  && median->get_active()) {
        listener->panelChanged (EvDPDNmetmed, methodmed->get_active_text ());
    }
}

void DirPyrDenoise::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    luma->setDefault    (defParams->dirpyrDenoise.luma);
    Ldetail->setDefault (defParams->dirpyrDenoise.Ldetail);
    chroma->setDefault  (defParams->dirpyrDenoise.chroma);
    redchro->setDefault  (defParams->dirpyrDenoise.redchro);
    bluechro->setDefault  (defParams->dirpyrDenoise.bluechro);
    gamma->setDefault   (defParams->dirpyrDenoise.gamma);
    passes->setDefault   (defParams->dirpyrDenoise.passes);

    if (pedited) {
        luma->setDefaultEditedState     (pedited->dirpyrDenoise.luma ? Edited : UnEdited);
        Ldetail->setDefaultEditedState  (pedited->dirpyrDenoise.Ldetail ? Edited : UnEdited);
        chroma->setDefaultEditedState   (pedited->dirpyrDenoise.chroma ? Edited : UnEdited);
        redchro->setDefaultEditedState   (pedited->dirpyrDenoise.redchro ? Edited : UnEdited);
        bluechro->setDefaultEditedState   (pedited->dirpyrDenoise.bluechro ? Edited : UnEdited);
        gamma->setDefaultEditedState    (pedited->dirpyrDenoise.gamma ? Edited : UnEdited);
        passes->setDefaultEditedState    (pedited->dirpyrDenoise.passes ? Edited : UnEdited);
    } else {
        luma->setDefaultEditedState     (Irrelevant);
        Ldetail->setDefaultEditedState  (Irrelevant);
        chroma->setDefaultEditedState   (Irrelevant);
        redchro->setDefaultEditedState   (Irrelevant);
        bluechro->setDefaultEditedState   (Irrelevant);
        gamma->setDefaultEditedState    (Irrelevant);
        passes->setDefaultEditedState    (Irrelevant);
    }
}

void DirPyrDenoise::adjusterChanged (Adjuster* a, double newval)
{

    Glib::ustring costr;
    costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());

    if (listener && getEnabled()) {
        if (a == Ldetail) {
            listener->panelChanged (EvDPDNLdetail, costr);
        } else if (a == luma) {
            listener->panelChanged (EvDPDNLuma, costr);
        } else if (a == chroma) {
            listener->panelChanged (EvDPDNChroma, costr);
        } else if (a == redchro) {
            listener->panelChanged (EvDPDNredchro, costr);
        } else if (a == bluechro) {
            listener->panelChanged (EvDPDNbluechro, costr);
        } else if (a == gamma) {
            listener->panelChanged (EvDPDNGamma, costr);
        } else if (a == passes  && median->get_active()) {
            listener->panelChanged (EvDPDNpasses, costr);
        }
    }
}

void DirPyrDenoise::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvDPDNEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvDPDNEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDPDNEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void DirPyrDenoise::enhanceChanged ()
{

    if (batchMode) {
        if (enhance->get_inconsistent()) {
            enhance->set_inconsistent (false);
            enhanConn.block (true);
            enhance->set_active (false);
            enhanConn.block (false);
        } else if (lastenhance) {
            enhance->set_inconsistent (true);
        }

        lastenhance = enhance->get_active ();
    }

    if (listener) {

        if (enhance->get_active ()) {
            listener->panelChanged (EvDPDNenhance, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDPDNenhance, M("GENERAL_DISABLED"));
        }
    }
}

void DirPyrDenoise::medianChanged ()
{

    if (batchMode) {
        if (median->get_inconsistent()) {
            median->set_inconsistent (false);
            medianConn.block (true);
            median->set_active (false);
            medianConn.block (false);
        } else if (lastmedian) {
            median->set_inconsistent (true);
        }

        lastmedian = median->get_active ();
    }

    if (listener) {
        if (median->get_active ()) {
            listener->panelChanged (EvDPDNmedian, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDPDNmedian, M("GENERAL_DISABLED"));
        }


    }
}
void DirPyrDenoise::autochromaChanged ()
{
//  printf("Autochroma\n");
    if (batchMode) {
        if (autochroma->get_inconsistent()) {
            autochroma->set_inconsistent (false);
            autochromaConn.block (true);
            autochroma->set_active (false);
            autochromaConn.block (false);
        } else if (lastautochroma) {
            autochroma->set_inconsistent (true);
        }

        lastautochroma = autochroma->get_active ();
    }

    if (autochroma->get_active ()) {
        chroma->set_sensitive(false);
        redchro->set_sensitive(false);
        bluechro->set_sensitive(false);
    } else {
        chroma->set_sensitive(true);
        redchro->set_sensitive(true);
        bluechro->set_sensitive(true);
    }

    if (listener) {
        if (autochroma->get_active ()) {
            listener->panelChanged (EvDPDNautochroma, M("GENERAL_ENABLED"));
            //  chroma->set_sensitive(false);
            //  redchro->set_sensitive(false);
            //  bluechro->set_sensitive(false);
        } else {
            listener->panelChanged (EvDPDNautochroma, M("GENERAL_DISABLED"));
            //chroma->set_sensitive(true);
            //redchro->set_sensitive(true);
            //bluechro->set_sensitive(true);
        }


    }
}


/*
void DirPyrDenoise::perform_toggled () {

    if (batchMode) {
        if (perform->get_inconsistent()) {
            perform->set_inconsistent (false);
            perfconn.block (true);
            perform->set_active (false);
            perfconn.block (false);
        }
        else if (lastperform)
            perform->set_inconsistent (true);

        lastperform = perform->get_active ();
    }

    if (listener) {
        if (perform->get_active ())
            listener->panelChanged (EvDPDNperform, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvDPDNperform, M("GENERAL_DISABLED"));
    }
}
*/

void DirPyrDenoise::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    luma->showEditedCB ();
    Ldetail->showEditedCB ();
    chroma->showEditedCB ();
    redchro->showEditedCB ();
    bluechro->showEditedCB ();
    gamma->showEditedCB ();
    passes->showEditedCB ();
    NoiscurveEditorG->setBatchMode (batchMode);
    CCcurveEditorG->setBatchMode (batchMode);

    dmethod->append (M("GENERAL_UNCHANGED"));
    Lmethod->append (M("GENERAL_UNCHANGED"));
    Cmethod->append (M("GENERAL_UNCHANGED"));
    C2method->append (M("GENERAL_UNCHANGED"));
    smethod->append (M("GENERAL_UNCHANGED"));
    medmethod->append (M("GENERAL_UNCHANGED"));
    methodmed->append (M("GENERAL_UNCHANGED"));
    rgbmethod->append (M("GENERAL_UNCHANGED"));

}

void DirPyrDenoise::setAdjusterBehavior (bool lumaadd, bool lumdetadd, bool chromaadd, bool chromaredadd, bool chromablueadd, bool gammaadd, bool passesadd)
{

    luma->setAddMode(lumaadd);
    Ldetail->setAddMode(lumdetadd);
    chroma->setAddMode(chromaadd);
    redchro->setAddMode(chromaredadd);
    bluechro->setAddMode(chromablueadd);
    gamma->setAddMode(gammaadd);
    passes->setAddMode(passesadd);

}
void DirPyrDenoise::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
    }

    else if (callerId == 4) {    // LH - bottom bar
        Color::hsv2rgb01(float(valX), 0.5f, float(valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}

void DirPyrDenoise::trimValues (rtengine::procparams::ProcParams* pp)
{

    luma->trimValue(pp->dirpyrDenoise.luma);
    Ldetail->trimValue(pp->dirpyrDenoise.Ldetail);
    chroma->trimValue(pp->dirpyrDenoise.chroma);
    redchro->trimValue(pp->dirpyrDenoise.redchro);
    bluechro->trimValue(pp->dirpyrDenoise.bluechro);
    gamma->trimValue(pp->dirpyrDenoise.gamma);
    passes->trimValue(pp->dirpyrDenoise.passes);
}
