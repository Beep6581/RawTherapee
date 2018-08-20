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
 */


#include "locallab.h"
#include "rtimage.h"
#include <iomanip>
#include "../rtengine/rt_math.h"
#include "options.h"
#include <cmath>
#include "edit.h"
#include "guiutils.h"
#include <string>
#include <unistd.h>
#include "../rtengine/improcfun.h"
#include "eventmapper.h"

#define MINCHRO 0.
#define MAXCHRO 150
#define CENTERCHRO 10
#define MAXCHROCC 100


using namespace rtengine;

extern Options options;


Locallab::Locallab():
    FoldableToolPanel(this, "locallab", M("TP_LOCALLAB_LABEL"), false, true),
    lastObject(-1),
    expcolor(new MyExpander(true, M("TP_LOCALLAB_COFR"))),
    expexpose(new MyExpander(true, M("TP_LOCALLAB_EXPOSE"))),
    expvibrance(new MyExpander(true, M("TP_LOCALLAB_VIBRANCE"))),
    expblur(new MyExpander(true, M("TP_LOCALLAB_BLUFR"))),
    exptonemap(new MyExpander(true, M("TP_LOCALLAB_TM"))),
    expreti(new MyExpander(true, M("TP_LOCALLAB_RETI"))),
    expsharp(new MyExpander(true, M("TP_LOCALLAB_SHARP"))),
    expcbdl(new MyExpander(true, M("TP_LOCALLAB_CBDL"))),
    expdenoi(new MyExpander(true, M("TP_LOCALLAB_DENOIS"))),
    expsettings(new ControlSpotPanel()),

    LocalcurveEditorgainT(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAIN"))),
    LocalcurveEditorgainTrab(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_TRANSMISSIONGAINRAB"))),
    llCurveEditorG(new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_LUM"))),


    anbspot(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ANBSPOT"), 0, 1, 1, 0))),
    locX(Gtk::manage(new Adjuster(M("TP_LOCAL_WIDTH"), 0, 2250, 1, 250))),
    locXL(Gtk::manage(new Adjuster(M("TP_LOCAL_WIDTH_L"), 0, 2250, 1, 250))),
    degree(Gtk::manage(new Adjuster(M("TP_LOCAL_DEGREE"), -180, 180, 1, 0))),
    locY(Gtk::manage(new Adjuster(M("TP_LOCAL_HEIGHT"), 0, 2250, 1, 250))),
    locYT(Gtk::manage(new Adjuster(M("TP_LOCAL_HEIGHT_T"), 0, 2250, 1, 250))),
    centerX(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTER_X"), -1000, 1000, 1, 0))),
    centerY(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTER_Y"), -1000, 1000, 1, 0))),
    circrad(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CIRCRADIUS"), 2, 150, 1, 18))),
    sensiexclu(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIEXCLU"), 0, 100, 1, 19))),
    struc(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRUC"), 0, 5, 1, 0))),
    thres(Gtk::manage(new Adjuster(M("TP_LOCALLAB_THRES"), 1, 35, 1, 18))),
    proxi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_PROXI"), 0, 60, 1, 0))),
    lightness(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LIGHTNESS"), -100, 100, 1, 0))),
    contrast(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CONTRAST"), -100, 100, 1, 0))),
    chroma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMA"), -100, 150, 1, 0))),
    sensi(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    expcomp(Gtk::manage(new Adjuster(M("TP_EXPOSURE_EXPCOMP"), -200, 200, 5, 0))),
    hlcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 20))),
    hlcomprthresh(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 33))),
    black(Gtk::manage(new Adjuster(M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0))),
    shcompr(Gtk::manage(new Adjuster(M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50))),
    sensiex(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    radius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RADIUS"), 1, 100, 1, 1))),
    strength(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STRENGTH"), 0, 100, 1, 0))),
    sensibn(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIBN"), 0, 100, 1, 40))),
    transit(Gtk::manage(new Adjuster(M("TP_LOCALLAB_TRANSIT"), 5, 95, 1, 60))),
    stren(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STREN"), -50, 100, 1, 0))),
    gamma(Gtk::manage(new Adjuster(M("TP_LOCALLAB_GAM"), 80, 150, 1, 100))),
    estop(Gtk::manage(new Adjuster(M("TP_LOCALLAB_ESTOP"), 10, 400, 1, 140))),
    scaltm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SCALTM"), 1, 100, 1, 10))),
    rewei(Gtk::manage(new Adjuster(M("TP_LOCALLAB_REWEI"), 0, 9, 1, 0))),
    sensitm(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19))),
    str(Gtk::manage(new Adjuster(M("TP_LOCALLAB_STR"), 0, 100, 1, 0))),
    neigh(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NEIGH"), 14, 150, 1, 50))),
    vart(Gtk::manage(new Adjuster(M("TP_LOCALLAB_VART"), 50, 500, 1, 200))),
    chrrt(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHRRT"), 0, 100, 1, 0))),
    sensih(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIH"), 0, 100, 1, 19))),
    retrab(Gtk::manage(new Adjuster(M("TP_LOCALLAB_RETRAB"), 0, 10000, 1, 500))),
    chromacbdl(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMACBDL"), 0, 300, 1, 0))),
    threshold(Gtk::manage(new Adjuster(M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 100, 1, 20))),
    sensicb(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSICB"), 0, 100, 1, 19))),
    sharradius(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARRADIUS"), 42, 500, 1, 4))),
    sharamount(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARAMOUNT"), 0, 100, 1, 75))),
    shardamping(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARDAMPING"), 0, 100, 1, 75))),
    shariter(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SHARITER"), 5, 100, 1, 30))),
    sensisha(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIS"), 0, 100, 1, 19))),
    noiselumdetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMDETAIL"), 0, 100, 1, 0))),
    noisechrodetail(Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHRODETAIL"), 0, 100, 1, 0))),
    bilateral(Gtk::manage(new Adjuster(M("TP_LOCALLAB_BILATERAL"), 0, 100, 1, 0))),
    sensiden(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSIDEN"), 0, 100, 1, 30))),
    hueref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_HUEREF"), -3.15, 3.15, 0.01, 0))),
    huerefblur(Gtk::manage(new Adjuster(M("TP_LOCALLAB_HUEREFBLUR"), -3.15, 3.15, 0.01, 0))),
    chromaref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CHROMAREF"), 0, 200, 0.01, 0))),
    lumaref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_LUMAMAREF"), 0, 100, 0.01, 0))),
    sobelref(Gtk::manage(new Adjuster(M("TP_LOCALLAB_SOBELREF"), 0, 100, 0.01, 0))),
    centerXbuf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTERBUF_X"), -1000, 1000, 1, 0))),
    centerYbuf(Gtk::manage(new Adjuster(M("TP_LOCALLAB_CENTERBUF_Y"), -1000, 1000, 1, 0))),

    shapemethod(Gtk::manage(new MyComboBoxText())),
    Smethod(Gtk::manage(new MyComboBoxText())),
    Exclumethod(Gtk::manage(new MyComboBoxText())),

    retinexMethod(Gtk::manage(new MyComboBoxText())),
    qualityMethod(Gtk::manage(new MyComboBoxText())),
    qualitycurveMethod(Gtk::manage(new MyComboBoxText())),
    blurMethod(Gtk::manage(new MyComboBoxText())),
    dustMethod(Gtk::manage(new MyComboBoxText())),

    shapeFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_SHFR")))),
    superFrame(Gtk::manage(new Gtk::Frame())),
    dustFrame(Gtk::manage(new Gtk::Frame(M("TP_LOCALLAB_DUST")))),
    wavFrame(Gtk::manage(new Gtk::Frame())),


    labmdh(Gtk::manage(new Gtk::Label(M("TP_LOCRETI_METHOD") + ":"))),
    labqualcurv(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_QUALCURV_METHOD") + ":"))),
    labmS(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_STYPE") + ":"))),
    labmEx(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_EXCLUTYPE") + ":"))),
    labmshape(Gtk::manage(new Gtk::Label(M("TP_LOCALLAB_SHAPETYPE") + ":"))),

    dhbox(Gtk::manage(new Gtk::HBox())),
    qualcurvbox(Gtk::manage(new Gtk::HBox())),

    avoid(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_AVOID")))),
    activlum(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_ACTIV")))),
    invers(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    curvactiv(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CURV")))),
    inversrad(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    inversret(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    inverssha(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_INVERS")))),
    cutpast(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_CUTPAST")))),
    lastdust(Gtk::manage(new Gtk::CheckButton(M("TP_LOCALLAB_LASTDUST"))))

{
    CurveListener::setMulti(true);
    ProcParams params;

    auto m = ProcEventMapper::getInstance();
    EvLocenacolor = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENACOLOR");//548
    EvLocenaexpose = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENAEXPOSE");//572
    EvLocenavibrance = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENAVIBR");//563
    EvLocenablur = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENABLUR");//549
    EvLocenatonemap = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENATM");//550
    EvLocenareti = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENARETI");//551
    EvLocenasharp = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENASHARP");//552
    EvLocenacbdl = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENACBDL");//553
    EvLocenadenoi = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENADENOI");//554

    EvlocallablocX = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLOCX"); //= 494,
    EvlocallabCenter = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCENTER"); //= 495,
    EvlocallabDegree = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCDEGRE"); //= 496,
    Evlocallablightness = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLIGHT"); //= 497,
    Evlocallabcontrast = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCONTRA"); //= 498,
    Evlocallabchroma = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCHROMA"); //= 499,
    Evlocallabtransit = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCTRANSIT"); //= 500,
    Evlocallabavoid = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCAVOID"); //= 501,
    EvlocallablocYT = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLOCYT"); // = 502,
    EvlocallablocXL = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCXL"); //= 503,
    EvlocallabSmet = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSMET"); //= 504,
    Evlocallabinvers = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCINVERS"); //= 505,
    Evlocallabradius = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRADIUS"); //= 506,
    Evlocallabinversrad = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCINVRAD"); //= 507,
    Evlocallabstrength = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSTRENGTH"); //= 508,
    Evlocallabsensi = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSI"); //= 509,
    EvlocallabretinexMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRETIMETH");//510
    Evlocallabstr = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRETISTR");//= 511,
    Evlocallabneigh = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRETINEIGH");//= 512,
    Evlocallabvart = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRETIVART");//= 513,
    EvlocallabCTgainCurve = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRETIGAINCURV");//= 514,
    Evlocallabchrrt = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCHRRT");//= 515,
    Evlocallabinversret = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCINVRET");//= 516,
    Evlocallabsensih = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIH");//= 517,
    Evlocallabnbspot = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNBSPOT");//= 518,
    Evlocallabactivlum = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCACTIVLUM");//= 519,
    Evlocallabanbspot = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCANBSPOT");//= 520,
    Evlocallabsharradius = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHARADIUS");//= 521,
    Evlocallabsharamount = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHAAMOUNT");//= 522,
    Evlocallabshardamping = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHADAMPING");//= 523,
    Evlocallabshariter = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHAITER");//= 524,
    Evlocallabsensis = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIS");//= 525,
    Evlocallabinverssha = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCINVSHA");//= 526,
    Evlocallabcircrad = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCIRCRAD");//= 527,
    Evlocallabthres = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCTHRES");//= 528,
    Evlocallabproxi = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCPROXI");//= 529,
    EvlocallabqualityMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCQUALMETH");//= 530,
    Evlocallabnoiselumf = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISLUMF");//= 531,
    Evlocallabnoiselumc = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISLUMC");//= 532,
    Evlocallabnoisechrof = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISCHROF");//= 533,
    Evlocallabnoisechroc = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISCHROC");//= 534,
    EvlocallabThresho = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCBDLTHRESHO");//= 535,
    EvlocallabEqualizer = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCBDLEQUALIZ");//= 536,
    Evlocallabsensicb = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSICB");//= 537,
    Evlocallabsensibn = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIBN");//= 538,
    Evlocallabstren = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSTREN");//= 539,
    Evlocallabgamma = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCGAMM");//= 540,
    Evlocallabestop = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCESTOP");//= 541,
    Evlocallabscaltm = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSCALTM");//= 542,
    Evlocallabrewei = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCREWEI");//= 543,
    Evlocallabsensitm = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSITM");//= 544,
    EvlocallabCTgainCurverab = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCGAINCURRAB");//= 545,
    Evlocallabretrab = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCRETRAB");//= 546,
    Evlocallabllshape = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLSHAPE");//= 547,
    EvlocallabLHshape = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLHSHAPE");// = 555,
    Evlocallabcurvactiv = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCURVACTIV");// = 556,
    Evlocallabccshape = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCCSHAPE");// = 557,
    EvlocallabqualitycurveMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCQUALCURVMETH");// = 558,
    Evlocallabhueref = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCHUEREF");// = 559,
    Evlocallabchromaref = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCHROMAREF");// = 560,
    Evlocallablumaref = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLUMAREF");// = 561,
    EvlocallabHHshape = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCHHSHAPE");// = 562,
    EvlocallabSkinTonesCurve = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSKINTONCURV");// = 564,
    EvlocallabProtectSkins = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCPROTSKIN");// = 565,
    EvlocallabAvoidColorShift = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCAVOIDCOLORSHIFT");// = 566,
    EvlocallabPastSatTog = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCPASTSATTOG");// = 567,
    EvlocallabPastels = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCPASTEL");// = 568,
    EvlocallabSaturated = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSATUR");// = 569,
    EvlocallabPastSatThreshold = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCPASTSATTHRES");// = 570,
    Evlocallabsensiv = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIV");// = 571,
    Evlocallabexpcomp = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCEXPCOMP");// = 573,
    Evlocallabhlcompr = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCHLCOMPR");// = 574,
    Evlocallabhlcomprthresh = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCHLCOMPRTHRESH");// = 575,
    Evlocallabblack = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCBLACK");// = 576,
    Evlocallabshcompr = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHCOMPR");// = 577,
    Evlocallabsensiex = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIEX");// = 578,
    Evlocallabshapeexpos = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHAPE");// = 579,
    EvlocallabCenterbuf = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCENTERBUF");// = 580,
    Evlocallabadjblur = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISEEQUALBLURED");// = 581,
    Evlocallabcutpast = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCUTPAST");// = 582,
    Evlocallabchromacbdl = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCCHROCBDL");// = 583,
    EvlocallabblurMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCBLURMETH"); //584
    EvlocallabdustMethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCDUSTMETH");// = 585,
    Evlocallablastdust = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLASTDUST");// = 586,
    Evlocallabsobelref = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSOBELREF");// = 587,
    Evlocallabexclumethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCEXCLUMETH");// = 588,
    Evlocallabsensiexclu = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIEXCL");// = 589,
    Evlocallabstruc = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSTRUC");// = 590,
    Evlocallabwarm = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCWARM");// = 591,
    Evlocallabnoiselumdetail = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISELUMDETAIL");// = 592,
    Evlocallabnoisechrodetail = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISECHRODETAIL");// = 593,
    Evlocallabsensiden = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSENSIDEN");// = 594,
    Evlocallabhuerefblur = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCHUEREFBLUR");// = 595,
    EvlocallabEnabled = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCENABLED");// = 596,
    EvlocallablocY = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCLOCY");// = 597,
    Evlocallabbilateral = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCBILATERAL");// = 598,
    Evlocallabnoiselequal = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCNOISELEQUAL");// = 599,
    Evlocallabshapemethod = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSHAPEMETH");// = 600,
    Evlocallabspotduplicated = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_LOCSPOTDUP");// = 601
    Evlocallabspotcreated = m->newEvent(LUMINANCECURVE, "Spot creation");// = 602

    int realnbspot;


    realnbspot = options.rtSettings.nspot;
    nbspot = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NBSPOT"), 1, realnbspot, 1, 1));

    if (options.rtSettings.locdelay) {

        if (nbspot->delay < 200) {
            nbspot->delay = 200;
        }
    }

    const LocallabParams default_params;

    shapeFrame->set_label_align(0.025, 0.5);

    expsettings->getExpander()->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsettings->getExpander()));

    expcolor->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcolor));
    enablecolorConn = expcolor->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcolor));

    expexpose->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expexpose));
    enableexposeConn = expexpose->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expexpose));

    expvibrance->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expvibrance));
    enablevibranceConn = expvibrance->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expvibrance));

    expblur->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expblur));
    enableblurConn = expblur->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expblur));

    exptonemap->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), exptonemap));
    enabletonemapConn = exptonemap->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), exptonemap));

    expreti->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expreti));
    enableretiConn = expreti->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expreti));

    expsharp->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expsharp));
    enablesharpConn = expsharp->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expsharp));

    expcbdl->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expcbdl));
    enablecbdlConn = expcbdl->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expcbdl));

    expdenoi->signal_button_release_event().connect_notify(sigc::bind(sigc::mem_fun(this, &Locallab::foldAllButMe), expdenoi));
    enabledenoiConn = expdenoi->signal_enabled_toggled().connect(sigc::bind(sigc::mem_fun(this, &Locallab::enableToggled), expdenoi));

    shapemethod->append(M("TP_LOCALLAB_ELI"));
    shapemethod->append(M("TP_LOCALLAB_RECT"));
    shapemethod->set_active(0);
    Exclumethod->append(M("TP_LOCALLAB_EXNORM"));
    Exclumethod->append(M("TP_LOCALLAB_EXECLU"));
    Exclumethod->set_active(0);

    struc->set_tooltip_text(M("TP_LOCALLAB_STRUC_TOOLTIP"));
    struc->setAdjusterListener(this);

    Smethod->append(M("TP_LOCALLAB_IND"));
    Smethod->append(M("TP_LOCALLAB_SYM"));
    Smethod->append(M("TP_LOCALLAB_INDSL"));
    Smethod->append(M("TP_LOCALLAB_SYMSL"));
    Smethod->set_active(0);
    qualityMethod->append(M("TP_LOCALLAB_STD"));
    qualityMethod->append(M("TP_LOCALLAB_ENH"));
    qualityMethod->append(M("TP_LOCALLAB_ENHDEN"));
    qualityMethod->set_active(0);

    std::vector<GradientMilestone> milestones;
    std::vector<double> defaultCurve;
    std::vector<double> defaultCurve2;
    std::vector<double> defaultCurve2rab;
    std::vector<double> defaultCurve3;
    std::vector<double> defaultCurve4;
    std::vector<double> defaultCurve5;

    irg   = Gtk::manage(new RTImage("Chanmixer-RG.png"));

    qualitycurveMethod->append(M("TP_LOCALLAB_CURVNONE"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVCURR"));
    qualitycurveMethod->append(M("TP_LOCALLAB_CURVENH"));
    qualitycurveMethod->set_active(0);
    qualitycurveMethodConn = qualitycurveMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::qualitycurveMethodChanged));
    qualitycurveMethod->set_tooltip_markup(M("TP_LOCALLAB_CURVEMETHOD_TOOLTIP"));

    llCurveEditorG->setCurveListener(this);

    llshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "L(L)"));
    llshape->setResetCurve(DCT_NURBS, {(double)DCT_NURBS, 0.0, 0.0, 1.0, 1.0});
    llshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    milestones.push_back(GradientMilestone(0., 0., 0., 0.));
    milestones.push_back(GradientMilestone(1., 1., 1., 1.));
    llshape->setBottomBarBgGradient(milestones);
    llshape->setLeftBarBgGradient(milestones);

    ccshape = static_cast<DiagonalCurveEditor*>(llCurveEditorG->addCurve(CT_Diagonal, "C(C)"));
    ccshape->setResetCurve(DCT_NURBS, {(double)DCT_NURBS, 0.0, 0.0, 1.0, 1.0});
    ccshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    milestones.push_back(GradientMilestone(0., 0., 0., 0.));
    milestones.push_back(GradientMilestone(1., 1., 1., 1.));
    ccshape->setBottomBarBgGradient(milestones);
    ccshape->setLeftBarBgGradient(milestones);


    LHshape = static_cast<FlatCurveEditor*>(llCurveEditorG->addCurve(CT_Flat, "L(H)", nullptr, false, true));

    LHshape->setIdentityValue(0.);
    LHshape->setResetCurve(FCT_MinMaxCPoints, {(double)FCT_MinMaxCPoints, 0.0, 0.50, 0.35, 0.35, 0.166, 0.50, 0.35, 0.35, 0.333, 0.50, 0.35, 0.35, 0.50, 0.50, 0.35, 0.35, 0.666, 0.50, 0.35, 0.35, 0.833, 0.50, 0.35, 0.35});
    LHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    LHshape->setCurveColorProvider(this, 1);
    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float (i) * (1.0f / 6.0);
        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    LHshape->setBottomBarBgGradient(milestones);


    HHshape = static_cast<FlatCurveEditor*>(llCurveEditorG->addCurve(CT_Flat, "H(H)", nullptr, false, true));

    HHshape->setIdentityValue(0.);
    HHshape->setResetCurve(FCT_MinMaxCPoints, {(double)FCT_MinMaxCPoints, 0.0, 0.50, 0.35, 0.35, 0.166, 0.50, 0.35, 0.35, 0.333, 0.50, 0.35, 0.35, 0.50, 0.50, 0.35, 0.35, 0.666, 0.50, 0.35, 0.35, 0.833, 0.50, 0.35, 0.35});
    HHshape->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_LL_TOOLTIP"));
    HHshape->setCurveColorProvider(this, 1);
    milestones.clear();

    for (int i = 0; i < 7; i++) {
        float R, G, B;
        float x = float (i) * (1.0f / 6.0);

        Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);
        milestones.push_back(GradientMilestone(double (x), double (R), double (G), double (B)));
    }

    HHshape->setBottomBarBgGradient(milestones);


    llCurveEditorG->curveListComplete();
    lightness->setAdjusterListener(this);

    contrast->setAdjusterListener(this);
    /*
        Gtk::Image* iblueredL = Gtk::manage(new RTImage("ajd-wb-bluered1.png"));
        Gtk::Image* iblueredR = Gtk::manage(new RTImage("ajd-wb-bluered2.png"));

        warm = Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., iblueredL, iblueredR));
        warm->setAdjusterListener(this);
    */
    chroma->setAdjusterListener(this);

    sensi->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensi->setAdjusterListener(this);

    centerXbuf->setAdjusterListener(this);;
    centerYbuf->setAdjusterListener(this);;
//    adjblur->setAdjusterListener(this);;

//exposure

    expcomp->setAdjusterListener(this);
    hlcomprthresh->setAdjusterListener(this);
    black->setAdjusterListener(this);
    hlcompr->setAdjusterListener(this);
    shcompr->setAdjusterListener(this);
    sensiex->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensiex->setAdjusterListener(this);
    Gtk::Image* iblueredL = Gtk::manage(new RTImage("ajd-wb-bluered1.png"));
    Gtk::Image* iblueredR = Gtk::manage(new RTImage("ajd-wb-bluered2.png"));

    warm = Gtk::manage(new Adjuster(M("TP_LOCALLAB_WARM"), -100., 100., 1., 0., iblueredL, iblueredR));
    warm->setAdjusterListener(this);
    warm->set_tooltip_text(M("TP_LOCALLAB_WARM_TOOLTIP"));

    radius->setAdjusterListener(this);
    strength->setAdjusterListener(this);


    sensibn->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensibn->setAdjusterListener(this);

    activlum->set_active(false);
    activlumConn  = activlum->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::activlumChanged));

    invers->set_active(false);
    inversConn  = invers->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversChanged));

    curvactiv->set_active(false);
    curvactivConn  = curvactiv->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::curvactivChanged));

    inversrad->set_active(false);
    inversradConn  = inversrad->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversradChanged));

    inversret->set_active(false);
    inversretConn  = inversret->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversretChanged));

    cutpast->set_active(false);
    cutpastConn  = cutpast->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::cutpastChanged));
    cutpast->set_tooltip_text(M("TP_LOCALLAB_CUTPAST_TOOLTIP"));
//tone mapping local

    lastdust->set_active(false);
    lastdustConn  = lastdust->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::lastdustChanged));
    lastdust->set_tooltip_text(M("TP_LOCALLAB_LASTDUST_TOOLTIP"));

    stren->setAdjusterListener(this);

    gamma->setAdjusterListener(this);

    estop->setAdjusterListener(this);

    scaltm->setAdjusterListener(this);

    rewei->setAdjusterListener(this);

    sensitm->set_tooltip_text(M("TP_LOCALLAB_SENSI_TOOLTIP"));
    sensitm->setAdjusterListener(this);

//end TM


//retinex local

    dhbox->pack_start(*labmdh, Gtk::PACK_SHRINK, 1);

    retinexMethod->append(M("TP_RETINEX_LOW"));
    retinexMethod->append(M("TP_RETINEX_UNIFORM"));
    retinexMethod->append(M("TP_RETINEX_HIGH"));
    retinexMethod->set_active(0);
    retinexMethodConn = retinexMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::retinexMethodChanged));
    retinexMethod->set_tooltip_markup(M("TP_LOCRETI_METHOD_TOOLTIP"));

    str->setAdjusterListener(this);
    neigh->setAdjusterListener(this);
    vart->setAdjusterListener(this);
    chrrt->setAdjusterListener(this);
    sensih->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensih->setAdjusterListener(this);
    retrab->setAdjusterListener(this);


    LocalcurveEditorgainT->setCurveListener(this);


    cTgainshape = static_cast<FlatCurveEditor*>(LocalcurveEditorgainT->addCurve(CT_Flat, "", nullptr, false, false));

    cTgainshape->setIdentityValue(0.);
    cTgainshape->setResetCurve(FCT_MinMaxCPoints, {(double)FCT_MinMaxCPoints, 0.0, 0.12, 0.35, 0.35, 0.70, 0.50, 0.35, 0.35, 1.00, 0.12, 0.35, 0.35});
    cTgainshape->setTooltip(M("TP_RETINEX_GAINTRANSMISSION_TOOLTIP"));

    LocalcurveEditorgainTrab->setCurveListener(this);



    cTgainshaperab = static_cast<FlatCurveEditor*>(LocalcurveEditorgainTrab->addCurve(CT_Flat, "", nullptr, false, false));


    cTgainshaperab->setIdentityValue(0.);
    // cTgainshaperab->setResetCurve(FlatCurveType(default_params.localTgaincurverab.at(0)), default_params.localTgaincurverab);
    cTgainshaperab->setTooltip(M("TP_RETINEX_GAINTRANSMISSIONRAB_TOOLTIP"));

    LocalcurveEditorgainT->curveListComplete();
    LocalcurveEditorgainT->show();
    LocalcurveEditorgainTrab->curveListComplete();
    LocalcurveEditorgainTrab->show();


// end reti
    avoid->set_active(false);
    avoidConn  = avoid->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidChanged));
    pack_start(*anbspot);

    hueref->setAdjusterListener(this);
    huerefblur->setAdjusterListener(this);
    chromaref->setAdjusterListener(this);
    lumaref->setAdjusterListener(this);
    sobelref->setAdjusterListener(this);

    pack_start(*hueref);
    pack_start(*huerefblur);
    pack_start(*chromaref);
    pack_start(*lumaref);
    pack_start(*sobelref);

    anbspot->hide();//keep anbspot  - i used it to test diffrent algo...
    hueref->hide();
    huerefblur->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();

    // expsettings->add(*spotPanel);
    expsettings->setLevel(2);
    pack_start(*expsettings->getExpander());



    Gtk::HBox * buttonBox1 = Gtk::manage(new Gtk::HBox(true, 10));

    Gtk::Button * lumacontrastMinusButton = Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")));
    buttonBox1->pack_start(*lumacontrastMinusButton);
    lumacontrastMinusPressedConn = lumacontrastMinusButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumacontrastMinusPressed));

    Gtk::Button * lumaneutralButton = Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")));
    buttonBox1->pack_start(*lumaneutralButton);
    lumaneutralPressedConn = lumaneutralButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumaneutralPressed));

    Gtk::Button * lumacontrastPlusButton = Gtk::manage(new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")));
    buttonBox1->pack_start(*lumacontrastPlusButton);
    lumacontrastPlusPressedConn = lumacontrastPlusButton->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::lumacontrastPlusPressed));
    ToolParamBlock* const cbdlBox = Gtk::manage(new ToolParamBlock());

    cbdlBox->pack_start(*buttonBox1);

    for (int i = 0; i < 5; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);

        if (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if (i == 4) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage(new Adjuster(ss, 0, 400, 1, 100));
        multiplier[i]->setAdjusterListener(this);
        cbdlBox->pack_start(*multiplier[i]);
    }

    Gtk::HSeparator *separator3 = Gtk::manage(new  Gtk::HSeparator());
    cbdlBox->pack_start(*separator3, Gtk::PACK_SHRINK, 2);

    chromacbdl->set_tooltip_text(M("TP_LOCALLAB_CHROMACB_TOOLTIP"));

    chromacbdl->setAdjusterListener(this);
    cbdlBox->pack_start(*chromacbdl);

    threshold->setAdjusterListener(this);
    cbdlBox->pack_start(*threshold);

    sensicb->set_tooltip_text(M("TP_LOCALLAB_SENSIH_TOOLTIP"));
    sensicb->setAdjusterListener(this);
    cbdlBox->pack_start(*sensicb);

    sharradius->setAdjusterListener(this);

    sharamount->setAdjusterListener(this);

    shardamping->setAdjusterListener(this);

    shariter->setAdjusterListener(this);


    sensisha->set_tooltip_text(M("TP_LOCALLAB_SENSIS_TOOLTIP"));
    sensisha->setAdjusterListener(this);

    inverssha->set_active(false);
    inversshaConn  = inverssha->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::inversshaChanged));
    ToolParamBlock* const sharpBox = Gtk::manage(new ToolParamBlock());

    sharpBox->pack_start(*sharradius);
    sharpBox->pack_start(*sharamount);
    sharpBox->pack_start(*shardamping);
    sharpBox->pack_start(*shariter);
    sharpBox->pack_start(*sensisha);
    sharpBox->pack_start(*inverssha);

    noiselumf = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMFINE"), MINCHRO, MAXCHRO, 1, 0));
    noiselumc = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELUMCOARSE"), MINCHRO, MAXCHROCC, 1, 0));

    noisechrof = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROFINE"), MINCHRO, MAXCHRO, 1, 0));
    noisechroc = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISECHROCOARSE"), MINCHRO, MAXCHROCC, 1, 0));

    Gtk::Image* iblueredL1 = Gtk::manage(new RTImage("ajd-wb-bluered1.png"));
    Gtk::Image* iblueredR1 = Gtk::manage(new RTImage("ajd-wb-bluered2.png"));

    adjblur = Gtk::manage(new Adjuster(M("TP_LOCALLAB_ADJ"), -100., 100., 1., 0., iblueredL1, iblueredR1));
    adjblur->setAdjusterListener(this);

    Gtk::Image *ar = Gtk::manage(new RTImage("adj-black.png"));
    Gtk::Image *al = Gtk::manage(new RTImage("adj-white.png"));

    noiselequal = Gtk::manage(new Adjuster(M("TP_LOCALLAB_NOISELEQUAL"), -2, 10, 1, 7, al, ar));
    noiselequal->setAdjusterListener(this);

    if (noiselequal->delay < 200) {
        noiselequal->delay = 200;
    }


    noiselumf->setAdjusterListener(this);

    if (noiselumf->delay < 200) {
        noiselumf->delay = 200;
    }


    noiselumc->setAdjusterListener(this);

    if (noiselumc->delay < 200) {
        noiselumc->delay = 200;
    }

    noiselumc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));


    noiselumdetail->setAdjusterListener(this);

    if (noiselumdetail->delay < 200) {
        noiselumdetail->delay = 200;
    }



    noisechrodetail->setAdjusterListener(this);

    if (noisechrodetail->delay < 200) {
        noisechrodetail->delay = 200;
    }


    bilateral->setAdjusterListener(this);

    if (bilateral->delay < 200) {
        bilateral->delay = 200;
    }


    sensiden->setAdjusterListener(this);

    if (sensiden->delay < 200) {
        sensiden->delay = 200;
    }


    noisechrof->setAdjusterListener(this);

    if (noisechrof->delay < 200) {
        noisechrof->delay = 200;
    }


    noisechroc->setAdjusterListener(this);

    if (noisechroc->delay < 200) {
        noisechroc->delay = 200;
    }


    noisechroc->set_tooltip_text(M("TP_LOCALLAB_NOISECHROC_TOOLTIP"));

    ToolParamBlock* const denoisBox = Gtk::manage(new ToolParamBlock());
    ToolParamBlock* const wavBox = Gtk::manage(new ToolParamBlock());

    wavBox->pack_start(*noiselumf);
    wavBox->pack_start(*noiselumc);
    wavBox->pack_start(*noiselumdetail);
    wavBox->pack_start(*noiselequal);
    wavBox->pack_start(*noisechrof);
    wavBox->pack_start(*noisechroc);
    //wavBox->pack_start(*noisechrodetail);
    wavBox->pack_start(*adjblur);


    wavFrame->add(*wavBox);
    denoisBox->pack_start(*wavFrame);

    denoisBox->pack_start(*bilateral);
    denoisBox->pack_start(*sensiden);

    neutrHBox1 = Gtk::manage(new Gtk::HBox());

    neutral1 = Gtk::manage(new Gtk::Button(M("TP_LOCALLAB_NEUTRAL")));
    RTImage *resetImg1 = Gtk::manage(new RTImage("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral1->set_image(*resetImg1);
    neutral1->set_tooltip_text(M("TP_LOCALLAB_NEUTRAL_TIP"));
    neutralconn1 = neutral1->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::neutral_pressed));
    neutral1->show();
    neutrHBox1->pack_start(*neutral1);
    pack_start(*neutrHBox1);

    superFrame->set_label_align(0.025, 0.5);
    //  Gtk::VBox *superVBox = Gtk::manage ( new Gtk::VBox());
    ToolParamBlock* const superBox = Gtk::manage(new ToolParamBlock());

    superFrame->set_label_widget(*curvactiv);
    ToolParamBlock* const colorBox = Gtk::manage(new ToolParamBlock());

//    ToolParamBlock* const dustBox = Gtk::manage (new ToolParamBlock());
    dustMethod->append(M("TP_LOCALLAB_DSCOP"));
    dustMethod->append(M("TP_LOCALLAB_DSMOV"));
    dustMethod->append(M("TP_LOCALLAB_DSPAS"));
    dustMethod->set_active(0);
    dustMethodConn = dustMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::dustMethodChanged));
    dustMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));

    superBox->pack_start(*lightness);
    superBox->pack_start(*contrast);
    superBox->pack_start(*chroma);

    superFrame->add(*superBox);
    colorBox->pack_start(*superFrame);

//    colorBox->pack_start(*warm);
    colorBox->pack_start(*sensi);
    /*
        dustFrame->set_label_align(0.025, 0.5);
        dustBox->pack_start (*dustMethod);
        dustBox->pack_start (*lastdust);
        dustBox->pack_start (*cutpast);
       dustBox->pack_start (*centerXbuf);
        dustBox->pack_start (*centerYbuf);
       dustBox->pack_start (*adjblur);
       dustFrame->add (*dustBox);
      colorBox->pack_start (*dustFrame);
      */
    centerXbuf->hide();
    centerYbuf->hide();

    qualcurvbox->pack_start(*labqualcurv, Gtk::PACK_SHRINK, 4);
    qualcurvbox->pack_start(*qualitycurveMethod);

    colorBox->pack_start(*qualcurvbox);


    colorBox->pack_start(*llCurveEditorG, Gtk::PACK_SHRINK, 2);
    colorBox->pack_start(*invers);

    expcolor->add(*colorBox);
    expcolor->setLevel(2);
    pack_start(*expcolor);

    ToolParamBlock* const exposeBox = Gtk::manage(new ToolParamBlock());

    curveEditorG = new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_LOCALLAB_CURVEEDITOR_TONES_LABEL"));
    curveEditorG->setCurveListener(this);

    shapeexpos = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, ""));
    shapeexpos->setTooltip(M("TP_LOCALLAB_CURVEEDITOR_TONES_TOOLTIP"));
    curveEditorG->curveListComplete();


    exposeBox->pack_start(*expcomp);
    exposeBox->pack_start(*hlcompr);
    exposeBox->pack_start(*hlcomprthresh);
    exposeBox->pack_start(*black);
    exposeBox->pack_start(*shcompr);
    exposeBox->pack_start(*warm);
    exposeBox->pack_start(*sensiex);
    exposeBox->pack_start(*curveEditorG);

    expexpose->add(*exposeBox);
    expexpose->setLevel(2);
    pack_start(*expexpose);

    ToolParamBlock* const vibranceBox = Gtk::manage(new ToolParamBlock());
    std::vector<GradientMilestone> milestonesvib;
    float R, G, B;
    // -0.1 rad < Hue < 1.6 rad
    Color::hsv2rgb01(0.92f, 0.45f, 0.6f, R, G, B);
    milestonesvib.push_back(GradientMilestone(0.0, double (R), double (G), double (B)));
    Color::hsv2rgb01(0.14056f, 0.45f, 0.6f, R, G, B);
    milestonesvib.push_back(GradientMilestone(1.0, double (R), double (G), double (B)));

    saturated = Gtk::manage(new Adjuster(M("TP_VIBRANCE_SATURATED"), -100., 100., 1., 0.));
    saturated->setAdjusterListener(this);
    saturated->set_sensitive(false);
    vibranceBox->pack_start(*saturated, Gtk::PACK_SHRINK, 0);

    pastels = Gtk::manage(new Adjuster(M("TP_VIBRANCE_PASTELS"), -100., 100., 1., 0.));
    pastels->setAdjusterListener(this);
    vibranceBox->pack_start(*pastels, Gtk::PACK_SHRINK, 0);

    psThreshold = Gtk::manage(new ThresholdAdjuster(M("TP_VIBRANCE_PSTHRESHOLD"), -100., 100., 0., M("TP_VIBRANCE_PSTHRESHOLD_WEIGTHING"), 0, 0., 100., 75., M("TP_VIBRANCE_PSTHRESHOLD_SATTHRESH"), 0, this, false));
    psThreshold->setAdjusterListener(this);
    psThreshold->set_tooltip_markup(M("TP_VIBRANCE_PSTHRESHOLD_TOOLTIP"));
    psThreshold->set_sensitive(false);
    vibranceBox->pack_start(*psThreshold, Gtk::PACK_SHRINK, 0);

    protectSkins = Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PROTECTSKINS")));
    protectSkins->set_active(true);
    vibranceBox->pack_start(*protectSkins, Gtk::PACK_SHRINK, 0);

    avoidColorShift = Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_AVOIDCOLORSHIFT")));
    avoidColorShift->set_active(true);
    vibranceBox->pack_start(*avoidColorShift, Gtk::PACK_SHRINK, 0);

    pastSatTog = Gtk::manage(new Gtk::CheckButton(M("TP_VIBRANCE_PASTSATTOG")));
    pastSatTog->set_active(true);
    vibranceBox->pack_start(*pastSatTog, Gtk::PACK_SHRINK, 0);

    sensiv = Gtk::manage(new Adjuster(M("TP_LOCALLAB_SENSI"), 0, 100, 1, 19));
    sensiv->setAdjusterListener(this);

    vibranceBox->pack_start(*sensiv, Gtk::PACK_SHRINK, 0);

    curveEditorGG = new CurveEditorGroup(options.lastlocalCurvesDir, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_LABEL"));
    curveEditorGG->setCurveListener(this);

    skinTonesCurve = static_cast<DiagonalCurveEditor*>(curveEditorGG->addCurve(CT_Diagonal, M("TP_VIBRANCE_CURVEEDITOR_SKINTONES")));
    skinTonesCurve->setTooltip(M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_TOOLTIP"));
    skinTonesCurve->setBottomBarBgGradient(milestonesvib);
    skinTonesCurve->setLeftBarBgGradient(milestonesvib);
    skinTonesCurve->setRangeLabels(
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE1"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE2"),
        M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE3"), M("TP_VIBRANCE_CURVEEDITOR_SKINTONES_RANGE4")
    );
    skinTonesCurve->setRangeDefaultMilestones(0.1, 0.4, 0.85);
    curveEditorGG->curveListComplete();

    vibranceBox->pack_start(*curveEditorGG, Gtk::PACK_SHRINK, 4);

    pskinsconn = protectSkins->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::protectskins_toggled));
    ashiftconn = avoidColorShift->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::avoidcolorshift_toggled));
    pastsattogconn = pastSatTog->signal_toggled().connect(sigc::mem_fun(*this, &Locallab::pastsattog_toggled));

    expvibrance->add(*vibranceBox);
    expvibrance->setLevel(2);
    pack_start(*expvibrance);

    ToolParamBlock* const blurrBox = Gtk::manage(new ToolParamBlock());
    blurMethod->append(M("TP_LOCALLAB_BLNORM"));
    blurMethod->append(M("TP_LOCALLAB_BLINV"));
    blurMethod->append(M("TP_LOCALLAB_BLSYM"));
    blurMethod->set_active(0);
    blurMethodConn = blurMethod->signal_changed().connect(sigc::mem_fun(*this, &Locallab::blurMethodChanged));
    blurMethod->set_tooltip_markup(M("TP_LOCALLAB_BLMETHOD_TOOLTIP"));

    blurrBox->pack_start(*radius);
    blurrBox->pack_start(*strength);
    blurrBox->pack_start(*sensibn);
    blurrBox->pack_start(*blurMethod);


    blurrBox->pack_start(*activlum);

    expblur->add(*blurrBox);
    expblur->setLevel(2);
    pack_start(*expblur);
    ToolParamBlock* const tmBox = Gtk::manage(new ToolParamBlock());

    tmBox->pack_start(*stren);
    tmBox->pack_start(*gamma);
    tmBox->pack_start(*estop);
    tmBox->pack_start(*scaltm);
    tmBox->pack_start(*rewei);
    tmBox->pack_start(*sensitm);

    exptonemap->add(*tmBox);
    exptonemap->setLevel(2);
    pack_start(*exptonemap);

    ToolParamBlock* const retiBox = Gtk::manage(new ToolParamBlock());

    retiBox->pack_start(*retinexMethod);
    retiBox->pack_start(*str);
    retiBox->pack_start(*chrrt);
    retiBox->pack_start(*neigh);
    retiBox->pack_start(*vart);
    retiBox->pack_start(*sensih);
    retiBox->pack_start(*retrab);

    retiBox->pack_start(*LocalcurveEditorgainTrab, Gtk::PACK_SHRINK, 4);

    retiBox->pack_start(*LocalcurveEditorgainT, Gtk::PACK_SHRINK, 4);
    retiBox->pack_start(*inversret);

    expreti->add(*retiBox);
    expreti->setLevel(2);
    pack_start(*expreti);


    expsharp->add(*sharpBox);
    expsharp->setLevel(2);
    pack_start(*expsharp);

    expcbdl->add(*cbdlBox);
    expcbdl->setLevel(2);
    pack_start(*expcbdl);

    expdenoi->add(*denoisBox);
    expdenoi->setLevel(2);
    pack_start(*expdenoi);


//    pack_start (*transit);
    pack_start(*avoid); //keep avoid clor shift in case of

    neutrHBox = Gtk::manage(new Gtk::HBox());

    neutral = Gtk::manage(new Gtk::Button(M("TP_LOCALLAB_NEUTRAL")));
    RTImage *resetImg = Gtk::manage(new RTImage("gtk-undo-ltr-small.png", "gtk-undo-rtl-small.png"));
    neutral->set_image(*resetImg);
    neutral->set_tooltip_text(M("TP_LOCALLAB_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect(sigc::mem_fun(*this, &Locallab::neutral_pressed));
    neutral->show();
    neutrHBox->pack_start(*neutral);
    pack_start(*neutrHBox);

    show_all();
}

Locallab::~Locallab()
{
    delete LocalcurveEditorgainT;
    delete LocalcurveEditorgainTrab;
    delete llCurveEditorG;

}
void Locallab::foldAllButMe(GdkEventButton* event, MyExpander *expander)
{
    if (event->button == 3) {
        expsettings->setExpanded(expsettings->getExpander() == expander);
        expcolor->set_expanded(expcolor == expander);
        expexpose->set_expanded(expexpose == expander);
        expvibrance->set_expanded(expvibrance == expander);
        expblur->set_expanded(expblur == expander);
        exptonemap->set_expanded(exptonemap == expander);
        expreti->set_expanded(expreti == expander);
        expsharp->set_expanded(expsharp == expander);
        expcbdl->set_expanded(expcbdl == expander);
        expdenoi->set_expanded(expdenoi == expander);

    }
}

void Locallab::enableToggled(MyExpander *expander)
{
    // TODO Locallab printf
    printf("enableToggled\n");

    if (getEnabled() && listener) {
        rtengine::ProcEvent event = NUMOFEVENTS;

        if (expander == expcolor) {
            event = EvLocenacolor;
        } else if (expander == expexpose) {
            event = EvLocenaexpose;
        } else if (expander == expvibrance) {
            event = EvLocenavibrance;
        } else if (expander == expblur) {
            event = EvLocenablur;
        } else if (expander == exptonemap) {
            event = EvLocenatonemap;
        } else if (expander == expreti) {
            event = EvLocenareti;
        } else if (expander == expsharp) {
            event = EvLocenasharp;
        } else if (expander == expcbdl) {
            event = EvLocenacbdl;
        } else if (expander == expdenoi) {
            event = EvLocenadenoi;
        } else {
            return;
        }

        if (expander->get_inconsistent()) {
            listener->panelChanged(event, M("GENERAL_UNCHANGED"));
        } else if (expander->getEnabled()) {
            listener->panelChanged(event, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(event, M("GENERAL_DISABLED"));
        }

    }
}

void Locallab::writeOptions(std::vector<int> &tpOpen)
{
    tpOpen.push_back(expsettings->getExpanded());
    tpOpen.push_back(expcolor->get_expanded());
    tpOpen.push_back(expexpose->get_expanded());
    tpOpen.push_back(expvibrance->get_expanded());
    tpOpen.push_back(expblur->get_expanded());
    tpOpen.push_back(exptonemap->get_expanded());
    tpOpen.push_back(expreti->get_expanded());
    tpOpen.push_back(expsharp->get_expanded());
    tpOpen.push_back(expcbdl->get_expanded());
    tpOpen.push_back(expdenoi->get_expanded());

}

void Locallab::updateToolState(std::vector<int> &tpOpen)
{
    if (tpOpen.size() >= 10) {
        expsettings->setExpanded(tpOpen.at(0));
        expcolor->set_expanded(tpOpen.at(1));
        expexpose->set_expanded(tpOpen.at(2));
        expvibrance->set_expanded(tpOpen.at(3));
        expblur->set_expanded(tpOpen.at(4));
        exptonemap->set_expanded(tpOpen.at(5));
        expreti->set_expanded(tpOpen.at(6));
        expsharp->set_expanded(tpOpen.at(7));
        expcbdl->set_expanded(tpOpen.at(8));
        expdenoi->set_expanded(tpOpen.at(9));
    }
}



void Locallab::neutral_pressed()
{
    // TODO Locallab printf
    printf("neutral_pressed\n");

    // TODO Locallab To be deleted
    /*
    Smethod->set_active(0);
    shapemethod->set_active(0);
    Exclumethod->set_active(0);
    locX->resetValue(false);
    locXL->resetValue(false);
    locY->resetValue(false);
    locYT->resetValue(false);
    centerX->resetValue(false);
    centerY->resetValue(false);
    circrad->resetValue(false);
    centerXbuf->resetValue(false);
    centerYbuf->resetValue(false);
    adjblur->resetValue(false);
    blurMethod->set_active(0);
    dustMethod->set_active(1);
    sensiexclu->resetValue(false);
    struc->resetValue(false);
    qualityMethod->set_active(1);
    qualitycurveMethod->set_active(0);
    thres->resetValue(false);
    proxi->resetValue(false);
    lightness->resetValue(false);
    chroma->resetValue(false);
    warm->resetValue(false);
    contrast->resetValue(false);
    sensi->resetValue(false);
    radius->resetValue(false);
    strength->resetValue(false);
    transit->resetValue(false);
    sensibn->resetValue(false);
    invers->set_active(false);
    cutpast->set_active(false);
    lastdust->set_active(false);
    curvactiv->set_active(false);
    inversrad->set_active(false);
    inversret->set_active(false);
    protectSkins->set_active(false);
    avoidColorShift->set_active(true);
    pastSatTog->set_active(true);
    expcolor->setEnabled(false);
    expexpose->setEnabled(false);
    expvibrance->setEnabled(false);
    expblur->setEnabled(false);
    exptonemap->setEnabled(false);
    expreti->setEnabled(false);
    expsharp->setEnabled(false);
    expcbdl->setEnabled(false);
    expdenoi->setEnabled(false);
    stren->resetValue(false);
    gamma->resetValue(false);
    estop->resetValue(false);
    scaltm->resetValue(false);
    rewei->resetValue(false);
    sensitm->resetValue(false);
    retinexMethod->set_active(2);
    str->resetValue(false);
    neigh->resetValue(false);
    vart->resetValue(false);
    chrrt->resetValue(false);
    sensih->resetValue(false);
    retrab->resetValue(false);
    expcomp->resetValue(false);
    hlcompr->resetValue(false);
    hlcomprthresh->resetValue(false);
    black->resetValue(false);
    shcompr->resetValue(false);
    sensiex->resetValue(false);
    cTgainshape->reset();
    llshape->reset();
    ccshape->reset();
    LHshape->reset();
    HHshape->reset();
    shapeexpos->reset();
    skinTonesCurve->reset();
    avoid->set_active(false);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->resetValue(false);
    }

    chromacbdl->resetValue(false);
    threshold->resetValue(false);
    sensicb->resetValue(false);
    sharradius->resetValue(false);
    sharamount->resetValue(false);
    shardamping->resetValue(false);
    shariter->resetValue(false);
    sensisha->resetValue(false);
    inverssha->set_active(false);
    noiselumf->resetValue(false);
    noiselumc->resetValue(false);
    noiselumdetail->resetValue(false);
    noiselequal->resetValue(false);
    noisechrof->resetValue(false);
    noisechroc->resetValue(false);
    noisechrodetail->resetValue(false);
    bilateral->resetValue(false);
    sensiden->resetValue(false);
    */
}


void Locallab::lumaneutralPressed()
{
    // TODO Locallab printf
    printf("lumaneutralPressed\n");

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setValue(100);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}


void Locallab::lumacontrastPlusPressed()
{
    // TODO Locallab printf
    printf("lumacontrastPlusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}


void Locallab::lumacontrastMinusPressed()
{
    // TODO Locallab printf
    printf("lumacontrastMinusPressed\n");

    for (int i = 0; i < 5; i++) {
        float inc = - (5 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
    }

    // Raise event (only for first multiplier because associated event concerns all multipliers)
    adjusterChanged(multiplier[0], multiplier[0]->getValue()); // Value isn't used
}

int localChangedUI(void* data)
{
    // TODO Locallab printf
    printf("localChangedUI\n");

    /*
    GThreadLock lock;
    (static_cast<Locallab*>(data))->localComputed_();
    */

    return 0;
}

int localretChangedUI(void* data)
{
    // TODO Locallab printf
    printf("localretChangedUI\n");

    /*
    GThreadLock lock;
    (static_cast<Locallab*>(data))->localretComputed_();
    */

    return 0;
}

bool Locallab::localretComputed_()
{
    // TODO Locallab printf
    printf("localretComputed_\n");

    /*
    disableListener();

    //Reticurv
    //update GUI and MIP specially for curve

    int *s_datc;
    s_datc = new int[70];
    int siz;
    //printf("nexts=%s\n", nextstr2.c_str());
    ImProcFunctions::strcurv_data(nextstr2, s_datc, siz);
    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back((double)(s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve(creti);

    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data(nextll_str2, s_datcl, sizl);
    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back((double)(s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;

    llshape->setCurve(cll);


    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data(nextcc_str2, s_datcc, sizc);
    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back((double)(s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;

    ccshape->setCurve(ccc);

    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data(nextlh_str2, s_datch, sizh);
    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back((double)(s_datch[j]) / 1000.);
    }

    delete [] s_datch;

    LHshape->setCurve(clh);


    int *s_datchh;
    s_datchh = new int[70];
    int sizhh;
    ImProcFunctions::strcurv_data(nexthh_str2, s_datchh, sizhh);
    std::vector<double>   chh;

    for (int j = 0; j < sizhh; j++) {
        chh.push_back((double)(s_datchh[j]) / 1000.);
    }

    delete [] s_datchh;

    HHshape->setCurve(chh);

    //skinTonesCurve
    int *s_datcsk;
    s_datcsk = new int[70];
    int sizsk;
    ImProcFunctions::strcurv_data(nextsk_str2, s_datcsk, sizsk);
    std::vector<double>   csk;

    for (int j = 0; j < sizsk; j++) {
        csk.push_back((double)(s_datcsk[j]) / 1000.);
    }

    delete [] s_datcsk;

    skinTonesCurve->setCurve(csk);

    //PSthreshold
    int sizps = 2;
    int s_datcps[sizps + 1];
    ImProcFunctions::strcurv_data(nextps_str2, s_datcps, sizps);
    psThreshold->setValue(s_datcps[0], s_datcps[1]);


    //exCurve
    int *s_datcex;
    s_datcex = new int[70];
    int sizex;
    ImProcFunctions::strcurv_data(nextex_str2, s_datcex, sizex);
    std::vector<double>   cex;

    for (int j = 0; j < sizex; j++) {
        cex.push_back((double)(s_datcex[j]) / 1000.);
    }

    delete [] s_datcex;

    shapeexpos->setCurve(cex);

    enableListener();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue(1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue(0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 0);

    }

    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve();

    if (cretirab.at(5) == 0.70) {
        cretirab.at(5) = 0.9;
        cTgainshaperab->setCurve(cretirab);

        curveChanged(cTgainshaperab);
    } else if (cretirab.at(5) == 0.90) {
        cretirab.at(5) = 0.7;
        cTgainshaperab->setCurve(cretirab);
        curveChanged(cTgainshaperab);

    }


    //    printf("G2 anbspot=%i\n", anbspot->getValue());

    if (listener) { //for all sliders
        listener->panelChanged(Evlocallabanbspot, ""); //anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabCTgainCurverab, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabCTgainCurve, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(Evlocallabllshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(Evlocallabccshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabLHshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabHHshape, M(""));
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabSkinTonesCurve, M(""));
    }

    if (listener) {//for PSthreshold
        listener->panelChanged(EvlocallabPastSatThreshold, M(""));
    }

    if (listener) {//for excurve
        listener->panelChanged(Evlocallabshapeexpos, M(""));
    }
    */

    return false;

}

bool Locallab::localComputed_()
{
    // TODO Locallab printf
    printf("localComputed_\n");

    /*
    //update GUI and MIP
    disableListener();

    //size spot
    circrad->setValue(nextdatasp[2]);
    //center and cursor
    locX->setValue(nextdatasp[3]);
    locY->setValue(nextdatasp[4]);
    locYT->setValue(nextdatasp[5]);
    locXL->setValue(nextdatasp[6]);
    centerX->setValue(nextdatasp[7]);
    centerY->setValue(nextdatasp[8]);

    //sliders
    lightness->setValue(nextdatasp[9]);
    contrast->setValue(nextdatasp[10]);
    chroma->setValue(nextdatasp[11]);
    sensi->setValue(nextdatasp[12]);
    transit->setValue(nextdatasp[13]);

    //inverse
    if (nextdatasp[14] == 0) {
        invers->set_active(false);
    } else {
        invers->set_active(true);
    }

    //method cursor
    if (nextdatasp[15] == 0) {
        Smethod->set_active(0);
    } else if (nextdatasp[15] == 1) {
        Smethod->set_active(1);
    } else if (nextdatasp[15] == 2) {
        Smethod->set_active(2);
    } else if (nextdatasp[15] == 3) {
        Smethod->set_active(3);
    }

    //nbspot->setValue(nextdatasp[16]);

    //sliders blurr
    radius->setValue(nextdatasp[17]);
    strength->setValue(nextdatasp[18]);
    sensibn->setValue(nextdatasp[19]);

    //inverse
    if (nextdatasp[20] == 0) {
        inversrad->set_active(false);
    } else {
        inversrad->set_active(true);
    }

    //sliders retinex
    str->setValue(nextdatasp[21]);
    chrrt->setValue(nextdatasp[22]);
    neigh->setValue(nextdatasp[23]);
    vart->setValue(nextdatasp[24]);
    sensih->setValue(nextdatasp[25]);

    //inverse
    if (nextdatasp[26] == 0) {
        inversret->set_active(false);
    } else {
        inversret->set_active(true);
    }

    //method retinex
    if (nextdatasp[27] == 0) {
        retinexMethod->set_active(0);
    } else if (nextdatasp[27] == 1) {
        retinexMethod->set_active(1);
    } else if (nextdatasp[27] == 2) {
        retinexMethod->set_active(2);
    }

    //sharpening
    sharradius->setValue(nextdatasp[28]);
    sharamount->setValue(nextdatasp[29]);
    shardamping->setValue(nextdatasp[30]);
    shariter->setValue(nextdatasp[31]);
    sensisha->setValue(nextdatasp[32]);

    if (nextdatasp[33] == 0) {
        inverssha->set_active(false);
    } else {
        inverssha->set_active(true);
    }

    if (nextdatasp[34] == 0) {
        qualityMethod->set_active(0);
    } else if (nextdatasp[34] == 1) {
        qualityMethod->set_active(1);
    } else if (nextdatasp[34] == 2) {
        qualityMethod->set_active(2);
    }

    thres->setValue(nextdatasp[35]);
    proxi->setValue(nextdatasp[36]);

    //denoise
    noiselumf->setValue(nextdatasp[37]);
    noiselumc->setValue(nextdatasp[38]);
    noisechrof->setValue(nextdatasp[39]);
    noisechroc->setValue(nextdatasp[40]);

    //cbdl
    multiplier[0]->setValue(nextdatasp[41]);
    multiplier[1]->setValue(nextdatasp[42]);
    multiplier[2]->setValue(nextdatasp[43]);
    multiplier[3]->setValue(nextdatasp[44]);
    multiplier[4]->setValue(nextdatasp[45]);
    threshold->setValue(nextdatasp[46]);
    sensicb->setValue(nextdatasp[47]);

    //blur luma
    if (nextdatasp[48] == 0) {
        activlum->set_active(false);
    } else {
        activlum->set_active(true);
    }

    //TM
    stren->setValue(nextdatasp[49]);
    gamma->setValue(nextdatasp[50]);
    estop->setValue(nextdatasp[51]);
    scaltm->setValue(nextdatasp[52]);
    rewei->setValue(nextdatasp[53]);
    sensitm->setValue(nextdatasp[54]);
    //  usleep(10000);

    //Reticurv
    retrab->setValue(nextdatasp[55]);

    //curvactiv
    if (nextdatasp[56] == 0) {
        curvactiv->set_active(false);
    } else {
        curvactiv->set_active(true);
    }

    if (nextdatasp[57] == 0) {
        qualitycurveMethod->set_active(0);
    } else if (nextdatasp[57] == 1) {
        qualitycurveMethod->set_active(1);
    } else if (nextdatasp[57] == 2) {
        qualitycurveMethod->set_active(2);
    }

    sensiv->setValue(nextdatasp[58]);
    pastels->setValue(nextdatasp[59]);
    saturated->setValue(nextdatasp[60]);

    //protectskin
    if (nextdatasp[61] == 0) {
        protectSkins->set_active(false);
    } else {
        protectSkins->set_active(true);
    }

    //avoidColorShift
    if (nextdatasp[62] == 0) {
        avoidColorShift->set_active(false);
    } else {
        avoidColorShift->set_active(true);
    }

    //pastSatTog
    if (nextdatasp[63] == 0) {
        pastSatTog->set_active(false);
        adjusterChanged(pastels, pastels->getValue());
        adjusterChanged(saturated, saturated->getValue());

    } else {
        pastSatTog->set_active(true);
        adjusterChanged(pastels, pastels->getValue());
        adjusterChanged(saturated, saturated->getValue());

    }

    expcomp->setValue(nextdatasp[64]);
    black->setValue(nextdatasp[65]);
    hlcompr->setValue(nextdatasp[66]);
    hlcomprthresh->setValue(nextdatasp[67]);
    shcompr->setValue(nextdatasp[68]);
    sensiex->setValue(nextdatasp[69]);

    centerXbuf->setValue(nextdatasp[70]);
    centerYbuf->setValue(nextdatasp[71]);
    adjblur->setValue(nextdatasp[72]);

    //protectskin
    if (nextdatasp[73] == 0) {
        cutpast->set_active(false);
    } else {
        cutpast->set_active(true);
    }

    chromacbdl->setValue(nextdatasp[74]);

    if (nextdatasp[75] == 0) {
        lastdust->set_active(false);
    } else {
        lastdust->set_active(true);
    }

    if (nextdatasp[76] == 0) {
        blurMethod->set_active(0);
    } else if (nextdatasp[76] == 1) {
        blurMethod->set_active(1);
    } else if (nextdatasp[76] == 2) {
        blurMethod->set_active(2);
    }

    if (nextdatasp[77] == 0) {
        dustMethod->set_active(0);
    } else if (nextdatasp[77] == 1) {
        dustMethod->set_active(1);
    } else if (nextdatasp[77] == 2) {
        dustMethod->set_active(2);
    }

    if (nextdatasp[78] == 0) {
        Exclumethod->set_active(0);
    } else if (nextdatasp[78] == 1) {
        Exclumethod->set_active(1);
    }

    sensiexclu->setValue(nextdatasp[79]);
    struc->setValue(nextdatasp[80]);
    warm->setValue(nextdatasp[81]);
    noiselumdetail->setValue(nextdatasp[82]);
    noisechrodetail->setValue(nextdatasp[83]);
    sensiden->setValue(nextdatasp[84]);

    //exp
    if (nextdatasp[85] == 0) {
        expdenoi->setEnabled(false);
    } else {
        expdenoi->setEnabled(true);
    }

    if (nextdatasp[86] == 0) {
        expcolor->setEnabled(false);
    } else {
        expcolor->setEnabled(true);
    }

    if (nextdatasp[87] == 0) {
        expvibrance->setEnabled(false);
    } else {
        expvibrance->setEnabled(true);
    }

    if (nextdatasp[88] == 0) {
        expblur->setEnabled(false);
    } else {
        expblur->setEnabled(true);
    }

    if (nextdatasp[89] == 0) {
        exptonemap->setEnabled(false);
    } else {
        exptonemap->setEnabled(true);
    }

    if (nextdatasp[90] == 0) {
        expreti->setEnabled(false);
    } else {
        expreti->setEnabled(true);
    }

    if (nextdatasp[91] == 0) {
        expsharp->setEnabled(false);
    } else {
        expsharp->setEnabled(true);
    }

    if (nextdatasp[92] == 0) {
        expcbdl->setEnabled(false);
    } else {
        expcbdl->setEnabled(true);
    }

    if (nextdatasp[93] == 0) {
        expexpose->setEnabled(false);
    } else {
        expexpose->setEnabled(true);
    }

    bilateral->setValue(nextdatasp[94]);
    noiselequal->setValue(nextdatasp[95]);

    if (nextdatasp[96] == 0) {
        shapemethod->set_active(0);
    } else if (nextdatasp[96] == 1) {
        shapemethod->set_active(1);
    }



    double intermedblur = 0.01 * (double) nextdatasp[nextlength - 5];
    huerefblur->setValue(intermedblur);
    double intermed = 0.01 * (double) nextdatasp[nextlength - 4];
    hueref->setValue(intermed);

    chromaref->setValue(nextdatasp[nextlength - 3]);
    lumaref->setValue(nextdatasp[nextlength - 2]);
    sobelref->setValue(nextdatasp[nextlength - 1]);

    int *s_datc;
    s_datc = new int[70];
    int siz;
    ImProcFunctions::strcurv_data(nextstr, s_datc, siz);


    std::vector<double>   creti;

    for (int j = 0; j < siz; j++) {
        creti.push_back((double)(s_datc[j]) / 1000.);
    }

    delete [] s_datc;

    cTgainshape->setCurve(creti);

    //LLcurv
    int *s_datcl;
    s_datcl = new int[70];
    int sizl;
    ImProcFunctions::strcurv_data(nextll_str, s_datcl, sizl);


    std::vector<double>   cll;

    for (int j = 0; j < sizl; j++) {
        cll.push_back((double)(s_datcl[j]) / 1000.);
    }

    delete [] s_datcl;
    llshape->setCurve(cll);

    //CCcurv
    int *s_datcc;
    s_datcc = new int[70];
    int sizc;
    ImProcFunctions::strcurv_data(nextcc_str, s_datcc, sizc);


    std::vector<double>   ccc;

    for (int j = 0; j < sizc; j++) {
        ccc.push_back((double)(s_datcc[j]) / 1000.);
    }

    delete [] s_datcc;
    ccshape->setCurve(ccc);


    //LHcurv
    int *s_datch;
    s_datch = new int[70];
    int sizh;
    ImProcFunctions::strcurv_data(nextlh_str, s_datch, sizh);


    std::vector<double>   clh;

    for (int j = 0; j < sizh; j++) {
        clh.push_back((double)(s_datch[j]) / 1000.);
    }

    delete [] s_datch;
    LHshape->setCurve(clh);

    //HHcurv
    int *s_datchh;
    s_datchh = new int[70];
    int sizhh;
    ImProcFunctions::strcurv_data(nexthh_str, s_datchh, sizhh);


    std::vector<double>   chh;

    for (int j = 0; j < sizhh; j++) {
        chh.push_back((double)(s_datchh[j]) / 1000.);
    }

    delete [] s_datchh;
    HHshape->setCurve(chh);

    //skinTonesCurve
    int *s_datcsk;
    s_datcsk = new int[70];
    int sizsk;
    ImProcFunctions::strcurv_data(nextsk_str, s_datcsk, sizsk);
    std::vector<double>   csk;

    for (int j = 0; j < sizsk; j++) {
        csk.push_back((double)(s_datcsk[j]) / 1000.);
    }

    delete [] s_datcsk;

    skinTonesCurve->setCurve(csk);

    //PSthreshold
    int sizps = 2;
    int s_datcps[sizps + 1];
    ImProcFunctions::strcurv_data(nextps_str, s_datcps, sizps);
    psThreshold->setValue(s_datcps[0], s_datcps[1]);


    //exCurve
    int *s_datcex;
    s_datcex = new int[70];
    int sizex;
    ImProcFunctions::strcurv_data(nextex_str, s_datcex, sizex);
    std::vector<double>   cex;

    for (int j = 0; j < sizex; j++) {
        cex.push_back((double)(s_datcex[j]) / 1000.);
    }

    delete [] s_datcex;
    shapeexpos->setCurve(cex);


    enableListener();

    //update all sliders by this strange process!
    if (anbspot->getValue() == 0) {
        anbspot->setValue(1);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 1);

    } else if (anbspot->getValue() == 1) {
        anbspot->setValue(0);

        if (options.rtSettings.locdelay) {
            if (anbspot->delay < 100) {
                anbspot->delay = 100;
            }
        }

        adjusterChanged(anbspot, 0);

    }


    //update all curves
    std::vector<double>   cretirab;
    cretirab = cTgainshaperab->getCurve();

    if (cretirab.at(5) == 0.70) {
        cretirab.at(5) = 0.9;
        cTgainshaperab->setCurve(cretirab);

        curveChanged(cTgainshaperab);
    } else if (cretirab.at(5) == 0.90) {
        cretirab.at(5) = 0.7;
        cTgainshaperab->setCurve(cretirab);
        curveChanged(cTgainshaperab);

    }

    //

    //   printf("G1 maj anbspot=%i  cretirab=%f\n", anbspot->getValue(), cretirab.at(5));


    //add events for each cases whitout that localalb does not work...
    //there is probably an other solution !

    if (listener) { //for all sliders
        listener->panelChanged(Evlocallabanbspot, ""); //anbspot->getTextValue());
    }

    if (listener) {//for curve
        listener->panelChanged(EvlocallabCTgainCurverab, M(""));
    }

    if (listener) {//for inverse color
        listener->panelChanged(Evlocallabinvers, M("GENERAL_ENABLED"));
    }

    if (listener) {//for cutpast
        //   listener->panelChanged (Evlocallabcutpast, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for cutpast
        //    listener->panelChanged (Evlocallablastdust, M ("GENERAL_ENABLED"));
    }

    if (listener) {//for blur method
        listener->panelChanged(EvlocallabblurMethod, blurMethod->get_active_text());
    }

    if (listener) {//for dust method
        //    listener->panelChanged (EvlocallabdustMethod, dustMethod->get_active_text ());
    }

    if (listener) {//for Exclude method
        listener->panelChanged(Evlocallabexclumethod, Exclumethod->get_active_text());
    }

    if (listener) {//for shape RT-spot method
        listener->panelChanged(Evlocallabshapemethod, shapemethod->get_active_text());
    }

    if (listener) {//for curvactiv
        listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse blurr
        listener->panelChanged(Evlocallabinversrad, M("GENERAL_ENABLED"));
    }

    if (listener) {//for quality method
        listener->panelChanged(EvlocallabqualityMethod, qualityMethod->get_active_text());
    }

    if (listener) {//for quality method
        listener->panelChanged(EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text());
    }

    if (listener) {//for inverse retinex
        listener->panelChanged(Evlocallabinversret, M("GENERAL_ENABLED"));
    }

    if (listener) {//for inverse sharpen
        listener->panelChanged(Evlocallabinverssha, M("GENERAL_ENABLED"));
    }

    if (listener) {//for Smethod : position of mouse cursor
        listener->panelChanged(EvlocallabSmet, Smethod->get_active_text());
    }

    if (listener) {//for retinex method
        listener->panelChanged(EvlocallabretinexMethod, retinexMethod->get_active_text());
    }

    if (listener) {//for curve reti
        listener->panelChanged(EvlocallabCTgainCurve, M(""));
    }

    if (listener) {//for curve LL
        listener->panelChanged(Evlocallabllshape, M(""));
    }

    if (listener) {//for curve LH
        listener->panelChanged(EvlocallabLHshape, M(""));
    }

    if (listener) {//for curve LH
        listener->panelChanged(Evlocallabccshape, M(""));
    }

    if (listener) {//for curve LH
        listener->panelChanged(EvlocallabHHshape, M(""));
    }

    if (listener) {//for curve Skin
        listener->panelChanged(EvlocallabSkinTonesCurve, M(""));
    }

    if (listener) {//for PSthreshold
        listener->panelChanged(EvlocallabPastSatThreshold, M(""));
    }

    //for checkbox
    if (listener) {//f
        listener->panelChanged(EvlocallabProtectSkins, M(""));
    }

    if (listener) {//fo
        listener->panelChanged(EvlocallabAvoidColorShift, M(""));
    }

    if (listener) {//for
        listener->panelChanged(EvlocallabPastSatTog, M(""));
    }

    if (listener) {//for expo curv
        listener->panelChanged(Evlocallabshapeexpos, M(""));
    }

    if (listener) {//for expander denoise
        listener->panelChanged(EvLocenadenoi, M(""));
    }

    if (listener) {//for expander color
        listener->panelChanged(EvLocenacolor, M(""));
    }

    if (listener) {//for expander vib
        listener->panelChanged(EvLocenavibrance, M(""));
    }

    if (listener) {//for expander tonem
        listener->panelChanged(EvLocenatonemap, M(""));
    }

    if (listener) {//for expander reti
        listener->panelChanged(EvLocenareti, M(""));
    }

    if (listener) {//for expander sharp
        listener->panelChanged(EvLocenasharp, M(""));
    }

    if (listener) {//for expander cbdl
        listener->panelChanged(EvLocenacbdl, M(""));
    }

    if (listener) {//for expander cbdl
        listener->panelChanged(EvLocenaexpose, M(""));
    }
    */


    return false;
}

void Locallab::localChanged(int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str, std::string ex_str, int sp, int maxdat)
{
    // TODO Locallab printf
    printf("localChanged\n");

    /*
    nextstr = datastr;
    nextll_str = ll_str;
    nextlh_str = lh_str;
    nextcc_str = cc_str;
    nexthh_str = hh_str;
    nextsk_str = sk_str;
    nextps_str = ps_str;
    nextex_str = ex_str;
    nextlength = maxdat;

    for (int i = 2; i < nextlength; i++) {//be carefull to this value 102 must be the same as above ==> sobelref->setValue(nextdatasp[101]);
        nextdatasp[i] = datasp[i][sp];
    }

    g_idle_add(localChangedUI, this);
    */
}

void Locallab::localretChanged(int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str, std::string ex_str, int sp, int maxdat)
{
    // TODO Locallab printf
    printf("localretChanged\n");

    /*
    nextlength = maxdat;
    nextstr2 = datastr;
    nextll_str2 = ll_str;
    nextlh_str2 = lh_str;
    nextcc_str2 = cc_str;
    nexthh_str2 = hh_str;
    nextsk_str2 = sk_str;
    nextps_str2 = ps_str;
    nextex_str2 = ex_str;

    g_idle_add(localretChangedUI, this);
    */
}


void Locallab::read(const ProcParams* pp, const ParamsEdited* pedited)
{
    printf("Locallab read\n");

    // Disable all listeners
    disableListener();

    if (pedited) {
        set_inconsistent(multiImage && !pedited->locallab.enabled);

        // Control spot settings
        ControlSpotPanel::SpotEdited* const se = new ControlSpotPanel::SpotEdited();

        if (pedited->locallab.nbspot) {
            se->addbutton = true;
            se->deletebutton = true;
        } else {
            se->addbutton = false;
            se->deletebutton = false;
        }

        se->treeview = pedited->locallab.nbspot || pedited->locallab.selspot;
        se->name = pedited->locallab.name;
        se->isvisible = pedited->locallab.isvisible;
        se->shape = pedited->locallab.shape;
        se->spotMethod = pedited->locallab.spotMethod;
        se->shapeMethod = pedited->locallab.shapeMethod;
        se->locX = pedited->locallab.locX;
        se->locXL = pedited->locallab.locXL;
        se->locY = pedited->locallab.locY;
        se->locYT = pedited->locallab.locYT;
        se->centerX = pedited->locallab.centerX;
        se->centerY = pedited->locallab.centerY;
        se->circrad = pedited->locallab.circrad;
        se->qualityMethod = pedited->locallab.qualityMethod;
        se->transit = pedited->locallab.transit;
        se->thresh = pedited->locallab.thresh;
        se->iter = pedited->locallab.iter;
        expsettings->setEditedStates(se);

        // Color & Light
        expcolor->set_inconsistent(!pedited->locallab.expcolor);
        curvactiv->set_inconsistent(multiImage && !pedited->locallab.curvactiv);
        lightness->setEditedState(pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setEditedState(pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setEditedState(pedited->locallab.chroma ? Edited : UnEdited);
        sensi->setEditedState(pedited->locallab.sensi ? Edited : UnEdited);

        if (!pedited->locallab.qualitycurveMethod) {
            qualitycurveMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        llshape->setUnChanged(!pedited->locallab.llcurve);
        ccshape->setUnChanged(!pedited->locallab.cccurve);
        LHshape->setUnChanged(!pedited->locallab.LHcurve);
        HHshape->setUnChanged(!pedited->locallab.HHcurve);
        invers->set_inconsistent(multiImage && !pedited->locallab.invers);

        // Exposure
        expexpose->set_inconsistent(!pedited->locallab.expexpose);
        expcomp->setEditedState(pedited->locallab.expcomp ? Edited : UnEdited);
        hlcompr->setEditedState(pedited->locallab.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setEditedState(pedited->locallab.hlcomprthresh ? Edited : UnEdited);
        black->setEditedState(pedited->locallab.black ? Edited : UnEdited);
        warm->setEditedState(pedited->locallab.warm ? Edited : UnEdited);
        shcompr->setEditedState(pedited->locallab.shcompr ? Edited : UnEdited);
        sensiex->setEditedState(pedited->locallab.sensiex ? Edited : UnEdited);
        shapeexpos->setUnChanged(!pedited->locallab.excurve);

        // Vibrance
        expvibrance->set_inconsistent(!pedited->locallab.expvibrance);
        saturated->setEditedState(pedited->locallab.saturated ? Edited : UnEdited);
        pastels->setEditedState(pedited->locallab.pastels ? Edited : UnEdited);
        psThreshold->setEditedState(pedited->locallab.psthreshold ? Edited : UnEdited);
        protectSkins->set_inconsistent(!pedited->locallab.protectskins);
        avoidColorShift->set_inconsistent(!pedited->locallab.avoidcolorshift);
        pastSatTog->set_inconsistent(!pedited->locallab.pastsattog);
        sensiv->setEditedState(pedited->locallab.sensiv ? Edited : UnEdited);
        skinTonesCurve->setUnChanged(!pedited->locallab.skintonescurve);

        // Blur & Noise
        expblur->set_inconsistent(!pedited->locallab.expblur);
        radius->setEditedState(pedited->locallab.radius ? Edited : UnEdited);
        strength->setEditedState(pedited->locallab.strength ? Edited : UnEdited);
        sensibn->setEditedState(pedited->locallab.sensibn ? Edited : UnEdited);

        if (!pedited->locallab.blurMethod) {
            blurMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        activlum->set_inconsistent(multiImage && !pedited->locallab.activlum);

        // Tone Mapping
        exptonemap->set_inconsistent(!pedited->locallab.exptonemap);
        stren->setEditedState(pedited->locallab.stren ? Edited : UnEdited);
        gamma->setEditedState(pedited->locallab.gamma ? Edited : UnEdited);
        estop->setEditedState(pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setEditedState(pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setEditedState(pedited->locallab.rewei ? Edited : UnEdited);
        sensitm->setEditedState(pedited->locallab.sensitm ? Edited : UnEdited);

        // Retinex
        expreti->set_inconsistent(!pedited->locallab.expreti);

        if (!pedited->locallab.retinexMethod) {
            retinexMethod->set_active_text(M("GENERAL_UNCHANGED"));
        }

        str->setEditedState(pedited->locallab.str ? Edited : UnEdited);
        chrrt->setEditedState(pedited->locallab.chrrt ? Edited : UnEdited);
        neigh->setEditedState(pedited->locallab.neigh ? Edited : UnEdited);
        vart->setEditedState(pedited->locallab.vart ? Edited : UnEdited);
        sensih->setEditedState(pedited->locallab.sensih ? Edited : UnEdited);
        cTgainshape->setUnChanged(!pedited->locallab.localTgaincurve);
        inversret->set_inconsistent(multiImage && !pedited->locallab.inversret);

        // Sharpening
        expsharp->set_inconsistent(!pedited->locallab.expsharp);
        sharradius->setEditedState(pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setEditedState(pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setEditedState(pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setEditedState(pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setEditedState(pedited->locallab.sensisha ? Edited : UnEdited);
        inverssha->set_inconsistent(multiImage && !pedited->locallab.inverssha);

        // Contrast by detail levels
        expcbdl->set_inconsistent(!pedited->locallab.expcbdl);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setEditedState(pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        chromacbdl->setEditedState(pedited->locallab.chromacbdl ? Edited : UnEdited);
        threshold->setEditedState(pedited->locallab.threshold ? Edited : UnEdited);
        sensicb->setEditedState(pedited->locallab.sensicb ? Edited : UnEdited);

        // Denoise
        expdenoi->set_inconsistent(!pedited->locallab.expdenoi);
        noiselumf->setEditedState(pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setEditedState(pedited->locallab.noiselumc ? Edited : UnEdited);
        noiselumdetail->setEditedState(pedited->locallab.noiselumdetail ? Edited : UnEdited);
        noiselequal->setEditedState(pedited->locallab.noiselequal ? Edited : UnEdited);
        noisechrof->setEditedState(pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setEditedState(pedited->locallab.noisechroc ? Edited : UnEdited);
        adjblur->setEditedState(pedited->locallab.adjblur ? Edited : UnEdited);
        bilateral->setEditedState(pedited->locallab.bilateral ? Edited : UnEdited);
        sensiden->setEditedState(pedited->locallab.sensiden ? Edited : UnEdited);
        avoid->set_inconsistent(multiImage && !pedited->locallab.avoid);
    }

    setEnabled(pp->locallab.enabled);

    // Add non existent spots and update existent ones
    ControlSpotPanel::SpotRow* const r = new ControlSpotPanel::SpotRow();

    for (int i = 0; i < pp->locallab.nbspot; i++) {
        r->id = pp->locallab.id.at(i);
        r->name = pp->locallab.name.at(i);
        r->isvisible = pp->locallab.isvisible.at(i);

        if (pp->locallab.shape.at(i) == "ELI") {
            r->shape = 0;
        } else {
            r->shape = 1;
        }

        if (pp->locallab.spotMethod.at(i) == "norm") {
            r->spotMethod = 0;
        } else {
            r->spotMethod = 1;
        }

        if (pp->locallab.shapeMethod.at(i) == "IND") {
            r->shapeMethod = 0;
        } else if (pp->locallab.shapeMethod.at(i) == "SYM") {
            r->shapeMethod = 1;
        } else if (pp->locallab.shapeMethod.at(i) == "INDSL") {
            r->shapeMethod = 2;
        } else {
            r->shapeMethod = 3;
        }

        r->locX = pp->locallab.locX.at(i);
        r->locXL = pp->locallab.locXL.at(i);
        r->locY = pp->locallab.locY.at(i);
        r->locYT = pp->locallab.locYT.at(i);
        r->centerX = pp->locallab.centerX.at(i);
        r->centerY = pp->locallab.centerY.at(i);
        r->circrad = pp->locallab.circrad.at(i);

        if (pp->locallab.qualityMethod.at(i) == "std") {
            r->qualityMethod = 0;
        } else if (pp->locallab.qualityMethod.at(i) == "enh") {
            r->qualityMethod = 1;
        } else {
            r->qualityMethod = 2;
        }

        r->transit = pp->locallab.transit.at(i);
        r->thresh = pp->locallab.thresh.at(i);
        r->iter = pp->locallab.iter.at(i);

        if (!expsettings->updateControlSpot(r)) {
            expsettings->addControlSpot(r);
        }
    }

    // Delete not anymore existent spots
    std::vector<int>* const list = expsettings->getSpotIdList();
    bool ispresent;

    for (int i = 0; i < (int)list->size(); i++) {
        ispresent = false;

        for (int j = 0; j < pp->locallab.nbspot; j++) {
            if (list->at(i) == pp->locallab.id.at(j)) {
                ispresent = true;
                break;
            }
        }

        if (!ispresent) {
            expsettings->deleteControlSpot(list->at(i));
        }
    }

    // Select active spot
    if (pp->locallab.nbspot > 0) {
        expsettings->setSelectedSpot(pp->locallab.id.at(pp->locallab.selspot));
    }

    // Update Locallab tools GUI
    updateLocallabGUI(pp, pp->locallab.selspot);

    // Enable all listeners
    enableListener();
}

void Locallab::write(ProcParams* pp, ParamsEdited* pedited)
{
    printf("Locallab write\n");

    pp->locallab.enabled = getEnabled();

    const int spotPanelEvent = expsettings->getEventType();
    int spotId;
    ControlSpotPanel::SpotRow* r;

    switch (spotPanelEvent) {
        case (1): // 1 = Spot creation event
            // Spot creation (default initialization)
            spotId = expsettings->getNewId();
            r = new ControlSpotPanel::SpotRow();
            r->id = spotId;
            r->name = "Control Spot #" + std::to_string(spotId);
            r->isvisible = true;
            r->shape = 0;
            r->spotMethod = 0;
            r->shapeMethod = 2;
            r->locX = 250;
            r->locXL = 250;
            r->locY = 250;
            r->locYT = 250;
            r->centerX = 0;
            r->centerY = 0;
            r->circrad = 18;
            r->qualityMethod = 0;
            r->transit = 60;
            r->thresh = 18;
            r->iter = 0;
            expsettings->addControlSpot(r);

            // ProcParams update
            pp->locallab.nbspot++;
            pp->locallab.selspot = pp->locallab.nbspot - 1;
            // Control spot settings
            pp->locallab.id.push_back(r->id);
            pp->locallab.name.push_back(r->name);
            pp->locallab.isvisible.push_back(r->isvisible);
            pp->locallab.shape.push_back("ELI");
            pp->locallab.spotMethod.push_back("norm");
            pp->locallab.shapeMethod.push_back("INDSL");
            pp->locallab.locX.push_back(r->locX);
            pp->locallab.locXL.push_back(r->locXL);
            pp->locallab.locY.push_back(r->locY);
            pp->locallab.locYT.push_back(r->locYT);
            pp->locallab.centerX.push_back(r->centerX);
            pp->locallab.centerY.push_back(r->centerY);
            pp->locallab.circrad.push_back(r->circrad);
            pp->locallab.qualityMethod.push_back("std");
            pp->locallab.transit.push_back(r->transit);
            pp->locallab.thresh.push_back(r->thresh);
            pp->locallab.iter.push_back(r->iter);
            // Color & Light
            pp->locallab.expcolor.push_back(0);
            pp->locallab.curvactiv.push_back(0);
            pp->locallab.lightness.push_back(0);
            pp->locallab.contrast.push_back(0);
            pp->locallab.chroma.push_back(0);
            pp->locallab.sensi.push_back(19);
            pp->locallab.qualitycurveMethod.push_back("none");
            pp->locallab.llcurve.push_back({(double)DCT_NURBS, 0.0, 0.0, 1.0, 1.0});
            pp->locallab.cccurve.push_back({(double)DCT_NURBS, 0.0, 0.0, 1.0, 1.0});
            pp->locallab.LHcurve.push_back({(double)FCT_MinMaxCPoints, 0.0, 0.50, 0.35, 0.35, 0.166, 0.50, 0.35, 0.35, 0.333, 0.50, 0.35, 0.35, 0.50, 0.50, 0.35, 0.35, 0.666, 0.50, 0.35, 0.35, 0.833, 0.50, 0.35, 0.35});
            pp->locallab.HHcurve.push_back({(double)FCT_MinMaxCPoints, 0.0, 0.50, 0.35, 0.35, 0.166, 0.50, 0.35, 0.35, 0.333, 0.50, 0.35, 0.35, 0.50, 0.50, 0.35, 0.35, 0.666, 0.50, 0.35, 0.35, 0.833, 0.50, 0.35, 0.35});
            pp->locallab.invers.push_back(0);
            // Exposure
            pp->locallab.expexpose.push_back(0);
            pp->locallab.expcomp.push_back(0);
            pp->locallab.hlcompr.push_back(20);
            pp->locallab.hlcomprthresh.push_back(33);
            pp->locallab.black.push_back(0);
            pp->locallab.shcompr.push_back(50);
            pp->locallab.warm.push_back(0);
            pp->locallab.sensiex.push_back(19);
            pp->locallab.excurve.push_back({(double)DCT_NURBS, 0.0, 0.0, 1.0, 1.0});
            // Vibrance
            pp->locallab.expvibrance.push_back(0);
            pp->locallab.saturated.push_back(0);
            pp->locallab.pastels.push_back(0);
            pp->locallab.psthreshold.push_back({0, 75, false});
            pp->locallab.protectskins.push_back(0);
            pp->locallab.avoidcolorshift.push_back(1);
            pp->locallab.pastsattog.push_back(1);
            pp->locallab.sensiv.push_back(19);
            pp->locallab.skintonescurve.push_back({(double)DCT_Linear});
            // Blur & Noise
            pp->locallab.expblur.push_back(0);
            pp->locallab.radius.push_back(1);
            pp->locallab.strength.push_back(0);
            pp->locallab.sensibn.push_back(40);
            pp->locallab.blurMethod.push_back("norm");
            pp->locallab.activlum.push_back(0);
            // Tone Mapping
            pp->locallab.exptonemap.push_back(0);
            pp->locallab.stren.push_back(1);
            pp->locallab.gamma.push_back(100);
            pp->locallab.estop.push_back(140);
            pp->locallab.scaltm.push_back(10);
            pp->locallab.rewei.push_back(0);
            pp->locallab.sensitm.push_back(19);
            // Retinex
            pp->locallab.expreti.push_back(0);
            pp->locallab.retinexMethod.push_back("high");
            pp->locallab.str.push_back(0);
            pp->locallab.chrrt.push_back(0);
            pp->locallab.neigh.push_back(50);
            pp->locallab.vart.push_back(200);
            pp->locallab.sensih.push_back(19);
            pp->locallab.localTgaincurve.push_back({(double)FCT_MinMaxCPoints, 0.0, 0.12, 0.35, 0.35, 0.70, 0.50, 0.35, 0.35, 1.00, 0.12, 0.35, 0.35});
            pp->locallab.inversret.push_back(0);
            // Sharpening
            pp->locallab.expsharp.push_back(0);
            pp->locallab.sharradius.push_back(40);
            pp->locallab.sharamount.push_back(75);
            pp->locallab.shardamping.push_back(75);
            pp->locallab.shariter.push_back(30);
            pp->locallab.sensisha.push_back(19);
            pp->locallab.inverssha.push_back(0);
            // Contrast by detail levels
            pp->locallab.expcbdl.push_back(0);

            for (int i = 0; i < 5; i++) {
                pp->locallab.mult[i].push_back(100.0);
            }

            pp->locallab.chromacbdl.push_back(0);
            pp->locallab.threshold.push_back(20.0);
            pp->locallab.sensicb.push_back(19);
            // Denoise
            pp->locallab.expdenoi.push_back(0);
            pp->locallab.noiselumf.push_back(0);
            pp->locallab.noiselumc.push_back(0);
            pp->locallab.noiselumdetail.push_back(0);
            pp->locallab.noiselequal.push_back(7);
            pp->locallab.noisechrof.push_back(0);
            pp->locallab.noisechroc.push_back(0);
            pp->locallab.adjblur.push_back(0);
            pp->locallab.bilateral.push_back(0);
            pp->locallab.sensiden.push_back(30);
            pp->locallab.avoid.push_back(0);

            // New created spot selection
            expsettings->setSelectedSpot(spotId);

            // Update Locallab tools GUI with new created spot
            disableListener();
            updateLocallabGUI(pp, pp->locallab.selspot);
            enableListener();

            break;

        case (2): // 2 = Spot deletion event
            // Get deleted spot index in ProcParams and update it
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot; i++) {
                if (pp->locallab.id.at(i) == spotId) {
                    // ProcParams update
                    pp->locallab.nbspot--;
                    pp->locallab.selspot = 0;
                    // Control spot settings
                    pp->locallab.id.erase(pp->locallab.id.begin() + i);
                    pp->locallab.name.erase(pp->locallab.name.begin() + i);
                    pp->locallab.isvisible.erase(pp->locallab.isvisible.begin() + i);
                    pp->locallab.shape.erase(pp->locallab.shape.begin() + i);
                    pp->locallab.spotMethod.erase(pp->locallab.spotMethod.begin() + i);
                    pp->locallab.shapeMethod.erase(pp->locallab.shapeMethod.begin() + i);
                    pp->locallab.locX.erase(pp->locallab.locX.begin() + i);
                    pp->locallab.locXL.erase(pp->locallab.locXL.begin() + i);
                    pp->locallab.locY.erase(pp->locallab.locY.begin() + i);
                    pp->locallab.locYT.erase(pp->locallab.locYT.begin() + i);
                    pp->locallab.centerX.erase(pp->locallab.centerX.begin() + i);
                    pp->locallab.centerY.erase(pp->locallab.centerY.begin() + i);
                    pp->locallab.circrad.erase(pp->locallab.circrad.begin() + i);
                    pp->locallab.qualityMethod.erase(pp->locallab.qualityMethod.begin() + i);
                    pp->locallab.transit.erase(pp->locallab.transit.begin() + i);
                    pp->locallab.thresh.erase(pp->locallab.thresh.begin() + i);
                    pp->locallab.iter.erase(pp->locallab.iter.begin() + i);
                    expsettings->deleteControlSpot(spotId);
                    // Color & Light
                    pp->locallab.expcolor.erase(pp->locallab.expcolor.begin() + i);
                    pp->locallab.curvactiv.erase(pp->locallab.curvactiv.begin() + i);
                    pp->locallab.lightness.erase(pp->locallab.lightness.begin() + i);
                    pp->locallab.contrast.erase(pp->locallab.contrast.begin() + i);
                    pp->locallab.chroma.erase(pp->locallab.chroma.begin() + i);
                    pp->locallab.sensi.erase(pp->locallab.sensi.begin() + i);
                    pp->locallab.qualitycurveMethod.erase(pp->locallab.qualitycurveMethod.begin() + i);
                    pp->locallab.llcurve.erase(pp->locallab.llcurve.begin() + i);
                    pp->locallab.cccurve.erase(pp->locallab.cccurve.begin() + i);
                    pp->locallab.LHcurve.erase(pp->locallab.LHcurve.begin() + i);
                    pp->locallab.HHcurve.erase(pp->locallab.HHcurve.begin() + i);
                    pp->locallab.invers.erase(pp->locallab.invers.begin() + i);
                    // Exposure
                    pp->locallab.expexpose.erase(pp->locallab.expexpose.begin() + i);
                    pp->locallab.expcomp.erase(pp->locallab.expcomp.begin() + i);
                    pp->locallab.hlcompr.erase(pp->locallab.hlcompr.begin() + i);
                    pp->locallab.hlcomprthresh.erase(pp->locallab.hlcomprthresh.begin() + i);
                    pp->locallab.black.erase(pp->locallab.black.begin() + i);
                    pp->locallab.shcompr.erase(pp->locallab.shcompr.begin() + i);
                    pp->locallab.warm.erase(pp->locallab.warm.begin() + i);
                    pp->locallab.sensiex.erase(pp->locallab.sensiex.begin() + i);
                    pp->locallab.excurve.erase(pp->locallab.excurve.begin() + i);
                    // Vibrance
                    pp->locallab.expvibrance.erase(pp->locallab.expvibrance.begin() + i);
                    pp->locallab.saturated.erase(pp->locallab.saturated.begin() + i);
                    pp->locallab.pastels.erase(pp->locallab.pastels.begin() + i);
                    pp->locallab.psthreshold.erase(pp->locallab.psthreshold.begin() + i);
                    pp->locallab.protectskins.erase(pp->locallab.protectskins.begin() + i);
                    pp->locallab.avoidcolorshift.erase(pp->locallab.avoidcolorshift.begin() + i);
                    pp->locallab.pastsattog.erase(pp->locallab.pastsattog.begin() + i);
                    pp->locallab.sensiv.erase(pp->locallab.sensiv.begin() + i);
                    pp->locallab.skintonescurve.erase(pp->locallab.skintonescurve.begin() + i);
                    // Blur & Noise
                    pp->locallab.expblur.erase(pp->locallab.expblur.begin() + i);
                    pp->locallab.radius.erase(pp->locallab.radius.begin() + i);
                    pp->locallab.strength.erase(pp->locallab.strength.begin() + i);
                    pp->locallab.sensibn.erase(pp->locallab.sensibn.begin() + i);
                    pp->locallab.blurMethod.erase(pp->locallab.blurMethod.begin() + i);
                    pp->locallab.activlum.erase(pp->locallab.activlum.begin() + i);
                    // Tone Mapping
                    pp->locallab.exptonemap.erase(pp->locallab.exptonemap.begin() + i);
                    pp->locallab.stren.erase(pp->locallab.stren.begin() + i);
                    pp->locallab.gamma.erase(pp->locallab.gamma.begin() + i);
                    pp->locallab.estop.erase(pp->locallab.estop.begin() + i);
                    pp->locallab.scaltm.erase(pp->locallab.scaltm.begin() + i);
                    pp->locallab.rewei.erase(pp->locallab.rewei.begin() + i);
                    pp->locallab.sensitm.erase(pp->locallab.sensitm.begin() + i);
                    // Retinex
                    pp->locallab.expreti.erase(pp->locallab.expreti.begin() + i);
                    pp->locallab.retinexMethod.erase(pp->locallab.retinexMethod.begin() + i);
                    pp->locallab.str.erase(pp->locallab.str.begin() + i);
                    pp->locallab.chrrt.erase(pp->locallab.chrrt.begin() + i);
                    pp->locallab.neigh.erase(pp->locallab.neigh.begin() + i);
                    pp->locallab.vart.erase(pp->locallab.vart.begin() + i);
                    pp->locallab.sensih.erase(pp->locallab.sensih.begin() + i);
                    pp->locallab.localTgaincurve.erase(pp->locallab.localTgaincurve.begin() + i);
                    pp->locallab.inversret.erase(pp->locallab.inversret.begin() + i);
                    // Sharpening
                    pp->locallab.expsharp.erase(pp->locallab.expsharp.begin() + i);
                    pp->locallab.sharradius.erase(pp->locallab.sharradius.begin() + i);
                    pp->locallab.sharamount.erase(pp->locallab.sharamount.begin() + i);
                    pp->locallab.shardamping.erase(pp->locallab.shardamping.begin() + i);
                    pp->locallab.shariter.erase(pp->locallab.shariter.begin() + i);
                    pp->locallab.sensisha.erase(pp->locallab.sensisha.begin() + i);
                    pp->locallab.inverssha.erase(pp->locallab.inverssha.begin() + i);
                    // Contrast by detail levels
                    pp->locallab.expcbdl.erase(pp->locallab.expcbdl.begin() + i);

                    for (int j = 0; j < 5; j++) {
                        pp->locallab.mult[j].erase(pp->locallab.mult[j].begin() + i);
                    }

                    pp->locallab.chromacbdl.erase(pp->locallab.chromacbdl.begin() + i);
                    pp->locallab.threshold.erase(pp->locallab.threshold.begin() + i);
                    pp->locallab.sensicb.erase(pp->locallab.sensicb.begin() + i);
                    // Denoise
                    pp->locallab.expdenoi.erase(pp->locallab.expdenoi.begin() + i);
                    pp->locallab.noiselumf.erase(pp->locallab.noiselumf.begin() + i);
                    pp->locallab.noiselumc.erase(pp->locallab.noiselumc.begin() + i);
                    pp->locallab.noiselumdetail.erase(pp->locallab.noiselumdetail.begin() + i);
                    pp->locallab.noiselequal.erase(pp->locallab.noiselequal.begin() + i);
                    pp->locallab.noisechrof.erase(pp->locallab.noisechrof.begin() + i);
                    pp->locallab.noisechroc.erase(pp->locallab.noisechroc.begin() + i);
                    pp->locallab.adjblur.erase(pp->locallab.adjblur.begin() + i);
                    pp->locallab.bilateral.erase(pp->locallab.bilateral.begin() + i);
                    pp->locallab.sensiden.erase(pp->locallab.sensiden.begin() + i);
                    pp->locallab.avoid.erase(pp->locallab.avoid.begin() + i);

                    // Select one remaining spot
                    if (pp->locallab.nbspot > 0) {
                        expsettings->setSelectedSpot(pp->locallab.id.at(pp->locallab.selspot));
                    }

                    // Update Locallab tools GUI with new created spot
                    disableListener();
                    updateLocallabGUI(pp, pp->locallab.selspot);
                    enableListener();

                    break;
                }
            }

            break;

        case (3):  // 3 = Spot selection event
            spotId = expsettings->getSelectedSpot();

            for (int i = 0; i < pp->locallab.nbspot; i++) {
                if (pp->locallab.id.at(i) == spotId) {
                    pp->locallab.selspot = i;
                    break;
                }
            }

            // Update control spots and Locallab tools GUI with selected spot
            expsettings->setSelectedSpot(spotId);
            disableListener();
            updateLocallabGUI(pp, pp->locallab.selspot);
            enableListener();

            break;

        default: // Spot or locallab GUI updated
            if (pp->locallab.nbspot > 0) {
                r = expsettings->getSpot(expsettings->getSelectedSpot());

                // ProcParams update
                // Control spot settings
                pp->locallab.name.at(pp->locallab.selspot) = r->name;
                pp->locallab.isvisible.at(pp->locallab.selspot) = r->isvisible;

                if (r->shape == 0) {
                    pp->locallab.shape.at(pp->locallab.selspot) = "ELI";
                } else {
                    pp->locallab.shape.at(pp->locallab.selspot) = "RECT";
                }

                if (r->spotMethod == 0) {
                    pp->locallab.spotMethod.at(pp->locallab.selspot) = "norm";
                } else {
                    pp->locallab.spotMethod.at(pp->locallab.selspot) = "exc";
                }

                if (r->shapeMethod == 0) {
                    pp->locallab.shapeMethod.at(pp->locallab.selspot) = "IND";
                } else if (r->shapeMethod == 1) {
                    pp->locallab.shapeMethod.at(pp->locallab.selspot) = "SYM";
                } else if (r->shapeMethod == 2) {
                    pp->locallab.shapeMethod.at(pp->locallab.selspot) = "INDSL";
                } else {
                    pp->locallab.shapeMethod.at(pp->locallab.selspot) = "SYMSL";
                }

                pp->locallab.locX.at(pp->locallab.selspot) = r->locX;
                pp->locallab.locXL.at(pp->locallab.selspot) = r->locXL;
                pp->locallab.locY.at(pp->locallab.selspot) = r->locY;
                pp->locallab.locYT.at(pp->locallab.selspot) = r->locYT;
                pp->locallab.centerX.at(pp->locallab.selspot) = r->centerX;
                pp->locallab.centerY.at(pp->locallab.selspot) = r->centerY;
                pp->locallab.circrad.at(pp->locallab.selspot) = r->circrad;

                if (r->qualityMethod == 0) {
                    pp->locallab.qualityMethod.at(pp->locallab.selspot) = "std";
                } else if (r->qualityMethod == 1) {
                    pp->locallab.qualityMethod.at(pp->locallab.selspot) = "enh";
                } else {
                    pp->locallab.qualityMethod.at(pp->locallab.selspot) = "enhden";
                }

                pp->locallab.transit.at(pp->locallab.selspot) = r->transit;
                pp->locallab.thresh.at(pp->locallab.selspot) = r->thresh;
                pp->locallab.iter.at(pp->locallab.selspot) = r->iter;
                // Color & Light
                pp->locallab.expcolor.at(pp->locallab.selspot) = (int)expcolor->getEnabled();
                pp->locallab.curvactiv.at(pp->locallab.selspot) = (int)curvactiv->get_active();
                pp->locallab.lightness.at(pp->locallab.selspot) = lightness->getIntValue();
                pp->locallab.contrast.at(pp->locallab.selspot) = contrast->getIntValue();
                pp->locallab.chroma.at(pp->locallab.selspot) = chroma->getIntValue();
                pp->locallab.sensi.at(pp->locallab.selspot) = sensi->getIntValue();

                if (qualitycurveMethod->get_active_row_number() == 0) {
                    pp->locallab.qualitycurveMethod.at(pp->locallab.selspot) = "none";
                } else if (qualitycurveMethod->get_active_row_number() == 1) {
                    pp->locallab.qualitycurveMethod.at(pp->locallab.selspot) = "std";
                } else if (qualitycurveMethod->get_active_row_number() == 2) {
                    pp->locallab.qualitycurveMethod.at(pp->locallab.selspot) = "enh";
                }

                pp->locallab.llcurve.at(pp->locallab.selspot) = llshape->getCurve();
                pp->locallab.cccurve.at(pp->locallab.selspot) = ccshape->getCurve();
                pp->locallab.LHcurve.at(pp->locallab.selspot) = LHshape->getCurve();
                pp->locallab.HHcurve.at(pp->locallab.selspot) = HHshape->getCurve();
                pp->locallab.invers.at(pp->locallab.selspot) = invers->get_active();
                // Exposure
                pp->locallab.expexpose.at(pp->locallab.selspot) = (int)expexpose->getEnabled();
                pp->locallab.expcomp.at(pp->locallab.selspot) = expcomp->getIntValue();
                pp->locallab.hlcompr.at(pp->locallab.selspot) = hlcompr->getIntValue();
                pp->locallab.hlcomprthresh.at(pp->locallab.selspot) = hlcomprthresh->getIntValue();
                pp->locallab.black.at(pp->locallab.selspot) = black->getIntValue();
                pp->locallab.shcompr.at(pp->locallab.selspot) = shcompr->getIntValue();
                pp->locallab.warm.at(pp->locallab.selspot) = warm->getIntValue();
                pp->locallab.sensiex.at(pp->locallab.selspot) = sensiex->getIntValue();
                pp->locallab.excurve.at(pp->locallab.selspot) = shapeexpos->getCurve();
                // Vibrance
                pp->locallab.expvibrance.at(pp->locallab.selspot) = (int)expvibrance->getEnabled();
                pp->locallab.saturated.at(pp->locallab.selspot) = saturated->getIntValue();
                pp->locallab.pastels.at(pp->locallab.selspot) = pastels->getIntValue();
                pp->locallab.psthreshold.at(pp->locallab.selspot) = psThreshold->getValue<int>();
                pp->locallab.protectskins.at(pp->locallab.selspot) = (int)protectSkins->get_active();
                pp->locallab.avoidcolorshift.at(pp->locallab.selspot) = (int)avoidColorShift->get_active();
                pp->locallab.pastsattog.at(pp->locallab.selspot) = (int)pastSatTog->get_active();
                pp->locallab.sensiv.at(pp->locallab.selspot) = sensiv->getIntValue();
                pp->locallab.skintonescurve.at(pp->locallab.selspot) = skinTonesCurve->getCurve();
                // Blur & Noise
                pp->locallab.expblur.at(pp->locallab.selspot) = (int)expblur->getEnabled();
                pp->locallab.radius.at(pp->locallab.selspot) = radius->getIntValue();
                pp->locallab.strength.at(pp->locallab.selspot) = strength->getIntValue();
                pp->locallab.sensibn.at(pp->locallab.selspot) = sensibn->getIntValue();

                if (blurMethod->get_active_row_number() == 0) {
                    pp->locallab.blurMethod.at(pp->locallab.selspot) = "norm";
                } else if (blurMethod->get_active_row_number() == 1) {
                    pp->locallab.blurMethod.at(pp->locallab.selspot) = "inv";
                } else if (blurMethod->get_active_row_number() == 2) {
                    pp->locallab.blurMethod.at(pp->locallab.selspot) = "sym";
                }

                pp->locallab.activlum.at(pp->locallab.selspot) = (int)activlum->get_active();
                // Tone Mapping
                pp->locallab.exptonemap.at(pp->locallab.selspot) = (int)exptonemap->getEnabled();
                pp->locallab.stren.at(pp->locallab.selspot) = stren->getIntValue();
                pp->locallab.gamma.at(pp->locallab.selspot) = gamma->getIntValue();
                pp->locallab.estop.at(pp->locallab.selspot) = estop->getIntValue();
                pp->locallab.scaltm.at(pp->locallab.selspot) = scaltm->getIntValue();
                pp->locallab.rewei.at(pp->locallab.selspot) = rewei->getIntValue();
                pp->locallab.sensitm.at(pp->locallab.selspot) = sensitm->getIntValue();
                // Retinex
                pp->locallab.expreti.at(pp->locallab.selspot) = (int)expreti->getEnabled();

                if (retinexMethod->get_active_row_number() == 0) {
                    pp->locallab.retinexMethod.at(pp->locallab.selspot) = "low";
                } else if (retinexMethod->get_active_row_number() == 1) {
                    pp->locallab.retinexMethod.at(pp->locallab.selspot) = "uni";
                } else if (retinexMethod->get_active_row_number() == 2) {
                    pp->locallab.retinexMethod.at(pp->locallab.selspot) = "high";
                }

                pp->locallab.str.at(pp->locallab.selspot) = str->getIntValue();
                pp->locallab.chrrt.at(pp->locallab.selspot) = chrrt->getIntValue();
                pp->locallab.neigh.at(pp->locallab.selspot) = neigh->getIntValue();
                pp->locallab.vart.at(pp->locallab.selspot) = vart->getIntValue();
                pp->locallab.sensih.at(pp->locallab.selspot) = sensih->getIntValue();
                pp->locallab.localTgaincurve.at(pp->locallab.selspot) = cTgainshape->getCurve();
                pp->locallab.inversret.at(pp->locallab.selspot) = (int)inversret->get_active();
                // Sharpening
                pp->locallab.expsharp.at(pp->locallab.selspot) = (int)expsharp->getEnabled();
                pp->locallab.sharradius.at(pp->locallab.selspot) = sharradius->getIntValue();
                pp->locallab.sharamount.at(pp->locallab.selspot) = sharamount->getIntValue();
                pp->locallab.shardamping.at(pp->locallab.selspot) = shardamping->getIntValue();
                pp->locallab.shariter.at(pp->locallab.selspot) = shariter->getIntValue();
                pp->locallab.sensisha.at(pp->locallab.selspot) = sensisha->getIntValue();
                pp->locallab.inverssha.at(pp->locallab.selspot) = (int)inverssha->get_active();
                // Contrast by detail levels
                pp->locallab.expcbdl.at(pp->locallab.selspot) = (int)expcbdl->getEnabled();

                for (int i = 0; i < 5; i++) {
                    pp->locallab.mult[i].at(pp->locallab.selspot) = multiplier[i]->getIntValue();
                }

                pp->locallab.chromacbdl.at(pp->locallab.selspot) = chromacbdl->getIntValue();
                pp->locallab.threshold.at(pp->locallab.selspot) = threshold->getValue();
                pp->locallab.sensicb.at(pp->locallab.selspot) = sensicb->getIntValue();
                // Denoise
                pp->locallab.expdenoi.at(pp->locallab.selspot) = (int)expdenoi->getEnabled();
                pp->locallab.noiselumf.at(pp->locallab.selspot) = noiselumf->getIntValue();
                pp->locallab.noiselumc.at(pp->locallab.selspot) = noiselumc->getIntValue();
                pp->locallab.noiselumdetail.at(pp->locallab.selspot) = noiselumdetail->getIntValue();
                pp->locallab.noiselequal.at(pp->locallab.selspot) = noiselequal->getIntValue();
                pp->locallab.noisechrof.at(pp->locallab.selspot) = noisechrof->getIntValue();
                pp->locallab.noisechroc.at(pp->locallab.selspot) = noisechroc->getIntValue();
                pp->locallab.adjblur.at(pp->locallab.selspot) = adjblur->getIntValue();
                pp->locallab.bilateral.at(pp->locallab.selspot) = bilateral->getIntValue();
                pp->locallab.sensiden.at(pp->locallab.selspot) = sensiden->getIntValue();
                pp->locallab.avoid.at(pp->locallab.selspot) = (int)avoid->get_active();
            }

            // Update Locallab tools GUI
            disableListener();
            updateLocallabGUI(pp, pp->locallab.selspot);
            enableListener();
    }

    if (pedited) {
        pedited->locallab.enabled = !get_inconsistent();
        // Control spot settings
        ControlSpotPanel::SpotEdited* const se = expsettings->getEditedStates();
        pedited->locallab.nbspot = se->addbutton || se->deletebutton;
        pedited->locallab.selspot = se->treeview;
        pedited->locallab.id = se->addbutton || se->deletebutton;
        pedited->locallab.name = se->name;
        pedited->locallab.isvisible = se->isvisible;
        pedited->locallab.shape = se->shape;
        pedited->locallab.spotMethod = se->spotMethod;
        pedited->locallab.shapeMethod = se->shapeMethod;
        pedited->locallab.locX = se->locX;
        pedited->locallab.locXL = se->locXL;
        pedited->locallab.locY = se->locY;
        pedited->locallab.locYT = se->locYT;
        pedited->locallab.centerX = se->centerX;
        pedited->locallab.centerY = se->centerY;
        pedited->locallab.circrad = se->circrad;
        pedited->locallab.qualityMethod = se->qualityMethod;
        pedited->locallab.transit = se->transit;
        pedited->locallab.thresh = se->thresh;
        pedited->locallab.iter = se->iter;
        // Color & Light
        pedited->locallab.expcolor = !expcolor->get_inconsistent();
        pedited->locallab.curvactiv = !curvactiv->get_inconsistent();
        pedited->locallab.lightness = lightness->getEditedState();
        pedited->locallab.contrast = contrast->getEditedState();
        pedited->locallab.chroma = chroma->getEditedState();
        pedited->locallab.sensi = sensi->getEditedState();
        pedited->locallab.qualitycurveMethod = qualitycurveMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.llcurve = !llshape->isUnChanged();
        pedited->locallab.cccurve = !ccshape->isUnChanged();
        pedited->locallab.LHcurve = !LHshape->isUnChanged();
        pedited->locallab.HHcurve = !HHshape->isUnChanged();
        pedited->locallab.invers = !invers->get_inconsistent();
        // Exposure
        pedited->locallab.expexpose = !expexpose->get_inconsistent();
        pedited->locallab.expcomp = expcomp->getEditedState();
        pedited->locallab.hlcompr = hlcompr->getEditedState();
        pedited->locallab.hlcomprthresh = hlcomprthresh->getEditedState();
        pedited->locallab.black = black->getEditedState();
        pedited->locallab.shcompr = shcompr->getEditedState();
        pedited->locallab.warm = warm->getEditedState();
        pedited->locallab.sensiex = sensiex->getEditedState();
        pedited->locallab.excurve = !shapeexpos->isUnChanged();
        // Vibrance
        pedited->locallab.expvibrance = !expvibrance->get_inconsistent();
        pedited->locallab.saturated = saturated->getEditedState();
        pedited->locallab.pastels = pastels->getEditedState();
        pedited->locallab.psthreshold = psThreshold->getEditedState();
        pedited->locallab.protectskins = !protectSkins->get_inconsistent();
        pedited->locallab.avoidcolorshift = !avoidColorShift->get_inconsistent();
        pedited->locallab.pastsattog = !pastSatTog->get_inconsistent();
        pedited->locallab.sensiv = sensiv->getEditedState();
        pedited->locallab.skintonescurve = !skinTonesCurve->isUnChanged();
        // Blur & Noise
        pedited->locallab.expblur = !expblur->get_inconsistent();
        pedited->locallab.radius = radius->getEditedState();
        pedited->locallab.strength = strength->getEditedState();
        pedited->locallab.sensibn = sensibn->getEditedState();
        pedited->locallab.blurMethod = blurMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.activlum = !activlum->get_inconsistent();
        // Tone Mapping
        pedited->locallab.exptonemap = !exptonemap->get_inconsistent();
        pedited->locallab.stren = stren->getEditedState();
        pedited->locallab.gamma = gamma->getEditedState();
        pedited->locallab.estop = estop->getEditedState();
        pedited->locallab.scaltm = scaltm->getEditedState();
        pedited->locallab.rewei = rewei->getEditedState();
        pedited->locallab.sensitm = sensitm->getEditedState();
        // Retinex
        pedited->locallab.expreti = !expreti->get_inconsistent();
        pedited->locallab.retinexMethod = retinexMethod->get_active_text() != M("GENERAL_UNCHANGED");
        pedited->locallab.str = str->getEditedState();
        pedited->locallab.chrrt = chrrt->getEditedState();
        pedited->locallab.neigh = neigh->getEditedState();
        pedited->locallab.vart = vart->getEditedState();
        pedited->locallab.sensih = sensih->getEditedState();
        pedited->locallab.localTgaincurve = !cTgainshape->isUnChanged();
        pedited->locallab.inversret = !inversret->get_inconsistent();
        // Sharpening
        pedited->locallab.expsharp = !expsharp->get_inconsistent();
        pedited->locallab.sharradius = sharradius->getEditedState();
        pedited->locallab.sharamount = sharamount->getEditedState();
        pedited->locallab.shardamping = shardamping->getEditedState();
        pedited->locallab.shariter = shariter->getEditedState();
        pedited->locallab.sensisha = sensisha->getEditedState();
        pedited->locallab.inverssha = !inverssha->get_inconsistent();
        // Contrast by detail levels
        pedited->locallab.expcbdl = !expcbdl->get_inconsistent();

        for (int i = 0; i < 5; i++) {
            pedited->locallab.mult[i] = multiplier[i]->getEditedState();
        }

        pedited->locallab.chromacbdl = chromacbdl->getEditedState();
        pedited->locallab.threshold = threshold->getEditedState();
        pedited->locallab.sensicb = sensicb->getEditedState();
        // Denoise
        pedited->locallab.expdenoi = !expdenoi->get_inconsistent();
        pedited->locallab.noiselumf = noiselumf->getEditedState();
        pedited->locallab.noiselumc = noiselumc->getEditedState();
        pedited->locallab.noiselumdetail = noiselumdetail->getEditedState();
        pedited->locallab.noiselequal = noiselequal->getEditedState();
        pedited->locallab.noisechrof = noisechrof->getEditedState();
        pedited->locallab.noisechroc = noisechroc->getEditedState();
        pedited->locallab.adjblur = adjblur->getEditedState();
        pedited->locallab.bilateral = bilateral->getEditedState();
        pedited->locallab.sensiden = sensiden->getEditedState();
        pedited->locallab.avoid = !avoid->get_inconsistent();
    }
}

void Locallab::protectskins_toggled()
{
    printf("protectskins_toggled\n");

    if (batchMode) {
        /*
        if (protectSkins->get_inconsistent()) {
            protectSkins->set_inconsistent(false);
            pskinsconn.block(true);
            protectSkins->set_active(false);
            pskinsconn.block(false);
        } else if (lastProtectSkins) {
            protectSkins->set_inconsistent(true);
        }

        lastProtectSkins = protectSkins->get_active();
        */
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            if (protectSkins->get_active()) {
                listener->panelChanged(EvlocallabProtectSkins, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvlocallabProtectSkins, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::avoidcolorshift_toggled()
{
    printf("avoidcolorshift_toggled\n");

    if (batchMode) {
        /*
        if (avoidColorShift->get_inconsistent()) {
            avoidColorShift->set_inconsistent(false);
            ashiftconn.block(true);
            avoidColorShift->set_active(false);
            ashiftconn.block(false);
        } else if (lastAvoidColorShift) {
            avoidColorShift->set_inconsistent(true);
        }

        lastAvoidColorShift = avoidColorShift->get_active();
        */
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            if (avoidColorShift->get_active()) {
                listener->panelChanged(EvlocallabAvoidColorShift, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvlocallabAvoidColorShift, M("GENERAL_DISABLED"));
            }
        }
    }
}


void Locallab::spotdupChanged(bool spotchan)
{
    // TODO Locallab printf
    printf("spotdupChanged\n");

    /*
    nextspotdup = spotchan;

    const auto func = [](gpointer data) -> gboolean {
        static_cast<Locallab*>(data)->spotdupComputed_();
        return FALSE;
    };

    idle_register.add(func, this);
    */
}

bool Locallab::spotdupComputed_()
{
    // TODO Locallab printf
    printf("spotdupComputed_\n");

    /*
    disableListener();
    // spotduplicated->set_active(nextspotdup);

    if (nextspotdup) {
        //labspotdup->show();
        //usleep(4000000);
    } else {
        //labspotdup->hide();
    }

    enableListener();
    */

    return false;
}


void Locallab::spotduplicatedChanged()
{
    printf("spotduplicatedChanged\n");

    /*
    if (batchMode) {
        if (spotduplicated->get_inconsistent()) {
            spotduplicated->set_inconsistent(false);
            spotduplicated->set_active(false);;
        } else if (lastspotduplicated) {
            spotduplicated->set_inconsistent(true);
        }

        lastspotduplicated = spotduplicated->get_active();
    }

    if (listener && getEnabled()) {

        if (spotduplicated->get_active()) {
            listener->panelChanged(Evlocallabspotduplicated, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabspotduplicated, M("GENERAL_DISABLED"));
        }
    }
    */
}

void Locallab::pastsattog_toggled()
{
    printf("pastsattog_toggled\n");

    if (batchMode) {
        /*
        if (pastSatTog->get_inconsistent()) {
            pastSatTog->set_inconsistent(false);
            pastsattogconn.block(true);
            pastSatTog->set_active(false);
            pastsattogconn.block(false);
        } else if (lastPastSatTog) {
            pastSatTog->set_inconsistent(true);
        }

        lastPastSatTog = pastSatTog->get_active();
        */
    }

    if (pastSatTog->get_active()) {
        saturated->setValue(pastels->getValue());
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            if (pastSatTog->get_active()) {
                listener->panelChanged(EvlocallabPastSatTog, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(EvlocallabPastSatTog, M("GENERAL_DISABLED"));
            }
        }
    }
}


void Locallab::curveChanged(CurveEditor* ce)
{
    // Color & Light
    if (getEnabled() && expcolor->getEnabled()) {
        if (ce == llshape) {
            if (listener) {
                listener->panelChanged(Evlocallabllshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == ccshape) {
            if (listener) {
                listener->panelChanged(Evlocallabccshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == LHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabLHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }

        if (ce == HHshape) {
            if (listener) {
                listener->panelChanged(EvlocallabHHshape, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    // Exposure
    if (getEnabled() && expexpose->getEnabled()) {
        if (ce == shapeexpos) {
            if (listener) {
                listener->panelChanged(Evlocallabshapeexpos, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    // Vibrance
    if (getEnabled() && expvibrance->getEnabled()) {
        if (ce == skinTonesCurve) {
            if (listener) {
                listener->panelChanged(EvlocallabSkinTonesCurve, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    // Retinex
    if (getEnabled() && expreti->getEnabled()) {
        if (ce == cTgainshape) {
            if (listener) {
                listener->panelChanged(EvlocallabCTgainCurve, M("HISTORY_CUSTOMCURVE"));
            }
        }
    }

    /*
    if (ce == skinTonesCurve) {
        listener->panelChanged(EvlocallabSkinTonesCurve, M("HISTORY_CUSTOMCURVE"));
        int strval = retrab->getValue();
        //update MIP
        retrab->setValue(strval + 1);
        adjusterChanged(retrab, strval + 1);
        usleep(10000);  //to test
        retrab->setValue(strval);
        adjusterChanged(retrab, strval);
    }
    */
}

void Locallab::retinexMethodChanged()
{
    printf("retinexMethodChanged\n");

    if (!batchMode) {
        /*
        retrab->hide();
        LocalcurveEditorgainTrab->hide();
        */
    }

    if (getEnabled() && expreti->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabretinexMethod, retinexMethod->get_active_text());
        }
    }
}

void Locallab::blurMethodChanged()
{
    printf("blurMethodChanged\n");

    if (getEnabled() && expblur->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabblurMethod, blurMethod->get_active_text());
        }
    }
}


void Locallab::dustMethodChanged()
{
    printf("dustMethodChanged\n");

    /*
    if (!batchMode) {

    }


    if (listener) {
        //   listener->panelChanged (EvlocallabdustMethod, dustMethod->get_active_text ());
    }
    */
}

void Locallab::qualityMethodChanged()
{
    /*
    if (!batchMode) {
        /*
        if (qualityMethod->get_active_row_number() == 0) { //STD
            proxi->hide();
            thres->hide();
        } else {//enh
            proxi->show();
            thres->show();
        }
        *//*
}

if (listener) {
listener->panelChanged(EvlocallabqualityMethod, qualityMethod->get_active_text());
}
*/
}

void Locallab::qualitycurveMethodChanged()
{
    printf("qualitycurveMethodChanged\n");

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabqualitycurveMethod, qualitycurveMethod->get_active_text());
        }
    }
}

void Locallab::ExclumethodChanged()
{
    /*
    if (!batchMode) {
        if (Exclumethod->get_active_row_number() == 0) {
        } else if (Exclumethod->get_active_row_number() == 1) {
        }
    }


    if (listener) {
        listener->panelChanged(Evlocallabexclumethod, Exclumethod->get_active_text());
    }
    */
}

void Locallab::shapemethodChanged()
{
    /*
    if (!batchMode) {
        //
    }


    if (listener) {
        listener->panelChanged(Evlocallabshapemethod, shapemethod->get_active_text());
    }
    */
}

void Locallab::SmethodChanged()
{
    /*
    if (!batchMode) {
        if (Smethod->get_active_row_number() == 0) { //IND 0
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();
        } else if (Smethod->get_active_row_number() == 1) {         // 1 SYM
            locX->hide();
            locXL->hide();
            locY->hide();
            locYT->hide();
            centerX->hide();
            centerY->hide();

        } else if (Smethod->get_active_row_number() == 2) {         //2 SYM
            locX->show();
            locXL->show();
            locY->show();
            locYT->show();
            centerX->show();
            centerY->show();

        } else if (Smethod->get_active_row_number() == 3) {         // 3 SYM
            locX->show();
            locXL->hide();
            locY->show();
            locYT->hide();
            centerX->show();
            centerY->show();

        }

        /*      else if(Smethod->get_active_row_number()==2) {              // LOC
                    locX->show();
                    locXL->hide();
                    locY->hide();
                    locYT->hide();
                }   *//*
}

if (listener && getEnabled()) {
if (Smethod->get_active_row_number() == 1  || Smethod->get_active_row_number() == 3) {
    listener->panelChanged(EvlocallabSmet, Smethod->get_active_text());
    locXL->setValue(locX->getValue());
    locYT->setValue(locY->getValue());
}
//   else if(Smethod->get_active_row_number()==2) {
//          listener->panelChanged (EvlocallabSmet, Smethod->get_active_text ());
//           locXL->setValue (locX->getValue());
//           locYT->setValue (locX->getValue());
//          locY->setValue (locX->getValue());
//     }
else

{
    listener->panelChanged(EvlocallabSmet, Smethod->get_active_text());

}
}
*/
}
void Locallab::inversChanged()
{
    printf("inversChanged\n");

    /*
    if (batchMode) {
        if (invers->get_inconsistent()) {
            invers->set_inconsistent(false);
            inversConn.block(true);
            invers->set_active(false);
            inversConn.block(false);
        } else if (lastinvers) {
            invers->set_inconsistent(true);
        }

        lastinvers = invers->get_active();
    }
    */

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            if (invers->get_active()) {
                listener->panelChanged(Evlocallabinvers, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinvers, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::cutpastChanged()
{
    /*
    if (batchMode) {
        if (cutpast->get_inconsistent()) {
            cutpast->set_inconsistent(false);
            cutpastConn.block(true);
            cutpast->set_active(false);
            cutpastConn.block(false);
        } else if (lastcutpast) {
            cutpast->set_inconsistent(true);
        }

        lastcutpast = cutpast->get_active();
    }


    if (listener) {
        if (getEnabled()) {
            //     listener->panelChanged (Evlocallabcutpast, M ("GENERAL_ENABLED"));
        } else {
            //     listener->panelChanged (Evlocallabcutpast, M ("GENERAL_DISABLED"));
        }
    }
    */
}

void Locallab::lastdustChanged()
{
    /*
    if (batchMode) {
        if (lastdust->get_inconsistent()) {
            lastdust->set_inconsistent(false);
            lastdustConn.block(true);
            lastdust->set_active(false);
            lastdustConn.block(false);
        } else if (lastlastdust) {
            lastdust->set_inconsistent(true);
        }

        lastlastdust = lastdust->get_active();
    }


    if (listener) {
        if (getEnabled()) {
            //    listener->panelChanged (Evlocallablastdust, M ("GENERAL_ENABLED"));
        } else {
            //    listener->panelChanged (Evlocallablastdust, M ("GENERAL_DISABLED"));
        }
    }
    */
}


void Locallab::curvactivChanged()
{
    printf("curvactivChanged\n");

    // TODO Batch mode
    /*
    if (batchMode) {
        if (curvactiv->get_inconsistent()) {
            curvactiv->set_inconsistent(false);
            curvactivConn.block(true);
            curvactiv->set_active(false);
            curvactivConn.block(false);
        } else if (lastcurvactiv) {
            curvactiv->set_inconsistent(true);
        }

        lastcurvactiv = curvactiv->get_active();
    }
    */

    if (getEnabled() && expcolor->getEnabled()) {
        if (listener) {
            if (curvactiv->get_active()) {
                listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabcurvactiv, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::activlumChanged()
{
    printf("activlumChanged\n");

    if (batchMode) {
        /*
        if (activlum->get_inconsistent()) {
            activlum->set_inconsistent(false);
            activlumConn.block(true);
            activlum->set_active(false);
            activlumConn.block(false);
        } else if (lastactivlum) {
            activlum->set_inconsistent(true);
        }

        lastactivlum = activlum->get_active();
        */
    }

    if (getEnabled() && expblur->getEnabled()) {
        if (listener) {
            if (activlum->get_active()) {
                listener->panelChanged(Evlocallabactivlum, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabactivlum, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::inversradChanged()
{
    /*
    if (batchMode) {
        if (inversrad->get_inconsistent()) {
            inversrad->set_inconsistent(false);
            inversradConn.block(true);
            inversrad->set_active(false);
            inversradConn.block(false);
        } else if (lastinversrad) {
            inversrad->set_inconsistent(true);
        }

        lastinversrad = inversrad->get_active();
    }

    if (inversrad->get_active()) {
        sensibn->hide();
    } else {
        sensibn->show();
    }


    if (listener) {
        if (getEnabled()) {
            listener->panelChanged(Evlocallabinversrad, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(Evlocallabinversrad, M("GENERAL_DISABLED"));
        }
    }
    */
}

void Locallab::inversshaChanged()
{
    printf("inversshaChanged\n");

    if (batchMode) {
        /*
        if (inverssha->get_inconsistent()) {
            inverssha->set_inconsistent(false);
            inversshaConn.block(true);
            inverssha->set_active(false);
            inversshaConn.block(false);
        } else if (lastinverssha) {
            inverssha->set_inconsistent(true);
        }

        lastinverssha = inverssha->get_active();
        */
    }

    if (getEnabled() && expsharp->getEnabled()) {
        if (listener) {
            if (inverssha->get_active()) {
                listener->panelChanged(Evlocallabinverssha, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinverssha, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::inversretChanged()
{
    printf("inversretChanged\n");

    if (batchMode) {
        /*
        if (inversret->get_inconsistent()) {
            inversret->set_inconsistent(false);
            inversretConn.block(true);
            inversret->set_active(false);
            inversretConn.block(false);
        } else if (lastinversret) {
            inversret->set_inconsistent(true);
        }

        lastinversret = inversret->get_active();
        */
    }

    if (getEnabled() && expsharp->getEnabled()) {
        if (listener) {
            if (inversret->get_active()) {
                listener->panelChanged(Evlocallabinversret, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabinversret, M("GENERAL_DISABLED"));
            }
        }
    }
}

// TODO
void Locallab::setDefaults(const ProcParams * defParams, const ParamsEdited * pedited)
{
    /*
    degree->setDefault(defParams->locallab.degree);
    locY->setDefault(defParams->locallab.locY);
    locX->setDefault(defParams->locallab.locX);
    locYT->setDefault(defParams->locallab.locYT);
    locXL->setDefault(defParams->locallab.locXL);
    centerX->setDefault(defParams->locallab.centerX);
    centerY->setDefault(defParams->locallab.centerY);
    circrad->setDefault(defParams->locallab.circrad);
    centerXbuf->setDefault(defParams->locallab.centerXbuf);
    centerYbuf->setDefault(defParams->locallab.centerYbuf);
    adjblur->setDefault(defParams->locallab.adjblur);
    thres->setDefault(defParams->locallab.thres);
    proxi->setDefault(defParams->locallab.proxi);
    lightness->setDefault(defParams->locallab.lightness);
    contrast->setDefault(defParams->locallab.contrast);
    chroma->setDefault(defParams->locallab.chroma);
    warm->setDefault(defParams->locallab.warm);
    expcomp->setDefault(defParams->locallab.expcomp);
    black->setDefault(defParams->locallab.black);
    hlcompr->setDefault(defParams->locallab.hlcompr);
    hlcomprthresh->setDefault(defParams->locallab.hlcomprthresh);
    shcompr->setDefault(defParams->locallab.shcompr);

    noiselumf->setDefault(defParams->locallab.noiselumf);
    noiselumc->setDefault(defParams->locallab.noiselumc);
    noiselumdetail->setDefault(defParams->locallab.noiselumdetail);
    noiselequal->setDefault(defParams->locallab.noiselequal);
    noisechrodetail->setDefault(defParams->locallab.noisechrodetail);
    bilateral->setDefault(defParams->locallab.bilateral);
    sensiden->setDefault(defParams->locallab.sensiden);
    noisechrof->setDefault(defParams->locallab.noisechrof);
    noisechroc->setDefault(defParams->locallab.noisechroc);
    sharradius->setDefault(defParams->locallab.sharradius);
    sharamount->setDefault(defParams->locallab.sharamount);
    shardamping->setDefault(defParams->locallab.shardamping);
    shariter->setDefault(defParams->locallab.shariter);
    sensisha->setDefault(defParams->locallab.sensisha);
    sensi->setDefault(defParams->locallab.sensi);
    sensiex->setDefault(defParams->locallab.sensiex);
    sensih->setDefault(defParams->locallab.sensih);
    retrab->setDefault(defParams->locallab.retrab);
    sensiexclu->setDefault(defParams->locallab.sensiexclu);
    struc->setDefault(defParams->locallab.struc);
    sensicb->setDefault(defParams->locallab.sensicb);
    sensibn->setDefault(defParams->locallab.sensibn);
    sensitm->setDefault(defParams->locallab.sensitm);
    transit->setDefault(defParams->locallab.transit);
    radius->setDefault(defParams->locallab.radius);
    strength->setDefault(defParams->locallab.strength);
    stren->setDefault(defParams->locallab.stren);
    gamma->setDefault(defParams->locallab.gamma);
    estop->setDefault(defParams->locallab.estop);
    gamma->setDefault(defParams->locallab.gamma);
    scaltm->setDefault(defParams->locallab.scaltm);
    rewei->setDefault(defParams->locallab.rewei);
    neigh->setDefault(defParams->locallab.neigh);
    nbspot->setDefault(defParams->locallab.nbspot);
    anbspot->setDefault(defParams->locallab.anbspot);
    hueref->setDefault(defParams->locallab.hueref);
    huerefblur->setDefault(defParams->locallab.huerefblur);
    chromaref->setDefault(defParams->locallab.chromaref);
    lumaref->setDefault(defParams->locallab.lumaref);
    sobelref->setDefault(defParams->locallab.sobelref);

    vart->setDefault(defParams->locallab.vart);
    chrrt->setDefault(defParams->locallab.chrrt);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->setDefault(defParams->locallab.mult[i]);
    }

    threshold->setDefault(defParams->locallab.threshold);
    chromacbdl->setDefault(defParams->locallab.chromacbdl);
    pastels->setDefault(defParams->locallab.pastels);
    saturated->setDefault(defParams->locallab.saturated);
    psThreshold->setDefault<int> (defParams->locallab.psthreshold);
    sensiv->setDefault(defParams->locallab.sensiv);


    if (pedited) {
        degree->setDefaultEditedState(pedited->locallab.degree ? Edited : UnEdited);
        locY->setDefaultEditedState(pedited->locallab.locY ? Edited : UnEdited);
        locX->setDefaultEditedState(pedited->locallab.locX ? Edited : UnEdited);
        locYT->setDefaultEditedState(pedited->locallab.locYT ? Edited : UnEdited);
        locXL->setDefaultEditedState(pedited->locallab.locXL ? Edited : UnEdited);
        centerX->setDefaultEditedState(pedited->locallab.centerX ? Edited : UnEdited);
        centerY->setDefaultEditedState(pedited->locallab.centerY ? Edited : UnEdited);
        circrad->setDefaultEditedState(pedited->locallab.circrad ? Edited : UnEdited);
        centerXbuf->setDefaultEditedState(pedited->locallab.centerXbuf ? Edited : UnEdited);
        centerYbuf->setDefaultEditedState(pedited->locallab.centerYbuf ? Edited : UnEdited);
        adjblur->setDefaultEditedState(pedited->locallab.adjblur ? Edited : UnEdited);
        thres->setDefaultEditedState(pedited->locallab.thres ? Edited : UnEdited);
        proxi->setDefaultEditedState(pedited->locallab.proxi ? Edited : UnEdited);
        lightness->setDefaultEditedState(pedited->locallab.lightness ? Edited : UnEdited);
        contrast->setDefaultEditedState(pedited->locallab.contrast ? Edited : UnEdited);
        chroma->setDefaultEditedState(pedited->locallab.chroma ? Edited : UnEdited);
        warm->setDefaultEditedState(pedited->locallab.warm ? Edited : UnEdited);
        expcomp->setDefaultEditedState(pedited->locallab.expcomp ? Edited : UnEdited);
        black->setDefaultEditedState(pedited->locallab.black ? Edited : UnEdited);
        hlcompr->setDefaultEditedState(pedited->locallab.hlcompr ? Edited : UnEdited);
        hlcomprthresh->setDefaultEditedState(pedited->locallab.hlcomprthresh ? Edited : UnEdited);
        shcompr->setDefaultEditedState(pedited->locallab.shcompr ? Edited : UnEdited);

        noiselumf->setDefaultEditedState(pedited->locallab.noiselumf ? Edited : UnEdited);
        noiselumc->setDefaultEditedState(pedited->locallab.noiselumc ? Edited : UnEdited);
        noiselumdetail->setDefaultEditedState(pedited->locallab.noiselumdetail ? Edited : UnEdited);
        noiselequal->setDefaultEditedState(pedited->locallab.noiselequal ? Edited : UnEdited);
        noisechrodetail->setDefaultEditedState(pedited->locallab.noisechrodetail ? Edited : UnEdited);
        bilateral->setDefaultEditedState(pedited->locallab.bilateral ? Edited : UnEdited);
        sensiden->setDefaultEditedState(pedited->locallab.sensiden ? Edited : UnEdited);
        noisechrof->setDefaultEditedState(pedited->locallab.noisechrof ? Edited : UnEdited);
        noisechroc->setDefaultEditedState(pedited->locallab.noisechroc ? Edited : UnEdited);
        sharradius->setDefaultEditedState(pedited->locallab.sharradius ? Edited : UnEdited);
        sharamount->setDefaultEditedState(pedited->locallab.sharamount ? Edited : UnEdited);
        shardamping->setDefaultEditedState(pedited->locallab.shardamping ? Edited : UnEdited);
        shariter->setDefaultEditedState(pedited->locallab.shariter ? Edited : UnEdited);
        sensisha->setDefaultEditedState(pedited->locallab.sensisha ? Edited : UnEdited);
        sensi->setDefaultEditedState(pedited->locallab.sensi ? Edited : UnEdited);
        sensiex->setDefaultEditedState(pedited->locallab.sensiex ? Edited : UnEdited);
        sensih->setDefaultEditedState(pedited->locallab.sensih ? Edited : UnEdited);
        retrab->setDefaultEditedState(pedited->locallab.retrab ? Edited : UnEdited);
        sensiexclu->setDefaultEditedState(pedited->locallab.sensiexclu ? Edited : UnEdited);
        struc->setDefaultEditedState(pedited->locallab.struc ? Edited : UnEdited);
        sensicb->setDefaultEditedState(pedited->locallab.sensicb ? Edited : UnEdited);
        sensibn->setDefaultEditedState(pedited->locallab.sensibn ? Edited : UnEdited);
        sensitm->setDefaultEditedState(pedited->locallab.sensitm ? Edited : UnEdited);
        radius->setDefaultEditedState(pedited->locallab.radius ? Edited : UnEdited);
        strength->setDefaultEditedState(pedited->locallab.strength ? Edited : UnEdited);
        stren->setDefaultEditedState(pedited->locallab.stren ? Edited : UnEdited);
        gamma->setDefaultEditedState(pedited->locallab.gamma ? Edited : UnEdited);
        estop->setDefaultEditedState(pedited->locallab.estop ? Edited : UnEdited);
        scaltm->setDefaultEditedState(pedited->locallab.scaltm ? Edited : UnEdited);
        rewei->setDefaultEditedState(pedited->locallab.rewei ? Edited : UnEdited);
        transit->setDefaultEditedState(pedited->locallab.transit ? Edited : UnEdited);
        str->setDefaultEditedState(pedited->locallab.str ? Edited : UnEdited);
        neigh->setDefaultEditedState(pedited->locallab.neigh ? Edited : UnEdited);
        nbspot->setDefaultEditedState(pedited->locallab.nbspot ? Edited : UnEdited);
        anbspot->setDefaultEditedState(pedited->locallab.anbspot ? Edited : UnEdited);
        hueref->setDefaultEditedState(pedited->locallab.hueref ? Edited : UnEdited);
        huerefblur->setDefaultEditedState(pedited->locallab.huerefblur ? Edited : UnEdited);
        chromaref->setDefaultEditedState(pedited->locallab.chromaref ? Edited : UnEdited);
        lumaref->setDefaultEditedState(pedited->locallab.lumaref ? Edited : UnEdited);
        sobelref->setDefaultEditedState(pedited->locallab.sobelref ? Edited : UnEdited);
        vart->setDefaultEditedState(pedited->locallab.vart ? Edited : UnEdited);
        chrrt->setDefaultEditedState(pedited->locallab.chrrt ? Edited : UnEdited);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(pedited->locallab.mult[i] ? Edited : UnEdited);
        }

        threshold->setDefaultEditedState(pedited->locallab.threshold ? Edited : UnEdited);
        chromacbdl->setDefaultEditedState(pedited->locallab.chromacbdl ? Edited : UnEdited);

        pastels->setDefaultEditedState(pedited->locallab.pastels ? Edited : UnEdited);
        saturated->setDefaultEditedState(pedited->locallab.saturated ? Edited : UnEdited);
        psThreshold->setDefaultEditedState(pedited->locallab.psthreshold ? Edited : UnEdited);
        sensiv->setDefaultEditedState(pedited->locallab.sensiv ? Edited : UnEdited);

    } else {
        degree->setDefaultEditedState(Irrelevant);
        locY->setDefaultEditedState(Irrelevant);
        locX->setDefaultEditedState(Irrelevant);
        locYT->setDefaultEditedState(Irrelevant);
        locXL->setDefaultEditedState(Irrelevant);
        centerX->setDefaultEditedState(Irrelevant);
        centerY->setDefaultEditedState(Irrelevant);
        circrad->setDefaultEditedState(Irrelevant);
        centerXbuf->setDefaultEditedState(Irrelevant);
        centerYbuf->setDefaultEditedState(Irrelevant);
        adjblur->setDefaultEditedState(Irrelevant);
        thres->setDefaultEditedState(Irrelevant);
        proxi->setDefaultEditedState(Irrelevant);
        lightness->setDefaultEditedState(Irrelevant);
        contrast->setDefaultEditedState(Irrelevant);
        chroma->setDefaultEditedState(Irrelevant);
        warm->setDefaultEditedState(Irrelevant);
        expcomp->setDefaultEditedState(Irrelevant);
        black->setDefaultEditedState(Irrelevant);
        hlcompr->setDefaultEditedState(Irrelevant);
        hlcomprthresh->setDefaultEditedState(Irrelevant);
        shcompr->setDefaultEditedState(Irrelevant);

        noiselumf->setDefaultEditedState(Irrelevant);
        noiselumc->setDefaultEditedState(Irrelevant);
        noiselumdetail->setDefaultEditedState(Irrelevant);
        noiselequal->setDefaultEditedState(Irrelevant);
        noisechrodetail->setDefaultEditedState(Irrelevant);
        bilateral->setDefaultEditedState(Irrelevant);
        sensiden->setDefaultEditedState(Irrelevant);
        noisechrof->setDefaultEditedState(Irrelevant);
        noisechroc->setDefaultEditedState(Irrelevant);
        sharradius->setDefaultEditedState(Irrelevant);
        sharamount->setDefaultEditedState(Irrelevant);
        shardamping->setDefaultEditedState(Irrelevant);
        shariter->setDefaultEditedState(Irrelevant);
        sensisha->setDefaultEditedState(Irrelevant);
        sensi->setDefaultEditedState(Irrelevant);
        sensiex->setDefaultEditedState(Irrelevant);
        sensih->setDefaultEditedState(Irrelevant);
        retrab->setDefaultEditedState(Irrelevant);
        sensicb->setDefaultEditedState(Irrelevant);
        sensiexclu->setDefaultEditedState(Irrelevant);
        sensicb->setDefaultEditedState(Irrelevant);
        struc->setDefaultEditedState(Irrelevant);
        sensitm->setDefaultEditedState(Irrelevant);
        radius->setDefaultEditedState(Irrelevant);
        strength->setDefaultEditedState(Irrelevant);
        stren->setDefaultEditedState(Irrelevant);
        gamma->setDefaultEditedState(Irrelevant);
        estop->setDefaultEditedState(Irrelevant);
        scaltm->setDefaultEditedState(Irrelevant);
        rewei->setDefaultEditedState(Irrelevant);
        transit->setDefaultEditedState(Irrelevant);
        str->setDefaultEditedState(Irrelevant);
        neigh->setDefaultEditedState(Irrelevant);
        nbspot->setDefaultEditedState(Irrelevant);
        anbspot->setDefaultEditedState(Irrelevant);
        hueref->setDefaultEditedState(Irrelevant);
        huerefblur->setDefaultEditedState(Irrelevant);
        chromaref->setDefaultEditedState(Irrelevant);
        lumaref->setDefaultEditedState(Irrelevant);
        sobelref->setDefaultEditedState(Irrelevant);
        vart->setDefaultEditedState(Irrelevant);
        chrrt->setDefaultEditedState(Irrelevant);

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setDefaultEditedState(Irrelevant);
        }

        threshold->setDefaultEditedState(Irrelevant);
        chromacbdl->setDefaultEditedState(Irrelevant);

        pastels->setDefaultEditedState(Irrelevant);
        saturated->setDefaultEditedState(Irrelevant);
        psThreshold->setDefaultEditedState(Irrelevant);
        sensiv->setDefaultEditedState(Irrelevant);

    }
    */
}

void Locallab::adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop)
{
    printf("adjusterChangedTS\n");

    if (getEnabled() && expvibrance->getEnabled()) {
        if (listener) {
            listener->panelChanged(EvlocallabPastSatThreshold, psThreshold->getHistoryString());
        }
    }
}

void Locallab::adjusterChanged(Adjuster * a, double newval)
{
    printf("adjusterChanged\n");

    // Color & Light
    if (getEnabled() && expcolor->getEnabled()) {
        if (a == lightness) {
            if (listener) {
                listener->panelChanged(Evlocallablightness, lightness->getTextValue());
            }
        }

        if (a == contrast) {
            if (listener) {
                listener->panelChanged(Evlocallabcontrast, contrast->getTextValue());
            }
        }

        if (a == chroma) {
            if (listener) {
                listener->panelChanged(Evlocallabchroma, chroma->getTextValue());
            }
        }

        if (a == sensi) {
            if (listener) {
                listener->panelChanged(Evlocallabsensi, sensi->getTextValue());
            }
        }
    }

    // Exposure
    if (getEnabled() && expexpose->getEnabled()) {
        if (a == expcomp) {
            if (listener) {
                listener->panelChanged(Evlocallabexpcomp, expcomp->getTextValue());
            }
        }

        if (a == hlcompr) {
            if (listener) {
                listener->panelChanged(Evlocallabhlcompr, hlcompr->getTextValue());
            }
        }

        if (a == hlcomprthresh) {
            if (listener) {
                listener->panelChanged(Evlocallabhlcomprthresh, hlcomprthresh->getTextValue());
            }
        }

        if (a == black) {
            if (listener) {
                listener->panelChanged(Evlocallabblack, black->getTextValue());
            }
        }

        if (a == shcompr) {
            if (listener) {
                listener->panelChanged(Evlocallabshcompr, shcompr->getTextValue());
            }
        }

        if (a == warm) {
            if (listener) {
                listener->panelChanged(Evlocallabwarm, warm->getTextValue());
            }
        }

        if (a == sensiex) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiex, sensiex->getTextValue());
            }
        }
    }

    // Vibrance
    if (a == pastels && pastSatTog->get_active()) {
        saturated->setValue(newval);
    }

    if (getEnabled() && expvibrance->getEnabled()) {
        if (a == saturated && !pastSatTog->get_active()) {
            if (listener) {
                listener->panelChanged(EvlocallabSaturated, saturated->getTextValue());
            }
        }

        if (a == pastels) {
            if (listener) {
                listener->panelChanged(EvlocallabPastels, pastels->getTextValue());
            }
        }

        if (a == sensiv) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiv, sensiv->getTextValue());
            }
        }
    }

    // Blur & Noise
    if (getEnabled() && expblur->getEnabled()) {
        if (a == radius) {
            if (listener) {
                listener->panelChanged(Evlocallabradius, radius->getTextValue());
            }
        }

        if (a == strength) {
            if (listener) {
                listener->panelChanged(Evlocallabstrength, strength->getTextValue());
            }
        }

        if (a == sensibn) {
            if (listener) {
                listener->panelChanged(Evlocallabsensibn, sensibn->getTextValue());
            }
        }
    }

    // Tone Mapping
    if (getEnabled() && exptonemap->getEnabled()) {
        if (a == stren) {
            if (listener) {
                listener->panelChanged(Evlocallabstren, stren->getTextValue());
            }
        }

        if (a == gamma) {
            if (listener) {
                listener->panelChanged(Evlocallabgamma, gamma->getTextValue());
            }
        }

        if (a == estop) {
            if (listener) {
                listener->panelChanged(Evlocallabestop, estop->getTextValue());
            }
        }

        if (a == scaltm) {
            if (listener) {
                listener->panelChanged(Evlocallabscaltm, scaltm->getTextValue());
            }
        }

        if (a == rewei) {
            if (listener) {
                listener->panelChanged(Evlocallabrewei, rewei->getTextValue());
            }
        }

        if (a == sensitm) {
            if (listener) {
                listener->panelChanged(Evlocallabsensitm, sensitm->getTextValue());
            }
        }
    }

    // Retinex
    if (getEnabled() && expreti->getEnabled()) {
        if (a == str) {
            if (listener) {
                listener->panelChanged(Evlocallabstr, str->getTextValue());
            }
        }

        if (a == chrrt) {
            if (listener) {
                listener->panelChanged(Evlocallabchrrt, chrrt->getTextValue());
            }
        }

        if (a == neigh) {
            if (listener) {
                listener->panelChanged(Evlocallabneigh, neigh->getTextValue());
            }
        }

        if (a == vart) {
            if (listener) {
                listener->panelChanged(Evlocallabvart, vart->getTextValue());
            }
        }

        if (a == sensih) {
            if (listener) {
                listener->panelChanged(Evlocallabsensih, sensih->getTextValue());
            }
        }
    }

    // Sharpening
    if (getEnabled() && expsharp->getEnabled()) {
        if (a == sharradius) {
            if (listener) {
                listener->panelChanged(Evlocallabsharradius, sharradius->getTextValue());
            }
        }

        if (a == sharamount) {
            if (listener) {
                listener->panelChanged(Evlocallabsharamount, sharamount->getTextValue());
            }
        }

        if (a == shardamping) {
            if (listener) {
                listener->panelChanged(Evlocallabshardamping, shardamping->getTextValue());
            }
        }

        if (a == shariter) {
            if (listener) {
                listener->panelChanged(Evlocallabshariter, shariter->getTextValue());
            }
        }

        if (a == sensisha) {
            if (listener) {
                listener->panelChanged(Evlocallabsensis, sensisha->getTextValue());
            }
        }
    }

    // Contrast by detail levels
    if (getEnabled() && expcbdl->getEnabled()) {
        if (a == multiplier[0] || a == multiplier[1] || a == multiplier[2] || a == multiplier[3] || a == multiplier[4]) {
            if (listener) {
                listener->panelChanged(EvlocallabEqualizer,
                                       Glib::ustring::compose("%1, %2, %3, %4, %5",
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[0]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[1]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[2]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[3]->getValue()),
                                               Glib::ustring::format(std::fixed, std::setprecision(0), multiplier[4]->getValue())));
            }
        }

        if (a == chromacbdl) {
            if (listener) {
                listener->panelChanged(Evlocallabchromacbdl, chromacbdl->getTextValue());
            }
        }

        if (a == threshold) {
            if (listener) {
                listener->panelChanged(EvlocallabThresho, threshold->getTextValue());
            }
        }

        if (a == sensicb) {
            if (listener) {
                listener->panelChanged(Evlocallabsensicb, sensicb->getTextValue());
            }
        }
    }

    // Denoise
    if (getEnabled() && expdenoi->getEnabled()) {
        if (a == noiselumf) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumf, noiselumf->getTextValue());
            }
        }

        if (a == noiselumc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumc, noiselumc->getTextValue());
            }
        }

        if (a == noiselumdetail) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselumdetail, noiselumdetail->getTextValue());
            }
        }

        if (a == noiselequal) {
            if (listener) {
                listener->panelChanged(Evlocallabnoiselequal, noiselequal->getTextValue());
            }
        }

        if (a == noisechrof) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechrof, noisechrof->getTextValue());
            }
        }

        if (a == noisechroc) {
            if (listener) {
                listener->panelChanged(Evlocallabnoisechroc, noisechroc->getTextValue());
            }
        }

        if (a == adjblur) {
            if (listener) {
                listener->panelChanged(Evlocallabadjblur, adjblur->getTextValue());
            }
        }

        if (a == bilateral) {
            if (listener) {
                listener->panelChanged(Evlocallabbilateral, bilateral->getTextValue());
            }
        }

        if (a == sensiden) {
            if (listener) {
                listener->panelChanged(Evlocallabsensiden, sensiden->getTextValue());
            }
        }
    }
}

void Locallab::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvlocallabEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void Locallab::avoidChanged()
{
    printf("avoidChanged\n");

    if (batchMode) {
        /*
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent(false);
            avoidConn.block(true);
            avoid->set_active(false);
            avoidConn.block(false);
        } else if (lastavoid) {
            avoid->set_inconsistent(true);
        }

        lastavoid = avoid->get_active();
        */
    }

    if (getEnabled() && expdenoi->getEnabled()) {
        if (listener) {
            if (avoid->get_active()) {
                listener->panelChanged(Evlocallabavoid, M("GENERAL_ENABLED"));
            } else {
                listener->panelChanged(Evlocallabavoid, M("GENERAL_DISABLED"));
            }
        }
    }
}

void Locallab::setAdjusterBehavior(bool degreeadd, bool locYadd, bool locXadd, bool locYTadd, bool locXLadd, bool centeradd, bool lightnessadd, bool contrastadd, bool chromaadd, bool sensiadd, bool transitadd, bool radiusadd,  bool strengthadd)
{
    // can I suppress ...no used
    degree->setAddMode(degreeadd);
    locY->setAddMode(locYadd);
    locX->setAddMode(locXadd);
    locYT->setAddMode(locYTadd);
    locXL->setAddMode(locXLadd);
    centerX->setAddMode(centeradd);
    centerY->setAddMode(centeradd);
    lightness->setAddMode(lightnessadd);
    contrast->setAddMode(contrastadd);
    chroma->setAddMode(chromaadd);
    sensi->setAddMode(sensiadd);
    transit->setAddMode(transitadd);
    radius->setAddMode(radiusadd);
    strength->setAddMode(strengthadd);

}

// TODO
void Locallab::trimValues(rtengine::procparams::ProcParams * pp)
{
    /*
    degree->trimValue(pp->locallab.degree);
    locY->trimValue(pp->locallab.locY);
    locX->trimValue(pp->locallab.locX);
    locYT->trimValue(pp->locallab.locYT);
    locXL->trimValue(pp->locallab.locXL);
    centerX->trimValue(pp->locallab.centerX);
    centerY->trimValue(pp->locallab.centerY);
    circrad->trimValue(pp->locallab.circrad);
    centerXbuf->trimValue(pp->locallab.centerXbuf);
    centerYbuf->trimValue(pp->locallab.centerYbuf);
    adjblur->trimValue(pp->locallab.adjblur);
    thres->trimValue(pp->locallab.thres);
    proxi->trimValue(pp->locallab.proxi);
    lightness->trimValue(pp->locallab.lightness);
    contrast->trimValue(pp->locallab.contrast);
    chroma->trimValue(pp->locallab.chroma);
    warm->trimValue(pp->locallab.warm);
    expcomp->trimValue(pp->locallab.expcomp);
    hlcompr->trimValue(pp->locallab.hlcompr);
    hlcomprthresh->trimValue(pp->locallab.hlcomprthresh);
    black->trimValue(pp->locallab.black);
    shcompr->trimValue(pp->locallab.shcompr);
    noiselumf->trimValue(pp->locallab.noiselumf);
    noiselumc->trimValue(pp->locallab.noiselumc);
    noiselumdetail->trimValue(pp->locallab.noiselumdetail);
    noiselequal->trimValue(pp->locallab.noiselequal);
    noisechrodetail->trimValue(pp->locallab.noisechrodetail);
    bilateral->trimValue(pp->locallab.bilateral);
    sensiden->trimValue(pp->locallab.sensiden);
    noisechrof->trimValue(pp->locallab.noisechrof);
    noisechroc->trimValue(pp->locallab.noisechroc);
    sharradius->trimValue(pp->locallab.sharradius);
    sharamount->trimValue(pp->locallab.sharamount);
    shardamping->trimValue(pp->locallab.shardamping);
    shariter->trimValue(pp->locallab.shariter);
    sensisha->trimValue(pp->locallab.sensisha);
    sensi->trimValue(pp->locallab.sensi);
    sensiex->trimValue(pp->locallab.sensiex);
    sensih->trimValue(pp->locallab.sensih);
    retrab->trimValue(pp->locallab.retrab);
    sensiexclu->trimValue(pp->locallab.sensiexclu);
    struc->trimValue(pp->locallab.struc);
    sensicb->trimValue(pp->locallab.sensicb);
    sensibn->trimValue(pp->locallab.sensibn);
    sensitm->trimValue(pp->locallab.sensitm);
    radius->trimValue(pp->locallab.radius);
    strength->trimValue(pp->locallab.strength);
    stren->trimValue(pp->locallab.stren);
    gamma->trimValue(pp->locallab.gamma);
    estop->trimValue(pp->locallab.estop);
    scaltm->trimValue(pp->locallab.scaltm);
    rewei->trimValue(pp->locallab.rewei);
    transit->trimValue(pp->locallab.transit);
    str->trimValue(pp->locallab.str);
    neigh->trimValue(pp->locallab.neigh);
    nbspot->trimValue(pp->locallab.nbspot);
    anbspot->trimValue(pp->locallab.anbspot);
    hueref->trimValue(pp->locallab.hueref);
    huerefblur->trimValue(pp->locallab.huerefblur);
    chromaref->trimValue(pp->locallab.chromaref);
    lumaref->trimValue(pp->locallab.lumaref);
    sobelref->trimValue(pp->locallab.sobelref);

    vart->trimValue(pp->locallab.vart);
    chrrt->trimValue(pp->locallab.chrrt);

    for (int i = 0; i < 5; i++) {
        multiplier[i]->trimValue(pp->locallab.mult[i]);
    }

    chromacbdl->trimValue(pp->locallab.chromacbdl);

    threshold->trimValue(pp->locallab.threshold);
    pastels->trimValue(pp->locallab.pastels);
    saturated->trimValue(pp->locallab.saturated);
    sensiv->trimValue(pp->locallab.sensiv);
     */
}

void Locallab::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    printf("BatchMode : %d\n", batchMode);

    hueref->hide();
    huerefblur->hide();
    chromaref->hide();
    lumaref->hide();
    sobelref->hide();
    degree->showEditedCB();
    locY->showEditedCB();
    locX->showEditedCB();
    locYT->showEditedCB();
    locXL->showEditedCB();
    centerX->showEditedCB();
    centerY->showEditedCB();
    circrad->showEditedCB();
    centerXbuf->showEditedCB();
    centerYbuf->showEditedCB();
    adjblur->showEditedCB();
    thres->showEditedCB();
    proxi->showEditedCB();
    lightness->showEditedCB();
    contrast->showEditedCB();
    chroma->showEditedCB();
    warm->showEditedCB();
    expcomp->showEditedCB();
    black->showEditedCB();
    hlcompr->showEditedCB();
    hlcomprthresh->showEditedCB();
    shcompr->showEditedCB();
    noiselumf->showEditedCB();
    noiselumc->showEditedCB();
    noiselumdetail->showEditedCB();
    noiselequal->showEditedCB();
    noisechrodetail->showEditedCB();
    bilateral->showEditedCB();
    sensiden->showEditedCB();
    noisechroc->showEditedCB();
    noiselumf->showEditedCB();
    sharradius->showEditedCB();
    sharamount->showEditedCB();
    shardamping->showEditedCB();
    shariter->showEditedCB();
    sensisha->showEditedCB();
    sensi->showEditedCB();
    sensiex->showEditedCB();
    sensih->showEditedCB();
    retrab->showEditedCB();
    sensiexclu->showEditedCB();
    struc->showEditedCB();
    sensicb->showEditedCB();
    sensibn->showEditedCB();
    sensitm->showEditedCB();
    radius->showEditedCB();
    strength->showEditedCB();
    stren->showEditedCB();
    gamma->showEditedCB();
    estop->showEditedCB();
    scaltm->showEditedCB();
    rewei->showEditedCB();
    transit->showEditedCB();
    Smethod->append(M("GENERAL_UNCHANGED"));
    str->showEditedCB();
    neigh->showEditedCB();
    nbspot->showEditedCB();
    anbspot->showEditedCB();
    hueref->showEditedCB();
    huerefblur->showEditedCB();
    chromaref->showEditedCB();
    lumaref->showEditedCB();
    sobelref->showEditedCB();
    vart->showEditedCB();
    LocalcurveEditorgainT->setBatchMode(batchMode);
    LocalcurveEditorgainTrab->setBatchMode(batchMode);
    llCurveEditorG->setBatchMode(batchMode);
//    llCurveEditorG2->setBatchMode (batchMode);
    chrrt->showEditedCB();

    for (int i = 0; i < 5; i++) {
        multiplier[i]->showEditedCB();
    }

    threshold->showEditedCB();
    chromacbdl->showEditedCB();
    pastels->showEditedCB();
    saturated->showEditedCB();
    psThreshold->showEditedCB();
    sensiv->showEditedCB();

    curveEditorGG->setBatchMode(batchMode);

}

std::vector<double> Locallab::getCurvePoints(ThresholdSelector* tAdjuster) const
{
    std::vector<double> points;
    double threshold, transitionWeighting;
    tAdjuster->getPositions<double> (transitionWeighting, threshold); // ( range -100;+100,   range 0;+100 )
    transitionWeighting /= 100.; // range -1., +1.
    threshold /= 100.;      // range  0., +1.

    // Initial point
    points.push_back(0.);
    points.push_back(0.);

    double p2 = 3.0 * threshold / 4.0;             // same one than in ipvibrance.cc
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
        chromaMean =      chromaMean  * transitionWeighting + chromaMean;
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




void Locallab::setEditProvider(EditDataProvider * provider)
{
    cTgainshape->setEditProvider(provider);
    cTgainshaperab->setEditProvider(provider);
    expsettings->setEditProvider(provider);
}

void Locallab::subscribe()
{
    expsettings->subscribe();
}

void Locallab::unsubscribe()
{
    expsettings->unsubscribe();
}

void Locallab::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R, G, B;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01(float (valX), float (valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float (valY), float (valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float (valX) + 0.7f;
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float (valY), float (valX), value, R, G, B);
    } else if (callerId == 4) {  // LH - bottom bar
        Color::hsv2rgb01(float (valX), 0.5f, float (valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float ((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
    }

    caller->ccRed = double (R);
    caller->ccGreen = double (G);
    caller->ccBlue = double (B);
}

void Locallab::setListener(ToolPanelListener* tpl)
{
    this->listener = tpl;
    expsettings->setListener(tpl);
}

void Locallab::enableListener()
{
    printf("enableListener\n");

    FoldableToolPanel::enableListener();
    // Color & Light
    enablecolorConn.block(false);
    curvactivConn.block(false);
    qualitycurveMethodConn.block(false);
    inversConn.block(false);
    // Exposure
    enableexposeConn.block(false);
    // Vibrance
    enablevibranceConn.block(false);
    pskinsconn.block(false);
    ashiftconn.block(false);
    pastsattogconn.block(false);
    // Blur & Noise
    enableblurConn.block(false);
    blurMethodConn.block(false);
    activlumConn.block(false);
    // Tone Mapping
    enabletonemapConn.block(false);
    // Retinex
    enableretiConn.block(false);
    retinexMethodConn.block(false);
    inversretConn.block(false);
    // Sharpening
    enablesharpConn.block(false);
    inversshaConn.block(false);
    // Contrast by detail levels
    enablecbdlConn.block(false);
    // Denoise
    enabledenoiConn.block(false);
    avoidConn.block(false);
}

void Locallab::disableListener()
{
    printf("disableListener\n");

    FoldableToolPanel::disableListener();
    // Color & Light
    enablecolorConn.block(true);
    curvactivConn.block(true);
    qualitycurveMethodConn.block(true);
    inversConn.block(true);
    // Exposure
    enableexposeConn.block(true);
    // Vibrance
    enablevibranceConn.block(true);
    pskinsconn.block(true);
    ashiftconn.block(true);
    pastsattogconn.block(true);
    // Blur & Noise
    enableblurConn.block(true);
    blurMethodConn.block(true);
    activlumConn.block(true);
    // Tone Mapping
    enabletonemapConn.block(true);
    // Retinex
    enableretiConn.block(true);
    retinexMethodConn.block(true);
    inversretConn.block(true);
    // Sharpening
    enablesharpConn.block(true);
    inversshaConn.block(true);
    // Contrast by detail levels
    enablecbdlConn.block(true);
    // Denoise
    enabledenoiConn.block(true);
    avoidConn.block(true);
}

void Locallab::updateLocallabGUI(const rtengine::procparams::ProcParams* pp, int index)
{
    printf("updateLocallabGUI\n");

    // Update GUI values
    if (index < pp->locallab.nbspot) {
        // Color & Light
        expcolor->setEnabled((bool)pp->locallab.expcolor.at(index));
        curvactiv->set_active((bool)pp->locallab.curvactiv.at(index));
        lightness->setValue(pp->locallab.lightness.at(index));
        contrast->setValue(pp->locallab.contrast.at(index));
        chroma->setValue(pp->locallab.chroma.at(index));
        sensi->setValue(pp->locallab.sensi.at(index));

        if (pp->locallab.qualitycurveMethod.at(index) == "none") {
            qualitycurveMethod->set_active(0);
        } else if (pp->locallab.qualitycurveMethod.at(index) == "std") {
            qualitycurveMethod->set_active(1);
        } else {
            qualitycurveMethod->set_active(2);
        }

        llshape->setCurve(pp->locallab.llcurve.at(index));
        ccshape->setCurve(pp->locallab.cccurve.at(index));
        LHshape->setCurve(pp->locallab.LHcurve.at(index));
        HHshape->setCurve(pp->locallab.HHcurve.at(index));
        invers->set_active((bool)pp->locallab.invers.at(index));

        // Exposure
        expexpose->setEnabled((bool)pp->locallab.expexpose.at(index));
        expcomp->setValue(pp->locallab.expcomp.at(index));
        hlcompr->setValue(pp->locallab.hlcompr.at(index));
        hlcomprthresh->setValue(pp->locallab.hlcomprthresh.at(index));
        black->setValue(pp->locallab.black.at(index));
        shcompr->setValue(pp->locallab.shcompr.at(index));
        warm->setValue(pp->locallab.warm.at(index));
        sensiex->setValue(pp->locallab.sensiex.at(index));
        shapeexpos->setCurve(pp->locallab.excurve.at(index));

        // Vibrance
        expvibrance->setEnabled((bool)pp->locallab.expvibrance.at(index));
        saturated->setValue(pp->locallab.saturated.at(index));
        pastels->setValue(pp->locallab.pastels.at(index));
        psThreshold->setValue<int>(pp->locallab.psthreshold.at(index));
        protectSkins->set_active((bool)pp->locallab.protectskins.at(index));
        avoidColorShift->set_active((bool)pp->locallab.avoidcolorshift.at(index));
        pastSatTog->set_active((bool)pp->locallab.pastsattog.at(index));
        sensiv->setValue(pp->locallab.sensiv.at(index));
        skinTonesCurve->setCurve(pp->locallab.skintonescurve.at(index));

        // Blur & Noise
        expblur->setEnabled((bool)pp->locallab.expblur.at(index));
        radius->setValue(pp->locallab.radius.at(index));
        strength->setValue(pp->locallab.strength.at(index));
        sensibn->setValue(pp->locallab.sensibn.at(index));

        if (pp->locallab.blurMethod.at(index) == "norm") {
            blurMethod->set_active(0);
        } else if (pp->locallab.blurMethod.at(index) == "inv") {
            blurMethod->set_active(1);
        } else {
            blurMethod->set_active(2);
        }

        activlum->set_active((bool)pp->locallab.activlum.at(index));

        // Tone Mapping
        exptonemap->setEnabled((bool)pp->locallab.exptonemap.at(index));
        stren->setValue(pp->locallab.stren.at(index));
        gamma->setValue(pp->locallab.gamma.at(index));
        estop->setValue(pp->locallab.estop.at(index));
        scaltm->setValue(pp->locallab.scaltm.at(index));
        rewei->setValue(pp->locallab.rewei.at(index));
        sensitm->setValue(pp->locallab.sensitm.at(index));

        // Retinex
        expreti->setEnabled((bool)pp->locallab.expreti.at(index));

        if (pp->locallab.retinexMethod.at(index) == "low") {
            retinexMethod->set_active(0);
        } else if (pp->locallab.retinexMethod.at(index) == "uni") {
            retinexMethod->set_active(1);
        } else {
            retinexMethod->set_active(2);
        }

        str->setValue(pp->locallab.str.at(index));
        chrrt->setValue(pp->locallab.chrrt.at(index));
        neigh->setValue(pp->locallab.neigh.at(index));
        vart->setValue(pp->locallab.vart.at(index));
        sensih->setValue(pp->locallab.sensih.at(index));
        cTgainshape->setCurve(pp->locallab.localTgaincurve.at(index));
        inversret->set_active((bool)pp->locallab.inversret.at(index));

        // Sharpening
        expsharp->setEnabled((bool)pp->locallab.expsharp.at(index));
        sharradius->setValue(pp->locallab.sharradius.at(index));
        sharamount->setValue(pp->locallab.sharamount.at(index));
        shardamping->setValue(pp->locallab.shardamping.at(index));
        shariter->setValue(pp->locallab.shariter.at(index));
        sensisha->setValue(pp->locallab.sensisha.at(index));
        inverssha->set_active((bool)pp->locallab.inverssha.at(index));

        // Contrast by detail levels
        expcbdl->setEnabled((bool)pp->locallab.expcbdl.at(index));

        for (int i = 0; i < 5; i++) {
            multiplier[i]->setValue(pp->locallab.mult[i].at(index));
        }

        chromacbdl->setValue(pp->locallab.chromacbdl.at(index));
        threshold->setValue(pp->locallab.threshold.at(index));
        sensicb->setValue(pp->locallab.sensicb.at(index));

        // Denoise
        expdenoi->setEnabled((bool)pp->locallab.expdenoi.at(index));
        noiselumf->setValue(pp->locallab.noiselumf.at(index));
        noiselumc->setValue(pp->locallab.noiselumc.at(index));
        noiselumdetail->setValue(pp->locallab.noiselumdetail.at(index));
        noiselequal->setValue(pp->locallab.noiselequal.at(index));
        noisechrof->setValue(pp->locallab.noisechrof.at(index));
        noisechroc->setValue(pp->locallab.noisechroc.at(index));
        adjblur->setValue(pp->locallab.adjblur.at(index));
        bilateral->setValue(pp->locallab.bilateral.at(index));
        sensiden->setValue(pp->locallab.sensiden.at(index));
        avoid->set_active((bool)pp->locallab.avoid.at(index));
    }

    // Update Color & Light GUI according to invers button state
    if (invers->get_active()) {
        sensi->show();
        llCurveEditorG->hide();
        curvactiv->hide();
        qualitycurveMethod->hide();
        labqualcurv->hide();

    } else {
        sensi->show();
        llCurveEditorG->show();
        curvactiv->show();
        qualitycurveMethod->show();
        labqualcurv->show();

    }

    // Update Exposure GUI according to black adjuster state
    shcompr->set_sensitive(!((int)black->getValue() == 0)); // At black = 0, shcompr value has no effect

    // Update Vibrance GUI according to pastsattog button state
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

    // Update Retinex GUI according to inversret button state
    if (inversret->get_active()) {
        sensih->hide();
    } else {
        sensih->show();
    }

    // Update Sharpening GUI according to inverssha button state
    if (inverssha->get_active()) {
        sensisha->hide();
    } else {
        sensisha->show();
    }

    /*
    lastactivlum = pp->locallab.activlum;
    lastProtectSkins = pp->locallab.protectskins;
    lastAvoidColorShift = pp->locallab.avoidcolorshift;
    lastPastSatTog = pp->locallab.pastsattog;
    lastavoid = pp->locallab.avoid;
    lastinvers = pp->locallab.invers;
    lastcutpast = pp->locallab.cutpast;
    lastlastdust = pp->locallab.lastdust;
    lastcurvactiv = pp->locallab.curvactiv;
    lastinversrad = pp->locallab.inversrad;
    lastinversret = pp->locallab.inversret;
    lastinverssha = pp->locallab.inverssha;
    */
}

void Locallab::autoOpenCurve()
{
    printf("autoOpenCurve\n");

    // TODO autoOpenCurve only considers linearity state of selected spot curve
    llshape->openIfNonlinear();
    ccshape->openIfNonlinear();
    LHshape->openIfNonlinear();
    HHshape->openIfNonlinear();
}
