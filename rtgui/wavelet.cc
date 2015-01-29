/*
 *  This file is part of RawTherapee.
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
 *
 *  2014 Jacques Desmis <jdesmis@gmail.com>
 */

#include "wavelet.h"
#include <iomanip>
#include <cmath>
#include "edit.h"
#include "guiutils.h"

using namespace rtengine;
using namespace rtengine::procparams;
extern Options options;

Wavelet::Wavelet () :  FoldableToolPanel(this) {
	std::vector<GradientMilestone> milestones;
	CurveListener::setMulti(true);
	nextnlevel=7.;
	float r, g, b;
	//from -PI to +PI (radians) convert to hsv and draw bottombar
	Color::hsv2rgb01(0.4199, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.0   , r, g, b) ); // hsv: 0.4199   rad: -3.14
	Color::hsv2rgb01(0.5000, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.054 , r, g, b) ); // hsv: 0.5   rad: -2.8
	Color::hsv2rgb01(0.6000, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.1336, r, g, b) ); // hsv: 0.60   rad: -2.3
	Color::hsv2rgb01(0.7500, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.3567, r, g, b) ); // hsv: 0.75   rad: -0.9
	Color::hsv2rgb01(0.8560, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.4363, r, g, b) ); // hsv: 0.856  rad: -0.4
	Color::hsv2rgb01(0.9200, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.4841, r, g, b) ); // hsv: 0.92   rad: -0.1
	Color::hsv2rgb01(0.9300, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.5000, r, g, b) ); // hsv: 0.93   rad:  0
	Color::hsv2rgb01(0.9600, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.5366, r, g, b) ); // hsv: 0.96   rad:  0.25
	Color::hsv2rgb01(1.0000, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.5955, r, g, b) ); // hsv: 1.     rad:  0.6
	Color::hsv2rgb01(0.0675, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.6911, r, g, b) ); // hsv: 0.0675 rad:  1.2
	Color::hsv2rgb01(0.0900, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.7229, r, g, b) ); // hsv: 0.09   rad:  1.4
	Color::hsv2rgb01(0.1700, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.7707, r, g, b) ); // hsv: 0.17   rad:  1.7
	Color::hsv2rgb01(0.2650, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.8503, r, g, b) ); // hsv: 0.265  rad:  2.1
	Color::hsv2rgb01(0.3240, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(0.8981, r, g, b) ); // hsv: 0.324  rad:  2.5
	Color::hsv2rgb01(0.4197, 0.5, 0.5, r, g, b);   milestones.push_back( GradientMilestone(1.    , r, g, b) ); // hsv: 0.419  rad:  3.14
	
	std::vector<GradientMilestone> milestones2;
	milestones2.push_back( GradientMilestone(0.0, 0.0, 0.0, 0.0) );
	milestones2.push_back( GradientMilestone(1.0, 1.0, 1.0, 1.0) );
   
    enabled = Gtk::manage (new Gtk::CheckButton (M("GENERAL_ENABLED")));
    enabled->set_active (true);
    pack_start(*enabled);
    enaConn = enabled->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::enabledToggled) );
	std::vector<double> defaultCurve;

	Gtk::Frame* residualFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_RESID")) );
	residualFrame->set_border_width(0);
	residualFrame->set_label_align(0.025, 0.5);

	Gtk::VBox * resBox = Gtk::manage (new Gtk::VBox());
	resBox->set_border_width(4);	
	resBox->set_spacing(2);

	rescon  = Gtk::manage (new Adjuster (M("TP_WAVELET_RESCON"), -100, 100, 1, 0));
    resBox->pack_start(*rescon, Gtk::PACK_SHRINK);
	rescon->setAdjusterListener (this);
	
	thr  = Gtk::manage (new Adjuster (M("TP_WAVELET_THR"), 0, 100, 1, 30));
    resBox->pack_start(*thr);
	thr->setAdjusterListener (this);

	resconH  = Gtk::manage (new Adjuster (M("TP_WAVELET_RESCONH"), -100, 100, 1, 0));
    resBox->pack_start(*resconH, Gtk::PACK_SHRINK);
	resconH->setAdjusterListener (this);
	
	thrH  = Gtk::manage (new Adjuster (M("TP_WAVELET_THRH"), 0, 100, 1, 70));
    resBox->pack_start(*thrH,Gtk::PACK_SHRINK);
	thrH->setAdjusterListener (this);
	
	reschro  = Gtk::manage (new Adjuster (M("TP_WAVELET_RESCHRO"), -100, 100, 1, 0));
    resBox->pack_start(*reschro);
	reschro->setAdjusterListener (this);


    hueskin2 = Gtk::manage (new ThresholdAdjuster (M("TP_WAVELET_HUESKY"), -314., 314., -260., -250, -130., -140., 0, false));  
    hueskin2->set_tooltip_markup (M("TP_WAVELET_HUESKY_TOOLTIP"));

    hueskin2->setBgGradient(milestones);
    resBox->pack_start(*hueskin2);
    hueskin2->setAdjusterListener (this); 
	
	sky  = Gtk::manage (new Adjuster (M("TP_WAVELET_SKY"), -100., 100.0, 1., 0.));
	sky->set_tooltip_text (M("TP_WAVELET_SKY_TOOLTIP"));
	sky->setAdjusterListener (this);
	resBox->pack_start(*sky);	

	
	Gtk::Frame* levelFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_LEVF")) );
	levelFrame->set_border_width(0);
	levelFrame->set_label_align(0.025, 0.5);

	Gtk::Frame* chromaFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_LEVCH")) );
	chromaFrame->set_border_width(0);
	chromaFrame->set_label_align(0.025, 0.5);

	toningFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_TON")) );
	toningFrame->set_border_width(0);
	toningFrame->set_label_align(0.025, 0.5);
	
	Gtk::Frame* controlFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_CONTR")) );
	controlFrame->set_border_width(0);
	controlFrame->set_label_align(0.025, 0.5);
	
    Gtk::HSeparator *separator1 = Gtk::manage (new  Gtk::HSeparator());
    pack_start(*separator1, Gtk::PACK_SHRINK, 2);
    
    Gtk::VBox * buBox = Gtk::manage (new Gtk::VBox());
	buBox->set_border_width(4);	
	buBox->set_spacing(2);
    
	Gtk::VBox * levBox = Gtk::manage (new Gtk::VBox());
	levBox->set_border_width(4);	
	levBox->set_spacing(2);

	Gtk::VBox * conBox = Gtk::manage (new Gtk::VBox());
	conBox->set_border_width(4);	
	conBox->set_spacing(2);
	
	Gtk::VBox * chBox = Gtk::manage (new Gtk::VBox());
	chBox->set_border_width(4);	
	chBox->set_spacing(2);

	Gtk::VBox * tonBox = Gtk::manage (new Gtk::VBox());
	tonBox->set_border_width(4);	
	tonBox->set_spacing(2);
	
    Gtk::HBox * buttonBox = Gtk::manage (new Gtk::HBox());
	wavLabels = Gtk::manage(new Gtk::Label("---", Gtk::ALIGN_CENTER));
	
    levBox->pack_start(*buttonBox, Gtk::PACK_SHRINK, 2);

        Gtk::Button * contrastMinusButton = Gtk::manage (new Gtk::Button(M("TP_WAVELET_CONTRAST_MINUS")));
        buttonBox->pack_start(*contrastMinusButton, Gtk::PACK_SHRINK, 2);
        contrastMinusPressedConn = contrastMinusButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::contrastMinusPressed));

        Gtk::Button * neutralButton = Gtk::manage (new Gtk::Button(M("TP_WAVELET_NEUTRAL")));
        buttonBox->pack_start(*neutralButton, Gtk::PACK_SHRINK, 2);
        neutralPressedConn = neutralButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::neutralPressed));
        
        Gtk::Button * contrastPlusButton = Gtk::manage (new Gtk::Button(M("TP_WAVELET_CONTRAST_PLUS")));
        buttonBox->pack_start(*contrastPlusButton, Gtk::PACK_SHRINK, 2);
        contrastPlusPressedConn = contrastPlusButton->signal_pressed().connect( sigc::mem_fun(*this, &Wavelet::contrastPlusPressed));

    buttonBox->show_all_children();
	
    Gtk::HSeparator *separator2 = Gtk::manage (new  Gtk::HSeparator());
    levBox->pack_start(*separator2, Gtk::PACK_SHRINK, 2);
	
	unif  = Gtk::manage (new Adjuster (M("TP_WAVELET_UNIF"), 0, 100, 1, 0));
	unif->set_tooltip_text (M("TP_WAVELET_UNIF_TOOLTIP"));
    //levBox->pack_start(*unif);  //keep the possibility to reinstall
	unif->setAdjusterListener (this);
	
	
    Gtk::HSeparator *separatorU = Gtk::manage (new  Gtk::HSeparator());
    levBox->pack_start(*separatorU, Gtk::PACK_SHRINK, 2);

    for(int i = 0; i < 9; i++)
    {
        Glib::ustring ss;
        switch( i ){
        case 0:
            ss =Glib::ustring::compose( "%1 (%2)",i, M("TP_WAVELET_FINEST"));break;
        case 8:
        	ss =Glib::ustring::compose( "%1 (%2)",i, M("TP_WAVELET_LARGEST"));break;
        default:
        	ss =Glib::ustring::compose( "%1",i);
        }
        
        correction[i] = Gtk::manage ( new Adjuster (ss, -100, 150, 1, 0) );
        correction[i]->setAdjusterListener(this);
        levBox->pack_start(*correction[i]);
    }
	sup  = Gtk::manage (new Adjuster (M("TP_WAVELET_SUPE"), -100, 150, 1, 0));
    levBox->pack_start(*sup);
	sup->setAdjusterListener (this);
	wavLabels->show();
	levBox->pack_start (*wavLabels);
	
	
    Gtk::HSeparator *separatorC = Gtk::manage (new  Gtk::HSeparator());
    levBox->pack_start(*separatorC, Gtk::PACK_SHRINK, 2);
	
	
	HSmethod = Gtk::manage (new MyComboBoxText ());
	HSmethod->append_text (M("TP_WAVELET_HS1"));
	HSmethod->append_text (M("TP_WAVELET_HS2"));
	HSmethodconn = HSmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::HSmethodChanged) );
    levBox->pack_start(*HSmethod);
	
    avoid = Gtk::manage (new Gtk::CheckButton (M("TP_WAVELET_AVOID")));
    avoid->set_active (true);
    avoidConn = avoid->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::avoidToggled) );
	
	
	hllev = Gtk::manage (new ThresholdAdjuster (M("TP_WAVELET_HIGHLIGHT"), 0., 100., 50., 75., 100., 98., 0, false));
	hllev->setAdjusterListener (this);
	hllev->setBgGradient(milestones2);
    levBox->pack_start(*hllev);
	
	threshold = Gtk::manage (new Adjuster (M("TP_WAVELET_THRESHOLD"), 1, 9, 1, 5));
    levBox->pack_start(*threshold);
    threshold->setAdjusterListener (this); 
	threshold->set_tooltip_text (M("TP_WAVELET_THRESHOLD_TOOLTIP"));

	bllev = Gtk::manage (new ThresholdAdjuster (M("TP_WAVELET_LOWLIGHT"), 0., 100., 0., 2., 50., 25., 0, false));
	bllev->setAdjusterListener (this);
	bllev->setBgGradient(milestones2);
    levBox->pack_start(*bllev);
	
	threshold2 = Gtk::manage (new Adjuster (M("TP_WAVELET_THRESHOLD2"), 1, 9, 1, 4));
    levBox->pack_start(*threshold2);
    threshold2->setAdjusterListener (this); 
	threshold2->set_tooltip_text (M("TP_WAVELET_THRESHOLD2_TOOLTIP"));


	CHmethod = Gtk::manage (new MyComboBoxText ());
	CHmethod->append_text (M("TP_WAVELET_CH1"));
	CHmethod->append_text (M("TP_WAVELET_CH2"));
	CHmethod->append_text (M("TP_WAVELET_CH3"));
	CHmethodconn = CHmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::CHmethodChanged) );
    chBox->pack_start(*CHmethod);
	
	chroma  = Gtk::manage (new Adjuster (M("TP_WAVELET_CHRO"), 1, 9, 1, 5));
	chroma->set_tooltip_text (M("TP_WAVELET_CHRO_TOOLTIP"));
    chBox->pack_start(*chroma);
	chroma->setAdjusterListener (this);

	satlev = Gtk::manage (new ThresholdAdjuster (M("TP_WAVELET_SAT"), 0., 130., 30., 45., 130., 100., 0, false));
	satlev->setAdjusterListener (this);
	satlev->setBgGradient(milestones2);

	pastlev = Gtk::manage (new ThresholdAdjuster (M("TP_WAVELET_PASTEL"), 0., 70., 0., 2., 30., 20., 0, false));
	pastlev->setAdjusterListener (this);
	pastlev->setBgGradient(milestones2);
    chBox->pack_start(*pastlev);
    chBox->pack_start(*satlev);
	
	chro  = Gtk::manage (new Adjuster (M("TP_WAVELET_CHR"), 0., 100., 1., 0.));
	chro->set_tooltip_text (M("TP_WAVELET_CHR_TOOLTIP"));
    chBox->pack_start(*chro);
	chro->setAdjusterListener (this);
	
	CLVcurveEditorG = new CurveEditorGroup (options.lastWaveletCurvesDir, M("TP_WAVELET_CLVCURVE"));
	CLVcurveEditorG->setCurveListener (this);
	rtengine::WaveletParams::getDefaultCLVCurve(defaultCurve);
	ccshape = static_cast<FlatCurveEditor*>(CLVcurveEditorG->addCurve(CT_Flat, "", NULL, false));
	ccshape->setIdentityValue(0.);
	ccshape->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
	ccshape->setTooltip(M("TP_WAVELET_CURVEEDITOR_CLV_TOOLTIP"));
//	ccshape->setBottomBarColorProvider(this, 2);

 	CLVcurveEditorG->curveListComplete();
    chBox->pack_start(*CLVcurveEditorG, Gtk::PACK_SHRINK, 4);
	
	
	//----------- Color Opacity curve RG ------------------------------


	opaCurveEditorG = new CurveEditorGroup (options.lastWaveletCurvesDir, M("TP_WAVELET_COLORT"));
	opaCurveEditorG->setCurveListener (this);
	std::vector<double> defaultCCurve;

	rtengine::WaveletParams::getDefaultOpacityCurveRG(defaultCCurve);
	opacityShapeRG = static_cast<FlatCurveEditor*>(opaCurveEditorG->addCurve(CT_Flat, "", NULL, false));
	opacityShapeRG->setIdentityValue(0.);
	opacityShapeRG->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCCurve);

	// This will add the reset button at the end of the curveType buttons
	opaCurveEditorG->curveListComplete();
	opaCurveEditorG->show();
	

	tonBox->pack_start( *opaCurveEditorG, Gtk::PACK_SHRINK, 2);
	
	//----------- Opacity curve BY------------------------------

	opacityCurveEditorG = new CurveEditorGroup (options.lastWaveletCurvesDir, M("TP_WAVELET_OPACITY"));
	opacityCurveEditorG->setCurveListener (this);
	std::vector<double> defaultCCBurve;

	rtengine::WaveletParams::getDefaultOpacityCurveBY(defaultCCBurve);
	opacityShapeBY = static_cast<FlatCurveEditor*>(opacityCurveEditorG->addCurve(CT_Flat, "", NULL, false));
	opacityShapeBY->setIdentityValue(0.);
	opacityShapeBY->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCCBurve);
//	opacityShape->setBottomBarBgGradient(milestones);
//	milestones.push_back( GradientMilestone(1., 1., 1., 0.) );
//	milestones.push_back( GradientMilestone(0., 0., 0., 1.) );
//	opacityShapeBY->setLeftBarBgGradient(milestones);

	// This will add the reset button at the end of the curveType buttons
	opacityCurveEditorG->curveListComplete();
	opacityCurveEditorG->show();

	tonBox->pack_start( *opacityCurveEditorG, Gtk::PACK_SHRINK, 2);

	
	

//-------------------------------------------------	
	
    median = Gtk::manage (new Gtk::CheckButton (M("TP_WAVELET_MEDI")));
    median->set_active (true);
    medianConn = median->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::medianToggled) );
    conBox->pack_start(*median);

    skinprotect	= Gtk::manage ( new Adjuster (M("TP_WAVELET_SKIN"), -100, 100, 1, 0.) );
    skinprotect->setAdjusterListener(this);
    conBox->pack_start(*skinprotect);
    skinprotect->set_tooltip_markup (M("TP_WAVELET_SKIN_TOOLTIP"));

    hueskin = Gtk::manage (new ThresholdAdjuster (M("TP_WAVELET_HUESKIN"), -314., 314., -5., 25., 170., 120., 0, false));  
    hueskin->set_tooltip_markup (M("TP_WAVELET_HUESKIN_TOOLTIP"));

    hueskin->setBgGradient(milestones);
    conBox->pack_start(*hueskin);
    hueskin->setAdjusterListener (this);
	
    Gtk::VBox * utiBox = Gtk::manage (new Gtk::VBox());
	utiBox->set_border_width(4);	
	utiBox->set_spacing(2);
	
    Gtk::VBox * diBox = Gtk::manage (new Gtk::VBox());
	diBox->set_border_width(4);	
	diBox->set_spacing(2);

	utilFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_UTIL")) );
	utilFrame->set_border_width(0);
	utilFrame->set_label_align(0.025, 0.5);

    display = Gtk::manage (new Gtk::CheckButton (M("TP_WAVELET_DISPLAY")));
    display->set_active (true);
    displayConn = display->signal_toggled().connect( sigc::mem_fun(*this, &Wavelet::displayToggled) );

	thres  = Gtk::manage (new Adjuster (M("TP_WAVELET_THRES"), 4, 9, 1, 7));
	thres->set_tooltip_text (M("TP_WAVELET_THRES_TOOLTIP"));
	Tilesmethod = Gtk::manage (new MyComboBoxText ());
	Tilesmethod->append_text (M("TP_WAVELET_TILESFULL"));
	Tilesmethod->append_text (M("TP_WAVELET_TILESBIG"));
	Tilesmethod->append_text (M("TP_WAVELET_TILESLIT"));
	Tilesmethodconn = Tilesmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::TilesmethodChanged) );
	Tilesmethod->set_tooltip_text (M("TP_WAVELET_TILES_TOOLTIP"));
	
	
	tiles  = Gtk::manage (new Adjuster (M("TP_WAVELET_TILES"), 12, 22, 1, 14));
	//tiles->set_tooltip_text (M("TP_WAVELET_TILES_TOOLTIP"));
	thres->setAdjusterListener (this);
	tiles->setAdjusterListener (this);
	utiBox->pack_start (*thres);
//	utiBox->pack_start (*tiles);
	utiBox->pack_start (*Tilesmethod);

	dispFrame = Gtk::manage (new Gtk::Frame (M("TP_WAVELET_DISP")) );
	dispFrame->set_border_width(0);
	dispFrame->set_label_align(0.025, 0.5);
	
	
	CLmethod = Gtk::manage (new MyComboBoxText ());
	CLmethod->append_text (M("TP_WAVELET_ONE"));
	CLmethod->append_text (M("TP_WAVELET_INF"));
	CLmethod->append_text (M("TP_WAVELET_SUP"));
	CLmethod->append_text (M("TP_WAVELET_ALL"));
	CLmethodconn = CLmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::CLmethodChanged) );
	diBox->pack_start (*CLmethod);

	Lmethod = Gtk::manage (new MyComboBoxText ());
	Lmethod->append_text (M("TP_WAVELET_0"));
	Lmethod->append_text (M("TP_WAVELET_1"));
	Lmethod->append_text (M("TP_WAVELET_2"));
	Lmethod->append_text (M("TP_WAVELET_3"));
	Lmethod->append_text (M("TP_WAVELET_4"));
	Lmethod->append_text (M("TP_WAVELET_5"));
	Lmethod->append_text (M("TP_WAVELET_6"));
	Lmethod->append_text (M("TP_WAVELET_7"));
	Lmethod->append_text (M("TP_WAVELET_8"));
	Lmethod->set_active(0);
	Lmethodconn = Lmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::LmethodChanged) );
	diBox->pack_start (*Lmethod);
	
	
	Dirmethod = Gtk::manage (new MyComboBoxText ());
	Dirmethod->append_text (M("TP_WAVELET_DONE"));
	Dirmethod->append_text (M("TP_WAVELET_DTWO"));
	Dirmethod->append_text (M("TP_WAVELET_DTHR"));
	Dirmethod->append_text (M("TP_WAVELET_DALL"));
	Dirmethodconn = Dirmethod->signal_changed().connect ( sigc::mem_fun(*this, &Wavelet::DirmethodChanged) );
	diBox->pack_start (*Dirmethod);
	
	pack_start(*display);
	utilFrame->add(*utiBox);
	pack_start (*utilFrame);
	dispFrame->add(*diBox);
	pack_start (*dispFrame);
	

	levelFrame->add(*levBox);
	pack_start (*levelFrame, Gtk::PACK_EXPAND_WIDGET, 4);
	chromaFrame->add(*chBox);
	pack_start (*chromaFrame, Gtk::PACK_EXPAND_WIDGET, 4);
	toningFrame->add(*tonBox);
	pack_start (*toningFrame, Gtk::PACK_EXPAND_WIDGET, 4);
	
	conBox->pack_start(*avoid);
	controlFrame->add(*conBox);
	pack_start (*controlFrame, Gtk::PACK_EXPAND_WIDGET, 4);
	
	residualFrame->add(*resBox);
	pack_start(*residualFrame, Gtk::PACK_EXPAND_WIDGET, 4);	
    show_all_children ();
}

Wavelet::~Wavelet () {
	delete CLVcurveEditorG;
	delete opaCurveEditorG;
	delete opacityCurveEditorG;
}
int wavChangedUI (void* data) {
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    (static_cast<Wavelet*>(data))->wavComputed_ ();
   return 0;
}

void Wavelet::wavChanged (double nlevel) 
{
	nextnlevel=nlevel;
    g_idle_add (wavChangedUI, this);
}
bool Wavelet::wavComputed_ () {

    disableListener ();
    enableListener ();
	updatewavLabel ();
    return false;
}
void Wavelet::updatewavLabel () {
	if (!batchMode) {
		float lv;
		lv=nextnlevel;
		{
		wavLabels->set_text(
				Glib::ustring::compose(M("TP_WAVELET_LEVLABEL"),
				Glib::ustring::format(std::fixed, std::setprecision(0), lv))
				
				);
			}
	}
}

void Wavelet::read (const ProcParams* pp, const ParamsEdited* pedited) {

    disableListener ();
    Lmethodconn.block(true);
    CLmethodconn.block(true);
    Tilesmethodconn.block(true);
    Dirmethodconn.block(true);
    CHmethodconn.block(true);
    HSmethodconn.block(true);
	HSmethod->set_active (1);
    if (pp->wavelet.HSmethod=="without")
        HSmethod->set_active (0);
    else if (pp->wavelet.HSmethod=="with")
       HSmethod->set_active (1);
	HSmethodChanged();
	
	CHmethod->set_active (1);
    if (pp->wavelet.CHmethod=="without")
        CHmethod->set_active (0);
    else if (pp->wavelet.CHmethod=="with")
       CHmethod->set_active (1);
    else if (pp->wavelet.CHmethod=="link")
       CHmethod->set_active (2);
	CHmethodChanged();
	
    CLmethod->set_active (3);
    if (pp->wavelet.CLmethod=="one")
        CLmethod->set_active (0);
    else if (pp->wavelet.CLmethod=="inf")
       CLmethod->set_active (1);
    else if (pp->wavelet.CLmethod=="sup")
       CLmethod->set_active (2);
    else if (pp->wavelet.CLmethod=="all")
       CLmethod->set_active (3);
	CLmethodChanged();

    Tilesmethod->set_active (2);
    if (pp->wavelet.Tilesmethod=="full")
        Tilesmethod->set_active (0);
    else if (pp->wavelet.Tilesmethod=="big")
       Tilesmethod->set_active (1);
    else if (pp->wavelet.Tilesmethod=="lit")
       Tilesmethod->set_active (2);
	TilesmethodChanged();
	
    Dirmethod->set_active (3);
    if (pp->wavelet.Dirmethod=="one")
        Dirmethod->set_active (0);
    else if (pp->wavelet.Dirmethod=="two")
       Dirmethod->set_active (1);
    else if (pp->wavelet.Dirmethod=="thr")
       Dirmethod->set_active (2);
    else if (pp->wavelet.Dirmethod=="all")
       Dirmethod->set_active (3);
	DirmethodChanged();
	
    Lmethod->set_active (4);
    if (pp->wavelet.Lmethod=="0_")
        Lmethod->set_active (0);
    else if (pp->wavelet.Lmethod=="1_")
       Lmethod->set_active (1);
    else if (pp->wavelet.Lmethod=="2_")
       Lmethod->set_active (2);
    else if (pp->wavelet.Lmethod=="3_")
       Lmethod->set_active (3);
    else if (pp->wavelet.Lmethod=="4_")
       Lmethod->set_active (4);
    else if (pp->wavelet.Lmethod=="5_")
       Lmethod->set_active (5);
    else if (pp->wavelet.Lmethod=="6_")
       Lmethod->set_active (6);
    else if (pp->wavelet.Lmethod=="7_")
       Lmethod->set_active (7);
    else if (pp->wavelet.Lmethod=="8_")
       Lmethod->set_active (8);
	   
	LmethodChanged();
	
    if (pedited) {
        if (!pedited->wavelet.Lmethod)
            Lmethod->set_active (8);
        if (!pedited->wavelet.CLmethod)
            CLmethod->set_active (3);
        if (!pedited->wavelet.Tilesmethod)
            Tilesmethod->set_active (2);
        if (!pedited->wavelet.Dirmethod)
            Dirmethod->set_active (3);
        if (!pedited->wavelet.CHmethod)
            CHmethod->set_active (1);
        if (!pedited->wavelet.HSmethod)
            HSmethod->set_active (1);

        enabled->set_inconsistent (!pedited->wavelet.enabled);
        ccshape->setUnChanged  (!pedited->wavelet.clvcurve);
		opacityShapeRG->setCurve (pp->wavelet.opacityCurveRG);
 		opacityShapeBY->setCurve (pp->wavelet.opacityCurveBY);
		avoid->set_inconsistent (!pedited->wavelet.avoid);
		tiles->setEditedState (pedited->wavelet.tiles ? Edited : UnEdited);
		rescon->setEditedState (pedited->wavelet.rescon ? Edited : UnEdited);
		resconH->setEditedState (pedited->wavelet.resconH ? Edited : UnEdited);
		reschro->setEditedState (pedited->wavelet.reschro ? Edited : UnEdited);
		sup->setEditedState (pedited->wavelet.sup ? Edited : UnEdited);
		sky->setEditedState (pedited->wavelet.sky ? Edited : UnEdited);
		thres->setEditedState (pedited->wavelet.thres ? Edited : UnEdited);
		threshold->setEditedState (pedited->wavelet.threshold ? Edited : UnEdited);
		threshold2->setEditedState (pedited->wavelet.threshold2 ? Edited : UnEdited);
        display->set_inconsistent (!pedited->wavelet.display);
		chroma->setEditedState (pedited->wavelet.chroma ? Edited : UnEdited);
		chro->setEditedState (pedited->wavelet.chro ? Edited : UnEdited);
        median->set_inconsistent (!pedited->wavelet.median);
		unif->setEditedState (pedited->wavelet.unif ? Edited : UnEdited);
		thr->setEditedState (pedited->wavelet.thr ? Edited : UnEdited);
		thrH->setEditedState (pedited->wavelet.thrH ? Edited : UnEdited);
        skinprotect->setEditedState (pedited->wavelet.skinprotect ? Edited : UnEdited);
        hueskin->setEditedState 	(pedited->wavelet.hueskin ? Edited : UnEdited);	
        hueskin2->setEditedState 	(pedited->wavelet.hueskin2 ? Edited : UnEdited);	
        hllev->setEditedState 	(pedited->wavelet.hllev ? Edited : UnEdited);	
        bllev->setEditedState 	(pedited->wavelet.bllev ? Edited : UnEdited);	
        pastlev->setEditedState 	(pedited->wavelet.pastlev ? Edited : UnEdited);	
        satlev->setEditedState 	(pedited->wavelet.satlev ? Edited : UnEdited);	
		
        for(int i = 0; i < 9; i++) {
            correction[i]->setEditedState (pedited->wavelet.c[i] ? Edited : UnEdited);
        }
    }
    ccshape->setCurve   (pp->wavelet.clvcurve);
    opacityShapeRG->setCurve   (pp->wavelet.opacityCurveRG);
    opacityShapeBY->setCurve   (pp->wavelet.opacityCurveBY);

    enaConn.block (true);
    enabled->set_active (pp->wavelet.enabled);
    enaConn.block (false);
    avoidConn.block (true);
    avoid->set_active (pp->wavelet.avoid);
    avoidConn.block (false);
    displayConn.block (true);
    display->set_active (pp->wavelet.display);
    displayConn.block (false);
    lastdisplay = pp->wavelet.display;
    medianConn.block (true);
    median->set_active (pp->wavelet.median);
    medianConn.block (false);
    lastmedian = pp->wavelet.median;
    lastEnabled = pp->wavelet.enabled;
    lastavoid = pp->wavelet.avoid;
    tiles->setValue (pp->wavelet.tiles);
    rescon->setValue (pp->wavelet.rescon);
    resconH->setValue (pp->wavelet.resconH);
    reschro->setValue (pp->wavelet.reschro);
    sup->setValue (pp->wavelet.sup);
    sky->setValue (pp->wavelet.sky);
    thres->setValue (pp->wavelet.thres);
    chroma->setValue (pp->wavelet.chroma);
    chro->setValue (pp->wavelet.chro);
    unif->setValue (pp->wavelet.unif);
    thr->setValue (pp->wavelet.thr);
    thrH->setValue (pp->wavelet.thrH);
    skinprotect->setValue(pp->wavelet.skinprotect);
    hueskin->setValue<int>(pp->wavelet.hueskin);
    hueskin2->setValue<int>(pp->wavelet.hueskin2);
    threshold->setValue(pp->wavelet.threshold);
    threshold2->setValue(pp->wavelet.threshold2);
    hllev->setValue<int>(pp->wavelet.hllev);
    bllev->setValue<int>(pp->wavelet.bllev);
    pastlev->setValue<int>(pp->wavelet.pastlev);
    satlev->setValue<int>(pp->wavelet.satlev);
   
    for (int i = 0; i < 9; i++) {
        correction[i]->setValue(pp->wavelet.c[i]);
    }
	int y;
	y=thres->getValue();
	int z;
//	for(z=y;z<9;z++) correction[z]->set_sensitive (false);
//	for(z=0;z<y;z++) correction[z]->set_sensitive (true);
	for(z=y;z<9;z++) correction[z]->hide();
	for(z=0;z<y;z++) correction[z]->show();
	if(z==9) sup->show(); else sup->hide();
    Lmethodconn.block(false);
    CLmethodconn.block(false);
    Tilesmethodconn.block(false);
    CHmethodconn.block(false);
    HSmethodconn.block(false);
    Dirmethodconn.block(false);
	
    enableListener ();
}
void Wavelet::setEditProvider  (EditDataProvider *provider) {
    ccshape->setEditProvider(provider);	
}
void Wavelet::autoOpenCurve () {
    ccshape->openIfNonlinear();	
	//opacityShapeRG->openIfNonlinear();
	//opacityShapeBY->openIfNonlinear();
}

void Wavelet::write (ProcParams* pp, ParamsEdited* pedited) {

    pp->wavelet.enabled 		= enabled->get_active ();
    pp->wavelet.avoid 		= avoid->get_active ();
	pp->wavelet.tiles 		=  tiles->getValue();
	pp->wavelet.rescon 		=  rescon->getValue();
	pp->wavelet.resconH 		=  resconH->getValue();
	pp->wavelet.reschro 		=  reschro->getValue();
	pp->wavelet.sup			=  sup->getValue();
	pp->wavelet.sky 			=  sky->getValue();
	pp->wavelet.thres			 =  thres->getValue();
    pp->wavelet.display 		= display->get_active ();
	pp->wavelet.chroma 		=  chroma->getValue();
	pp->wavelet.chro 			=  chro->getValue();
    pp->wavelet.median 		= median->get_active ();
	pp->wavelet.unif 			=  unif->getValue();
	pp->wavelet.thr 			=  thr->getValue();
	pp->wavelet.thrH 			=  thrH->getValue();
    pp->wavelet.hueskin      	= hueskin->getValue<int> ();
    pp->wavelet.hueskin2      	= hueskin2->getValue<int> ();
	pp->wavelet.skinprotect	=  skinprotect->getValue();
    pp->wavelet.threshold  	= threshold->getValue();
    pp->wavelet.threshold2    = threshold2->getValue();
    pp->wavelet.hllev      	= hllev->getValue<int> ();
    pp->wavelet.bllev       	= bllev->getValue<int> ();
    pp->wavelet.clvcurve  	= ccshape->getCurve ();
    pp->wavelet.opacityCurveRG  	= opacityShapeRG->getCurve ();
    pp->wavelet.opacityCurveBY  = opacityShapeBY->getCurve ();
    pp->wavelet.pastlev       = pastlev->getValue<int> ();
    pp->wavelet.satlev        = satlev->getValue<int> ();

    for (int i = 0; i < 9; i++) {
        pp->wavelet.c[i] = (int) correction[i]->getValue();
    }

    if (pedited) {
        pedited->wavelet.enabled 			=  !enabled->get_inconsistent();
        pedited->wavelet.avoid 			=  !avoid->get_inconsistent();
        pedited->wavelet.display			=  !display->get_inconsistent();
        pedited->wavelet.median 			=  !median->get_inconsistent();
        pedited->wavelet.Lmethod 		 	= Lmethod->get_active_row_number() != 8;
        pedited->wavelet.CLmethod 	 	= CLmethod->get_active_row_number() != 3;
        pedited->wavelet.Tilesmethod 	 	= Tilesmethod->get_active_row_number() != 2;
		pedited->wavelet.CHmethod 	 	= CHmethod->get_active_row_number() != 2;
        pedited->wavelet.HSmethod  		= HSmethod->get_active_row_number() != 1;
        pedited->wavelet.Dirmethod 		= Dirmethod->get_active_row_number() != 3;
        pedited->wavelet.tiles 			= tiles->getEditedState();
        pedited->wavelet.rescon 			= rescon->getEditedState();
        pedited->wavelet.resconH 			= resconH->getEditedState();
        pedited->wavelet.reschro 			= reschro->getEditedState();
        pedited->wavelet.sup 				= sup->getEditedState();
        pedited->wavelet.sky				= sky->getEditedState();
        pedited->wavelet.thres 			= thres->getEditedState();
        pedited->wavelet.threshold 		= threshold->getEditedState();
        pedited->wavelet.threshold2		= threshold2->getEditedState();
        pedited->wavelet.chroma 			= chroma->getEditedState();
        pedited->wavelet.chro 			= chro->getEditedState();
        pedited->wavelet.unif				= unif->getEditedState();
        pedited->wavelet.thr 				= thr->getEditedState();
        pedited->wavelet.thrH 			= thrH->getEditedState();
        pedited->wavelet.hueskin 			= hueskin->getEditedState ();
        pedited->wavelet.hueskin2 		= hueskin2->getEditedState ();
        pedited->wavelet.skinprotect 		= skinprotect->getEditedState();
        pedited->wavelet.hllev 			= hllev->getEditedState ();
        pedited->wavelet.clvcurve    		= !ccshape->isUnChanged ();
        pedited->wavelet.opacityCurveRG    	= !opacityShapeRG->isUnChanged ();
        pedited->wavelet.opacityCurveBY    	= !opacityShapeBY->isUnChanged ();
		pedited->wavelet.bllev 			= bllev->getEditedState ();
		pedited->wavelet.pastlev 			= pastlev->getEditedState ();
		pedited->wavelet.satlev 			= satlev->getEditedState ();

        for(int i = 0; i < 9; i++) {
            pedited->wavelet.c[i] = correction[i]->getEditedState();
        }
    }
		if (CHmethod->get_active_row_number()==0)
			pp->wavelet.CHmethod = "without";
		else if (CHmethod->get_active_row_number()==1)
			pp->wavelet.CHmethod = "with";
		else if (CHmethod->get_active_row_number()==2)
			pp->wavelet.CHmethod = "link";
		
		if (HSmethod->get_active_row_number()==0)
			pp->wavelet.HSmethod = "without";
		else if (HSmethod->get_active_row_number()==1)
			pp->wavelet.HSmethod = "with";
	
		if (CLmethod->get_active_row_number()==0)
			pp->wavelet.CLmethod = "one";
		else if (CLmethod->get_active_row_number()==1)
			pp->wavelet.CLmethod = "inf";
		else if (CLmethod->get_active_row_number()==2)
			pp->wavelet.CLmethod = "sup";
		else if (CLmethod->get_active_row_number()==3)
			pp->wavelet.CLmethod = "all";

		if (Tilesmethod->get_active_row_number()==0)
			pp->wavelet.Tilesmethod = "full";
		else if (Tilesmethod->get_active_row_number()==1)
			pp->wavelet.Tilesmethod = "big";
		else if (Tilesmethod->get_active_row_number()==2)
			pp->wavelet.Tilesmethod = "lit";
		
		if (Dirmethod->get_active_row_number()==0)
			pp->wavelet.Dirmethod = "one";
		else if (Dirmethod->get_active_row_number()==1)
			pp->wavelet.Dirmethod = "two";
		else if (Dirmethod->get_active_row_number()==2)
			pp->wavelet.Dirmethod = "thr";
		else if (Dirmethod->get_active_row_number()==3)
			pp->wavelet.Dirmethod = "all";
	
		if (Lmethod->get_active_row_number()==0)
			pp->wavelet.Lmethod = "0_";
		else if (Lmethod->get_active_row_number()==1)
			pp->wavelet.Lmethod = "1_";
		else if (Lmethod->get_active_row_number()==2)
			pp->wavelet.Lmethod = "2_";
		else if (Lmethod->get_active_row_number()==3)
			pp->wavelet.Lmethod = "3_";
		else if (Lmethod->get_active_row_number()==4)
			pp->wavelet.Lmethod = "4_";
		else if (Lmethod->get_active_row_number()==5)
			pp->wavelet.Lmethod = "5_";
		else if (Lmethod->get_active_row_number()==6)
			pp->wavelet.Lmethod = "6_";
		else if (Lmethod->get_active_row_number()==7)
			pp->wavelet.Lmethod = "7_";
		else if (Lmethod->get_active_row_number()==8)
			pp->wavelet.Lmethod = "8_";
	
	
	
}
void Wavelet::curveChanged (CurveEditor* ce) {

    if (listener && enabled->get_active()) {
	    if (ce == ccshape)
            listener->panelChanged (EvWavCLVCurve, M("HISTORY_CUSTOMCURVE"));	
		else if (ce == opacityShapeRG)
			listener->panelChanged (EvWavColor, M("HISTORY_CUSTOMCURVE"));
		else if (ce == opacityShapeBY)
			listener->panelChanged (EvWavOpac, M("HISTORY_CUSTOMCURVE"));
		}
}

void Wavelet::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited) {

    for (int i = 0; i < 9; i++) {
        correction[i]->setDefault(defParams->wavelet.c[i]);
    }
	tiles->setDefault (defParams->wavelet.tiles);
	rescon->setDefault (defParams->wavelet.rescon);
	resconH->setDefault (defParams->wavelet.resconH);
	reschro->setDefault (defParams->wavelet.reschro);
	sup->setDefault (defParams->wavelet.sup);
	sky->setDefault (defParams->wavelet.sky);
	thres->setDefault (defParams->wavelet.thres);
	threshold->setDefault (defParams->wavelet.threshold);
	threshold2->setDefault (defParams->wavelet.threshold2);
	chroma->setDefault (defParams->wavelet.chroma);
	chro->setDefault (defParams->wavelet.chro);
	unif->setDefault (defParams->wavelet.unif);
	thr->setDefault (defParams->wavelet.thr);
	thrH->setDefault (defParams->wavelet.thrH);
    hueskin->setDefault<int> (defParams->wavelet.hueskin);
    hueskin2->setDefault<int> (defParams->wavelet.hueskin2);
    hllev->setDefault<int> (defParams->wavelet.hllev);
    bllev->setDefault<int> (defParams->wavelet.bllev);
    pastlev->setDefault<int> (defParams->wavelet.pastlev);
    satlev->setDefault<int> (defParams->wavelet.satlev);
    
    if (pedited) {
	tiles->setDefault (defParams->wavelet.tiles);
	rescon->setDefault (defParams->wavelet.rescon);
	resconH->setDefault (defParams->wavelet.resconH);
	reschro->setDefault (defParams->wavelet.reschro);
	sup->setDefault (defParams->wavelet.sup);
    sky->setDefaultEditedState(pedited->wavelet.sky ? Edited : UnEdited);
    thres->setDefaultEditedState(pedited->wavelet.thres ? Edited : UnEdited);
    threshold->setDefaultEditedState(pedited->wavelet.threshold ? Edited : UnEdited);
    threshold2->setDefaultEditedState(pedited->wavelet.threshold2 ? Edited : UnEdited);
    chroma->setDefaultEditedState(pedited->wavelet.chroma ? Edited : UnEdited);
    chro->setDefaultEditedState(pedited->wavelet.chro ? Edited : UnEdited);
    unif->setDefaultEditedState(pedited->wavelet.unif ? Edited : UnEdited);
    thr->setDefaultEditedState(pedited->wavelet.thr ? Edited : UnEdited);
    thrH->setDefaultEditedState(pedited->wavelet.thrH ? Edited : UnEdited);
    skinprotect->setDefaultEditedState(pedited->wavelet.skinprotect ? Edited : UnEdited);
    hueskin->setDefaultEditedState	(pedited->wavelet.hueskin ? Edited : UnEdited);
    hueskin2->setDefaultEditedState	(pedited->wavelet.hueskin2 ? Edited : UnEdited);
    hllev->setDefaultEditedState	(pedited->wavelet.hllev ? Edited : UnEdited);
    bllev->setDefaultEditedState	(pedited->wavelet.bllev ? Edited : UnEdited);
    pastlev->setDefaultEditedState	(pedited->wavelet.pastlev ? Edited : UnEdited);
    satlev->setDefaultEditedState	(pedited->wavelet.satlev ? Edited : UnEdited);
	
        for (int i = 0; i < 9; i++) {
            correction[i]->setDefaultEditedState(pedited->wavelet.c[i] ? Edited : UnEdited);
        }
    }
    else {
            tiles->setDefaultEditedState(Irrelevant);
            rescon->setDefaultEditedState(Irrelevant);
            resconH->setDefaultEditedState(Irrelevant);
            reschro->setDefaultEditedState(Irrelevant);
            sup->setDefaultEditedState(Irrelevant);
            sky->setDefaultEditedState(Irrelevant);
            thres->setDefaultEditedState(Irrelevant);
            threshold->setDefaultEditedState(Irrelevant);
            threshold2->setDefaultEditedState(Irrelevant);
            chroma->setDefaultEditedState(Irrelevant);
            chro->setDefaultEditedState(Irrelevant);
            unif->setDefaultEditedState(Irrelevant);
            thr->setDefaultEditedState(Irrelevant);
            thrH->setDefaultEditedState(Irrelevant);
			skinprotect->setDefaultEditedState(Irrelevant);
			hueskin->setDefaultEditedState (Irrelevant);
			hueskin2->setDefaultEditedState (Irrelevant);
			hllev->setDefaultEditedState (Irrelevant);
			bllev->setDefaultEditedState (Irrelevant);
			pastlev->setDefaultEditedState (Irrelevant);
			satlev->setDefaultEditedState (Irrelevant);
	
        for (int i = 0; i < 9; i++) {
            correction[i]->setDefaultEditedState(Irrelevant);
        }
    }
}


void Wavelet::adjusterChanged2 (ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) {
    if (listener && (multiImage||enabled->get_active()) ) {
		if(a==hueskin) { 
		listener->panelChanged (EvWavHueskin,hueskin->getHistoryString());
		}
		else if(a==hueskin2) { 
		listener->panelChanged (EvWavHueskin2,hueskin2->getHistoryString());
		}
		else if(a==hllev) {
		listener->panelChanged (EvWavlhl,hllev->getHistoryString());
		}	
		else if(a==bllev) {
		listener->panelChanged (EvWavlbl,bllev->getHistoryString());
		}
		else if(a==pastlev) {
		listener->panelChanged (EvWavpast,pastlev->getHistoryString());
		}
		else if(a==satlev) {
		listener->panelChanged (EvWavsat,satlev->getHistoryString());
		}
		
    }
}
void Wavelet::HSmethodChanged() {
	if (!batchMode) {
		if(HSmethod->get_active_row_number()==0) {//without
		hllev->hide();
		bllev->hide();
		threshold->hide();
		threshold2->hide();
		}
		else {//with
		hllev->show();
		bllev->show();
		threshold->show();
		threshold2->show();
		}		
	}
	
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvWavHSmet, HSmethod->get_active_text ());
	}	
}

void Wavelet::CHmethodChanged() {
	if (!batchMode) {
		if(CHmethod->get_active_row_number()==0) {//without
			pastlev->hide();
			satlev->hide();
			chroma->hide();
			chro->hide();
			CLVcurveEditorG->show();	
		}
		else if(CHmethod->get_active_row_number()==1) {//with
			pastlev->show();
			satlev->show();
			chroma->show();
			chro->hide();
			CLVcurveEditorG->show();			
		}
		else {//link
			chro->show();
			pastlev->hide();
			satlev->hide();
			chroma->hide();
			CLVcurveEditorG->hide();	
		}	
	}
	
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvWavCHmet, CHmethod->get_active_text ());
	}	
}

void Wavelet::CLmethodChanged() {
	if (!batchMode) {
	}
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvWavCLmet, CLmethod->get_active_text ());
	}	
}

void Wavelet::TilesmethodChanged() {
	if (!batchMode) {
	}
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvWavTilesmet, Tilesmethod->get_active_text ());
	}	
}

void Wavelet::DirmethodChanged() {
	if (!batchMode) {	
	}
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvWavDirmeto, Dirmethod->get_active_text ());
	}	
}

void Wavelet::LmethodChanged() {
	if (!batchMode) {
	}
	if (listener && (multiImage||enabled->get_active()) ) {
		listener->panelChanged (EvWavLmet, Lmethod->get_active_text ());
	}	
}

void Wavelet::setBatchMode (bool batchMode) {
    Lmethod->append_text (M("GENERAL_UNCHANGED"));
    CLmethod->append_text (M("GENERAL_UNCHANGED"));
    Tilesmethod->append_text (M("GENERAL_UNCHANGED"));
    CHmethod->append_text (M("GENERAL_UNCHANGED"));
    HSmethod->append_text (M("GENERAL_UNCHANGED"));
    Dirmethod->append_text (M("GENERAL_UNCHANGED"));
	CLVcurveEditorG->setBatchMode (batchMode);
 	opaCurveEditorG->setBatchMode (batchMode);
 	opacityCurveEditorG->setBatchMode (batchMode);
	tiles->showEditedCB ();
	rescon->showEditedCB ();
	resconH->showEditedCB ();
	reschro->showEditedCB ();
	sup->showEditedCB ();
    sky->showEditedCB ();
    thres->showEditedCB ();
    threshold->showEditedCB ();
    threshold2->showEditedCB ();
    chroma->showEditedCB ();
    chro->showEditedCB ();
    unif->showEditedCB ();
    thr->showEditedCB ();
    thrH->showEditedCB ();
    skinprotect->showEditedCB();
    hueskin->showEditedCB ();
    hueskin2->showEditedCB ();
    hllev->showEditedCB ();
    bllev->showEditedCB ();
    pastlev->showEditedCB ();
    satlev->showEditedCB ();
    ToolPanel::setBatchMode (batchMode);
    
    for (int i = 0; i < 9; i++) {
        correction[i]->showEditedCB();
    }
}

void Wavelet::adjusterChanged (Adjuster* a, double newval) {
    if (listener && enabled->get_active()) {
        if (a == tiles) {
            listener->panelChanged (EvWavtiles,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), tiles->getValue()))
            );
        }
        else if (a == rescon ) {
            listener->panelChanged (EvWavrescon,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), rescon->getValue()))
            );
        }
        else if (a == resconH ) {
            listener->panelChanged (EvWavresconH,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), resconH->getValue()))
            );
        }
		
        else if (a == reschro ) {
            listener->panelChanged (EvWavreschro,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), reschro->getValue()))
            );
        }
		
        else if (a == sky ) {
            listener->panelChanged (EvWavsky,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(2), sky->getValue()))
            );
        }
        else if (a == sup ) {
            listener->panelChanged (EvWavsup,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), sup->getValue()))
            );
        }
		
        else if (a == chroma ) {
            listener->panelChanged (EvWavchroma,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), chroma->getValue()))
            );
        }
        else if (a == chro ) {
            listener->panelChanged (EvWavchro,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), chro->getValue()))
            );
        }
        else if (a == unif ) {
            listener->panelChanged (EvWavunif,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), unif->getValue()))
            );
        }
        else if (a == thr ) {
            listener->panelChanged (EvWavthr,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), thr->getValue()))
            );
        }
        else if (a == thrH ) {
            listener->panelChanged (EvWavthrH,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), thrH->getValue()))
            );
        }	
        else if (a == threshold ) {
            listener->panelChanged (EvWavThreshold,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), threshold->getValue()))
            );
        }
        else if (a == threshold2 ) {
            listener->panelChanged (EvWavThreshold2,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), threshold2->getValue()))
            );			
        }
		
        else if (a == thres ) {
			int y;
			y=thres->getValue();
			int z;
		//	for(z=y;z<9;z++) correction[z]->set_sensitive (false);;
		//	for(z=0;z<y;z++) correction[z]->set_sensitive (true);;
			for(z=y;z<9;z++) correction[z]->hide();
			for(z=0;z<y;z++) correction[z]->show();
			if(z==9) sup->show(); else sup->hide();
			
            listener->panelChanged (EvWavthres,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(0), thres->getValue()))
            );
        }
        else if (a == skinprotect) {
            listener->panelChanged (EvWavSkin,
                         Glib::ustring::compose("%1",
                         Glib::ustring::format(std::fixed, std::setprecision(2), skinprotect->getValue()))
            );
        }
		
	
        else {
            listener->panelChanged (EvWavelet,
                         Glib::ustring::compose("%1, %2, %3, %4, %5, %6, %7, %8, %9",
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[0]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[1]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[2]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[3]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[4]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[5]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[6]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[7]->getValue()),
                         Glib::ustring::format(std::fixed, std::setprecision(0), correction[8]->getValue()))
            );
		}
    }    
}

void Wavelet::enabledToggled () {

    if (batchMode) {
        if (enabled->get_inconsistent()) {
            enabled->set_inconsistent (false);
            enaConn.block (true);
            enabled->set_active (false);
            enaConn.block (false);
        }
        else if (lastEnabled)
            enabled->set_inconsistent (true);

        lastEnabled = enabled->get_active ();
    }
	else {
		int y=thres->getValue();
		int z;
		for(z=y;z<9;z++) correction[z]->hide();
		for(z=0;z<y;z++) correction[z]->show();
		if(z==9) sup->show(); else sup->hide();
	}
    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvWavEnabled, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvWavEnabled, M("GENERAL_DISABLED"));
    }
}

void Wavelet::medianToggled () {

    if (batchMode) {
        if (median->get_inconsistent()) {
            median->set_inconsistent (false);
            medianConn.block (true);
            median->set_active (false);
            medianConn.block (false);
        }
        else if (lastmedian)
            median->set_inconsistent (true);

        lastmedian = median->get_active ();
    }

    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvWavmedian, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvWavmedian, M("GENERAL_DISABLED"));
    }
}

void Wavelet::avoidToggled () {

    if (batchMode) {
        if (avoid->get_inconsistent()) {
            avoid->set_inconsistent (false);
            avoidConn.block (true);
            avoid->set_active (false);
            avoidConn.block (false);
        }
        else if (lastavoid)
            avoid->set_inconsistent (true);

        lastavoid = avoid->get_active ();
    }
    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvWavavoid, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvWavavoid, M("GENERAL_DISABLED"));
    }
}

void Wavelet::displayToggled () {

    if (batchMode) {
        if (display->get_inconsistent()) {
            display->set_inconsistent (false);
            displayConn.block (true);
            display->set_active (false);
            displayConn.block (false);
        }
        else if (lastdisplay)
            display->set_inconsistent (true);

        lastdisplay = display->get_active ();
    }
	if (display->get_active ()) {
		utilFrame->show();
		dispFrame->show();
	}
	else {
		utilFrame->hide();
		dispFrame->hide();
	
	}
	/*
    if (listener) {
        if (enabled->get_active ())
            listener->panelChanged (EvEqldisplay, M("GENERAL_ENABLED"));
        else
            listener->panelChanged (EvEqldisplay, M("GENERAL_DISABLED"));
    }
	*/
}
void Wavelet::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller) {

	float R, G, B;

	if (elemType==ColorCaller::CCET_VERTICAL_BAR)
		valY = 0.5;

	if (callerId == 1) {         // ch - main curve

		Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
	}
/*	else if (callerId == 2) {    // cc - bottom bar

	//	float value = (1.f - 0.7f) * float(valX) + 0.7f;
		float value = (1.f - 0.7f) * float(valX) + 0.7f;
		// whole hue range
		// Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
		Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
	}
	*/
	else if (callerId == 4) {    // LH - bottom bar
		Color::hsv2rgb01(float(valX), 0.5f, float(valY), R, G, B);
	}
	else if (callerId == 5) {    // HH - bottom bar
		float h = float((valY - 0.5) * 0.3 + valX);
		if (h > 1.0f)
			h -= 1.0f;
		else if (h < 0.0f)
			h += 1.0f;
		Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
	}
	caller->ccRed = double(R);
	caller->ccGreen = double(G);
	caller->ccBlue = double(B);
}
void Wavelet::setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool threshold2add, bool thresadd, bool chroadd,bool chromaadd, bool unifadd, bool skinadd, bool reschroadd, bool resconadd, bool resconHadd, bool thradd, bool thrHadd, bool skyadd ) {

	for (int i=0; i<9; i++)
		correction[i]->setAddMode(multiplieradd);
		threshold->setAddMode(thresholdadd);
		skinprotect->setAddMode(skinadd);
		threshold2->setAddMode(threshold2add);
		thres->setAddMode(thresadd);
		chro->setAddMode(chroadd);
		chroma->setAddMode(chromaadd);
		unif->setAddMode(unifadd);
		rescon->setAddMode(resconadd);		
		resconH->setAddMode(resconHadd);		
		reschro->setAddMode(reschroadd);		
		thr->setAddMode(thradd);		
		thrH->setAddMode(thrHadd);		
		sky->setAddMode(skyadd);		
		
}


void Wavelet::neutralPressed () {

    for (int i = 0; i < 9; i++) {
        correction[i]->setValue(0);
       adjusterChanged(correction[i], 0);
    }
}


void Wavelet::contrastPlusPressed () {

    for (int i = 0; i < 9; i++) {
        int inc = 1 * (9 - i);
        correction[i]->setValue(correction[i]->getValue() + inc);
        adjusterChanged(correction[i], correction[i]->getValue());
    }
}


void Wavelet::contrastMinusPressed () {

    for (int i = 0; i < 9; i++) {
        int inc = -1 * (9 - i);
        correction[i]->setValue(correction[i]->getValue() + inc);
        adjusterChanged(correction[i], correction[i]->getValue());
    }
}

