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
#ifndef _PREPROCESS_H_
#define _PREPROCESS_H_

#include <gtkmm.h>
#include <adjuster.h>
#include <toolpanel.h>
#include <rawimage.h>

class DFProvider {
  public:
    virtual rtengine::RawImage* getDF() {}
    // add other info here
};

class FFProvider {
  public:
    virtual rtengine::RawImage* getFF() {}
    // add other info here
};

class PreProcess : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel{

  protected:

    Gtk::ComboBoxText* darkFrameMethod;
    Gtk::FileChooserButton *darkFrameFile;
    Gtk::HBox *hbdf;
    Gtk::Button *btnReset;
    Gtk::Label *dfLabel;
	Gtk::Label *dfInfo;
    bool dfChanged;
	
	Gtk::FileChooserButton *flatFieldFile;
	Gtk::Label *ffLabel;
	Gtk::Label *ffInfo;
	Gtk::Button *flatFieldFileReset; 
	Gtk::CheckButton* flatFieldAutoSelect;
	Adjuster* flatFieldBlurRadius; 
	Gtk::ComboBoxText* flatFieldBlurType; 
	Gtk::HBox *hbff; 
	bool ffChanged; 

	Adjuster* caRed;
    Adjuster* caBlue;
	Adjuster* PexPos;
	Adjuster* PexPreser;
	
    Adjuster* lineDenoise;
    Adjuster* greenEqThreshold;
    Gtk::CheckButton* caAutocorrect;
    Gtk::CheckButton* hotDeadPixel;
    Gtk::CheckButton* dfAuto;
	bool lastCA,lastHot,lastDFauto, lastFFAutoSelect;
	
	DFProvider *dfp;
	FFProvider *ffp;

	sigc::connection caacsconn,dfautoconn,hdpixelconn,dfFile,flatFieldFileconn,flatFieldAutoSelectconn,flatFieldBlurTypeconn;

  public:

    PreProcess ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);

    void adjusterChanged     (Adjuster* a, double newval);
    void caCorrectionChanged();
    void hotDeadPixelChanged();
    void darkFrameChanged();
    void darkFrameReset();
    void dfAutoChanged();
    
    void flatFieldFileChanged();      
    void flatFieldFile_Reset();       
    void flatFieldAutoSelectChanged();
    void flatFieldBlurRadiusChanged();
    void flatFieldBlurTypeChanged();

    void setDFProvider (DFProvider* p) { dfp = p; };
    void setFFProvider (FFProvider* p) { ffp = p; };
};

#endif
