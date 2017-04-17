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
#ifndef _ICMPANEL_
#define _ICMPANEL_

#include <memory>
#include <gtkmm.h>
#include "adjuster.h"
#include "guiutils.h"

#include "toolpanel.h"
#include "popupbutton.h"
#include "../rtengine/imagedata.h"

class ICMPanelListener
{

public:
    virtual ~ICMPanelListener() {}
    virtual void saveInputICCReference (Glib::ustring fname, bool apply_wb) {}
};

class ICMPanel : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{

protected:
    Gtk::Frame*        dcpFrame;
    Adjuster* gampos;
    Adjuster* slpos;
    bool lastgamfree;
    sigc::connection  gamcsconn;
    //bool freegamma;
    bool lastToneCurve;
    sigc::connection tcurveconn;
    bool lastApplyLookTable;
    sigc::connection ltableconn;
    bool lastApplyBaselineExposureOffset;
    sigc::connection beoconn;
    bool lastApplyHueSatMap;
    sigc::connection hsmconn;
    bool lastobpc;
    sigc::connection obpcconn;
    bool lastBlendCMSMatrix;
    bool isBatchMode;
    sigc::connection blendcmsconn;

private:
    Gtk::VBox       *  iVBox;

    Gtk::CheckButton*  obpc;
    Gtk::CheckButton*  freegamma;
    Gtk::RadioButton*  inone;

    Gtk::RadioButton*  iembedded;
    Gtk::RadioButton*  icamera;
    Gtk::RadioButton*  icameraICC;
    Gtk::RadioButton*  ifromfile;
    Gtk::Label*        dcpIllLabel;
    MyComboBoxText*    dcpIll;
    sigc::connection   dcpillconn;
    Gtk::CheckButton*  ckbToneCurve;
    Gtk::CheckButton*  ckbApplyLookTable;
    Gtk::CheckButton*  ckbApplyBaselineExposureOffset;
    Gtk::CheckButton*  ckbApplyHueSatMap;
    Gtk::CheckButton*  ckbBlendCMSMatrix;
    MyComboBoxText*    wnames;
    sigc::connection   wnamesconn;
    MyComboBoxText*    wgamma;
    sigc::connection   wgammaconn;

    MyComboBoxText*    onames;
    sigc::connection   onamesconn;
    PopUpButton*       ointent;
    sigc::connection   ointentconn;
    Gtk::RadioButton*  ofromdir;
    Gtk::RadioButton*  ofromfile;
    Gtk::RadioButton*  iunchanged;
    MyFileChooserButton* ipDialog;
    Gtk::RadioButton::Group opts;
    Gtk::Button*        saveRef;
    sigc::connection   ipc;
    Glib::ustring      oldip;
    ICMPanelListener*  icmplistener;

    double dcpTemperatures[2];
    bool enableLastICCWorkDirChange;
    Glib::ustring lastRefFilename;
    Glib::ustring camName;
    void updateDCP(int dcpIlluminant, Glib::ustring dcp_name);
    void updateRenderingIntent (const Glib::ustring &profile);
public:
    ICMPanel ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode   (bool batchMode);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void adjusterChanged (Adjuster* a, double newval);
    void setAdjusterBehavior (bool gammaadd, bool slopeadd);

    void wpChanged ();
    void opChanged ();
    void oiChanged (int n);
    void oBPCChanged ();
    void ipChanged ();
    void gpChanged ();
    void GamChanged ();
    void ipSelectionChanged ();
    void blendCMSMatrixChanged();
    void dcpIlluminantChanged();
    void toneCurveChanged();
    void applyLookTableChanged();
    void applyBaselineExposureOffsetChanged();
    void applyHueSatMapChanged();

    void setRawMeta (bool raw, const rtengine::ImageData* pMeta);
    void saveReferencePressed ();

    void setICMPanelListener (ICMPanelListener* ipl)
    {
        icmplistener = ipl;
    }
};

#endif
