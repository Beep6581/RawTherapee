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
#ifndef _CROP_H_
#define _CROP_H_

#include <gtkmm.h>
#include <cropguilistener.h>
#include <toolpanel.h>

class CropPanelListener {

  public:
    virtual void cropSelectRequested () {}
};


class Crop : public Gtk::VBox, public CropGUIListener, public ToolPanel, public rtengine::SizeListener {

  protected:
    Gtk::CheckButton* enabled;
    Gtk::CheckButton* fixr;
    Gtk::ComboBoxText* ratio;
    Gtk::ComboBoxText* orientation;
    Gtk::ComboBoxText* guide;
    Gtk::Button* selectCrop;
    CropPanelListener* clistener;
    int opt;
    Gtk::SpinButton* x;
    Gtk::SpinButton* y;
    Gtk::SpinButton* w;
    Gtk::SpinButton* h;
    Gtk::SpinButton* dpi;
    Gtk::Label* sizecm;
    Gtk::Label* sizein;
    Gtk::VBox* dpibox;
    int maxw, maxh;
    int nx, ny, nw, nh;
    double nsx, nsy, nsw, nsh, lastScale;
    int lastRotationDeg;
    sigc::connection xconn, yconn, wconn, hconn, econn, fconn, rconn, oconn, gconn;
    bool wDirty, hDirty, xDirty, yDirty, lastEnabled, lastAspect;

  public:

    Crop ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);
    
    void ratioChanged   ();
    void refreshSize    ();
    void selectPressed  ();
    void setDimensions  (int mw, int mh);
    void enabledChanged ();
    void positionChanged ();
    void widthChanged   ();
    void heightChanged  ();
    bool refreshSpins   (bool notify=false);
    void notifyListener ();
    void sizeChanged    (int w, int h, int ow, int oh);
    
    void readOptions    ();
    void writeOptions   ();

    void cropMoved          (int &x, int &y, int &w, int &h);
    void cropWidth1Resized  (int &x, int &y, int &w, int &h);
    void cropWidth2Resized  (int &x, int &y, int &w, int &h);
    void cropHeight1Resized (int &x, int &y, int &w, int &h);
    void cropHeight2Resized (int &x, int &y, int &w, int &h);
    void cropInit           (int &x, int &y, int &w, int &h);
    void cropResized        (int &x, int &y, int& x2, int& y2);
    void cropManipReady     ();
    double getRatio         ();

    void setCropPanelListener (CropPanelListener* cl) { clistener = cl; }

    void resizeScaleChanged (double rsc);
    void hFlipCrop          ();
    void vFlipCrop          ();
    void rotateCrop         (int deg);
};

#endif
