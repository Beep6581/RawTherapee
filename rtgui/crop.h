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
#include "cropguilistener.h"
#include "toolpanel.h"
#include "guiutils.h"
#include <vector>

class CropPanelListener
{

public:
    virtual void cropSelectRequested () {}
};

class CropRatio
{

public:
    Glib::ustring label;
    double value;
};

class Crop : public ToolParamBlock, public CropGUIListener, public FoldableToolPanel, public rtengine::SizeListener
{

protected:
    Gtk::CheckButton* fixr;
    MyComboBoxText* ratio;
    MyComboBoxText* orientation;
    MyComboBoxText* guide;
    Gtk::Button* selectCrop;
    CropPanelListener* clistener;
    int opt;
    MySpinButton* x;
    MySpinButton* y;
    MySpinButton* w;
    MySpinButton* h;
    MySpinButton* ppi;
    Gtk::Label* sizecm;
    Gtk::Label* sizein;
    Gtk::VBox* ppibox;
    Gtk::VBox* sizebox;
    int maxw, maxh;
    double nx, ny;
    int nw, nh;
    int lastRotationDeg;
    sigc::connection xconn, yconn, wconn, hconn, fconn, rconn, oconn, gconn;
    bool wDirty, hDirty, xDirty, yDirty, lastFixRatio;
    void adjustCropToRatio();
    std::vector<CropRatio>   cropratio;

public:

    Crop ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode   (bool batchMode);

    void ratioChanged   ();
    void ratioFixedChanged ();  // The toggle button
    void refreshSize    ();
    void selectPressed  ();
    void setDimensions   (int mw, int mh);
    void enabledChanged ();
    void positionChanged ();
    void widthChanged   ();
    void heightChanged  ();
    bool refreshSpins   (bool notify = false);
    void notifyListener ();
    void sizeChanged    (int w, int h, int ow, int oh);
    void trim           (rtengine::procparams::ProcParams* pp, int ow, int oh);
    void readOptions    ();
    void writeOptions   ();

    void cropMoved          (int &x, int &y, int &w, int &h);
    void cropWidth1Resized  (int &x, int &y, int &w, int &h);
    void cropWidth2Resized  (int &x, int &y, int &w, int &h);
    void cropHeight1Resized (int &x, int &y, int &w, int &h);
    void cropHeight2Resized (int &x, int &y, int &w, int &h);
    void cropTopLeftResized     (int &x, int &y, int &w, int &h);
    void cropTopRightResized    (int &x, int &y, int &w, int &h);
    void cropBottomLeftResized  (int &x, int &y, int &w, int &h);
    void cropBottomRightResized (int &x, int &y, int &w, int &h);
    void cropInit           (int &x, int &y, int &w, int &h);
    void cropResized        (int &x, int &y, int& x2, int& y2);
    void cropManipReady     ();
    bool inImageArea        (int x, int y);
    double getRatio         ();

    void setCropPanelListener (CropPanelListener* cl)
    {
        clistener = cl;
    }

    void resizeScaleChanged (double rsc);
    void hFlipCrop          ();
    void vFlipCrop          ();
    void rotateCrop         (int deg, bool hflip, bool vflip);
};

#endif
