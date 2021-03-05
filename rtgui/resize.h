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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "guiutils.h"
#include "guiutils.h"
#include "toolpanel.h"

class Resize final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::SizeListener
{
public:
    Resize ();
    ~Resize () override;

    Gtk::Box* getPackBox ()
    {
        return packBox;
    }

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged  (Adjuster* a, double newval) override;
    void entryWChanged    ();
    void entryHChanged    ();
    void appliesToChanged ();
    void methodChanged    ();
    void specChanged      ();
    void update           (bool isCropped, int cw, int ch, int ow = 0, int oh = 0);
    void setGUIFromCrop   (bool isCropped, int cw, int ch);
    void sizeChanged      (int w, int h, int ow, int oh) override;
    void setDimensions    ();
    void enabledChanged   () override;

    void setAdjusterBehavior (bool scaleadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;

private:
    void fitBoxScale ();
    int getComputedWidth ();
    int getComputedHeight ();
    void notifyBBox ();
    void updateGUI ();
    void allowUpscalingChanged();

    rtengine::ProcEvent EvResizeAllowUpscaling;
    Adjuster*          scale;
    Gtk::Box*         sizeBox;
    MyComboBoxText*    appliesTo;
    MyComboBoxText*    method;
    MyComboBoxText*    spec;
    MySpinButton*      w;
    MySpinButton*      h;
    Gtk::CheckButton *allowUpscaling;
    int                maxw, maxh;
    int                cropw, croph;
    sigc::connection   sconn, aconn, wconn, hconn;
    bool               wDirty, hDirty;
    ToolParamBlock*    packBox;
    IdleRegister       idle_register;

    static constexpr int MAX_SCALE = 16; // 16 to match the main preview max scale of 1600%
};
