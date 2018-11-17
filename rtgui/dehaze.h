/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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
#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class Dehaze: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{
private:
    Adjuster *strength;
    Adjuster *depth;
    Gtk::CheckButton *showDepthMap;    

    rtengine::ProcEvent EvDehazeEnabled;
    rtengine::ProcEvent EvDehazeStrength;
    rtengine::ProcEvent EvDehazeDepth;
    rtengine::ProcEvent EvDehazeShowDepthMap;
    
public:

    Dehaze();

    void read(const rtengine::procparams::ProcParams *pp, const ParamsEdited *pedited=nullptr);
    void write(rtengine::procparams::ProcParams *pp, ParamsEdited *pedited=nullptr);
    void setDefaults(const rtengine::procparams::ProcParams *defParams, const ParamsEdited *pedited=nullptr);
    void setBatchMode(bool batchMode);

    void adjusterChanged(Adjuster *a, double newval);
    void enabledChanged();
    void showDepthMapChanged();
    void setAdjusterBehavior(bool strengthAdd);
    void adjusterAutoToggled(Adjuster* a, bool newval) {}
};

