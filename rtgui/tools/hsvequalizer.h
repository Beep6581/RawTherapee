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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  2010 Ilya Popov <ilia_popov@rambler.ru>
 */
#pragma once

#include <gtkmm.h>

#include "colorprovider.h"
#include "curvelistener.h"
#include "guiutils.h"
#include "toolpanel.h"

class CurveEditor;
class CurveEditorGroup;
class FlatCurveEditor;

class HSVEqualizer final :
    public ToolParamBlock,
    public FoldableToolPanel,
    public CurveListener,
    public ColorProvider
{

protected:

    CurveEditorGroup*  curveEditorG;
    FlatCurveEditor*   hshape;
    FlatCurveEditor*   sshape;
    FlatCurveEditor*   vshape;

public:
    static const Glib::ustring TOOL_NAME;

    HSVEqualizer ();
    ~HSVEqualizer () override;

    void read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void curveChanged    (CurveEditor* ce) override;
    //void setDefaults     (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode    (bool batchMode) override;
    void setEditProvider (EditDataProvider *provider) override;
    void autoOpenCurve   () override;
    void colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller* caller) override;

    void enabledChanged() override;
};
