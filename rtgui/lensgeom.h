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
#ifndef _LENSGEOM_H_
#define _LENSGEOM_H_

#include <gtkmm.h>
#include "toolpanel.h"
#include "lensgeomlistener.h"

class LensGeometry : public ToolParamBlock, public FoldableToolPanel
{

protected:
    Gtk::Button*        autoCrop;
    LensGeomListener*   rlistener;
    Gtk::CheckButton*   fill;
    bool                lastFill;
    sigc::connection    fillConn;
    ToolParamBlock*           packBox;

public:

    LensGeometry ();

    Gtk::Box* getPackBox ()
    {
        return packBox;
    }

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL);
    void setBatchMode   (bool batchMode);

    void fillPressed            ();
    void autoCropPressed        ();
    void setLensGeomListener    (LensGeomListener* l)
    {
        rlistener = l;
    }
};

#endif
