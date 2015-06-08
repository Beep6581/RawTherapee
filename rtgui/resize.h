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
#ifndef _RESIZE_H_
#define _RESIZE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "guiutils.h"
#include "toolpanel.h"
#include "guiutils.h"

class Resize : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public rtengine::SizeListener {

  protected:
    Adjuster*          scale;
    Gtk::VBox*         sizeBox;
    MyComboBoxText*    appliesTo;
    MyComboBoxText*    method;
    MyComboBoxText*    spec;
    MySpinButton*      w;
    MySpinButton*      h;
    int                maxw, maxh;
    int                cropw, croph;
    sigc::connection   sconn, aconn, wconn, hconn;
    bool               wDirty, hDirty;

  public:

    Resize ();
    ~Resize ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
    void setBatchMode   (bool batchMode);

    void adjusterChanged  (Adjuster* a, double newval);
    void entryWChanged    ();
    void entryHChanged    ();
    void appliesToChanged ();
    void methodChanged    ();
    void specChanged      ();
    void update           (bool isCropped, int cw, int ch, int ow=0, int oh=0);
    void setGUIFromCrop   (bool isCropped, int cw, int ch);
    void sizeChanged      (int w, int h, int ow, int oh);
    void setDimensions    ();
    void enabledChanged   ();

  private:
    void fitBoxScale ();
    int getComputedWidth ();
   	int getComputedHeight ();
   	void notifyBBox ();
    void updateGUI ();
};

#endif
