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
#ifndef _RAWCROP_H_
#define _RAWCROP_H_

#include <gtkmm.h>
#include "toolpanel.h"
#include "guiutils.h"
#include "edit.h"
#include <vector>

class RawCrop final :
    public ToolParamBlock,
    public FoldableToolPanel,
    public rtedit::EditSubscriber,
    public rtengine::SizeListener
{
public:
    RawCrop();
    ~RawCrop();

    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setBatchMode(bool batchMode);

    void enabledChanged();
    void spinChanged(MySpinButton *spinButton);
    void trim(rtengine::procparams::ProcParams* pp, int ow, int oh);

    bool inImageArea(int x, int y);

    void sizeChanged(int w, int h, int ow, int oh);
    void resizeScaleChanged (double rsc);

    void setEditProvider (rtedit::EditDataProvider* provider);

    // EditSubscriber interface
    CursorShape getCursor(const int objectID);
    bool mouseOver(const int modifierKey);
    bool button1Pressed(const int modifierKey);
    bool button1Released();
    bool drag1(const int modifierKey);
    void switchOffEditMode ();

private:

    enum class LitArea {
        none,
        topLeft,
        top,
        topRight,
        left,
        center,
        right,
        bottomLeft,
        bottom,
        bottomRight
    };

    LitArea litArea;

    Gtk::ToggleButton* edit;
    MySpinButton* x;
    MySpinButton* y;
    MySpinButton* w;
    MySpinButton* h;

    int maxw, maxh;
    sigc::connection xconn, yconn, wconn, hconn, editConn;
    bool xSet, ySet, wSet, hSet;

    rtedit::Circle centerCircle;
    rtedit::Polyline topLeft;
    rtedit::Polyline topRight;
    rtedit::Polyline bottomRight;
    rtedit::Polyline bottomLeft;
    rtedit::Polyline side;

    int oldX[2];  // Left and right absolute position of the side
    int oldY[2];  // Top and bottom absolute position of the side

    IdleRegister idle_register;

    void notifyListener(int nx, int ny, int nw, int nh);
    void updateState();
    void updateGeometry  (const int x, const int y, const int width, const int height, const int fullWidth=-1, const int fullHeight=-1);
    bool inArea(const int x1, const int y1, const int x2, const int y2, const rtengine::Coord &mousePos);
    void editToggled();
    void setDimensions(int mw, int mh);
};

#endif
