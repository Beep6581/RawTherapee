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
#ifndef __TOOLBAR_H__
#define __TOOLBAR_H__

#include <gtkmm.h>
#include "toolenum.h"
#include "rtimage.h"
#include "lockablecolorpicker.h"

class ToolBarListener
{

public:
    virtual ~ToolBarListener() {}
    /// Callback when a tool is selected
    virtual void toolSelected (ToolMode tool) {}

    /// Callback when the Edit mode is stopped
    virtual void editModeSwitchedOff () {}
};

class ToolBar : public Gtk::HBox
{
private:
    RTImage* handimg;
    RTImage* editinghandimg;
    RTImage* showcolpickersimg;
    RTImage* hidecolpickersimg;
    bool showColPickers;

    void hand_pressed ();
    void wb_pressed ();
    void colPicker_pressed (GdkEventButton* event);
    void crop_pressed ();
    void stra_pressed ();
    bool showColorPickers(bool showCP);
    void switchColorPickersVisibility();

protected:
    Gtk::ToggleButton* handTool;
    Gtk::ToggleButton* wbTool;
    Gtk::ToggleButton* colPickerTool;
    Gtk::ToggleButton* cropTool;
    Gtk::ToggleButton* straTool;
    ToolBarListener* listener;
    LockablePickerToolListener* pickerListener;
    ToolMode current;
    bool allowNoTool;
    bool editingMode;  // true if the cursor is being used to remotely edit tool's values
    sigc::connection  handConn;
    sigc::connection  wbConn;
    sigc::connection  cpConn;
    sigc::connection  cropConn;
    sigc::connection  straConn;

public:
    ToolBar ();
    ~ToolBar ();

    void     setTool (ToolMode tool);
    ToolMode getTool ()
    {
        return current;
    }

    bool showColorPickers() {
        return showColPickers;
    }

    void setToolBarListener (ToolBarListener* tpl)
    {
        listener = tpl;
    }

    void setLockablePickerToolListener (LockablePickerToolListener* lptl)
    {
        pickerListener = lptl;
    }

    void startEditMode();
    void stopEditMode();

    bool handleShortcutKey (GdkEventKey* event);
    void setBatchMode();
};

#endif
