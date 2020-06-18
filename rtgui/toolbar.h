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

#include "toolenum.h"

class RTImage;
class LockablePickerToolListener;

class ToolBarListener
{

public:
    virtual ~ToolBarListener() = default;
    /// Callback when a tool is selected
    virtual void toolSelected(ToolMode tool) = 0;

    /// Callback when the Edit mode is stopped
    virtual void editModeSwitchedOff() = 0;
};

class ToolBar final : public Gtk::HBox
{
private:
    std::unique_ptr<RTImage> handimg;
    std::unique_ptr<RTImage> editinghandimg;
    std::unique_ptr<RTImage> showcolpickersimg;
    std::unique_ptr<RTImage> hidecolpickersimg;
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
    bool blockEdit; // true if edit tool shouldn't be disabled when pressing hand button or h/H key
    sigc::connection  handConn;
    sigc::connection  wbConn;
    sigc::connection  cpConn;
    sigc::connection  cropConn;
    sigc::connection  straConn;

public:
    ToolBar ();

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

    void blockEditDeactivation(bool cond = true)
    {
        blockEdit = cond;
    }
};
