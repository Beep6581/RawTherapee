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
#include "toolbar.h"
#include "multilangmgr.h"
#include "guiutils.h"

ToolBar::ToolBar () : showColPickers(true), listener (NULL)
{

    editingMode = false;

    handimg = Gtk::manage (new RTImage ("openhand.png"));
    handimg->reference();
    editinghandimg = Gtk::manage (new RTImage ("editmodehand.png"));
    editinghandimg->reference();

    handTool = Gtk::manage (new Gtk::ToggleButton ());
    handTool->add (*handimg);
    handimg->show ();
    handTool->set_relief(Gtk::RELIEF_NONE);
    handTool->show ();

    pack_start (*handTool);

    wbTool = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* wbimg = Gtk::manage (new RTImage ("gtk-color-picker.png"));
    wbTool->add (*wbimg);
    wbimg->show ();
    wbTool->set_relief(Gtk::RELIEF_NONE);
    wbTool->show ();

    pack_start (*wbTool);

    showcolpickersimg = Gtk::manage (new RTImage ("colorPickers-show.png"));
    showcolpickersimg->reference();
    hidecolpickersimg = Gtk::manage (new RTImage ("colorPickers-hide.png"));
    hidecolpickersimg->reference();

    colPickerTool = Gtk::manage (new Gtk::ToggleButton ());
    colPickerTool->add (*showcolpickersimg);
    showcolpickersimg->show ();
    colPickerTool->set_relief(Gtk::RELIEF_NONE);
    colPickerTool->show ();

    pack_start (*colPickerTool);

    cropTool = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* cropimg = Gtk::manage (new RTImage ("crop.png"));
    cropTool->add (*cropimg);
    cropimg->show ();
    cropTool->set_relief(Gtk::RELIEF_NONE);
    cropTool->show ();

    pack_start (*cropTool);

    straTool = Gtk::manage (new Gtk::ToggleButton ());
    Gtk::Image* straimg = Gtk::manage (new RTImage ("straighten.png"));
    straTool->add (*straimg);
    straimg->show ();
    straTool->set_relief(Gtk::RELIEF_NONE);
    straTool->show ();

    pack_start (*straTool);


    handTool->set_active (true);
    current = TMHand;

    handConn = handTool->signal_toggled().connect( sigc::mem_fun(*this, &ToolBar::hand_pressed));
    wbConn   = wbTool->signal_toggled().connect( sigc::mem_fun(*this, &ToolBar::wb_pressed));
    cpConn   = colPickerTool->signal_button_press_event().connect_notify( sigc::mem_fun(*this, &ToolBar::colPicker_pressed));
    cropConn = cropTool->signal_toggled().connect( sigc::mem_fun(*this, &ToolBar::crop_pressed));
    straConn = straTool->signal_toggled().connect( sigc::mem_fun(*this, &ToolBar::stra_pressed));

    handTool->set_tooltip_markup (M("TOOLBAR_TOOLTIP_HAND"));
    wbTool->set_tooltip_markup (M("TOOLBAR_TOOLTIP_WB"));
    colPickerTool->set_tooltip_markup (M("TOOLBAR_TOOLTIP_COLORPICKER"));
    cropTool->set_tooltip_markup (M("TOOLBAR_TOOLTIP_CROP"));
    straTool->set_tooltip_markup (M("TOOLBAR_TOOLTIP_STRAIGHTEN"));
}

ToolBar::~ToolBar ()
{
    handimg->unreference();
    editinghandimg->unreference();
    showcolpickersimg->unreference();
    hidecolpickersimg->unreference();

}
//
// Selects the desired tool without notifying the listener
//
void ToolBar::setTool (ToolMode tool)
{

    handConn.block (true);
    cropConn.block (true);

    if (wbTool) {
        wbConn.block (true);
    }
    if (colPickerTool) {
        cpConn.block (true);
    }

    straConn.block (true);

    bool stopEdit = tool == TMHand && handTool->get_active() && editingMode;

    handTool->set_active (false);

    if (wbTool) {
        wbTool->set_active (false);
    }

    cropTool->set_active (false);
    straTool->set_active (false);
    if (colPickerTool) {
        colPickerTool->set_active (false);
    }

    if (tool == TMHand) {
        handTool->set_active (true);
        handTool->grab_focus(); // switch focus to the handTool button
    } else if (tool == TMSpotWB) {
        if (wbTool) {
            wbTool->set_active (true);
        }
    } else if (tool == TMCropSelect) {
        cropTool->set_active (true);
    } else if (tool == TMStraighten) {
        straTool->set_active (true);
    } else if (tool == TMColorPicker) {
        if (colPickerTool) {
            colPickerTool->set_active (true);
        }
    }

    current = tool;

    handConn.block (false);
    cropConn.block (false);

    if (wbTool) {
        wbConn.block (false);
    }
    if (colPickerTool) {
        cpConn.block (false);
    }

    straConn.block (false);

    if (stopEdit) {
        stopEditMode();

        if (listener) {
            listener->editModeSwitchedOff();
        }
    }
}

void ToolBar::startEditMode()
{
    if (!editingMode) {
        handTool->set_active(true);  // will call hand_pressed, with editingMode=false
        editingMode = true;
        handTool->set_image(*editinghandimg);
    }

#ifndef NDEBUG
    else {
        printf("Editing mode already active!\n");
    }

#endif
}

void ToolBar::stopEditMode()
{
    if (editingMode) {
        editingMode = false;
        /* WARNING: Should we toggle the Hand button on?
         * This method can be called while another tool is active, e.g. if the user toggle off
         * the Subscriber's Edit button. For now, we keep that other tool active. If one want to
         * switch to the Hand tool, uncommenting the following line should suffice (not tested).
         *
         * handTool->set_active(true);
         */
        handTool->set_image(*handimg);
    }
}

void ToolBar::hand_pressed ()
{

    handConn.block (true);
    cropConn.block (true);
    if (colPickerTool) {
        cpConn.block (true);
    }

    if (wbTool) {
        wbConn.block (true);
    }

    straConn.block (true);

    if (current != TMHand) {
        if (colPickerTool) {
            colPickerTool->set_active(false);
        }
        if (wbTool) {
            wbTool->set_active (false);
        }

        cropTool->set_active (false);
        straTool->set_active (false);
        current = TMHand;
    } else {
        if (editingMode) {
            stopEditMode();
        }

        if (listener) {
            listener->editModeSwitchedOff ();
        }
    }

    handTool->set_active (true);
    handConn.block (false);
    cropConn.block (false);
    if (colPickerTool) {
        cpConn.block (false);
    }
    if (wbTool) {
        wbConn.block (false);
    }

    straConn.block (false);

    if (listener) {
        listener->toolSelected (TMHand);
    }
}

void ToolBar::wb_pressed ()
{

    handConn.block (true);
    cropConn.block (true);
    if (colPickerTool) {
        cpConn.block (true);
    }
    if (wbTool) {
        wbConn.block (true);
    }

    straConn.block (true);

    if (current != TMSpotWB) {
        handTool->set_active (false);
        cropTool->set_active (false);
        straTool->set_active (false);
        if (colPickerTool) {
            colPickerTool->set_active(false);
        }
        current = TMSpotWB;
    }

    if (wbTool) {
        wbTool->set_active (true);
    }

    handConn.block (false);
    cropConn.block (false);
    if (colPickerTool) {
        cpConn.block (false);
    }
    if (wbTool) {
        wbConn.block (false);
    }

    straConn.block (false);

    if (listener) {
        listener->toolSelected (TMSpotWB);
    }
}

void ToolBar::colPicker_pressed (GdkEventButton* event)
{

    if (event->button == 1) {
        handConn.block (true);
        cropConn.block (true);
        cpConn.block (true);
        if (wbTool) {
            wbConn.block (true);
        }

        straConn.block (true);

        cropTool->set_active (false);
        if (wbTool) {
            wbTool->set_active (false);
        }
        straTool->set_active (false);

        if (current != TMColorPicker) {
            // Disabling all other tools, enabling the Picker tool and entering the "visible pickers" mode
            handTool->set_active (false);
            showColorPickers(true);
            current = TMColorPicker;
        } else {
            // Disabling the picker tool, enabling the Hand tool and keeping the "visible pickers" mode
            handTool->set_active (true);
            colPickerTool->set_active (false);
            current = TMHand;
        }

        handConn.block (false);
        cropConn.block (false);
        cpConn.block (false);
        wbConn.block (false);
        straConn.block (false);

        if (listener) {
            listener->toolSelected (current);
        }
    } else if (event->button == 3) {
        if (current == TMColorPicker) {
            // Disabling the Picker tool and entering into the "invisible pickers" mode
            cpConn.block (true);
            handConn.block (true);
            handTool->set_active (true);
            colPickerTool->set_active (false);
            current = TMHand;
            showColorPickers(false);
            cpConn.block (false);
            handConn.block (false);
        } else {
            // The Picker tool is already disabled, entering into the "invisible pickers" mode
            switchColorPickersVisibility();
        }
        if (pickerListener) {
            pickerListener->switchPickerVisibility (showColPickers);
        }
    }
}

bool ToolBar::showColorPickers(bool showCP)
{
    if (showColPickers != showCP) {
        // Inverting the state
        colPickerTool->set_image(showCP ? *showcolpickersimg : *hidecolpickersimg);
        showColPickers = showCP;
        return true;
    }

    return false;
}

void ToolBar::switchColorPickersVisibility()
{
    // Inverting the state
    showColPickers = !showColPickers;
    colPickerTool->set_image(showColPickers ? *showcolpickersimg : *hidecolpickersimg);
}

void ToolBar::crop_pressed ()
{

    handConn.block (true);
    cropConn.block (true);
    cpConn.block(true);

    if (wbTool) {
        wbConn.block (true);
    }

    straConn.block (true);

    if (current != TMCropSelect) {
        handTool->set_active (false);
        if (colPickerTool) {
            colPickerTool->set_active(false);
        }
        if (wbTool) {
            wbTool->set_active (false);
        }

        straTool->set_active (false);
        current = TMCropSelect;
    }

    cropTool->set_active (true);
    handConn.block (false);
    cropConn.block (false);
    cpConn.block(false);
    wbConn.block (false);
    straConn.block (false);

    if (listener) {
        listener->toolSelected (TMCropSelect);
    }
}

void ToolBar::stra_pressed ()
{

    handConn.block (true);
    cropConn.block (true);
    if (colPickerTool) {
        cpConn.block (true);
    }

    if (wbTool) {
        wbConn.block (true);
    }

    straConn.block (true);

    if (current != TMStraighten) {
        handTool->set_active (false);
        if (colPickerTool) {
            colPickerTool->set_active(false);
        }
        if (wbTool) {
            wbTool->set_active (false);
        }

        cropTool->set_active (false);
        current = TMStraighten;
    }

    straTool->set_active (true);
    handConn.block (false);
    cropConn.block (false);
    cpConn.block (false);

    if (wbTool) {
        wbConn.block (false);
    }

    straConn.block (false);

    if (listener) {
        listener->toolSelected (TMStraighten);
    }
}

bool ToolBar::handleShortcutKey (GdkEventKey* event)
{

    bool ctrl = event->state & GDK_CONTROL_MASK;
    //bool shift = event->state & GDK_SHIFT_MASK;
    bool alt = event->state & GDK_MOD1_MASK;

    if (!ctrl && !alt) {
        switch(event->keyval) {
        case GDK_w:
        case GDK_W:
            if(wbTool) {
                wb_pressed ();
                return true;
            }

            return false;

        case GDK_c:
        case GDK_C:
            crop_pressed ();
            return true;

        case GDK_s:
        case GDK_S:
            stra_pressed ();
            return true;

        case GDK_h:
        case GDK_H:
            hand_pressed ();
            return true;
        }
    } else {
        switch (event->keyval) {
        }
    }

    return false;
}

void ToolBar::setBatchMode()
{
    if (wbTool) {
        wbConn.disconnect();
        removeIfThere(this, wbTool, false);
        wbTool = nullptr;
    }
    if (colPickerTool) {
        cpConn.disconnect();
        removeIfThere(this, colPickerTool, false);
        colPickerTool = nullptr;
    }
}

