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

class ToolBarListener {

  public:
    virtual void toolSelected (ToolMode tool) {}

};

class ToolBar : public Gtk::HBox {

  protected:
    Gtk::ToggleButton* handTool;
    Gtk::ToggleButton* wbTool;
    Gtk::ToggleButton* cropTool;
    Gtk::ToggleButton* straTool;
    ToolBarListener* listener;
    ToolMode current;
    sigc::connection  handConn;
    sigc::connection  wbConn;
    sigc::connection  cropConn;
    sigc::connection  straConn;

  public:
    ToolBar ();

    void     setTool (ToolMode tool);
    ToolMode getTool () { return current; }

    void setToolBarListener (ToolBarListener* tpl) { listener = tpl; }

    void hand_pressed ();
    void wb_pressed ();
    void crop_pressed ();
    void stra_pressed ();

    bool handleShortcutKey (GdkEventKey* event);
    void removeWbTool();
};

#endif
