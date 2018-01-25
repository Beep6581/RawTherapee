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
#ifndef __TOOLPANEL__
#define __TOOLPANEL__

#include <gtkmm.h>
#include <glibmm.h>
#include "../rtengine/rtengine.h"
#include "../rtengine/procparams.h"
#include "guiutils.h"
#include "multilangmgr.h"
#include "paramsedited.h"
#include "edit.h"

class ToolPanel;
class FoldableToolPanel;

class ToolPanelListener
{

public:

    virtual ~ToolPanelListener() {}
    virtual void panelChanged   (rtengine::ProcEvent event, const Glib::ustring& descr) {}
};

/// @brief This class control the space around the group of tools inside a tab, as well as the space separating each tool. */
class ToolVBox : public Gtk::VBox
{
public:
    ToolVBox();
};

/// @brief This class control the space around a tool's block of parameter. */
class ToolParamBlock : public Gtk::VBox
{
public:
    ToolParamBlock();
};

class ToolPanel
{

protected:
    Glib::ustring toolName;
    ToolPanelListener* listener;
    ToolPanelListener* tmp;
    bool batchMode;  // True if the ToolPanel is used in Batch mode
    bool multiImage; // True if more than one image are being edited at the same time (also imply that batchMode=true), false otherwise
    bool need100Percent;

public:

    ToolPanel (Glib::ustring toolName = "", bool need11 = false) : toolName(toolName), listener(nullptr), tmp(nullptr), batchMode(false), multiImage(false), need100Percent(need11) {}
    virtual ~ToolPanel() {}

    virtual void           setParent       (Gtk::Box* parent) {}
    virtual Gtk::Box*      getParent       ()
    {
        return nullptr;
    }
    virtual MyExpander*    getExpander     ()
    {
        return nullptr;
    }
    virtual void           setExpanded     (bool expanded) {}
    virtual bool           getExpanded     ()
    {
        return false;
    }
    void           setMultiImage   (bool m)
    {
        multiImage = m;
    }
    virtual void           setListener     (ToolPanelListener* tpl)
    {
        listener = tpl;
    }
    virtual void           setEditProvider (EditDataProvider *provider) {}
    virtual void           read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) {}
    virtual void           write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) {}
    virtual void           trimValues      (rtengine::procparams::ProcParams* pp)
    {
        return;
    }
    virtual void           setDefaults     (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) {}
    virtual void           autoOpenCurve   () {}

    /** @brief Disable the event broadcasting mechanism
     *
     * @return Return the previous state of the broadcast (true: enabled ; false: disabled)
     */
    bool disableListener ()
    {
        if (tmp == nullptr) {
            tmp = listener;
        }

        bool prevState = listener != nullptr;
        listener = nullptr;
        return prevState;
    }

    /** @brief Enable the event broadcasting mechanism
     */
    void enableListener  ()
    {
        if (tmp != nullptr) {
            listener = tmp;
        }

        tmp = nullptr;
    }

    virtual void setBatchMode    (bool batchMode)
    {
        this->batchMode = batchMode;
    }

};

class FoldableToolPanel : public ToolPanel
{

protected:
    Gtk::Box* parentContainer;
    MyExpander* exp;
    bool lastEnabled;
    sigc::connection enaConn;
    void foldThemAll (GdkEventButton* event);
    void enabled_toggled();

public:

    FoldableToolPanel(Gtk::Box* content, Glib::ustring toolName, Glib::ustring UILabel, bool need11 = false, bool useEnabled = false);

    MyExpander* getExpander()
    {
        return exp;
    }
    void setExpanded (bool expanded)
    {
        if (exp) {
            exp->set_expanded( expanded );
        }
    }

    void hide() {
        if (exp && !batchMode) {  // conditional hide
            exp->hide();
        }
    }

    void show() {
        if (exp) {                // always show
            exp->show();
        }
    }
    bool getExpanded ()
    {
        if (exp) {
            return exp->get_expanded();
        }

        return false;
    }
    void setParent (Gtk::Box* parent)
    {
        parentContainer = parent;
    }
    Gtk::Box* getParent ()
    {
        return parentContainer;
    }

    virtual void enabledChanged  () {}

    bool getUseEnabled ()
    {
        if (exp) {
            return exp->getUseEnabled();
        } else {
            return true;
        }
    }
    bool getEnabled();  // related to the enabled/disabled state
    void setEnabled(bool isActive);  // related to the enabled/disabled state
    void setEnabledTooltipMarkup(Glib::ustring tooltipMarkup);
    void setEnabledTooltipText(Glib::ustring tooltipText);
    bool get_inconsistent();  // related to the enabled/disabled state
    void set_inconsistent(bool isInconsistent);  // related to the enabled/disabled state

    void setLevel (int level);

    // Functions that want to receive an enabled/disabled event from this class
    // will have to receive it from MyExpander directly, we do not create
    // a relaying event
    MyExpander::type_signal_enabled_toggled signal_enabled_toggled()
    {
        return exp->signal_enabled_toggled();
    }
};

#endif
