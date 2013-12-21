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
#include "multilangmgr.h"
#include "paramsedited.h"

class ToolPanel;
class FoldableToolPanel;

class ToolPanelListener {

    public:
    
        virtual ~ToolPanelListener() {}
        virtual void panelChanged   (rtengine::ProcEvent event, const Glib::ustring& descr) {}
};

class ToolPanel {

    protected:
        ToolPanelListener* listener;
        ToolPanelListener* tmp;
        bool batchMode;  // True if the ToolPanel is used in Batch mode
        bool multiImage; // True if more than one image are being edited at the same time (also imply that batchMode=true), false otherwise

    public:

        ToolPanel () : listener(NULL), tmp(NULL), batchMode(false), multiImage(false) {}
        virtual ~ToolPanel() {}

                void setParent       (Gtk::Box* parent) {}
           Gtk::Box* getParent       () { return NULL; }
                void setMultiImage   (bool m) { multiImage = m; }
                void setListener     (ToolPanelListener* tpl) { listener = tpl; }
        virtual void read            (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL) {}
        virtual void write           (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL) {}
        virtual void trimValues      (rtengine::procparams::ProcParams* pp) { return; }
        virtual void setDefaults     (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL) {}
        virtual void autoOpenCurve   () {}

        /** @brief Disable the event broadcasting mechanism
         *
         * @return Return the previous state of the broadcast (true: enabled ; false: disabled)
         */
                bool disableListener () { if (tmp==NULL) tmp = listener; bool prevState = listener!=NULL; listener = NULL; return prevState; }

        /** @brief Enable the event broadcasting mechanism
         */
                void enableListener  () { if (tmp!=NULL) listener = tmp; tmp = NULL; }

        virtual void setBatchMode    (bool batchMode) { this->batchMode = batchMode; }

};

class FoldableToolPanel : public ToolPanel {

	protected:
        Gtk::Box* parentContainer;
        void foldThemAll (GdkEventButton* event);

    public:
        Gtk::Expander* exp;

        FoldableToolPanel(Gtk::Box* content);

        void setParent (Gtk::Box* parent) { parentContainer = parent; }
        Gtk::Box* getParent () { return parentContainer; }
        void setLabel (Glib::ustring label, bool need100Percent=false);
};

#endif
