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

#include <glibmm.h>
#include <rtengine.h>
#include <procparams.h>
#include <multilangmgr.h>
#include <paramsedited.h>

class ToolPanelListener {

    public:
    
        virtual void panelChanged   (rtengine::ProcEvent event, const Glib::ustring& descr) {}
};

class ToolPanel {

    protected:
        ToolPanelListener* listener;
        ToolPanelListener* tmp;
        bool batchMode;

    public:

        ToolPanel () : listener(NULL), tmp(NULL), batchMode(false) {}

                void setListener (ToolPanelListener* tpl) { listener = tpl; }
        virtual void read        (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL) {}
        virtual void write       (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL) {}
        virtual void setDefaults (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL) {}
        
                void disableListener () { if (tmp==NULL) tmp = listener; listener = NULL; }
                void enableListener  () { if (tmp!=NULL) listener = tmp; tmp = NULL; }
                
        virtual void setBatchMode    (bool batchMode) { this->batchMode = batchMode; }
                
};  

#endif
