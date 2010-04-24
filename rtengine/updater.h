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
#ifndef _UPDATER_
#define _UPDATER_

#include <procparams.h>
#include <glibmm.h>
#include <progresslistener.h>
#include <improccoordinator.h>

class Updater {

    protected:
        int change;
        ProcParams params;
        Glib::Mutex mutex;
        Glib::Mutex tstart;
        Glib::Mutex processing;
        ImProcCoordinator* ipc;        
        ProgressListener* pl;
        bool running;
                
    public:
    
        Updater ();
        void        setProgressListener (ProgressListener* l) { pl = l; }
        void        setIPC              (ImProcCoordinator* ipc) { this->ipc = ipc; }

        ProcParams* changing        (int what);
        void        changed         ();
        void        clearState      ();
        int         getClear        ();
        ProcParams* getParams       ();
        int         getChange       ();
        void        startProcessing ();
        void        process         ();

        void        stop                ();
};

#endif
