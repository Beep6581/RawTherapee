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
#ifndef _PROGRESSCONNECTOR_
#define _PROGRESSCONNECTOR_

#include <sigc++/sigc++.h>
#include <gtkmm.h>
#include "../rtengine/rtengine.h"
#include "guiutils.h"

#undef THREAD_PRIORITY_NORMAL

class PLDBridge : public rtengine::ProgressListener {

	    rtengine::ProgressListener* pl;

    public:
        PLDBridge ( rtengine::ProgressListener* pb)
            : pl(pb) {}

    // ProgressListener interface
    void setProgress (double p) {
        GThreadLock lock;
        pl->setProgress(p);
    }
    void setProgressStr (Glib::ustring str) {
        GThreadLock lock;
        Glib::ustring progrstr;
        progrstr = M(str);
        pl->setProgressStr(progrstr);
    }

    void setProgressState (bool inProcessing){
        GThreadLock lock;
        pl->setProgressState(inProcessing);
    }

    void error (Glib::ustring descr){
        GThreadLock lock;
        pl->error(descr);
    }
};

template<class T>
class ProgressConnector {

        sigc::signal0<T> opStart;
        sigc::signal0<bool> opEnd;
        T retval;
        Glib::Thread *workThread;

		static int emitEndSignalUI (void* data) {

			sigc::signal0<bool>* opEnd = (sigc::signal0<bool>*) data;
			int r = opEnd->emit ();
			delete opEnd;

			return r;
		}
        
        void workingThread () {
            retval = opStart.emit ();
            g_idle_add (ProgressConnector<T>::emitEndSignalUI, new sigc::signal0<bool> (opEnd));
            workThread = 0;
        }
        
    public:
    
        ProgressConnector ():workThread( 0 ) { }

        void startFunc (const sigc::slot0<T>& startHandler, const sigc::slot0<bool>& endHandler ) {
        	if( !workThread ){
				opStart.connect (startHandler);
				opEnd.connect (endHandler);
				workThread = Glib::Thread::create(sigc::mem_fun(*this, &ProgressConnector<T>::workingThread), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        	}
        }

        T returnValue(){
        	return retval;
        }
};
#endif
