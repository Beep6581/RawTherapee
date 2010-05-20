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
#include <rtengine.h>

#undef THREAD_PRIORITY_NORMAL

class PLDBridge : public rtengine::ProgressListener {

        Gtk::Label* label;
        Gtk::ProgressBar* progBar;

    public:
        PLDBridge ( Gtk::Label* l, Gtk::ProgressBar* pb)
            : label(l), progBar(pb) {}

    // ProgressListener interface
    void setProgress (double p) {
        gdk_threads_enter ();
        progBar->set_fraction (p);
        gdk_threads_leave ();
    }
    void setProgressStr (Glib::ustring str) {
        gdk_threads_enter ();
        Glib::ustring progrstr;
        if (str=="Decoding...")
            progrstr = M("PROGRESSBAR_DECODING");
        else if (str=="Ready.")
            progrstr = M("PROGRESSBAR_READY");
        else if (str=="Demosaicing...")
            progrstr = M("PROGRESSBAR_DEMOSAICING");
        else if (str=="Loading...")
            progrstr = M("PROGRESSBAR_LOADING");
        else if (str=="Loading PNG file...")
            progrstr = M("PROGRESSBAR_LOADPNG");
        else if (str=="Loading JPEG file...")
            progrstr = M("PROGRESSBAR_LOADJPEG");
        else if (str=="Loading TIFF file...")
            progrstr = M("PROGRESSBAR_LOADTIFF");
        else if (str=="Saving PNG file...")
            progrstr = M("PROGRESSBAR_SAVEPNG");
        else if (str=="Saving JPEG file...")
            progrstr = M("PROGRESSBAR_SAVEJPEG");
        else if (str=="Saving TIFF file...")
            progrstr = M("PROGRESSBAR_SAVETIFF");
        else if (str=="Processing...")
            progrstr = M("PROGRESSBAR_PROCESSING");
        else 
            progrstr = str;

        label->set_text (progrstr);
        gdk_threads_leave ();
    }
};

template<class T>
class ProgressConnector {

        sigc::signal0<T> opStart;
        sigc::signal0<bool> opEnd;
        T retval;
        Glib::Thread *workThread;

		static int emitEndSignal (void* data) {
			gdk_threads_enter ();
			sigc::signal0<bool>* opEnd = (sigc::signal0<bool>*) data;
			int r = opEnd->emit ();
			delete opEnd;
			gdk_threads_leave ();
			return r;
		}
        
        void workingThread () {
            retval = opStart.emit ();
            g_idle_add (ProgressConnector<T>::emitEndSignal, new sigc::signal0<bool> (opEnd));
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
