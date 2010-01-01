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
#ifndef _PROGRESSDIALOG_
#define _PROGRESSDIALOG_

#include <sigc++/sigc++.h>
#include <gtkmm.h>
#include <rtengine.h>

class PLDBridge : public rtengine::ProgressListener {

        Gtk::Dialog* dialog;
        Gtk::Label* label;
        Gtk::ProgressBar* progBar;

    public:
        PLDBridge (Gtk::Dialog* d, Gtk::Label* l, Gtk::ProgressBar* pb) 
            : dialog(d), label(l), progBar(pb) {}

    // progresslistener interface
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
    void setProgressState (int state) {}
    void error (Glib::ustring descr) {}
};

template<class T>
class ProgressDialog : public Gtk::Dialog {

        sigc::signal0<T> operation;
        T* retval;
        Gtk::Label prLabel;
        Gtk::ProgressBar prProgBar;
        
        PLDBridge* pldBridge;
        
        
        void workingThread () {
            *retval = operation.emit ();
            gdk_threads_enter ();
            response (1);
            gdk_threads_leave ();
        }
        
    public:
    
        ProgressDialog (Glib::ustring label) : Gtk::Dialog (label, true) {
            pldBridge = new PLDBridge (this, &prLabel, &prProgBar); 
            get_vbox()->pack_start (prLabel, Gtk::PACK_SHRINK, 4);
            get_vbox()->pack_start (prProgBar, Gtk::PACK_SHRINK, 4);
            set_size_request (300, -1);
            show_all_children ();
        }
        
        ~ProgressDialog () {
            delete pldBridge;
        }
        
        rtengine::ProgressListener* getProgressListener () { return pldBridge; }
        
        void setFunc (const sigc::slot0<T>& slot, T* rv) {
            retval = rv;
            operation.connect (slot);
        }
        
        void start () {
            Glib::Thread *thread = Glib::Thread::create(sigc::mem_fun(*this, &ProgressDialog<T>::workingThread), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
            int x = run ();
            if (x<0) {
                gdk_threads_leave ();
                thread->join ();
                gdk_threads_enter ();
            }
        }
};
#endif
