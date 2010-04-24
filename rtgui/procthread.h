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
#ifndef _PROCTHREAD_
#define _PROCTHREAD_

#include <gtkmm.h>

template <class T>
class ProcessingThread {

  protected:
    bool tostop;
    bool stopped;
    std::list<T> jqueue;
    Glib::Thread* thread;
    bool fifo;

  public:
    ProcessingThread () : tostop(false), stopped(true), fifo(true) {}
    void init () { tostop = false; stopped = true; }
    void add (T de) { jqueue.push_back (de); }
    void removejobs () { while (!jqueue.empty()) jqueue.pop_front (); }
    void stop () { if (stopped) { tostop = true; return; } gdk_threads_leave(); tostop = true; Glib::Thread::self()->yield(); if (!stopped) thread->join (); gdk_threads_enter();}
    
    virtual void start () {}
    virtual void process (T& e) {}
    virtual void processCustomOrder () {}
    virtual void end () {}

    void process () { 
        if (stopped)
            #undef THREAD_PRIORITY_NORMAL
            thread = Glib::Thread::create(sigc::mem_fun(*this, &ProcessingThread::process_), (unsigned long int)0, true, true, Glib::THREAD_PRIORITY_NORMAL);
    }
    void process_ () { 
        stopped = false;
        tostop = false;
    
        start ();     //   jqueue.sort ();


        while (!tostop && !jqueue.empty ()) {
            if (fifo) {
                T current = jqueue.front ();
                jqueue.pop_front ();
                process (current);
            }
            else
                processCustomOrder ();
        }
        stopped = true;
        end ();
    }
    
    bool runs () { return !stopped; }
    
    void terminate () { stop (); removejobs (); }
    
    int numOfJobs () { return jqueue.size(); }
};

#endif
