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
#include <thumbimageupdater.h>
#include <gtkmm.h>

#define threadNum 8
ThumbImageUpdater thumbImageUpdater;

ThumbImageUpdater::ThumbImageUpdater () 
    : tostop(false), stopped(true), qMutex(NULL), startMutex(NULL), threadPool(NULL) {
        
}

ThumbImageUpdater::~ThumbImageUpdater ()
{
    delete threadPool;
}

void ThumbImageUpdater::add (Thumbnail* t, const rtengine::procparams::ProcParams& params, int height, bool* priority, ThumbImageUpdateListener* l) {

    if (!qMutex)
        qMutex = new Glib::Mutex ();
    if (!startMutex)
        startMutex = new Glib::Mutex ();

    qMutex->lock ();
    // look up if an older version is in the queue
    std::list<Job>::iterator i;
    for (i=jqueue.begin(); i!=jqueue.end(); i++)
        if (i->thumbnail==t && i->listener==l) {
            i->pparams = params;
            i->height = height;
            i->priority = priority;
            break;
        }
    // not found, create and append new job
    if (i==jqueue.end ()) {
        Job j;
        j.thumbnail = t;
        j.pparams = params;
        j.height = height;
        j.listener = l;
        j.priority = priority;
        jqueue.push_back (j);
    }
    qMutex->unlock ();
}

void ThumbImageUpdater::process () {

    if (stopped) {
        stopped = false;

        if(!threadPool)
            threadPool = new Glib::ThreadPool(threadNum,0);
     
        //thread = Glib::Thread::create (sigc::mem_fun(*this, &ThumbImageUpdater::process_), (unsigned long int)0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        process_();
    }
}

void ThumbImageUpdater::process_ () { 

    stopped = false;
    tostop = false;

    while (!tostop && !jqueue.empty ()) {

            std::list<Job>::iterator i;
            for (i=jqueue.begin (); i!=jqueue.end(); i++)
                if (*(i->priority))
                    break;

            if (i==jqueue.end())
                i = jqueue.begin();

            Job current = *i;         
            if (current.listener)                    
                    threadPool->push(sigc::bind(sigc::mem_fun(*this, &ThumbImageUpdater::processJob), current));

            jqueue.erase (i);
    }
    
    stopped = true;
    //printf("Threads # %d \n", threadPool->get_num_threads());

}

void ThumbImageUpdater::processJob (Job current) {

    if (current.listener) {
        double scale = 1.0;
        rtengine::IImage8* img = current.thumbnail->processThumbImage (current.pparams, current.height, scale);
        if (img)
            current.listener->updateImage (img, scale, current.pparams.crop);
    }

}

void ThumbImageUpdater::stop () {
        
    gdk_threads_leave(); 
    tostop = true;

    if (threadPool)    {
        threadPool->shutdown(TRUE);
        threadPool = NULL;
    }
    gdk_threads_enter();
}

void ThumbImageUpdater::removeJobs () {

    if (!qMutex)
        return;

    qMutex->lock ();
    while (!jqueue.empty()) 
        jqueue.pop_front (); 
    qMutex->unlock ();
}

void ThumbImageUpdater::removeJobs (ThumbImageUpdateListener* listener) {

    if (!qMutex)
        return;

    qMutex->lock ();
    bool ready = false;
    while (!ready) {
        ready = true;
        std::list<Job>::iterator i;
        for (i=jqueue.begin(); i!=jqueue.end(); i++)
            if (i->listener == listener) {
                jqueue.erase (i);
                ready = false;
                break;
            }
    }
    qMutex->unlock ();
}

void ThumbImageUpdater::terminate  () { 

    stop (); 
    removeJobs (); 
}


