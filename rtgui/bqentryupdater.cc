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
#include <bqentryupdater.h>
#include <gtkmm.h>
#include <guiutils.h>

BatchQueueEntryUpdater batchQueueEntryUpdater;

BatchQueueEntryUpdater::BatchQueueEntryUpdater () 
    : tostop(false), stopped(true), qMutex(NULL) {
}

void BatchQueueEntryUpdater::add (guint8* oimg, int ow, int oh, int newh, BQEntryUpdateListener* listener) {

    if (!qMutex)
        qMutex = new Glib::Mutex ();    

    qMutex->lock ();
    // look up if an older version is in the queue
    std::list<Job>::iterator i;
    for (i=jqueue.begin(); i!=jqueue.end(); i++)
        if (i->oimg==oimg && i->listener==listener) {
            i->ow = ow;
            i->oh = oh;
            i->newh = newh;
            i->listener = listener;
            break;
        }
    // not found, create and append new job
    if (i==jqueue.end ()) {
        Job j;
        j.oimg = oimg;
        j.ow = ow;
        j.oh = oh;
        j.newh = newh;
        j.listener = listener;
        jqueue.push_back (j);
    }
    qMutex->unlock ();
}

void BatchQueueEntryUpdater::process () {

    if (stopped)
        #undef THREAD_PRIORITY_NORMAL
        thread = Glib::Thread::create(sigc::mem_fun(*this, &BatchQueueEntryUpdater::process_), (unsigned long int)0, true, true, Glib::THREAD_PRIORITY_NORMAL);
}

void BatchQueueEntryUpdater::process_ () { 

    stopped = false;
    tostop = false;

// TODO: process visible jobs first
    while (!tostop && !jqueue.empty ()) {
        qMutex->lock ();
        Job current = jqueue.front ();
        jqueue.pop_front ();
        qMutex->unlock ();
        if (current.listener) {
            int neww = current.newh * current.ow / current.oh;
            guint8* img = new guint8 [current.newh*neww*3];
            thumbInterp (current.oimg, current.ow, current.oh, img, neww, current.newh);
            current.listener->updateImage (img, neww, current.newh);
        }
    }
    stopped = true;
}

void BatchQueueEntryUpdater::stop () {

    if (stopped) {
        tostop = true; 
        return; }
        
    gdk_threads_leave(); 
    tostop = true; 
    Glib::Thread::self()->yield(); 
    if (!stopped) 
        thread->join (); 
    gdk_threads_enter();
}

void BatchQueueEntryUpdater::removeJobs () {

    if (!qMutex)
        return;

    qMutex->lock ();
    while (!jqueue.empty()) 
        jqueue.pop_front (); 
    qMutex->unlock ();
}

void BatchQueueEntryUpdater::removeJobs (BQEntryUpdateListener* listener) {

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

void BatchQueueEntryUpdater::terminate  () { 

    stop (); 
    removeJobs (); 
}


