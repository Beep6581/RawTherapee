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
#include "bqentryupdater.h"
#include <gtkmm.h>
#include "guiutils.h"

BatchQueueEntryUpdater batchQueueEntryUpdater;

BatchQueueEntryUpdater::BatchQueueEntryUpdater ()
    : tostop(false), stopped(true), thread(NULL), qMutex(NULL)
{
}

void BatchQueueEntryUpdater::process (guint8* oimg, int ow, int oh, int newh, BQEntryUpdateListener* listener, rtengine::ProcParams* pparams, Thumbnail* thumbnail)
{
    if (!oimg && (!pparams || !thumbnail)) {
        //printf("WARNING! !oimg && (!pparams || !thumbnail)\n");
        return;
    }

    if (!qMutex) {
        qMutex = new MyMutex ();
    }

    qMutex->lock ();
    // look up if an older version is in the queue
    std::list<Job>::iterator i;

    for (i = jqueue.begin(); i != jqueue.end(); ++i)
        if (i->oimg == oimg && i->listener == listener) {
            i->ow = ow;
            i->oh = oh;
            i->newh = newh;
            i->listener = listener;
            i->pparams = pparams;
            i->thumbnail = thumbnail;
            break;
        }

    // not found, create and append new job
    if (i == jqueue.end ()) {
        Job j;
        j.oimg = oimg;
        j.ow = ow;
        j.oh = oh;
        j.newh = newh;
        j.listener = listener;
        j.pparams = pparams;
        j.thumbnail = thumbnail;
        jqueue.push_back (j);
    }

    qMutex->unlock ();

    // Start thread if not running yet
    if (stopped) {
        stopped = false;
        tostop  = false;

#if __GNUC__ == 4 && __GNUC_MINOR__ == 8 && defined( WIN32 ) && defined(__x86_64__)
#undef THREAD_PRIORITY_NORMAL
        // See Issue 2384 comment #3
        thread = Glib::Thread::create(sigc::mem_fun(*this, &BatchQueueEntryUpdater::processThread), (unsigned long int)0, true, true, Glib::THREAD_PRIORITY_NORMAL);
#else
#undef THREAD_PRIORITY_LOW
        thread = Glib::Thread::create(sigc::mem_fun(*this, &BatchQueueEntryUpdater::processThread), (unsigned long int)0, true, true, Glib::THREAD_PRIORITY_LOW);
#endif
    }
}

void BatchQueueEntryUpdater::processThread ()
{
    // TODO: process visible jobs first
    bool isEmpty = false;

    while (!tostop && !isEmpty) {

        qMutex->lock ();
        isEmpty = jqueue.empty (); // do NOT put into while() since it must be within mutex section
        Job current;

        if (!isEmpty) {
            current = jqueue.front ();
            jqueue.pop_front ();
        }

        qMutex->unlock ();

        if(isEmpty) {
            break;
        }

        rtengine::IImage8* img = NULL;
        bool newBuffer = false;

        if (current.thumbnail && current.pparams) {
            // the thumbnail and the pparams are provided, it means that we have to build the original preview image
            double tmpscale;
            img = current.thumbnail->processThumbImage (*current.pparams, current.oh, tmpscale);

            //current.thumbnail->decreaseRef (); // WARNING: decreasing refcount (and maybe deleting) thumbnail, with or without processed image
            if (img) {
                int prevw = img->getWidth();
                int prevh = img->getHeight();
#ifndef NDEBUG

                if (current.ow != img->getW() || current.oh != img->getH()) {
                    printf("WARNING!  Expected image size: %dx%d ; image size is: %dx%d\n", current.ow, current.oh, img->getW(), img->getH());
                }

                assert ((current.ow + 1)*current.oh >= img->getW()*img->getH());
#endif
                current.ow = prevw;
                current.oh = prevh;

                if (!current.oimg) {
                    current.oimg = new guint8[prevw * prevh * 3];
                    newBuffer = true;
                }

                memcpy(current.oimg, img->getData(), prevw * prevh * 3);
                img->free();
            }
        }

        if (current.oimg && !isEmpty && current.listener) {
            int neww = current.newh * current.ow / current.oh;
            guint8* img = new guint8 [current.newh * neww * 3];
            thumbInterp (current.oimg, current.ow, current.oh, img, neww, current.newh);
            current.listener->updateImage (img, neww, current.newh, current.ow, current.oh, newBuffer ? current.oimg : NULL);
        }

        if(current.oimg) {
            delete[] current.oimg;
            current.oimg = NULL;
        }
    }

    stopped = true;
}


void BatchQueueEntryUpdater::removeJobs (BQEntryUpdateListener* listener)
{
    if (!qMutex) {
        return;
    }

    qMutex->lock ();
    bool ready = false;

    while (!ready) {
        ready = true;
        std::list<Job>::iterator i;

        for (i = jqueue.begin(); i != jqueue.end(); ++i)
            if (i->listener == listener) {
                jqueue.erase (i);
                ready = false;
                break;
            }
    }

    qMutex->unlock ();
}

void BatchQueueEntryUpdater::terminate  ()
{
    // never started or currently not running?
    if (!qMutex || stopped) {
        return;
    }

    if (!stopped) {
        // Yield to currently running thread and wait till it's finished
        GThreadUnLock lock;
        tostop = true;
        Glib::Thread::self()->yield();

        if (!stopped) {
            thread->join ();
        }
    }

    // Remove remaining jobs
    qMutex->lock ();

    while (!jqueue.empty()) {
        jqueue.pop_front ();
    }

    qMutex->unlock ();
}


