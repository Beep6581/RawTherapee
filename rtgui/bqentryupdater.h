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
#ifndef _BQENTRYUPDATER_
#define _BQENTRYUPDATER_

#include <glibmm.h>
#include <rtengine.h>
#include <thumbnail.h>

class BQEntryUpdateListener {

    public:
        virtual void updateImage (guint8* img, int w, int h) {}
};

class BatchQueueEntryUpdater {

    struct Job {
        guint8* oimg;
        int ow, oh, newh;
        BQEntryUpdateListener* listener;
    };

  protected:
    bool tostop;
    bool stopped;
    std::list<Job> jqueue;
    Glib::Thread* thread;
    Glib::Mutex* qMutex;

  public:
    BatchQueueEntryUpdater ();

    void add        (guint8* oimg, int ow, int oh, int newh, BQEntryUpdateListener* listener);
    void process    ();
    void stop       ();
    void removeJobs ();
    void removeJobs (BQEntryUpdateListener* listener);
    void terminate  ();

    void process_   ();
};

extern BatchQueueEntryUpdater batchQueueEntryUpdater;

#endif
