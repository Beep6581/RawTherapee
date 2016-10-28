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
#include "../rtengine/rtengine.h"
#include "threadutils.h"
#include "thumbnail.h"

class BQEntryUpdateListener
{

public:
    virtual ~BQEntryUpdateListener () {}
    virtual void updateImage (guint8* img, int w, int h, int origw, int origh, guint8* newOPreview) {}
};

class BatchQueueEntryUpdater
{

    struct Job {
        guint8* oimg;
        int ow, oh, newh;
        BQEntryUpdateListener* listener;
        rtengine::ProcParams* pparams;
        Thumbnail* thumbnail;
    };

protected:
    bool tostop;
    bool stopped;
    std::list<Job> jqueue;
    Glib::Thread* thread;
    MyMutex* qMutex;

public:
    BatchQueueEntryUpdater ();

    void process    (guint8* oimg, int ow, int oh, int newh, BQEntryUpdateListener* listener, rtengine::ProcParams* pparams = nullptr, Thumbnail* thumbnail = nullptr);
    void removeJobs (BQEntryUpdateListener* listener);
    void terminate  ();

    void processThread ();
};

extern BatchQueueEntryUpdater batchQueueEntryUpdater;

#endif
