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
#ifndef _THUMBIMAGEUPDATER_
#define _THUMBIMAGEUPDATER_

#include <glibmm.h>
#include <rtengine.h>
#include <thumbnail.h>

class ThumbImageUpdateListener {

    public:
        virtual void updateImage (rtengine::IImage8* img, double scale, rtengine::procparams::CropParams cropParams) {}
};

class ThumbImageUpdater {

    struct Job {
        Thumbnail* thumbnail;
        rtengine::procparams::ProcParams pparams;
        int height;
        bool* priority;
        ThumbImageUpdateListener* listener;
    };

  protected:
    bool tostop;
    bool stopped;
    std::list<Job> jqueue;
    Glib::Thread* thread;
    Glib::Mutex* qMutex;
    Glib::Mutex* startMutex;

  public:
    ThumbImageUpdater ();

    void add        (Thumbnail* t, const rtengine::procparams::ProcParams& params, int height, bool* priority, ThumbImageUpdateListener* l);
    void process    ();
    void stop       ();
    void removeJobs ();
    void removeJobs (ThumbImageUpdateListener* listener);
    void terminate  ();

    void process_   ();
    void processJob (Job current);
};

extern ThumbImageUpdater thumbImageUpdater;

#endif
