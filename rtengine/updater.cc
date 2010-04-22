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
#include <updater.h>
#include <gdkmm.h>

Updater::Updater () :
    change (0),
    ipc (NULL),
    pl (NULL),
    running (false)
{
}

ProcParams* Updater::changing (int what) {

    mutex.lock ();
    change |= what;
    return &params;
}

void Updater::changed () {

    mutex.unlock ();
    startProcessing ();
}

void Updater::clearState  () {

    mutex.lock ();
    change = 0;
    mutex.unlock ();
}

ProcParams* Updater::getParams () {

    return &params;
}

void Updater::startProcessing () {

    #undef THREAD_PRIORITY_NORMAL

    tstart.lock ();
    if (ipc && !running) {
        running = true;
        tstart.unlock ();
        Glib::Thread::create(sigc::mem_fun(*this, &Updater::process), 0, false, true, Glib::THREAD_PRIORITY_NORMAL);    
    }
    else
        tstart.unlock ();
}

void Updater::process () {

    processing.lock ();
    if (pl)
        pl->setProcessingState (true);

    int ch;
 
    mutex.lock ();
    while (change && ipc) {
        ipc->params.copy (&params);
        ch = change;
        change = 0;
        mutex.unlock ();
        if (ch<=16384)
            ipc->update (ch);
        mutex.lock ();
    }    
    mutex.unlock ();
    tstart.lock ();
    running = false;
    tstart.unlock ();

    if (pl)
        pl->setProcessingState (false);
    processing.unlock ();
}

void Updater::stop () {

    if (running) {
        change = 0;
        processing.lock ();
        processing.unlock ();
        
    }
}

