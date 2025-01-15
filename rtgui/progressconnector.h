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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>

#include <sigc++/sigc++.h>

#include "guiutils.h"
#include "multilangmgr.h"
#include "rtengine/rtengine.h"

#undef THREAD_PRIORITY_NORMAL

class PLDBridge final :
    public rtengine::ProgressListener
{
public:
    explicit PLDBridge(rtengine::ProgressListener* pb) :
        pl(pb)
    {
    }

    // ProgressListener interface
    void setProgress(double p) override
    {
        GThreadLock lock;
        pl->setProgress(p);
    }
    void setProgressStr(const Glib::ustring& str) override
    {
        GThreadLock lock;
        Glib::ustring progrstr;
        progrstr = M(str);
        pl->setProgressStr(progrstr);
    }

    void setProgressState(bool inProcessing) override
    {
        GThreadLock lock;
        pl->setProgressState(inProcessing);
    }

    void error(const Glib::ustring& descr) override
    {
        GThreadLock lock;
        pl->error(descr);
    }

private:
    rtengine::ProgressListener* const pl;
};

template<class T>
class ProgressConnector
{

    sigc::signal0<T> opStart;
    sigc::signal0<bool> opEnd;
    T retval;
    Glib::Thread *workThread;

    static int emitEndSignalUI (void* data)
    {

        const sigc::signal0<bool>* lopEnd = reinterpret_cast<sigc::signal0<bool>*>(data);
        const int r = lopEnd->emit ();
        delete lopEnd;

        return r;
    }

    void workingThread ()
    {
        retval = opStart.emit ();
        gdk_threads_add_idle(ProgressConnector<T>::emitEndSignalUI, new sigc::signal0<bool>(opEnd));
        workThread = nullptr;
    }

public:

    ProgressConnector (): retval( 0 ), workThread( nullptr ) { }

    void startFunc (const sigc::slot0<T>& startHandler, const sigc::slot0<bool>& endHandler )
    {
        if( !workThread ) {
            opStart.connect (startHandler);
            opEnd.connect (endHandler);
            workThread = Glib::Thread::create(sigc::mem_fun(*this, &ProgressConnector<T>::workingThread), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);
        }
    }

    T returnValue()
    {
        return retval;
    }
};
