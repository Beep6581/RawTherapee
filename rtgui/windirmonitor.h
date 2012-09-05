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
#ifndef _WINDIRMONITOR_
#define _WINDIRMONITOR_

#include <glibmm.h>
#include <windows.h>

class WinDirChangeListener {

    public:
        virtual void winDirChanged () {}
};

class WinDirMonitor : public Glib::Object {

    public:
        struct MonitorData {
    		OVERLAPPED overlapped;
    		DWORD buffer_allocated_bytes;
    		char *file_notify_buffer;
    		DWORD buffer_filled_bytes;    
            HANDLE hDirectory;
            WinDirChangeListener* listener;
            int bigyo;
            DWORD lastTimeUpdateTick;  // for filtering multiple updates events
        };
    
    private:
        MonitorData* monData;
    
    public:
        WinDirMonitor (Glib::ustring dirName, WinDirChangeListener* listener);
        ~WinDirMonitor ();
};

#endif

