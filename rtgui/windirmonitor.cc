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
#include <windirmonitor.h>

static void CALLBACK current_directory_monitor_callback (DWORD error, DWORD nBytes, LPOVERLAPPED lpOverlapped) {

//    printf ("hivas: %d %d %d\n", error, nBytes, lpOverlapped);
  
    WinDirMonitor::MonitorData* monData = (WinDirMonitor::MonitorData*)lpOverlapped;
    if (!nBytes) {
        delete monData;
        return;
    }

    if (monData->listener)
        monData->listener->winDirChanged ();
	
    ReadDirectoryChangesW (monData->hDirectory,
			 monData->file_notify_buffer,
			 monData->buffer_allocated_bytes,
			 FALSE, 
			 FILE_NOTIFY_CHANGE_FILE_NAME |
			 FILE_NOTIFY_CHANGE_DIR_NAME |
			 FILE_NOTIFY_CHANGE_ATTRIBUTES |
			 FILE_NOTIFY_CHANGE_SIZE,
			 &monData->buffer_filled_bytes,
			 &monData->overlapped,
			 current_directory_monitor_callback);
}

WinDirMonitor::WinDirMonitor (Glib::ustring dirName, WinDirChangeListener* listener) : monData(NULL) {

    wchar_t* wdirname = (wchar_t*)g_utf8_to_utf16 (dirName.raw().c_str(), -1, NULL, NULL, NULL);
    HANDLE hDirectory = CreateFileW (wdirname, FILE_LIST_DIRECTORY,FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OVERLAPPED, NULL); 
    g_free (wdirname);
  
    if (hDirectory != INVALID_HANDLE_VALUE) {

        monData = new MonitorData ();
        monData->listener = listener;
        monData->buffer_allocated_bytes = 32768;
        monData->file_notify_buffer = new char [monData->buffer_allocated_bytes];
        monData->hDirectory = hDirectory;
        
//        printf ("mondata=%d\n", monData);
        ReadDirectoryChangesW (monData->hDirectory,
				  monData->file_notify_buffer,
				  monData->buffer_allocated_bytes,
				  FALSE, 
				  FILE_NOTIFY_CHANGE_FILE_NAME |
				  FILE_NOTIFY_CHANGE_DIR_NAME |
				  FILE_NOTIFY_CHANGE_ATTRIBUTES |
				  FILE_NOTIFY_CHANGE_SIZE,
				  &monData->buffer_filled_bytes,
				  &monData->overlapped,
				  current_directory_monitor_callback);
    }
}

WinDirMonitor::~WinDirMonitor () {

    if (monData && monData->hDirectory != INVALID_HANDLE_VALUE)
        CloseHandle (monData->hDirectory);
}
