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
#include <options.h>

static void CALLBACK current_directory_monitor_callback (DWORD error, DWORD nBytes, LPOVERLAPPED lpOverlapped) {

    DWORD dwOffset = 0;
    FILE_NOTIFY_INFORMATION* pInfo = NULL;

    WinDirMonitor::MonitorData* monData = (WinDirMonitor::MonitorData*)lpOverlapped;
    if (!nBytes) {
        delete monData;
        return;
    }

    bool notify = false;
    // Analysis of the modifications
    do {
        Glib::ustring fname = "";
        Glib::ustring action = "";
        int strLen = 0;

        // Get a pointer to the first change record...
        pInfo = (FILE_NOTIFY_INFORMATION*) &monData->file_notify_buffer[dwOffset];

        char fnameC[(MAX_PATH+1)*2] = {0};
        strLen = WideCharToMultiByte(CP_UTF8,0,pInfo->FileName,pInfo->FileNameLength/sizeof(WCHAR),fnameC,sizeof(fnameC),0,0);
        fnameC[strLen] = 0;
        fname = fnameC;

        if (options.has_retained_extention(fname))
        {
        	switch (pInfo->Action)
        	{
			case (FILE_ACTION_ADDED):
				action = "FILE_ACTION_ADDED";
				break;
			case (FILE_ACTION_REMOVED):
				action = "FILE_ACTION_REMOVED";
				break;
			case (FILE_ACTION_MODIFIED):
				action = "FILE_ACTION_MODIFIED";
				break;
			case (FILE_ACTION_RENAMED_OLD_NAME):
				action = "FILE_ACTION_RENAMED_OLD_NAME";
				break;
			case (FILE_ACTION_RENAMED_NEW_NAME):
				action = "FILE_ACTION_RENAMED_NEW_NAME";
				break;
			case (FILE_ACTION_REMOVED_BY_DELETE):
				action = "FILE_ACTION_REMOVED_BY_DELETE";
			/*	break;
			case (FILE_ACTION_ADDED_STREAM):
				action = "FILE_ACTION_ADDED_STREAM";
				break;
			case (FILE_ACTION_REMOVED_STREAM):
				action = "FILE_ACTION_REMOVED_STREAM";
				break;
			case (FILE_ACTION_MODIFIED_STREAM):
				action = "FILE_ACTION_MODIFIED_STREAM";
				break;
			case (FILE_ACTION_ID_NOT_TUNNELLED):
				action = "FILE_ACTION_ID_NOT_TUNNELLED";
				break;
			case (FILE_ACTION_TUNNELLED_ID_COLLISION):
				action = "FILE_ACTION_TUNNELLED_ID_COLLISION";*/
			default:
				break;
			}
			notify = true;
        }

        // More than one change may happen at the same time.  Load the next change and continue...
        dwOffset += pInfo->NextEntryOffset;
    }
    while (pInfo->NextEntryOffset != 0);

	// ReadDirectoryChangesW sometimes emits multiple events per change (one for each change type)
	// To make sure it's not flooding update, this gets filtered.
	time_t curTime= ::time(NULL);
    if (notify && monData->listener && ::difftime(curTime, monData->lastTimeUpdateDir)>1.0) {
    	printf("----- Appel de WinDirChanged -----\n");
        monData->listener->winDirChanged ();
		monData->lastTimeUpdateDir = curTime;
	}
	
    ReadDirectoryChangesW (monData->hDirectory,
			 monData->file_notify_buffer,
			 monData->buffer_allocated_bytes,
			 FALSE, 
			 FILE_NOTIFY_CHANGE_FILE_NAME |
			 FILE_NOTIFY_CHANGE_DIR_NAME |
			 FILE_NOTIFY_CHANGE_LAST_WRITE,
			 &monData->buffer_filled_bytes,
			 &monData->overlapped,
			 current_directory_monitor_callback);
}

WinDirMonitor::WinDirMonitor (Glib::ustring dirName, WinDirChangeListener* listener) : monData(NULL) {
    wchar_t* wdirname = (wchar_t*)g_utf8_to_utf16 (dirName.c_str(), -1, NULL, NULL, NULL);
    HANDLE hDirectory = CreateFileW (wdirname, FILE_LIST_DIRECTORY,FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS | FILE_FLAG_OVERLAPPED, NULL); 
    g_free (wdirname);
  
    if (hDirectory != INVALID_HANDLE_VALUE) {

        monData = new MonitorData ();
        monData->listener = listener;
        monData->buffer_allocated_bytes = 32768;
        monData->file_notify_buffer = new char [monData->buffer_allocated_bytes];
        monData->hDirectory = hDirectory;
		monData->lastTimeUpdateDir = ::time(NULL);
        
//        printf ("mondata=%d\n", monData);
        ReadDirectoryChangesW (monData->hDirectory,
				  monData->file_notify_buffer,
				  monData->buffer_allocated_bytes,
				  FALSE, 
				  FILE_NOTIFY_CHANGE_FILE_NAME |
				  FILE_NOTIFY_CHANGE_DIR_NAME |
				  FILE_NOTIFY_CHANGE_LAST_WRITE,
				  &monData->buffer_filled_bytes,
				  &monData->overlapped,
				  current_directory_monitor_callback);
    }
}

WinDirMonitor::~WinDirMonitor () {

    if (monData && monData->hDirectory != INVALID_HANDLE_VALUE)
        CloseHandle (monData->hDirectory);
}
