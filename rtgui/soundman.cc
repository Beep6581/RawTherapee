/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2010 Oliver Duis <www.oliverduis.de>
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
*
*/

#include "soundman.h"
#include "options.h"

#ifdef WIN32
#include <windows.h>
#include <mmsystem.h>
#endif

#ifdef __linux__
#include <canberra-gtk.h>
#endif


void SoundManager::init()
{
#ifdef WIN32
// TODO: On Windows Vista/7 RT should register with the OS sound system, so it can enjoy application specific
// volume, safed, process independent etc. from the start.
// Function call is IAudioClient::Initialize
// Unfortunately MinGW does not support this yet. If audioclient.h is available, add an Init
// called once on program start.
    //
    // This mitigation plays an empty file on start, so RT is immidiately avaible in the Windows mixer at least
    playSoundAsync(Glib::ustring("sounds\\Empty.wav"));
#endif
}

// Plays a sound in async mode to not block the main thread
// param is either file name or name of the system event on Windows (e.g. "SystemAsterisk" or "SystemDefault").
void SoundManager::playSoundAsync(const Glib::ustring &sound)
{
    if (sound.empty() || !options.sndEnable) {
        return;
    }

#ifdef WIN32
    DWORD sndParam = SND_ASYNC | SND_NODEFAULT;

    if (sound.find('.') != Glib::ustring::npos) {
        // contain dot, so it's a filename
        sndParam |= SND_FILENAME;
    } else {
        // no dot, so it's a system event
        sndParam |= SND_ALIAS;
    }

    wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (sound.c_str(), -1, NULL, NULL, NULL);
    PlaySoundW(wfilename, NULL, sndParam);
    g_free( wfilename );
#elif defined(__linux__)
    ca_context_play(ca_gtk_context_get(), 0, CA_PROP_EVENT_ID, sound.c_str(), CA_PROP_MEDIA_FILENAME, sound.c_str(), NULL);
#endif
}
