/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <www.oliverduis.de>
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
#include <cstring>

#include "extprog.h"
#include "multilangmgr.h"
#include "../rtengine/safegtk.h"
#ifdef WIN32
#include <windows.h>
// for GCC32
#ifndef _WIN32_IE
#define _WIN32_IE 0x0600
#endif
#include <shlobj.h>
#endif
using namespace std;


ExtProgAction::ExtProgAction() {}

ExtProgAction::ExtProgAction(const ExtProgAction* other, int target)
    : target(target), filePathEXE(other->filePathEXE), preparams(other->preparams), name(other->name) { }

Glib::ustring ExtProgAction::GetFullName()
{
    return name + " [" + M(Glib::ustring::compose("EXTPROGTARGET_%1", target)) + "]";
}

bool ExtProgAction::Execute(std::vector<Glib::ustring> fileNames)
{
    if (fileNames.empty()) {
        return false;
    }

    // Check if they all exists (maybe not precessed yet)
    for (int i = 0; i < fileNames.size(); i++) {
        if (!safe_file_test(fileNames[i], Glib::FILE_TEST_EXISTS)) {
            Gtk::MessageDialog msgd (M("MAIN_MSG_IMAGEUNPROCESSED"), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
            msgd.run ();
            return false;
        }
    }

    Glib::ustring cmdLine = "\"" + filePathEXE + "\"";

    if (preparams.length() > 0) {
        cmdLine += " " + preparams;
    }

    for (int i = 0; i < fileNames.size(); i++) {
        cmdLine += " \"" + fileNames[i] + "\"";
    }

    return safe_spawn_command_line_async (cmdLine);
}


// Generates as singleton
ExtProgStore* ExtProgStore::getInstance()
{
    static ExtProgStore instance_;
    return &instance_;
}

ExtProgStore::~ExtProgStore()
{
    for (list<ExtProgAction*>::iterator it = lActions.begin(); it != lActions.end(); it++) {
        delete *it;
    }
}

// Reads all profiles from the given profiles dir
void ExtProgStore::init ()
{
    MyMutex::MyLock lock(mtx);

    lActions.clear();

#ifdef WIN32

    SearchProg("Photoshop", "Adobe\\Adobe Photoshop CS%1 (64 Bit)\\Photoshop.exe", "Adobe\\Adobe Photoshop CS%1\\Photoshop.exe", 9, false, true);
    SearchProg("Photomatix Pro", "PhotomatixPro%1\\PhotomatixPro.exe", "", 9, true, true);
    SearchProg("Paint.NET", "Paint.NET\\PaintDotNet.exe", "", 0, false, true);
    SearchProg("MS Image Composition Editor", "Microsoft Research\\Image Composite Editor\\ICE.exe", "", 0, false, true);
    SearchProg("PTGui", "PTGui\\PTGui.exe", "", 0, false, true);
    SearchProg("GeoSetter", "GeoSetter\\GeoSetter.exe", "", 0, true, true);
    SearchProg("FastStone Image Viewer", "FastStone Image Viewer\\FSViewer.exe", "", 0, true, true);
    SearchProg("FastPictureViewer", "FastPictureViewer\\FastPictureViewer.exe", "", 0, true, true);

    if (!SearchProg("Autopano Giga 3", "Kolor\\Autopano Giga 3.%1\\AutopanoGiga_x64.exe", "Kolor\\Autopano Giga 3.%1\\AutopanoGiga.exe", 15, true, true)) {
        if ( !SearchProg("Autopano Pro 3", "Kolor\\Autopano Pro 3.%1\\AutopanoPro_x64.exe", "Kolor\\Autopano Pro 3.%1\\AutopanoPro.exe", 15, true, true))   {
            if (!SearchProg("Autopano Giga 2", "Kolor\\Autopano Giga 2.%1\\AutopanoGiga_x64.exe", "Kolor\\Autopano Giga 2.%1\\AutopanoGiga.exe", 6, true, true)) {
                SearchProg("Autopano Pro 2", "Kolor\\Autopano Pro 2.%1\\AutopanoPro_x64.exe", "Kolor\\Autopano Pro 2.%1\\AutopanoPro.exe", 6, true, true);
            }
        }
    }

    // DO NOT add obscure little tools here, only widely used programs with proper setup program to have a standard path
#endif

}

bool ExtProgStore::SearchProg(Glib::ustring name, Glib::ustring exePath, Glib::ustring exePath86, int maxVer, bool allowRaw, bool allowQueueProcess)
{
    bool found = false;

#ifdef WIN32
    // get_user_special_dir crashes on some Windows configurations.
    // so we use the safe native functions here
    static Glib::ustring progFilesDir, progFilesDirx86;

    if (progFilesDir.empty()) {
        WCHAR pathW[MAX_PATH] = {0};
        char pathA[MAX_PATH];

        // First prio folder (64bit, otherwise 32bit)
        if (SHGetSpecialFolderPathW(NULL, pathW, CSIDL_PROGRAM_FILES, false)) {
            char pathA[MAX_PATH];
            WideCharToMultiByte(CP_UTF8, 0, pathW, -1, pathA, MAX_PATH, 0, 0);
            progFilesDir = Glib::ustring(pathA);
        }

        if (SHGetSpecialFolderPathW(NULL, pathW, CSIDL_PROGRAM_FILESX86, false)) {
            WideCharToMultiByte(CP_UTF8, 0, pathW, -1, pathA, MAX_PATH, 0, 0);
            progFilesDirx86 = Glib::ustring(pathA);
        }
    }

    if (exePath86.empty()) {
        exePath86 = exePath;
    }

    ExtProgAction *pAct = new ExtProgAction();
    pAct->name = name;
    pAct->target = (allowRaw ? 1 : 2);

    if (maxVer > 0) {
        for (int verNo = maxVer; verNo >= 0; verNo--) {
            pAct->filePathEXE = progFilesDir + "\\" + Glib::ustring::compose(exePath, verNo);

            if (safe_file_test(pAct->filePathEXE, Glib::FILE_TEST_EXISTS)) {
                break;
            }

            pAct->filePathEXE = progFilesDirx86 + "\\" + Glib::ustring::compose(exePath86, verNo);

            if (safe_file_test(pAct->filePathEXE, Glib::FILE_TEST_EXISTS)) {
                break;
            }

            pAct->filePathEXE = "";
        }
    } else {
        pAct->filePathEXE = progFilesDir + "\\" + exePath;

        if (!safe_file_test(pAct->filePathEXE, Glib::FILE_TEST_EXISTS)) {

            pAct->filePathEXE = progFilesDirx86 + "\\" + exePath86;

            if (!safe_file_test(pAct->filePathEXE, Glib::FILE_TEST_EXISTS)) {
                pAct->filePathEXE = "";
            }
        }
    }

    if (pAct->filePathEXE.length() > 0) {
        lActions.push_back(pAct);

        // Copy for second target
        if (allowRaw && allowQueueProcess) {
            lActions.push_back(new ExtProgAction(pAct, 2));
        }

        found = true;
    } else {
        delete pAct;
    }

#endif

    return found;
}
