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
*  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#include <vector>

#include <glibmm/refptr.h>
#include <glibmm/ustring.h>

#include "threadutils.h"

namespace Gio
{
    class AppInfo;
}

struct ExtProgAction
{
    Glib::ustring filePathEXE;
    Glib::ustring preparams;  // after EXE and before file names
    Glib::ustring name;  // already localized if necessary
    int target;  // 1=RAW files, 2=batch converted files

    Glib::ustring getFullName () const;  // e.g. "Photoshop (RAW)"

    bool execute (const std::vector<Glib::ustring>& fileNames) const;
};

struct EditorInfo
{
    Glib::ustring name;
    Glib::ustring commandline;
    bool isNativeCommand;
};

// Stores all external programs that could be called by the user
class ExtProgStore
{
    MyMutex mtx;  // covers actions
    std::vector<ExtProgAction> actions;

#ifdef _WIN32
    bool searchProgram (const Glib::ustring& name,
                        const Glib::ustring& exePath,
                        const Glib::ustring& exePath86,
                        int maxVer,
                        bool allowRaw,
                        bool allowQueueProcess);
#endif
public:
    static ExtProgStore* getInstance();

    // searches computer for installed standard programs
    void init();

    const std::vector<ExtProgAction>& getActions () const;

    static bool spawnCommandAsync (const Glib::ustring& cmd);
    static bool spawnCommandSync (const Glib::ustring& cmd);

    static bool openInGimp (const Glib::ustring& fileName);
    static bool openInPhotoshop (const Glib::ustring& fileName);
    static bool openInCustomEditor (const Glib::ustring& fileName, const Glib::ustring* command = nullptr);
    static bool openInExternalEditor(const Glib::ustring &fileName, const EditorInfo &editorInfo);
};

#define extProgStore ExtProgStore::getInstance()

inline const std::vector<ExtProgAction>& ExtProgStore::getActions () const
{
    return actions;
}
