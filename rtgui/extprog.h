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

#ifndef _EXTPROG_
#define _EXTPROG_

#include <glibmm/ustring.h>

#include <vector>

#include "threadutils.h"

struct ExtProgAction
{
    Glib::ustring filePathEXE;
    Glib::ustring preparams;  // after EXE and before file names
    Glib::ustring name;  // already localized if necessary
    int target;  // 1=RAW files, 2=batch converted files

    Glib::ustring getFullName () const;  // e.g. "Photoshop (RAW)"

    bool execute (const std::vector<Glib::ustring>& fileNames) const;
};

// Stores all external programs that could be called by the user
class ExtProgStore
{
    MyMutex mtx;  // covers actions
    std::vector<ExtProgAction> actions;

    bool searchProgram (const Glib::ustring& name,
                        const Glib::ustring& exePath,
                        const Glib::ustring& exePath86,
                        int maxVer,
                        bool allowRaw,
                        bool allowQueueProcess);

public:
    static ExtProgStore* getInstance();

    // searches computer for installed standard programs
    void init();

    const std::vector<ExtProgAction>& getActions () const;

    static bool spawnCommandAsync (const Glib::ustring& cmd);
    static bool spawnCommandSync (const Glib::ustring& cmd);

    static bool openInGimp (const Glib::ustring& fileName);
    static bool openInPhotoshop (const Glib::ustring& fileName);
    static bool openInCustomEditor (const Glib::ustring& fileName);
};

#define extProgStore ExtProgStore::getInstance()

inline const std::vector<ExtProgAction>& ExtProgStore::getActions () const
{
    return actions;
}

#endif
