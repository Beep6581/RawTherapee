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

#include <glibmm.h>
#include <list>
#include "threadutils.h"

class ExtProgAction {
public:
    ExtProgAction();
    ExtProgAction(const ExtProgAction* other, int target);

    Glib::ustring filePathEXE;
    Glib::ustring preparams;  // after EXE and before file names
    Glib::ustring name;  // already localized if necessary
    int target;  // 1=RAW files, 2=batch converted files

    Glib::ustring GetFullName();  // e.g. "Photoshop (RAW)"

    virtual bool Execute(std::vector<Glib::ustring> fileNames);
};

// Stores all external programs that could be called by the user
class ExtProgStore {
    MyMutex mtx;  // covers actions

    bool SearchProg(Glib::ustring name, Glib::ustring exePath, Glib::ustring exePath86, int maxVer, bool allowRaw, bool allowQueueProcess);

public:
    ~ExtProgStore();

    void init();  // searches computer for installed standard programs
    static ExtProgStore* getInstance();

    std::list<ExtProgAction*> lActions;
};

#define extProgStore ExtProgStore::getInstance()

#endif
