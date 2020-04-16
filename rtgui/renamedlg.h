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

#include "guiutils.h"

#define RESPONSE_ALL 100

class CacheImageData;

class RenameDialog final :
    public Gtk::Dialog
{

protected:


    Gtk::Window* p;
    Gtk::Label* oldName;
    Gtk::Entry* newName;
    const CacheImageData* imageData;

public:
    explicit RenameDialog (Gtk::Window* parent);

    void initName (const Glib::ustring& iname, const CacheImageData* cid);
    Glib::ustring getNewName ();

};
