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

#include <memory>
#include <string>
#include <vector>

#include <glibmm/ustring.h>

namespace rtengine
{

struct badPix;

class RawImage;

class DFManager final
{
public:
    static DFManager& getInstance();

    void init(const Glib::ustring& pathname);
    Glib::ustring getPathname() const;
    void getStat(int& totFiles, int& totTemplates) const;
    const RawImage* searchDarkFrame(const std::string& mak, const std::string& mod, int iso, double shut, time_t t);
    const RawImage* searchDarkFrame(const Glib::ustring& filename);
    const std::vector<badPix>* getHotPixels(const std::string& mak, const std::string& mod, int iso, double shut, time_t t);
    const std::vector<badPix>* getHotPixels(const Glib::ustring& filename);
    const std::vector<badPix>* getBadPixels(const std::string& mak, const std::string& mod, const std::string& serial) const;

private:
    DFManager();
    ~DFManager();

    class Implementation;

    const std::unique_ptr<Implementation> implementation;
};

}
