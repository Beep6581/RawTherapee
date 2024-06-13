/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 Rawtherapee developers
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
#include <array>
#include <iostream>

#include "lensmetadata.h"

namespace rtengine
{

std::unique_ptr<MetadataLensCorrection> MetadataLensCorrectionFinder::findCorrection(const FramesMetaData *meta)
{
    static const std::unordered_set<std::string> makers = {};

    std::string make = Glib::ustring(meta->getMake()).uppercase();

    if (makers.find(make) == makers.end()) {
        return nullptr;
    }

    std::unique_ptr<MetadataLensCorrection> correction;

    return correction;
}

} // namespace rtengine
