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
 *
 *  2024-2024 Daniel Gao <daniel.gao.work@gmail.com>
 */

#include "aspectratios.h"
#include "options.h"

void fillAspectRatios(std::vector<AspectRatio>& ratios) {
    const auto& options = App::get().options();
    ratios.reserve(ratios.size() + options.sysAspectRatiosLabels.size() + options.userAspectRatiosLabels.size());
    int sysEnabledCount = 0;
    for (size_t i = 0; i < options.sysAspectRatiosLabels.size(); i++) {
        if (options.sysAspectRatiosEnabled[i]) {
            ratios.push_back({options.sysAspectRatiosLabels[i], options.sysAspectRatiosValues[i]});
            sysEnabledCount++;
        }
    }
    int userEnabledCount = 0;
    for (size_t i = 0; i < options.userAspectRatiosLabels.size(); i++) {
        if (options.userAspectRatiosEnabled[i]) {
            ratios.push_back({options.userAspectRatiosLabels[i], options.userAspectRatiosValues[i]});
            userEnabledCount++;
        }
    }
}
