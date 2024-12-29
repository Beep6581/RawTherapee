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

namespace {

static const std::vector<AspectRatio> ASPECT_RATIOS {
    {"3:2", 3.0 / 2.0},                 // L1.5,        P0.666...
    {"4:3", 4.0 / 3.0},                 // L1.333...,   P0.75
    {"16:9", 16.0 / 9.0},               // L1.777...,   P0.5625
    {"16:10", 16.0 / 10.0},             // L1.6,        P0.625
    {"1:1", 1.0 / 1.0},                 // L1,          P1
    {"2:1", 2.0 / 1.0},                 // L2,          P0.5
    {"3:1", 3.0 / 1.0},                 // L3,          P0.333...
    {"4:1", 4.0 / 1.0},                 // L4,          P0.25
    {"5:1", 5.0 / 1.0},                 // L5,          P0.2
    {"6:1", 6.0 / 1.0},                 // L6,          P0.1666...
    {"7:1", 7.0 / 1.0},                 // L7,          P0.142...
    {"4:5", 4.0 / 5.0},                 // L1.25,       P0.8
    {"5:7", 5.0 / 7.0},                 // L1.4,        P0.714...
    {"6:7", 6.0 / 7.0},                 // L1.166...,   P0.857...
    {"6:17", 6.0 / 17.0},               // L2.833...,   P0.352...
    {"24:65 - XPAN", 24.0 / 65.0},      // L2.708...,   P0.369...
    {"1.414 - DIN EN ISO 216", 1.414},  // L1.414,      P0.707...
    {"3.5:5", 3.5 / 5.0},               // L1.428...,   P0.7
    {"8.5:11 - US Letter", 8.5 / 11.0}, // L1.294...,   P0.772...
    {"9.5:12", 9.5 / 12.0},             // L1.263...,   P0.791...
    {"10:12", 10.0 / 12.0},             // L1.2,        P0.833...
    {"11:14", 11.0 / 14.0},             // L1.272...,   P0.785...
    {"11:17 - Tabloid", 11.0 / 17.0},   // L1.545...,   P0.647...
    {"13:19", 13.0 / 19.0},             // L1.461...,   P0.684...
    {"17:22", 17.0 / 22.0},             // L1.294...,   P0.772...
    {"45:35 - ePassport", 45.0 / 35.0}, // L1.285,...   P0.777...
    {"64:27", 64.0 / 27.0},             // L2.370...,   P0.421...
    {"13:18", 13.0 / 18.0},             // L1.384...,   P0.722...
};

}  // namespace

void fillAspectRatios(std::vector<AspectRatio>& ratios) {
    ratios.reserve(ratios.size() + ASPECT_RATIOS.size());

    for (const auto& ratio : ASPECT_RATIOS) {
        ratios.push_back(ratio);
    }
}
