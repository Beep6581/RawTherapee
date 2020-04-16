/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2019 Gabor Horvath <hgabor@rawtherapee.com>
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

// For compatibility and simplicity reason, order shouldn't change, and must be identical to the order specified in the curveType widget
enum DiagonalCurveType {
    DCT_Empty = -1,     // Also used for identity curves
    DCT_Linear,         // 0
    DCT_Spline,         // 1
    DCT_Parametric,     // 2
    DCT_NURBS,          // 3
    DCT_CatumullRom,    // 4
    // Insert new curve type above this line
    DCT_Unchanged       // Must remain the last of the enum
};
