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

/// @brief List of pipette editing operation
enum EditUniqueID : int {
    EUID_None,  /// special value (default)

    EUID_ToneCurve1,
    EUID_ToneCurve2,
    EUID_Lab_LCurve,
    EUID_Lab_CCurve,
    EUID_Lab_LCCurve,
    EUID_Lab_CLCurve,
    EUID_Lab_LHCurve,
    EUID_Lab_CHCurve,
    EUID_Lab_HHCurve,
    EUID_Lab_aCurve,
    EUID_Lab_bCurve,
    EUID_RGB_R,
    EUID_RGB_G,
    EUID_RGB_B,
    EUID_HSV_H,
    EUID_HSV_S,
    EUID_HSV_V,
    EUID_BlackWhiteLuminance,
    EUID_BlackWhiteBeforeCurve,
    EUID_BlackWhiteAfterCurve,
    EUID_WW_HHCurve,

};

/// @brief Editing mechanisms
//  TODO: Look out if it has to be a bitfield to allow both mechanisms at a time
enum EditType {
    ET_PIPETTE,  /// Will trigger dedicated methods; can have a geometry list to be displayed, but without "mouse over" capabilities
    ET_OBJECTS   /// The objects are geometrical widgets with "mouse over" capabilities
};

/// @brief Buffer type for ET_PIPETTE type editing
enum BufferType {
    BT_IMAGEFLOAT,          /// Imagefloat buffer type (3 channels of float values)
    BT_LABIMAGE,            /// LabImage buffer type (3 channels of float values)
    BT_SINGLEPLANE_FLOAT    /// All purpose, 1 channel buffer of float values
};

/// @brief Number of object to be handled (for optimization purpose)
enum ObjectMode {
    OM_255,   /// less or equal than 255 objects
    OM_65535  /// less or equal than 65535 objects
};
