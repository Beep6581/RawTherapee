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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _EDITID_H_
#define _EDITID_H_


/// @brief List of pipette editing operation
enum EditUniqueID {
	EUID_None,  /// special value (default)

	EUID_ToneCurve1,
	EUID_ToneCurve2,
	EUID_Lab_LCurve,
	EUID_RGB_R,
	EUID_RGB_G,
	EUID_RGB_B,
	EUID_HSV_H,
	EUID_HSV_S,
	EUID_HSV_V,
	EUID_BlackWhiteLuminance,
	EUID_BlackWhiteBeforeCurve,
	EUID_BlackWhiteAfterCurve,
};

/// @brief Editing mechanisms
//  TODO: Look out if it has to be a bitfield to allow both mechanisms at a time
enum EditType {
	ET_PIPETTE,  /// Will trigger dedicated methods; can have a geometry list to be displayed, but without "mouse over" capabilities
	ET_OBJECTS   /// The objects are geometrical widgets with "mouse over" capabilities
};

/// @brief Buffer type for ET_PIPETTE type editing
enum BufferType {
	BT_IMAGEFLOAT,
	BT_LABIMAGE,
	BT_SINGLEPLANE_FLOAT
};

/// @brief Number of object to be handled (for optimization purpose)
enum ObjectMode {
	OM_255,   /// less or equal than 255 objects
	OM_65535  /// less or equal than 65535 objects
};

#endif
