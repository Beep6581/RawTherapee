/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (C) 2025 Daniel Gao <daniel.gao.work@gmail.com>
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

#ifdef _WIN32

// Smaller Windows headers with less pollution.
// If some needed API is missing, include the specific header instead.
#define WIN32_LEAN_AND_MEAN

// Minimum supported version is Windows 7 if not already set by MinGW.
#ifndef WINVER
#define WINVER 0x0601
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif

// Include before Windows.h to fix compiler warning on MinGW
#include <winsock2.h>

#include <windows.h>

#endif // _WIN32
