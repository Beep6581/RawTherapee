/*
 *  This file is part of RawTherapee.
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

#ifdef RT_SIMDE
#include <simde/simde-common.h>
#endif

// Check if the current version of SIMDe is at least major.minor.patch.
#if defined(SIMDE_VERSION_MAJOR) && defined(SIMDE_VERSION_MINOR) && defined(SIMDE_VERSION_MICRO)
#define SIMDE_VERSION_CHECK(major, minor, patch) (                                                  \
    (SIMDE_VERSION_MAJOR > (major)) ||                                                             \
        (SIMDE_VERSION_MAJOR == (major) && SIMDE_VERSION_MINOR > (minor)) ||                        \
        (SIMDE_VERSION_MAJOR == (major) && SIMDE_VERSION_MINOR == (minor) && SIMDE_VERSION_MICRO >= (patch)))
#else
#define SIMDE_VERSION_CHECK(major, minor, patch) (0)
#endif

