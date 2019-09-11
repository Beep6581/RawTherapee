# This file is part of RawTherapee.
#
# Copyright (C) 2018 Flössie <floessie.mail@gmail.com>
#
# RawTherapee is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RawTherapee is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_QUIET_COPY "${CMAKE_REQUIRED_QUIET}")
set(CMAKE_REQUIRED_QUIET ON)

set(TEST_SOURCE
"
#include <cstddef>
#include <type_traits>

int main()
{
#if defined(__SSE2__) && (defined(__i386) || defined(_M_IX86)) && defined(__linux)
    static_assert(std::alignment_of<std::max_align_t>::value >= 16, \"Unaligned heap objects possible\");
#endif

    return 0;
}
")

CHECK_CXX_SOURCE_COMPILES("${TEST_SOURCE}" HAVE_ALIGNED_MALLOC)

if(NOT HAVE_ALIGNED_MALLOC)
    set(HAVE_UNALIGNED_MALLOC 1)
else()
    unset(HAVE_ALIGNED_MALLOC)
endif()

unset(TEST_SOURCE)

set(CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET_COPY}")
unset(CMAKE_REQUIRED_QUIET_COPY)
