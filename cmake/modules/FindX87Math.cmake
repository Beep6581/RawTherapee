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
# along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_QUIET_COPY "${CMAKE_REQUIRED_QUIET}")
set(CMAKE_REQUIRED_QUIET ON)

set(TEST_SOURCE
"
#if !defined(__i386) && !defined(_M_IX86)
#error
#endif

#if defined(__SSE2__)
#error
#endif

int main()
{
}
")

CHECK_CXX_SOURCE_COMPILES("${TEST_SOURCE}" HAVE_X87_MATH)

set(TEST_SOURCE
"
#if !defined(__i386) && !defined(_M_IX86)
#error
#endif

#if !defined(__SSE2__)
#error
#endif

int main()
{
}
")

CHECK_CXX_SOURCE_COMPILES("${TEST_SOURCE}" HAVE_X86_SSE_MATH)

unset(TEST_SOURCE)

set(CMAKE_REQUIRED_QUIET "${CMAKE_REQUIRED_QUIET_COPY}")
unset(CMAKE_REQUIRED_QUIET_COPY)
