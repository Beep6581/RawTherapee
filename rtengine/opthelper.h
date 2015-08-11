////////////////////////////////////////////////////////////////
//
//  opthelper.h includes some #defines which help to make optimizations easier and better readable
//
//  copyright (c) 2013  Ingo Weyrich
//
//  this is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#ifndef OPTHELPER_H
#define OPTHELPER_H

#ifdef __SSE2__
#include "sleefsseavx.c"
#ifdef __GNUC__
#if defined(WIN32) && !defined( __x86_64__ )
// needed for actual versions of GCC with 32-Bit Windows
#define SSEFUNCTION __attribute__((force_align_arg_pointer))
#else
#define SSEFUNCTION
#endif
#else
#define SSEFUNCTION
#endif
#else
#ifdef __SSE__
#ifdef __GNUC__
#if defined(WIN32) && !defined( __x86_64__ )
// needed for actual versions of GCC with 32-Bit Windows
#define SSEFUNCTION __attribute__((force_align_arg_pointer))
#else
#define SSEFUNCTION
#endif
#else
#define SSEFUNCTION
#endif
#else
#define SSEFUNCTION
#endif
#endif

#ifdef __GNUC__
#define RESTRICT    __restrict__
#define LIKELY(x)   __builtin_expect (!!(x), 1)
#define UNLIKELY(x) __builtin_expect (!!(x), 0)
#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 9) || __GNUC__ > 4
#define ALIGNED64 __attribute__ ((aligned (64)))
#define ALIGNED16 __attribute__ ((aligned (16)))
#else // there is a bug in gcc 4.7.x when using openmp and aligned memory and -O3
#define ALIGNED64
#define ALIGNED16
#endif
#else
#define RESTRICT
#define LIKELY(x)    (x)
#define UNLIKELY(x)  (x)
#define ALIGNED64
#define ALIGNED16
#endif
#ifndef __clang__
#define _RT_NESTED_OPENMP _OPENMP
#endif
#endif
