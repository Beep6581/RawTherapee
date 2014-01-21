////////////////////////////////////////////////////////////////
//
//  opthelper.h includes some #defines which help to make optimizations easier and better readable
// 
//	copyright (c) 2013  Ingo Weyrich
//
//	this is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////

#ifndef OPTHELPER_H
	#define	OPTHELPER_H

	#ifdef __SSE2__
		#include "sleefsseavx.c"
		#ifdef __GNUC__
			#ifdef WIN32
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
				#ifdef WIN32
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
		#define RESTRICT 	__restrict__
		#define LIKELY(x)   __builtin_expect (!!(x), 1)
		#define UNLIKELY(x) __builtin_expect (!!(x), 0)
	#else
		#define RESTRICT
		#define LIKELY(x)    (x)
		#define UNLIKELY(x)  (x)
	#endif
#endif
