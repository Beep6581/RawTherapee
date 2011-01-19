/*
 * util.h
 *
 *  Created on: Jan 18, 2011
 *      Author: gabor
 */

#ifndef UTIL_H_
#define UTIL_H_

#undef CLIPTO

#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))

namespace rtengine {

template<class T, int ArraySize>
inline float lutInterp (T *array, float f)
{
		int index = CLIPTO(floor(f),0,ArraySize-2);
		float part = ((f)-(float)index)*(T)(array[index+1]-array[index]);
		return (T)array[index]+part;
}

}
#endif /* UTIL_H_ */
