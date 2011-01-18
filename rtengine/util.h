/*
 * util.h
 *
 *  Created on: Jan 18, 2011
 *      Author: gabor
 */

#ifndef UTIL_H_
#define UTIL_H_

#undef CLIPI
#undef CLIPTO

#define CLIPI(a) ((a)>0?((a)<65534?(a):65534):0)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))

namespace rtengine {

inline float lutInterp (int *array, float f)
{
		int index = CLIPI(floor(f));
		float part = (float)((f)-index)*(float)(array[index+1]-array[index]);
		return (float)array[index]+part;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// accurately determine value from float array with float as index
//linearly interpolate from ends of range if arg is out of bounds
inline float lutInterp (float *array, float f)
{
		int index = CLIPI(floor(f));
		float part = ((f)-(float)index)*(array[index+1]-array[index]);
		return array[index]+part;
}

inline float lutInterp (float *array, float f, int arraySize)
{
		int index = CLIPTO(floor(f),0,arraySize-1);
		float part = ((f)-(float)index)*(array[index+1]-array[index]);
		return array[index]+part;
}

}
#endif /* UTIL_H_ */
