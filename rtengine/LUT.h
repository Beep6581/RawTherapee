/*
 * LUT.h
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2011 Jan Rinze Peterzon (janrinze@gmail.com)
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

/*
 *  Declaration of flexible 2D arrays
 *
 *  Usage:
 *
 *  	LUT<type> name (size);
 *		LUT<type> name (size, flags);
 *
 *		creates an array which is valid within the normal C/C++ scope "{ ... }"
 *
 *      access to elements is a simple as:
 *
 *      	LUT<float> my_lut (10);
 *      	float value = my_lut[3];
 *          float value = my_lut[2.5]; // this will interpolate
 *
 *      when using a float type index it will interpolate the lookup values
 *
 *      extra setting in flags: (clipping is set by default)
 *      LUT_CLIP_ABOVE
 *      LUT_CLIP_BELOW
 *
 *      example:
 *      	LUT<float> my_lut (10,LUT_CLIP_BELOW);
 *          float value = my_lut[22.5];  // this will extrapolate
 *          float value = my_lut[-22.5]; // this will not extrapolate
 *
 *          LUT<float> my_lut (10,0); // this will extrapolate on either side
 *
 *
 */

#ifndef LUT_H_
#define LUT_H_

#define LUT_CLIP_BELOW 1
#define LUT_CLIP_ABOVE 	2

template <typename T>
class LUT
{
private:
  int size,owner;
  unsigned int clip;
  T * data;
public:
  LUT(int s,int flags=0xfffffff)
  {
	  clip=flags;
	  data = new T [s];
      owner=1;
	  size=s;
  }

  LUT(int s,T * source)
  {
	  data = new T [s];
      owner=1;
	  size=s;
	  for (int i=0;i<s;i++)
	  {
		  data[i]=source[i];
	  }
  }

  ~LUT()
  {
	  if (owner) delete [] data;
  }

  // use with integer indices
  T& operator[](int index)
  {
	  if(index<0) return data[0];
	  if(index>=size) return data[size-1];
	  return data[index];
  }
  // use with float indices
  T operator[](float index)
  {
	  int idx = floor(index);
	  if(idx<=0)
	  {
		  idx=0;
		  if (clip&LUT_CLIP_BOTTOM)
			  return data[0];
	  }
	  if(idx>=size-1)
	  {
		  idx=size-2;
		  if (clip&LUT_CLIP_ABOVE)
			  return data[size-1];
	  }
	  float diff=index - (float)idx;
	  T p1 = data[idx]* (1.0-diff);
	  T p2 = data[idx+1]*diff;
	  return p1+p2;
  }
};

#endif /* LUT_H_ */
