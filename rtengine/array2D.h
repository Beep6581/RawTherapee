/*
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
 *  	array2D<type> name (Y-size,X-size);
 *		array2D<type> name (Y-size,X-size, type ** data);
 *
 *		creates an array which is valid within the normal C/C++ scope "{ ... }"
 *
 *      access to elements is a simple as:
 *
 *      	array2D<float> my_array (10,10); // creates 10x10 array of floats
 *      	value =  my_array[3][5];
 *      	my_array[4][6]=value;
 *
 *      or copy an existing 2D array
 *
 *      	float ** mydata;
 *      	array2D<float> my_array (10,10,mydata);
 *
 *
 *		Useful extra pointers
 *
 *			<type> ** my_array		gives access to the pointer for access with [][]
 *			<type> *  my_array		gives access to the flat stored data.
 *
 *		Advanced usage:
 *			array2D<float> my_array				; // empty container.
 *			my_array(10,10) 					; // resize to 10x10 array
 *			my_array(10,10,ARRAY2D_CLEAR_DATA)  ; // resize to 10x10 and clear data
 *			my_array(10,10,ARRAY2D_CLEAR_DATA|ARRAY2D_LOCK_DATA)  ; same but set a lock on changes
 *
 *			!! locked arrays cannot be resized and cannot be unlocked again !!
 */
#ifndef ARRAY2D_H_
#define ARRAY2D_H_
#include <unistd.h>  // for sleep()
// flags for use
#define ARRAY2D_LOCK_DATA	1
#define ARRAY2D_CLEAR_DATA	2

template<typename T>
class array2D {
private:
	int x, y, owner;
	T ** ptr;
	T * data;
	bool lock; // useful lock to ensure data is not changed anymore.
public:
	array2D(int h, int w, unsigned int flags = 0) {
		lock = flags & ARRAY2D_LOCK_DATA;
		data = new T[h * w];
		owner = 1;
		x = w;
		y = h;
		ptr = new T*[h];
		for (int i = 0; i < h; i++)
			ptr[i] = data + i * w;
		if (flags & ARRAY2D_CLEAR_DATA)
			memset(data, 0, w * h * sizeof(T));
	}

	array2D(int h, int w, T ** source, unsigned int flags = 0) {
		lock = flags & ARRAY2D_LOCK_DATA;
		data = new T[h * w];
		owner = 1;
		x = w;
		y = h;
		ptr = new T*[h];
		for (int i = 0; i < h; i++) {
			ptr[i] = data + i * w;
			for (int j = 0; j < w; j++)
				ptr[i][j] = source[i][j];
		}
	}

	~array2D() {
		if ((owner) && (data))
			delete[] data;
		if (ptr)
			delete[] ptr;
	}

	// use with indices
	T * operator[](size_t index) {
		return ptr[index];
	}

	// use as pointer to T**
	operator T**() {
		return ptr;
	}

	// use as pointer to data
	operator T*() {
		return data;
	}

	// use as empty declaration, resize before use!
	// very useful as a member object
	array2D() :
		x(0), y(0), owner(0), data(NULL), ptr(NULL), lock(0) {
	}

	// useful within init of parent object
	// or use as resize of 2D array
	void operator()(int h, int w, unsigned int flags = 0) {
		if (lock) // our object was locked so don't allow a change.
		{
			// we do a 'hang' here.
			// in a multithreaded app it should be possible to save
			// what is left.
			while (1) {
				sleep(1);
			};
		}
		lock = flags & ARRAY2D_LOCK_DATA;

		// can we reuse the current allocated data?
		// if less than a quarter of our data then reallocate.
		if ((ptr) && ((h > y)||(4*h<y))) {
			delete[] ptr;
			ptr = NULL;
		}
		if ((data) && (((h * w) > (x * y)) || ((h * w) < ((x * y) / 4)))) {
			delete[] data;
			data = NULL;
		}
		if (ptr == NULL)
			ptr = new T*[h];
		if (data == NULL)
			data = new T[h * w];
		if (flags & ARRAY2D_CLEAR_DATA)
			memset(data, 0, w * h * sizeof(T));
		x = w;
		y = h;
		for (int i = 0; i < h; i++)
			ptr[i] = data + w * i;
	}

	// import from flat data
	void operator()(int h, int w, T* copy,unsigned int flags = 0) {
		if (lock) // our object was locked so don't allow a change.
		{
			// we do a 'hang' here.
			// in a multithreaded app it should be possible to save
			// what is left.
			while (1) {
				sleep(1);
			};
		}
		lock = flags & ARRAY2D_LOCK_DATA;

		// can we reuse the current allocated data?
		// if less than a quarter of our data then reallocate.
		if (data) delete[] data;
		data = new T[h * w];
		if ((ptr) && ((h > y)||(4*h<y))) {
			delete[] ptr;
			ptr = NULL;
		}
		if (ptr == NULL)
			ptr = new T*[h];
		memcpy(data, copy, w * h * sizeof(T));
		x = w;
		y = h;
		for (int i = 0; i < h; i++)
			ptr[i] = data + w * i;
	}
};
#endif /* array2D_H_ */
