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
 *      	array2D<float> my_array (10,10);
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
 *			<type> ** my_array		gives access to the pointer
 *			<type> *  my_array		gives access to the stored data.
 *
 */

template<typename T>
class array2D {
private:
	int x, y, owner;
	T ** ptr;
	T * data;
public:
	array2D(int h, int w) {
		data = new T[h * w];
		owner = 1;
		x = w;
		y = h;
		ptr = new T*[h];
		for (int i = 0; i < h; i++)
			ptr[i] = data + i * w;
	}

	array2D(int h, int w, T ** source) {
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
		if (owner)
			delete[] data;
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
};
