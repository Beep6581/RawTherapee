/*
 * matrix33.h
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#ifndef MATRIX33_H_
#define MATRIX33_H_

#include <string.h>

namespace rtengine {

// This class represents a 3x3 transformation matrix used for converting
// color (r,g,b) triplets between different color spaces
class Matrix33 {

		float data[3][3];

	public:
		// initializes the matrix with the array passed
        Matrix33 (float (*values)[4]);
        Matrix33 ();

		// applies transformation on the given (r,g,b) (column) vector
		// result is written back to the variables passed
        void transform (float& r, float& g, float& b);
        void transform (unsigned short& r, unsigned short& g, unsigned short& b);

		// the same, result is stored separately
        void transform (float r, float g, float b, float& nr, float& ng, float& nb);
        void transform (unsigned short r, unsigned short g, unsigned short b, unsigned short& nr, unsigned short& ng, unsigned short& nb);

		// returns inverse of the transformation matrix
		Matrix33 inverse ();

		// multiplies this matrix from the given one from the right
		void multiply (const float (*m)[3]);
};

}

#endif /* MATRIX33_H_ */
