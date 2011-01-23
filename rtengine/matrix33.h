/*
 * matrix33.h
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#ifndef MATRIX33_H_
#define MATRIX33_H_

#include <string.h>
#include <iostream>

namespace rtengine {

// This class represents a 3x3 transformation matrix used for converting
// color (r,g,b) triplets between different color spaces
class Matrix33 {

	public:

		float data[3][3];

		// initializes the matrix with the array passed
        Matrix33 (float (*values)[4]);
        Matrix33 (float (*values)[3]);
        Matrix33 (float d00, float d01, float d02, float d10, float d11, float d12, float d20, float d21, float d22);
        Matrix33 ();
        Matrix33 (const Matrix33& other);

		// applies transformation on the given (r,g,b) (column) vector
		// result is written back to the variables passed
        void transform (float& r, float& g, float& b) const;
        void transform (unsigned short& r, unsigned short& g, unsigned short& b) const;

		// the same, result is stored separately
        void transform (float r, float g, float b, float& nr, float& ng, float& nb) const;
        void transform (unsigned short r, unsigned short g, unsigned short b, unsigned short& nr, unsigned short& ng, unsigned short& nb) const;

		// returns inverse of the transformation matrix
		Matrix33 inverse () const;

		// multiplies this matrix from the given one from the right
		void multiply (const Matrix33& other);

		// gives back the sum of the elements of the ith row
		float rowsum (int i) const;

};

std::ostream& operator<< (std::ostream& os, const Matrix33& m);

}

#endif /* MATRIX33_H_ */
