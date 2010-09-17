/*
 * matrix33.h
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#ifndef MATRIX33_H_
#define MATRIX33_H_

namespace rtengine {

// This class represents a 3x3 transformation matrix used for converting
// color (r,g,b) triplets between different color spaces
class Matrix33 {

		double data[3][3];

	public:
		// initializes the matrix with the array passed
		Matrix33 (double (*values)[3] = NULL);

		// applies transformation on the given (r,g,b) (column) vector
		// result is written back to the variables passed
		void transform (double& r, double& g, double& b);

		// the same, result is stored separately
		void transform (double r, double g, double b, double& nr, double& ng, double& nb);

		// returns inverse of the transformation matrix
		Matrix33 inverse ();
};

}

#endif /* MATRIX33_H_ */
