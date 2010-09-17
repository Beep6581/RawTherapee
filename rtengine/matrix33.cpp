/*
 * matrix33.cpp
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#include "matrix33.h"

namespace rtengine {

Matrix33::Matrix33 (double (*values)[3]) {

	if (values)
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				data[i][j] = values[i][j];
}

// applies transformation on the given (r,g,b) (column) vector
// result is written back to the variables passed
void Matrix33::transform (double& r, double& g, double& b) {

	double nr, ng, nb;
	nr = data[0][0]*r + data[0][1]*g + data[0][2]*b;
    ng = data[1][0]*r + data[1][1]*g + data[1][2]*b;
    nb = data[2][0]*r + data[2][1]*g + data[2][2]*b;

}

// returns invers of the transformation matrix
Matrix33 Matrix33::inverse () {

}

}
