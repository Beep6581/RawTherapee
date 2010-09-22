/*
 * matrix33.cc
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#include "matrix33.h"

namespace rtengine {

Matrix33::Matrix33 (float (*values)[4]) {

	if (values)
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				data[i][j] = values[i][j];
}

Matrix33::Matrix33 () {

    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            data[i][j] = i==j ? 1.0 : 0.0;
}

// applies transformation on the given (r,g,b) (column) vector
// result is written back to the variables passed
void Matrix33::transform (double& r, double& g, double& b) {

	double nr, ng, nb;
	transform (r, g, b, nr, ng, nb);
	r = nr;
	g = ng;
	b = nb;
}
void Matrix33::transform (double r, double g, double b, double& nr, double& ng, double& nb) {

	nr = data[0][0]*r + data[0][1]*g + data[0][2]*b;
    ng = data[1][0]*r + data[1][1]*g + data[1][2]*b;
    nb = data[2][0]*r + data[2][1]*g + data[2][2]*b;
}

// returns inverse of the transformation matrix
Matrix33 Matrix33::inverse () {

	Matrix33 res;
    double nom = data[0][2]*data[1][1]*data[2][0] - data[0][1]*data[1][2]*data[2][0] - data[0][2]*data[1][0]*data[2][1] + data[0][0]*data[1][2]*data[2][1] + data[0][1]*data[1][0]*data[2][2] - data[0][0]*data[1][1]*data[2][2];
    res.data[0][0] = (data[1][2]*data[2][1]-data[1][1]*data[2][2]) / nom;
    res.data[0][1] = -(data[0][2]*data[2][1]-data[0][1]*data[2][2]) / nom;
    res.data[0][2] = (data[0][2]*data[1][1]-data[0][1]*data[1][2]) / nom;
    res.data[1][0] = -(data[1][2]*data[2][0]-data[1][0]*data[2][2]) / nom;
    res.data[1][1] = (data[0][2]*data[2][0]-data[0][0]*data[2][2]) / nom;
    res.data[1][2] = -(data[0][2]*data[1][0]-data[0][0]*data[1][2]) / nom;
    res.data[2][0] = (data[1][1]*data[2][0]-data[1][0]*data[2][1]) / nom;
    res.data[2][1] = -(data[0][1]*data[2][0]-data[0][0]*data[2][1]) / nom;
    res.data[2][2] = (data[0][1]*data[1][0]-data[0][0]*data[1][1]) / nom;

    return res;
}

}
