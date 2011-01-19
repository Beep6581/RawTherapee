/*
 * matrix33.cc
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#include "matrix33.h"
#include "macros.h"

namespace rtengine {

Matrix33::Matrix33 (float (*values)[4]) {

	if (values)
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				data[i][j] = values[i][j];
}

Matrix33::Matrix33 (float (*values)[3]) {

	if (values)
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				data[i][j] = values[i][j];
}

Matrix33::Matrix33 (float d00, float d01, float d02, float d10, float d11, float d12, float d20, float d21, float d22) {

	data[0][0] = d00;
	data[0][1] = d01;
	data[0][2] = d02;
	data[1][0] = d10;
	data[1][1] = d11;
	data[1][2] = d12;
	data[2][0] = d20;
	data[2][1] = d21;
	data[2][2] = d22;
}


Matrix33::Matrix33 () {

    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            data[i][j] = i==j ? 1.0 : 0.0;
}

// applies transformation on the given (r,g,b) (column) vector
// result is written back to the variables passed
void Matrix33::transform (float& r, float& g, float& b) const {

	float nr, ng, nb;
	transform (r, g, b, nr, ng, nb);
	r = nr;
	g = ng;
	b = nb;
}

void Matrix33::transform (unsigned short& r, unsigned short& g, unsigned short& b) const {

    unsigned short nr, ng, nb;
    transform (r, g, b, nr, ng, nb);
    r = nr;
    g = ng;
    b = nb;
}

void Matrix33::transform (float r, float g, float b, float& nr, float& ng, float& nb) const {

    nr = data[0][0]*r + data[0][1]*g + data[0][2]*b;
    ng = data[1][0]*r + data[1][1]*g + data[1][2]*b;
    nb = data[2][0]*r + data[2][1]*g + data[2][2]*b;
}

void Matrix33::transform (unsigned short r, unsigned short g, unsigned short b, unsigned short& nr, unsigned short& ng, unsigned short& nb) const {

	float dnr = data[0][0]*(float)r + data[0][1]*(float)g + data[0][2]*(float)b;
	float dng = data[1][0]*(float)r + data[1][1]*(float)g + data[1][2]*(float)b;
	float dnb = data[2][0]*(float)r + data[2][1]*(float)g + data[2][2]*(float)b;
    nr = CLIP((unsigned short)dnr);
    ng = CLIP((unsigned short)dng);
    nb = CLIP((unsigned short)dnb);
}

// returns inverse of the transformation matrix
Matrix33 Matrix33::inverse () const {

	Matrix33 res;
	float nom = data[0][2]*data[1][1]*data[2][0] - data[0][1]*data[1][2]*data[2][0] - data[0][2]*data[1][0]*data[2][1] + data[0][0]*data[1][2]*data[2][1] + data[0][1]*data[1][0]*data[2][2] - data[0][0]*data[1][1]*data[2][2];
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

void Matrix33::multiply (const Matrix33& m) {

    Matrix33 r;
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            r.data[i][j] = 0.0;
            for (int k=0; k<3; k++)
                r.data[i][j] += data[i][k] * m.data[k][j];
        }
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            data[i][j] = r.data[i][j];
}

float Matrix33::rowsum (int i) const {

	float r = 0.0;
	for (int j=0; j<3; j++)
		r += data[i][j];

	return r;
}

}
