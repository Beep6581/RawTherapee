////////////////////////////////////////////////////////////////
//
//		Gamut mapping algorithm
//
//		copyright (c) 2010-2011  Emil Martinec <ejmartin@uchicago.edu>
//
//
// code dated: February 2, 2011
//
//	sRGBgamutbdy.cc is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//#include "color.h"

#define SQR(x) ((x)*(x))

#define u0 4.0*D50x/(D50x+15+3*D50z)
#define v0 9.0/(D50x+15+3*D50z)

#define sgn(x) (x<0 ? -1 : 1)
#define sqrt2 (sqrt(2))
#define sqrt3 (sqrt(3))
#define sqrt6 (sqrt(6))
#define K 20
#define rmax	(3*K + 1)
#define rctr	(2*K)
#define cctr	(2*K + 1)
#define maxindx	(2*(3*K+1) + SQR(3*K+1))



// solutions to scaling u and v to XYZ paralleliped boundaries
/* some equations:
 fu(X,Y,Z) = 4 X/(X + 15 Y + 3 Z);
 fv(X,Y,Z) = 9 Y/(X + 15 Y + 3 Z);
 
 take the plane spanned by X=a*Xr+b*Xg+c*Xb etc with one of a,b,c equal to 0 or 1,
 and itersect with the line u0+lam*u, or in other words solve
 
 u0+lam*u=fu(X,Y,Z)
 v0+lam*v=fv(X,Y,Z)
 
 the value of lam is the scale factor that takes the color to the gamut boundary
*/
// columns of the matrix p=xyz_rgb are RGB tristimulus primaries in XYZ
// c is the color fixed on the boundary; and m=0 for c=0, m=1 for c=255

void Color::gamutmap(float &X, float &Y, float &Z, const double p[3][3])
{
	float u = 4*X/(X+15*Y+3*Z)-u0;
	float v = 9*Y/(X+15*Y+3*Z)-v0;
	
	float lam[3][2];
	float lam_min = 1.0;
	
	for (int c=0; c<3; c++)
		for (int m=0; m<2; m++) {
			
			int c1=(c+1)%3;
			int c2=(c+2)%3;
			
			lam[c][m] = (-(p[0][c1]*p[1][c]*((-12 + 3*u0 + 20*v0)*Y + 4*m*65535*v0*p[2][c2])) + \
						 p[0][c]*p[1][c1]*((-12 + 3*u0 + 20*v0)*Y + 4*m*65535*v0*p[2][c2]) - \
						 4*v0*p[0][c1]*(Y - m*65535*p[1][c2])*p[2][c] + 4*v0*p[0][c]*(Y - m*65535*p[1][c2])*p[2][c1] - \
						 (4*m*65535*v0*p[0][c2] - 9*u0*Y)*(p[1][c1]*p[2][c] - p[1][c]*p[2][c1]));
			
			lam[c][m] /= (3*u*Y*(p[0][c1]*p[1][c] - p[1][c1]*(p[0][c] + 3*p[2][c]) + 3*p[1][c]*p[2][c1]) + \
					4*v*(p[0][c1]*(5*Y*p[1][c] + m*65535*p[1][c]*p[2][c2] + Y*p[2][c] - m*65535*p[1][c2]*p[2][c]) - \
						 p[0][c]*(5*Y*p[1][c1] + m*65535*p[1][c1]*p[2][c2] + Y*p[2][c1] - m*65535*p[1][c2]*p[2][c1]) + \
						 m*65535*p[0][c2]*(p[1][c1]*p[2][c] - p[1][c]*p[2][c1])));
			
			if (lam[c][m]<lam_min && lam[c][m]>0) lam_min=lam[c][m];

		}
	
	u = u*lam_min + u0;
	v = v*lam_min + v0;
	
	X = (9*u*Y)/(4*v); 
	Z = (12 - 3*u - 20*v)*Y/(4*v);

}

	
//};//namespace
