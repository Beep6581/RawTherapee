/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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

inline void cubintch (unsigned short** src, int xs, int ys, double Dx, double Dy, unsigned short *r, double mul) {

  register double w[4];

  { register double t1, t2;
  t1 = -A*(Dx-1.0)*Dx;
  t2 = (3.0-2.0*Dx)*Dx*Dx;
  w[3] = t1*Dx;
  w[2] = t1*(Dx-1.0) + t2;
  w[1] = -t1*Dx + 1.0 - t2;
  w[0] = -t1*(Dx-1.0);
  }

  register double rd;  
  double yr[4];

  for (int k=ys, kx=0; k<ys+4; k++, kx++) {                                                               
    rd = 0.0;                                         
    for (int i=xs, ix=0; i<xs+4; i++, ix++) {                                                           
      rd += src[k][i] * w[ix];
    }                                                           
    yr[kx] = rd;
  }                                                               

                                                                    
  { register double t1, t2;
  t1 = -A*(Dy-1.0)*Dy;
  t2 = (3.0-2.0*Dy)*Dy*Dy;
  w[3] = t1*Dy;
  w[2] = t1*(Dy-1.0) + t2;
  w[1] = -t1*Dy + 1.0 - t2;
  w[0] = -t1*(Dy-1.0);
  }

  rd = 0.0;                                             
  for (int i=0; i<4; i++)
    rd += yr[i] * w[i];

  rd*=mul;

  *r = (int)CLIP(rd); 
}
