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
#ifndef _CIECAM02_
#define _CIECAM02_
#include <cmath>
#include "LUT.h"
#include "opthelper.h"

namespace rtengine
{

class Ciecam02
{
private:
    static double d_factor( double f, double la );
    static float d_factorfloat( float f, float la );
    static double calculate_fl_from_la_ciecam02( double la );
    static float calculate_fl_from_la_ciecam02float( float la );
    static double achromatic_response_to_white( double x, double y, double z, double d, double fl, double nbb, int gamu );
    static float achromatic_response_to_whitefloat( float x, float y, float z, float d, float fl, float nbb, int gamu );
    static void xyz_to_cat02 ( double &r,  double &g,  double &b,  double x, double y, double z, int gamu );
    static void cat02_to_hpe ( double &rh, double &gh, double &bh, double r, double g, double b, int gamu );
    static void cat02_to_xyz ( double &x,  double &y,  double &z,  double r, double g, double b, int gamu );
    static void hpe_to_xyz   ( double &x,  double &y,  double &z,  double r, double g, double b );

    static void xyz_to_cat02float ( float &r,  float &g,  float &b,  float x, float y, float z, int gamu );
    static void cat02_to_hpefloat ( float &rh, float &gh, float &bh, float r, float g, float b, int gamu );

#ifdef __SSE2__
    static void xyz_to_cat02float ( vfloat &r,  vfloat &g,  vfloat &b,  vfloat x, vfloat y, vfloat z );
    static void cat02_to_hpefloat ( vfloat &rh, vfloat &gh, vfloat &bh, vfloat r, vfloat g, vfloat b );
    static vfloat nonlinear_adaptationfloat( vfloat c, vfloat fl );
#endif

    static void Aab_to_rgb( double &r, double &g, double &b, double A, double aa, double bb, double nbb );
    static void calculate_ab( double &aa, double &bb, double h, double e, double t, double nbb, double a );

    static double nonlinear_adaptation( double c, double fl );
    static float nonlinear_adaptationfloat( float c, float fl );
    static double inverse_nonlinear_adaptation( double c, double fl );


    static float inverse_nonlinear_adaptationfloat( float c, float fl );
    static void calculate_abfloat( float &aa, float &bb, float h, float e, float t, float nbb, float a );
    static void Aab_to_rgbfloat( float &r, float &g, float &b, float A, float aa, float bb, float nbb );
    static void hpe_to_xyzfloat   ( float &x,  float &y,  float &z,  float r, float g, float b );
    static void cat02_to_xyzfloat ( float &x,  float &y,  float &z,  float r, float g, float b, int gamu );
#ifdef __SSE2__
    static vfloat inverse_nonlinear_adaptationfloat( vfloat c, vfloat fl );
    static void calculate_abfloat( vfloat &aa, vfloat &bb, vfloat h, vfloat e, vfloat t, vfloat nbb, vfloat a );
    static void Aab_to_rgbfloat( vfloat &r, vfloat &g, vfloat &b, vfloat A, vfloat aa, vfloat bb, vfloat nbb );
    static void hpe_to_xyzfloat   ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b );
    static void cat02_to_xyzfloat ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b );
#endif

public:
    Ciecam02 () {}
    static void curvecolor(double satind, double satval, double &sres, double parsat);
    static void curvecolorfloat(float satind, float satval, float &sres, float parsat);
    static void curveJ (double br, double contr, int db, LUTf & outCurve , LUTu & histogram ) ;
    static void curveJfloat (float br, float contr, const LUTu & histogram, LUTf & outCurve ) ;

    /**
     * Inverse transform from CIECAM02 JCh to XYZ.
     */
    static void jch2xyz_ciecam02( double &x, double &y, double &z,
                                  double J, double C, double h,
                                  double xw, double yw, double zw,
                                  double yb, double la,
                                  double f, double c, double nc, int gamu, double n, double nbb, double ncb, double fl, double cz, double d, double aw);


    static void jch2xyz_ciecam02float( float &x, float &y, float &z,
                                       float J, float C, float h,
                                       float xw, float yw, float zw,
                                       float f, float c, float nc, int gamu, float n, float nbb, float ncb, float fl, float cz, float d, float aw );
#ifdef __SSE2__
    static void jch2xyz_ciecam02float( vfloat &x, vfloat &y, vfloat &z,
                                       vfloat J, vfloat C, vfloat h,
                                       vfloat xw, vfloat yw, vfloat zw,
                                       vfloat f, vfloat nc, vfloat n, vfloat nbb, vfloat ncb, vfloat fl, vfloat d, vfloat aw, vfloat reccmcz );
#endif
    /**
     * Forward transform from XYZ to CIECAM02 JCh.
     */
    static void initcam1(double gamu, double yb, double pilotd, double f, double la, double xw, double yw, double zw, double &n, double &d, double &nbb, double &ncb,
                         double &cz, double &aw, double &wh, double &pfl, double &fl, double &c);

    static void initcam2(double gamu, double yb, double f, double la, double xw, double yw, double zw, double &n, double &d, double &nbb, double &ncb,
                         double &cz, double &aw, double &fl);

    static void initcam1float(float gamu, float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                              float &cz, float &aw, float &wh, float &pfl, float &fl, float &c);

    static void initcam2float(float gamu, float yb, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                              float &cz, float &aw, float &fl);

    static void xyz2jchqms_ciecam02( double &J, double &C, double &h,
                                     double &Q, double &M, double &s, double &aw, double &fl, double &wh,
                                     double x, double y, double z,
                                     double xw, double yw, double zw,
                                     double yb, double la,
                                     double f, double c, double nc,  double pilotd, int gamu , double n, double nbb, double ncb, double pfl, double cz, double d );

    static void xyz2jch_ciecam02float( float &J, float &C, float &h,
                                       float aw, float fl,
                                       float x, float y, float z,
                                       float xw, float yw, float zw,
                                       float c, float nc, float n, float nbb, float ncb, float cz, float d  );

    static void xyz2jchqms_ciecam02float( float &J, float &C, float &h,
                                          float &Q, float &M, float &s, float aw, float fl, float wh,
                                          float x, float y, float z,
                                          float xw, float yw, float zw,
                                          float c, float nc, int gamu, float n, float nbb, float ncb, float pfl, float cz, float d  );

#ifdef __SSE2__
    static void xyz2jchqms_ciecam02float( vfloat &J, vfloat &C, vfloat &h,
                                          vfloat &Q, vfloat &M, vfloat &s, vfloat aw, vfloat fl, vfloat wh,
                                          vfloat x, vfloat y, vfloat z,
                                          vfloat xw, vfloat yw, vfloat zw,
                                          vfloat c, vfloat nc, vfloat n, vfloat nbb, vfloat ncb, vfloat pfl, vfloat cz, vfloat d  );


#endif

};
}
#endif
