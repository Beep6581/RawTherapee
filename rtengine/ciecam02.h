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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include <cmath>
#include <cstdint>

#include "opthelper.h"

template<typename T>
class LUT;

using LUTu = LUT<uint32_t>;
using LUTf = LUT<float>;

namespace rtengine
{

class Ciecam02
{
private:
    static float d_factorfloat ( float f, float la );
    static float calculate_fl_from_la_ciecam02float ( float la );
    static float achromatic_response_to_whitefloat ( float x, float y, float z, float d, float fl, float nbb);
    static void xyz_to_cat02float ( float &r,  float &g,  float &b,  float x, float y, float z);
    static void cat02_to_hpefloat ( float &rh, float &gh, float &bh, float r, float g, float b);

#ifdef __SSE2__
    static void xyz_to_cat02float ( vfloat &r,  vfloat &g,  vfloat &b,  vfloat x, vfloat y, vfloat z );
    static void cat02_to_hpefloat ( vfloat &rh, vfloat &gh, vfloat &bh, vfloat r, vfloat g, vfloat b );
    static vfloat nonlinear_adaptationfloat ( vfloat c, vfloat fl );
#endif

    static float nonlinear_adaptationfloat ( float c, float fl );

    static float inverse_nonlinear_adaptationfloat ( float c, float fl );
    static void calculate_abfloat ( float &aa, float &bb, float h, float e, float t, float nbb, float a );
    static void Aab_to_rgbfloat ( float &r, float &g, float &b, float A, float aa, float bb, float nbb );
    static void hpe_to_xyzfloat   ( float &x,  float &y,  float &z,  float r, float g, float b );
    static void cat02_to_xyzfloat ( float &x,  float &y,  float &z,  float r, float g, float b);
#ifdef __SSE2__
    static vfloat inverse_nonlinear_adaptationfloat ( vfloat c, vfloat fl );
    static void calculate_abfloat ( vfloat &aa, vfloat &bb, vfloat h, vfloat e, vfloat t, vfloat nbb, vfloat a );
    static void Aab_to_rgbfloat ( vfloat &r, vfloat &g, vfloat &b, vfloat A, vfloat aa, vfloat bb, vfloat nbb );
    static void hpe_to_xyzfloat   ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b );
    static void cat02_to_xyzfloat ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b );
#endif

public:
    Ciecam02 () {}
    static void curvecolorfloat (float satind, float satval, float &sres, float parsat);
    static void curveJfloat (float br, float contr, const LUTu & histogram, LUTf & outCurve ) ;

    /**
     * Inverse transform from CIECAM02 JCh to XYZ.
     */
    static void jch2xyz_ciecam02float ( float &x, float &y, float &z,
                                        float J, float C, float h,
                                        float xw, float yw, float zw,
                                        float c, float nc, float n, float nbb, float ncb, float fl, float cz, float d, float aw );
#ifdef __SSE2__
    static void jch2xyz_ciecam02float ( vfloat &x, vfloat &y, vfloat &z,
                                        vfloat J, vfloat C, vfloat h,
                                        vfloat xw, vfloat yw, vfloat zw,
                                        vfloat nc, vfloat n, vfloat nbb, vfloat ncb, vfloat fl, vfloat d, vfloat aw, vfloat reccmcz );
#endif
    /**
     * Forward transform from XYZ to CIECAM02 JCh.
     */
    static void initcam1float (float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                               float &cz, float &aw, float &wh, float &pfl, float &fl, float c);

    static void initcam2float (float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                               float &cz, float &aw, float &fl);

    static void xyz2jch_ciecam02float ( float &J, float &C, float &h,
                                        float aw, float fl,
                                        float x, float y, float z,
                                        float xw, float yw, float zw,
                                        float c, float nc, float n, float nbb, float ncb, float cz, float d  );

    static void xyz2jchqms_ciecam02float ( float &J, float &C, float &h,
                                           float &Q, float &M, float &s, float aw, float fl, float wh,
                                           float x, float y, float z,
                                           float xw, float yw, float zw,
                                           float c, float nc, float n, float nbb, float ncb, float pfl, float cz, float d  );

#ifdef __SSE2__
    static void xyz2jchqms_ciecam02float ( vfloat &J, vfloat &C, vfloat &h,
                                           vfloat &Q, vfloat &M, vfloat &s, vfloat aw, vfloat fl, vfloat wh,
                                           vfloat x, vfloat y, vfloat z,
                                           vfloat xw, vfloat yw, vfloat zw,
                                           vfloat c, vfloat nc, vfloat n, vfloat nbb, vfloat ncb, vfloat pfl, vfloat cz, vfloat d  );


#endif

};
}
