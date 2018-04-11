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
#include "ciecam02.h"
#include "rtengine.h"
#include "curves.h"
#include <math.h>
#include "sleef.c"

#ifdef _DEBUG
#include "settings.h"
#include <stdio.h>
#endif

#undef CLIPD
#define CLIPD(a) ((a)>0.0?((a)<1.0?(a):1.0):0.0)
#define MAXR(a,b) ((a) > (b) ? (a) : (b))

namespace rtengine
{

#ifdef _DEBUG
extern const Settings* settings;
#endif

void Ciecam02::curvecolorfloat (float satind, float satval, float &sres, float parsat)
{
    if (satind > 0.f) {
        if (satval >= 1.f) { // The calculation below goes wrong direction when satval > 1
            sres = satval;
        } else {
            sres = (1.f - (satind) / 100.f) * satval + (satind) / 100.f * (1.f - SQR (SQR (1.f - min (satval, 1.0f))));
        }

        if (sres > parsat) {
            sres = max (parsat, satval);
        }
    } else if (satind < 0.f) {
        sres = satval * (1.f + (satind) / 100.f);
    } else { // satind == 0 means we don't want to change the value at all
        sres = satval;
    }
}

void Ciecam02::curveJfloat (float br, float contr, const LUTu & histogram, LUTf & outCurve)
{

    // check if brightness curve is needed
    if (br > 0.00001f || br < -0.00001f) {

        std::vector<double> brightcurvePoints (9);
        brightcurvePoints[0] = double (DCT_NURBS);

        brightcurvePoints[1] = 0.f; // black point.  Value in [0 ; 1] range
        brightcurvePoints[2] = 0.f; // black point.  Value in [0 ; 1] range

        if (br > 0) {
            brightcurvePoints[3] = 0.1f; // toe point
            brightcurvePoints[4] = 0.1f + br / 150.0f; //value at toe point

            brightcurvePoints[5] = 0.7f; // shoulder point
            brightcurvePoints[6] = min (1.0f, 0.7f + br / 300.0f); //value at shoulder point
        } else {
            brightcurvePoints[3] = 0.1f - br / 150.0f; // toe point
            brightcurvePoints[4] = 0.1f; // value at toe point

            brightcurvePoints[5] = min (1.0f, 0.7f - br / 300.0f); // shoulder point
            brightcurvePoints[6] = 0.7f; // value at shoulder point
        }

        brightcurvePoints[7] = 1.f; // white point
        brightcurvePoints[8] = 1.f; // value at white point

        DiagonalCurve brightcurve (brightcurvePoints, CURVES_MIN_POLY_POINTS);

        // Applying brightness curve
        for (int i = 0; i < 32768; i++) {

            // change to [0,1] range
            float val = (float)i / 32767.0f;

            // apply brightness curve
            val = brightcurve.getVal (val);

            // store result
            outCurve[i] = CLIPD (val);
        }

    } else {
        // set the identity curve
        outCurve.makeIdentity (32767.f);
    }


    if (contr > 0.00001f || contr < -0.00001f) {

        // compute mean luminance of the image with the curve applied
        float sum = 0.f;
        float avg = 0.f;

        for (int i = 0; i < 32768; i++) {
            avg += outCurve[i] * histogram[i];//approximation for average : usage of L (lab) instead of J
            sum += histogram[i];
        }

        avg /= sum;
        std::vector<double> contrastcurvePoints (9);

        contrastcurvePoints[0] = double (DCT_NURBS);

        contrastcurvePoints[1] = 0.f; // black point.  Value in [0 ; 1] range
        contrastcurvePoints[2] = 0.f; // black point.  Value in [0 ; 1] range

        contrastcurvePoints[3] = avg - avg * (0.6f - contr / 250.0f); // toe point
        contrastcurvePoints[4] = avg - avg * (0.6f + contr / 250.0f); // value at toe point

        contrastcurvePoints[5] = avg + (1 - avg) * (0.6f - contr / 250.0f); // shoulder point
        contrastcurvePoints[6] = avg + (1 - avg) * (0.6f + contr / 250.0f); // value at shoulder point

        contrastcurvePoints[7] = 1.f; // white point
        contrastcurvePoints[8] = 1.f; // value at white point

        DiagonalCurve contrastcurve (contrastcurvePoints, CURVES_MIN_POLY_POINTS);

        // apply contrast enhancement
        for (int i = 0; i < 32768; i++) {
            outCurve[i] = contrastcurve.getVal (outCurve[i]);
        }

    }

    outCurve *= 32767.f;
	//printf("out500=%f out15000=%f\n", outCurve[500], outCurve[15000]);
	//outCurve.dump("brig");
}

/**
 * Copyright (c) 2003 Billy Biggs <vektor@dumbterm.net>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

float Ciecam02::d_factorfloat ( float f, float la )
{
    return f * (1.0f - ((1.0f / 3.6f) * xexpf ((-la - 42.0f) / 92.0f)));
}

float Ciecam02::calculate_fl_from_la_ciecam02float ( float la )
{
    float la5 = la * 5.0f;
    float k = 1.0f / (la5 + 1.0f);

    /* Calculate k^4. */
    k = k * k;
    k = k * k;

    return (0.2f * k * la5) + (0.1f * (1.0f - k) * (1.0f - k) * std::cbrt (la5));
}

float Ciecam02::achromatic_response_to_whitefloat ( float x, float y, float z, float d, float fl, float nbb, int gamu )
{
    float r, g, b;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    gamu = 1;
    xyz_to_cat02float ( r, g, b, x, y, z, gamu );

    rc = r * (((y * d) / r) + (1.0f - d));
    gc = g * (((y * d) / g) + (1.0f - d));
    bc = b * (((y * d) / b) + (1.0f - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, gamu );

    if (gamu == 1) { //gamut correction M.H.Brill S.Susstrunk
        rp = MAXR (rp, 0.0f);
        gp = MAXR (gp, 0.0f);
        bp = MAXR (bp, 0.0f);
    }

    rpa = nonlinear_adaptationfloat ( rp, fl );
    gpa = nonlinear_adaptationfloat ( gp, fl );
    bpa = nonlinear_adaptationfloat ( bp, fl );

    return ((2.0f * rpa) + gpa + ((1.0f / 20.0f) * bpa) - 0.305f) * nbb;
}

void Ciecam02::xyz_to_cat02float ( float &r, float &g, float &b, float x, float y, float z, int gamu )
{
    gamu = 1;

    if (gamu == 0) {
        r = ( 0.7328f * x) + (0.4296f * y) - (0.1624f * z);
        g = (-0.7036f * x) + (1.6975f * y) + (0.0061f * z);
        b = ( 0.0030f * x) + (0.0136f * y) + (0.9834f * z);
    } else if (gamu == 1) { //gamut correction M.H.Brill S.Susstrunk
        //r = ( 0.7328 * x) + (0.4296 * y) - (0.1624 * z);
        //g = (-0.7036 * x) + (1.6975 * y) + (0.0061 * z);
        //b = ( 0.0000 * x) + (0.0000 * y) + (1.0000 * z);
        r = ( 1.007245f * x) + (0.011136f * y) - (0.018381f * z); //Changjun Li
        g = (-0.318061f * x) + (1.314589f * y) + (0.003471f * z);
        b = ( 0.0000f * x) + (0.0000f * y) + (1.0000f * z);
    }
}
#ifdef __SSE2__
void Ciecam02::xyz_to_cat02float ( vfloat &r, vfloat &g, vfloat &b, vfloat x, vfloat y, vfloat z )
{
    //gamut correction M.H.Brill S.Susstrunk
    r = ( F2V (1.007245f) * x) + (F2V (0.011136f) * y) - (F2V (0.018381f) * z); //Changjun Li
    g = (F2V (-0.318061f) * x) + (F2V (1.314589f) * y) + (F2V (0.003471f) * z);
    b = z;
}
#endif

void Ciecam02::cat02_to_xyzfloat ( float &x, float &y, float &z, float r, float g, float b, int gamu )
{
    gamu = 1;

    if (gamu == 0) {
        x = ( 1.096124f * r) - (0.278869f * g) + (0.182745f * b);
        y = ( 0.454369f * r) + (0.473533f * g) + (0.072098f * b);
        z = (-0.009628f * r) - (0.005698f * g) + (1.015326f * b);
    } else if (gamu == 1) { //gamut correction M.H.Brill S.Susstrunk
        //x = ( 1.0978566 * r) - (0.277843 * g) + (0.179987 * b);
        //y = ( 0.455053 * r) + (0.473938 * g) + (0.0710096* b);
        //z = ( 0.000000 * r) - (0.000000 * g) + (1.000000 * b);
        x = ( 0.99015849f * r) - (0.00838772f * g) + (0.018229217f * b); //Changjun Li
        y = ( 0.239565979f * r) + (0.758664642f * g) + (0.001770137f * b);
        z = ( 0.000000f * r) - (0.000000f * g) + (1.000000f * b);
    }
}
#ifdef __SSE2__
void Ciecam02::cat02_to_xyzfloat ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b )
{
    //gamut correction M.H.Brill S.Susstrunk
    x = ( F2V (0.99015849f) * r) - (F2V (0.00838772f) * g) + (F2V (0.018229217f) * b); //Changjun Li
    y = ( F2V (0.239565979f) * r) + (F2V (0.758664642f) * g) + (F2V (0.001770137f) * b);
    z = b;
}
#endif

void Ciecam02::hpe_to_xyzfloat ( float &x, float &y, float &z, float r, float g, float b )
{
    x = (1.910197f * r) - (1.112124f * g) + (0.201908f * b);
    y = (0.370950f * r) + (0.629054f * g) - (0.000008f * b);
    z = b;
}
#ifdef __SSE2__
void Ciecam02::hpe_to_xyzfloat ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b )
{
    x = (F2V (1.910197f) * r) - (F2V (1.112124f) * g) + (F2V (0.201908f) * b);
    y = (F2V (0.370950f) * r) + (F2V (0.629054f) * g) - (F2V (0.000008f) * b);
    z = b;
}
#endif

void Ciecam02::cat02_to_hpefloat ( float &rh, float &gh, float &bh, float r, float g, float b, int gamu )
{
    gamu = 1;

    if (gamu == 0) {
        rh = ( 0.7409792f * r) + (0.2180250f * g) + (0.0410058f * b);
        gh = ( 0.2853532f * r) + (0.6242014f * g) + (0.0904454f * b);
        bh = (-0.0096280f * r) - (0.0056980f * g) + (1.0153260f * b);
    } else if (gamu == 1) { //Changjun Li
        rh = ( 0.550930835f * r) + (0.519435987f * g) - ( 0.070356303f * b);
        gh = ( 0.055954056f * r) + (0.89973132f * g) + (0.044315524f * b);
        bh = (0.0f * r) - (0.0f * g) + (1.0f * b);
    }
}

#ifdef __SSE2__
void Ciecam02::cat02_to_hpefloat ( vfloat &rh, vfloat &gh, vfloat &bh, vfloat r, vfloat g, vfloat b)
{
    //Changjun Li
    rh = ( F2V (0.550930835f) * r) + (F2V (0.519435987f) * g) - ( F2V (0.070356303f) * b);
    gh = ( F2V (0.055954056f) * r) + (F2V (0.89973132f) * g) + (F2V (0.044315524f) * b);
    bh = b;
}
#endif

void Ciecam02::Aab_to_rgbfloat ( float &r, float &g, float &b, float A, float aa, float bb, float nbb )
{
    float x = (A / nbb) + 0.305f;

    /*       c1              c2               c3       */
    r = (0.32787f * x) + (0.32145f * aa) + (0.20527f * bb);
    /*       c1              c4               c5       */
    g = (0.32787f * x) - (0.63507f * aa) - (0.18603f * bb);
    /*       c1              c6               c7       */
    b = (0.32787f * x) - (0.15681f * aa) - (4.49038f * bb);
}
#ifdef __SSE2__
void Ciecam02::Aab_to_rgbfloat ( vfloat &r, vfloat &g, vfloat &b, vfloat A, vfloat aa, vfloat bb, vfloat nbb )
{
    vfloat c1 = F2V (0.32787f) * ((A / nbb) + F2V (0.305f));

    /*       c1              c2               c3       */
    r = c1 + (F2V (0.32145f) * aa) + (F2V (0.20527f) * bb);
    /*       c1              c4               c5       */
    g = c1 - (F2V (0.63507f) * aa) - (F2V (0.18603f) * bb);
    /*       c1              c6               c7       */
    b = c1 - (F2V (0.15681f) * aa) - (F2V (4.49038f) * bb);
}
#endif

void Ciecam02::calculate_abfloat ( float &aa, float &bb, float h, float e, float t, float nbb, float a )
{
    float2 sincosval = xsincosf(h * rtengine::RT_PI_F_180);
    float sinh = sincosval.x;
    float cosh = sincosval.y;
    float x = (a / nbb) + 0.305f;
    constexpr float p3 = 1.05f;
    const bool swapValues = fabs(sinh) > fabs(cosh);

    if (swapValues) {
        std::swap(sinh, cosh);
    }

    float c1 = 1.f;
    float c2 = sinh / cosh;

    if (swapValues) {
        std::swap(c1, c2);
    }

    float div = ((e / (t * cosh)) - (-0.31362f - (p3 * 0.15681f)) * c1 - ((0.01924f - (p3 * 4.49038f)) * c2));
    // for large values of t the above calculation can change its sign which results in a hue shift of 180 degree
    // so we have to check the sign to avoid this shift.
    // Additionally it seems useful to limit the minimum value of div
    // I limited it, but I'm sure the actual limit is not the best one

    if (signf(div) != signf(cosh) || fabsf(div) <= fabsf(cosh) * 2.f) {
        div = cosh * 2.f;
    }

    aa = ((0.32787f * x) * (2.0f + p3)) / div;
    bb = (aa * sinh) / cosh;

    if (swapValues) {
        std::swap(aa, bb);
    }
}
#ifdef __SSE2__
void Ciecam02::calculate_abfloat ( vfloat &aa, vfloat &bb, vfloat h, vfloat e, vfloat t, vfloat nbb, vfloat a )
{
    vfloat2 sincosval = xsincosf ((h * F2V (rtengine::RT_PI)) / F2V (180.0f));
    vfloat sinh = sincosval.x;
    vfloat cosh = sincosval.y;
    vfloat x = (a / nbb) + F2V (0.305f);
    vfloat p3 = F2V (1.05f);
    vmask swapMask = vmaskf_gt (vabsf (sinh), vabsf (cosh));
    vswap (swapMask, sinh, cosh);
    vfloat c1 = F2V (1.f);
    vfloat c2 = sinh / cosh;
    vswap (swapMask, c1, c2);

    vfloat div = ((e / (t * cosh)) - (F2V (-0.31362f) - (p3 * F2V (0.15681f))) * c1 - ((F2V (0.01924f) - (p3 * F2V (4.49038f))) * (c2)));
    // for large values of t the above calculation can change its sign which results in a hue shift of 180 degree
    // so we have to check the sign to avoid this shift.
    // Additionally it seems useful to limit the minimum value of div
    // I limited it, but I'm sure the actual limit is not the best one

    vmask limitMask = vmaskf_neq (vsignf (div), vsignf (cosh));
    limitMask = vorm (limitMask, vmaskf_le (vabsf (div), vabsf (cosh) * F2V (2.f)));
    div = vself (limitMask, cosh * F2V (2.f), div);

    aa = ((F2V (0.32787f) * x) * (F2V (2.0f) + p3)) / div;
    bb = (aa * sinh) / cosh;

    vswap (swapMask, aa, bb);
}

#endif

void Ciecam02::initcam1float (float gamu, float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                              float &cz, float &aw, float &wh, float &pfl, float &fl, float &c)
{
    n = yb / yw;

    if (pilotd == 2.0) {
        d = d_factorfloat ( f, la );
    } else {
        d = pilotd;
    }

    fl = calculate_fl_from_la_ciecam02float ( la );
    nbb = ncb = 0.725f * pow_F ( 1.0f / n, 0.2f );
    cz = 1.48f + sqrt ( n );
    aw = achromatic_response_to_whitefloat ( xw, yw, zw, d, fl, nbb, gamu );
    wh = ( 4.0f / c ) * ( aw + 4.0f ) * pow_F ( fl, 0.25f );
    pfl = pow_F ( fl, 0.25f );
#ifdef _DEBUG

    if (settings->verbose) {
        printf ("Source float d=%f aw=%f fl=%f wh=%f c=%f  awc=%f\n", d, aw, fl, wh, c, (4.f / c) * (aw + 4.f));
    }

#endif
}

void Ciecam02::initcam2float (float gamu, float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                              float &cz, float &aw, float &fl)
{
    n = yb / yw;

    if (pilotd == 2.0) {
        d = d_factorfloat ( f, la );
    } else {
        d = pilotd;
    }

//   d = d_factorfloat( f, la );
    fl = calculate_fl_from_la_ciecam02float ( la );
    nbb = ncb = 0.725f * pow_F ( 1.0f / n, 0.2f );
    cz = 1.48f + sqrt ( n );
    aw = achromatic_response_to_whitefloat ( xw, yw, zw, d, fl, nbb, gamu );
#ifdef _DEBUG

    if (settings->verbose) {
        printf ("Viewing float d=%f aw=%f fl=%f n=%f\n", d, aw, fl, n);
    }

#endif
}

void Ciecam02::xyz2jchqms_ciecam02float ( float &J, float &C, float &h, float &Q, float &M, float &s, float aw, float fl, float wh,
        float x, float y, float z, float xw, float yw, float zw,
        float c, float nc, int gamu, float pow1, float nbb, float ncb, float pfl, float cz, float d)

{
    float r, g, b;
    float rw, gw, bw;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    float a, ca, cb;
    float e, t;
    float myh;
    gamu = 1;
    xyz_to_cat02float ( r, g, b, x, y, z, gamu );
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw, gamu );
    rc = r * (((yw * d) / rw) + (1.f - d));
    gc = g * (((yw * d) / gw) + (1.f - d));
    bc = b * (((yw * d) / bw) + (1.f - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, gamu );

    if (gamu == 1) { //gamut correction M.H.Brill S.Susstrunk
        rp = MAXR (rp, 0.0f);
        gp = MAXR (gp, 0.0f);
        bp = MAXR (bp, 0.0f);
    }

    rpa = nonlinear_adaptationfloat ( rp, fl );
    gpa = nonlinear_adaptationfloat ( gp, fl );
    bpa = nonlinear_adaptationfloat ( bp, fl );

    ca = rpa - ((12.0f * gpa) - bpa) / 11.0f;
    cb = (0.11111111f) * (rpa + gpa - (2.0f * bpa));

    myh = xatan2f ( cb, ca );

    if ( myh < 0.0f ) {
        myh += (2.f * rtengine::RT_PI);
    }

    a = ((2.0f * rpa) + gpa + (0.05f * bpa) - 0.305f) * nbb;

    if (gamu == 1) {
        a = MAXR (a, 0.0f); //gamut correction M.H.Brill S.Susstrunk
    }

    J = pow_F ( a / aw, c * cz * 0.5f);

    e = ((961.53846f) * nc * ncb) * (xcosf ( myh + 2.0f ) + 3.8f);
    t = (e * sqrtf ( (ca * ca) + (cb * cb) )) / (rpa + gpa + (1.05f * bpa));

    C = pow_F ( t, 0.9f ) * J * pow1;

    Q = wh * J;
    J *= J * 100.0f;
    M = C * pfl;
    Q = (Q == 0.f ? 0.0001f : Q); // avoid division by zero
    s = 100.0f * sqrtf ( M / Q );
    h = (myh * 180.f) / (float)rtengine::RT_PI;
}
#ifdef __SSE2__
void Ciecam02::xyz2jchqms_ciecam02float ( vfloat &J, vfloat &C, vfloat &h, vfloat &Q, vfloat &M, vfloat &s, vfloat aw, vfloat fl, vfloat wh,
        vfloat x, vfloat y, vfloat z, vfloat xw, vfloat yw, vfloat zw,
        vfloat c, vfloat nc, vfloat pow1, vfloat nbb, vfloat ncb, vfloat pfl, vfloat cz, vfloat d)

{
    vfloat r, g, b;
    vfloat rw, gw, bw;
    vfloat rc, gc, bc;
    vfloat rp, gp, bp;
    vfloat rpa, gpa, bpa;
    vfloat a, ca, cb;
    vfloat e, t;

    xyz_to_cat02float ( r, g, b, x, y, z);
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw);
    vfloat onev = F2V (1.f);
    rc = r * (((yw * d) / rw) + (onev - d));
    gc = g * (((yw * d) / gw) + (onev - d));
    bc = b * (((yw * d) / bw) + (onev - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc);
    //gamut correction M.H.Brill S.Susstrunk
    rp = _mm_max_ps (rp, ZEROV);
    gp = _mm_max_ps (gp, ZEROV);
    bp = _mm_max_ps (bp, ZEROV);
    rpa = nonlinear_adaptationfloat ( rp, fl );
    gpa = nonlinear_adaptationfloat ( gp, fl );
    bpa = nonlinear_adaptationfloat ( bp, fl );

    ca = rpa - ((F2V (12.0f) * gpa) - bpa) / F2V (11.0f);
    cb = F2V (0.11111111f) * (rpa + gpa - (bpa + bpa));

    vfloat myh = xatan2f ( cb, ca );
    vfloat temp = F2V (rtengine::RT_PI);
    temp += temp;
    temp += myh;
    myh = vself (vmaskf_lt (myh, ZEROV), temp, myh);

    a = ((rpa + rpa) + gpa + (F2V (0.05f) * bpa) - F2V (0.305f)) * nbb;
    a = _mm_max_ps (a, ZEROV);  //gamut correction M.H.Brill S.Susstrunk

    J = pow_F ( a / aw, c * cz * F2V (0.5f));

    e = ((F2V (961.53846f)) * nc * ncb) * (xcosf ( myh + F2V (2.0f) ) + F2V (3.8f));
    t = (e * _mm_sqrt_ps ( (ca * ca) + (cb * cb) )) / (rpa + gpa + (F2V (1.05f) * bpa));

    C = pow_F ( t, F2V (0.9f) ) * J * pow1;

    Q = wh * J;
    J *= J * F2V (100.0f);
    M = C * pfl;
    Q = _mm_max_ps (Q, F2V (0.0001f)); // avoid division by zero
    s = F2V (100.0f) * _mm_sqrt_ps ( M / Q );
    h = (myh * F2V (180.f)) / F2V (rtengine::RT_PI);
}
#endif

void Ciecam02::xyz2jch_ciecam02float ( float &J, float &C, float &h, float aw, float fl,
                                       float x, float y, float z, float xw, float yw, float zw,
                                       float c, float nc, float pow1, float nbb, float ncb, float cz, float d)

{
    float r, g, b;
    float rw, gw, bw;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    float a, ca, cb;
    float e, t;
    float myh;
    int gamu = 1;
    xyz_to_cat02float ( r, g, b, x, y, z, gamu );
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw, gamu );
    rc = r * (((yw * d) / rw) + (1.f - d));
    gc = g * (((yw * d) / gw) + (1.f - d));
    bc = b * (((yw * d) / bw) + (1.f - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, gamu );

    if (gamu == 1) { //gamut correction M.H.Brill S.Susstrunk
        rp = MAXR (rp, 0.0f);
        gp = MAXR (gp, 0.0f);
        bp = MAXR (bp, 0.0f);
    }

#ifdef __SSE2__
    vfloat pv = _mm_setr_ps(rp, gp, bp, 1.f);
    vfloat fv = F2V(fl);
    vfloat outv = nonlinear_adaptationfloat(pv, fv);
    rpa = outv[0];
    gpa = outv[1];
    bpa = outv[2];
#else
    rpa = nonlinear_adaptationfloat(rp, fl);
    gpa = nonlinear_adaptationfloat(gp, fl);
    bpa = nonlinear_adaptationfloat(bp, fl);
#endif

    ca = rpa - ((12.0f * gpa) - bpa) / 11.0f;
    cb = (0.11111111f) * (rpa + gpa - (2.0f * bpa));

    myh = xatan2f ( cb, ca );

    if ( myh < 0.0f ) {
        myh += (2.f * rtengine::RT_PI);
    }

    a = ((2.0f * rpa) + gpa + (0.05f * bpa) - 0.305f) * nbb;

    if (gamu == 1) {
        a = MAXR (a, 0.0f); //gamut correction M.H.Brill S.Susstrunk
    }

    J = pow_F ( a / aw, c * cz * 0.5f);

    e = ((961.53846f) * nc * ncb) * (xcosf ( myh + 2.0f ) + 3.8f);
    t = (e * sqrtf ( (ca * ca) + (cb * cb) )) / (rpa + gpa + (1.05f * bpa));

    C = pow_F ( t, 0.9f ) * J * pow1;

    J *= J * 100.0f;
    h = (myh * 180.f) / (float)rtengine::RT_PI;
}

void Ciecam02::jch2xyz_ciecam02float ( float &x, float &y, float &z, float J, float C, float h,
                                       float xw, float yw, float zw,
                                       float c, float nc, int gamu, float pow1, float nbb, float ncb, float fl, float cz, float d, float aw)
{
    float r, g, b;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    float rw, gw, bw;
    float a, ca, cb;
    float e, t;
    gamu = 1;
    xyz_to_cat02float(rw, gw, bw, xw, yw, zw, gamu);
    e = ((961.53846f) * nc * ncb) * (xcosf(h * rtengine::RT_PI_F_180 + 2.0f) + 3.8f);

#ifdef __SSE2__
    vfloat powinv1 = _mm_setr_ps(J / 100.0f, 10.f * C / (sqrtf(J) * pow1), 1.f, 1.f);
    vfloat powinv2 = _mm_setr_ps(1.0f / (c * cz), 1.1111111f, 1.f, 1.f);
    vfloat powoutv = pow_F(powinv1, powinv2);
    a = powoutv[0] * aw;
    t = powoutv[1];
#else
    a = pow_F(J / 100.0f, 1.0f / (c * cz)) * aw;
    t = pow_F(10.f * C / (sqrtf(J) * pow1), 1.1111111f);
#endif

    calculate_abfloat(ca, cb, h, e, t, nbb, a);
    Aab_to_rgbfloat(rpa, gpa, bpa, a, ca, cb, nbb);

#ifdef __SSE2__
    vfloat pav = _mm_setr_ps(rpa, gpa, bpa, 1.f);
    vfloat fv = F2V(fl);
    vfloat outv = inverse_nonlinear_adaptationfloat(pav, fv);
    rp = outv[0];
    gp = outv[1];
    bp = outv[2];
#else
    rp = inverse_nonlinear_adaptationfloat(rpa, fl);
    gp = inverse_nonlinear_adaptationfloat(gpa, fl);
    bp = inverse_nonlinear_adaptationfloat(bpa, fl);
#endif
    hpe_to_xyzfloat(x, y, z, rp, gp, bp);
    xyz_to_cat02float(rc, gc, bc, x, y, z, gamu);

    r = rc / (((yw * d) / rw) + (1.0f - d));
    g = gc / (((yw * d) / gw) + (1.0f - d));
    b = bc / (((yw * d) / bw) + (1.0f - d));

    cat02_to_xyzfloat(x, y, z, r, g, b, gamu);
}

#ifdef __SSE2__
void Ciecam02::jch2xyz_ciecam02float ( vfloat &x, vfloat &y, vfloat &z, vfloat J, vfloat C, vfloat h,
                                       vfloat xw, vfloat yw, vfloat zw,
                                       vfloat nc, vfloat pow1, vfloat nbb, vfloat ncb, vfloat fl, vfloat d, vfloat aw, vfloat reccmcz)
{
    vfloat r, g, b;
    vfloat rc, gc, bc;
    vfloat rp, gp, bp;
    vfloat rpa, gpa, bpa;
    vfloat rw, gw, bw;
    vfloat a, ca, cb;
    vfloat e, t;
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw);
    e = ((F2V (961.53846f)) * nc * ncb) * (xcosf ( ((h * F2V (rtengine::RT_PI)) / F2V (180.0f)) + F2V (2.0f) ) + F2V (3.8f));
    a = pow_F ( J / F2V (100.0f), reccmcz ) * aw;
    t = pow_F ( F2V (10.f) * C / (_mm_sqrt_ps ( J ) * pow1), F2V (1.1111111f) );

    calculate_abfloat ( ca, cb, h, e, t, nbb, a );
    Aab_to_rgbfloat ( rpa, gpa, bpa, a, ca, cb, nbb );

    rp = inverse_nonlinear_adaptationfloat ( rpa, fl );
    gp = inverse_nonlinear_adaptationfloat ( gpa, fl );
    bp = inverse_nonlinear_adaptationfloat ( bpa, fl );

    hpe_to_xyzfloat ( x, y, z, rp, gp, bp );
    xyz_to_cat02float ( rc, gc, bc, x, y, z );

    r = rc / (((yw * d) / rw) + (F2V (1.0f) - d));
    g = gc / (((yw * d) / gw) + (F2V (1.0f) - d));
    b = bc / (((yw * d) / bw) + (F2V (1.0f) - d));

    cat02_to_xyzfloat ( x, y, z, r, g, b );
}
#endif

float Ciecam02::nonlinear_adaptationfloat ( float c, float fl )
{
    float p;

    if (c < 0.0f) {
        p = pow_F ( (-1.0f * fl * c) / 100.0f, 0.42f );
        return ((-1.0f * 400.0f * p) / (27.13f + p)) + 0.1f;
    } else {
        p = pow_F ( (fl * c) / 100.0f, 0.42f );
        return ((400.0f * p) / (27.13f + p)) + 0.1f;
    }
}

#ifdef __SSE2__
vfloat Ciecam02::nonlinear_adaptationfloat ( vfloat c, vfloat fl )
{
    vfloat c100 = F2V (100.f);
    vfloat czd42 = F2V (0.42f);
    vfloat c400 = vmulsignf (F2V (400.f), c);
    fl = vmulsignf (fl, c);
    vfloat p = pow_F ( (fl * c) / c100, czd42 );
    vfloat c27d13 = F2V (27.13);
    vfloat czd1 = F2V (0.1f);
    return ((c400 * p) / (c27d13 + p)) + czd1;
}
#endif

float Ciecam02::inverse_nonlinear_adaptationfloat ( float c, float fl )
{
    c -= 0.1f;

    if (c < 0.f) {
        fl *= -1.f;

        if (c < -399.99f) { // avoid nan values
            c = -399.99f;
        }
    } else if (c > 399.99f) { // avoid nan values
        c = 399.99f;
    }

    return (100.0f / fl) * pow_F ( (27.13f * fabsf ( c )) / (400.0f - fabsf ( c )), 2.38095238f );
}

#ifdef __SSE2__
vfloat Ciecam02::inverse_nonlinear_adaptationfloat ( vfloat c, vfloat fl )
{
    c -= F2V (0.1f);
    fl = vmulsignf (fl, c);
    c = vabsf (c);
    c = _mm_min_ps ( c, F2V (399.99f));
    return (F2V (100.0f) / fl) * pow_F ( (F2V (27.13f) * c) / (F2V (400.0f) - c), F2V (2.38095238f) );
}
#endif
//end CIECAM Billy Bigg

}
