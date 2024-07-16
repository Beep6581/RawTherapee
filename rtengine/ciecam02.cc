/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Changes in Ciecam02 with Ciecam16 Jacques Desmis jdesmis@gmail.com   12/2020 
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
#include "ciecam02.h"
#include "rt_math.h"
#include "curves.h"
#include <math.h>
#include "sleef.h"

#undef CLIPD
#define CLIPD(a) ((a)>0.f?((a)<1.f?(a):1.f):0.f)
#define MAXR(a,b) ((a) > (b) ? (a) : (b))
#define Jzazbz_b 1.15
#define Jzazbz_g 0.66
#define Jzazbz_c1 (3424/4096.0)
#define Jzazbz_c2 (2413/128.0)
#define Jzazbz_c3 (2392/128.0)
#define Jzazbz_n (2610/16384.0)
#define Jzazbz_p (1.7*2523/32.0)
#define Jzazbz_d (-0.56)
#define Jzazbz_d0 (1.6295499532821566e-11)
#define Jzazbz_ni (16384.0/2610.0)
#define Jzazbz_pi (32.0/4289.1)  //4289.1 = 2523 * 1.7



namespace rtengine
{

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

void Ciecam02::curveJfloat (float br, float contr, float thr, const LUTu & histogram, LUTf & outCurve)
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
            brightcurvePoints[3] = max(0.0, 0.1 - (double) br / 150.0); // toe point
        //    brightcurvePoints[3] = 0.1f - br / 150.0f; // toe point
            brightcurvePoints[4] = 0.1f; // value at toe point

        //    brightcurvePoints[5] = min (1.0f, 0.7f - br / 300.0f); // shoulder point
            brightcurvePoints[5] = 0.7f - br / 300.0f; // shoulder point
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
       // printf("avg=%f \n", (double) avg);
        float thrmin = (thr - contr / 250.0f);
        float thrmax = (thr + contr / 250.0f);
        
        std::vector<double> contrastcurvePoints (9);

        contrastcurvePoints[0] = double (DCT_NURBS);

        contrastcurvePoints[1] = 0.f; // black point.  Value in [0 ; 1] range
        contrastcurvePoints[2] = 0.f; // black point.  Value in [0 ; 1] range

        contrastcurvePoints[3] = avg - avg * thrmin; // toe point
        contrastcurvePoints[4] = avg - avg * thrmax;// value at toe point

        contrastcurvePoints[5] = avg + (1.f - avg) * thrmin; // shoulder point
        contrastcurvePoints[6] = avg + (1.f - avg) * thrmax; // value at shoulder point

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

float Ciecam02::achromatic_response_to_whitefloat ( float x, float y, float z, float d, float fl, float nbb, int c16, float plum)
{
    float r, g, b;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
//    gamu = 1;
    xyz_to_cat02float ( r, g, b, x, y, z, c16, plum);

    rc = r * (((y * d) / r) + (1.0f - d));
    gc = g * (((y * d) / g) + (1.0f - d));
    bc = b * (((y * d) / b) + (1.0f - d));


    if(c16 == 1) {
        cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, c16);
        rp = MAXR (rp, 0.0f);
        gp = MAXR (gp, 0.0f);
        bp = MAXR (bp, 0.0f);
    } else {
        rp = MAXR (rc, 0.0f);
        gp = MAXR (gc, 0.0f);
        bp = MAXR (bc, 0.0f);
    }

    rpa = nonlinear_adaptationfloat ( rp, fl );
    gpa = nonlinear_adaptationfloat ( gp, fl );
    bpa = nonlinear_adaptationfloat ( bp, fl );

    return ((2.0f * rpa) + gpa + (0.05f * bpa) - 0.305f) * nbb;
}

void Ciecam02::xyz_to_cat02float ( float &r, float &g, float &b, float x, float y, float z, int c16, float plum)
{   //I use isnan() because I have tested others solutions with std::max(xxx,0) and in some cases crash
    //original cat02
    //r = ( 0.7328 * x) + (0.4296 * y) - (0.1624 * z);
    //g = (-0.7036 * x) + (1.6975 * y) + (0.0061 * z);
    //b = ( 0.0000 * x) + (0.0000 * y) + (1.0000 * z);
    float peakLum = 1.f/ plum;
    if(c16 == 1) {//cat02
        r = ( 1.007245f * x) + (0.011136f * y) - (0.018381f * z); //Changjun Li
        g = (-0.318061f * x) + (1.314589f * y) + (0.003471f * z);
        b = ( 0.0000f * x) + (0.0000f * y) + (1.0000f * z);
    } else if (c16 == 16) {//cat16
        r = ( 0.401288f * x) + (0.650173f * y) - (0.051461f * z); //cat16
        g = (-0.250268f * x) + (1.204414f * y) + (0.045854f * z);
        b = ( -0.002079f * x) + (0.048952f * y) + (0.953127f * z);
    } else if (c16 == 21) {//cam16 PQ
        float rp = ( 0.401288f * x) + (0.650173f * y) - (0.051461f * z); //cat16
        float gp = (-0.250268f * x) + (1.204414f * y) + (0.045854f * z);
        float bp = ( -0.002079f * x) + (0.048952f * y) + (0.953127f * z);
        rp *= 0.01f;
        gp *= 0.01f;
        bp *= 0.01f;
        float tmp = pow_F(rp * peakLum, Jzazbz_n);
        if(std::isnan(tmp)) {//to avoid crash
            tmp = 0.f;
        }
        r = 100.f * pow((Jzazbz_c1 + Jzazbz_c2 * tmp) / (1. + Jzazbz_c3 * tmp), Jzazbz_p);
        if(std::isnan(r) || r < 0.f) {
            r = 0.f;
        }

        tmp = pow_F(gp * peakLum, Jzazbz_n);
        if(std::isnan(tmp)) {//to avoid crash
            tmp = 0.f;
        }
        g = 100.f * pow((Jzazbz_c1 + Jzazbz_c2 * tmp) / (1. + Jzazbz_c3 * tmp), Jzazbz_p);
        if(std::isnan(g) || g < 0.f) {
            g = 0.f;
        }
        tmp = pow_F(bp * peakLum, Jzazbz_n);
        if(std::isnan(tmp)) {//to avoid crash
            tmp = 0.f;
        }
        
        b = 100.f * pow((Jzazbz_c1 + Jzazbz_c2 * tmp) / (1. + Jzazbz_c3 * tmp), Jzazbz_p);
        if(std::isnan(b) || b < 0.f) {
            b = 0.f;
        }
    }

}

#ifdef __SSE2__
void Ciecam02::xyz_to_cat02float ( vfloat &r, vfloat &g, vfloat &b, vfloat x, vfloat y, vfloat z, int c16, vfloat plum)
{   //I use isnan() because I have tested others solutions with std::max(xxx,0) and in some cases crash
    //gamut correction M.H.Brill S.Susstrunk
    if(c16 == 1) {
        r = ( F2V (1.007245f) * x) + (F2V (0.011136f) * y) - (F2V (0.018381f) * z); //Changjun Li
        g = (F2V (-0.318061f) * x) + (F2V (1.314589f) * y) + (F2V (0.003471f) * z);
        b = z;
    } else  if (c16 == 16) {
    //cat16
        r = ( F2V (0.401288f) * x) + (F2V (0.650173f) * y) - (F2V (0.051461f) * z);
        g = -(F2V (0.250268f) * x) + (F2V (1.204414f) * y) + (F2V (0.045854f) * z);
        b = -(F2V(0.002079f) * x) + (F2V(0.048952f) * y) + (F2V(0.953127f) * z);
    } else  if (c16 == 21) {
        vfloat rp = ( F2V (0.401288f) * x) + (F2V (0.650173f) * y) - (F2V (0.051461f) * z);
        vfloat gp = -(F2V (0.250268f) * x) + (F2V (1.204414f) * y) + (F2V (0.045854f) * z);
        vfloat bp = -(F2V(0.002079f) * x) + (F2V(0.048952f) * y) + (F2V(0.953127f) * z);
        vfloat Jzazbz_c1v = F2V(Jzazbz_c1);
        vfloat Jzazbz_c2v = F2V(Jzazbz_c2);
        vfloat Jzazbz_nv = F2V(Jzazbz_n);
        vfloat Jzazbz_c3v = F2V(Jzazbz_c3);
        vfloat Jzazbz_pv = F2V(Jzazbz_p);
        vfloat mulone = F2V(0.01f);
        vfloat mulhund = F2V(100.f);
        float RR, GG, BB;
        vfloat one = F2V(1.);
        vfloat peakLumv = one / plum;
        rp *= mulone;
        gp *= mulone;
        bp *= mulone;
        vfloat tmp = pow_F(rp * peakLumv, Jzazbz_nv );
        STVF(RR, tmp);
        if(std::isnan(RR)) {//to avoid crash
            tmp = F2V(0.f);;
        }
        r = mulhund * pow_F((Jzazbz_c1v + Jzazbz_c2v * tmp) / (one + Jzazbz_c3v * tmp), Jzazbz_pv);
            STVF(RR, r);
            if(std::isnan(RR) || RR < 0.f) {//to avoid crash
                r = F2V(0.f);;
            }
        tmp = pow_F(gp * peakLumv, Jzazbz_nv );
        STVF(RR, tmp);
        if(std::isnan(RR)) {//to avoid crash
            tmp = F2V(0.f);;
        }

        g = mulhund * pow_F((Jzazbz_c1v + Jzazbz_c2v * tmp) / (one + Jzazbz_c3v * tmp), Jzazbz_pv);
            STVF(GG, g);
            if(std::isnan(GG) || GG < 0.f) {//to avoid crash
                g = F2V(0.f);;
            }

        tmp = pow_F(bp * peakLumv, Jzazbz_nv );
        STVF(RR, tmp);
        if(std::isnan(RR)) {//to avoid crash
            tmp = F2V(0.f);;
        }

        b = mulhund * pow_F((Jzazbz_c1v + Jzazbz_c2v * tmp) / (one + Jzazbz_c3v * tmp), Jzazbz_pv);
            STVF(BB, b);
            if(std::isnan(BB) || BB < 0.f) {//to avoid crash
                b = F2V(0.f);;
            }
        
    }
}
#endif

void Ciecam02::cat02_to_xyzfloat ( float &x, float &y, float &z, float r, float g, float b, int c16, float plum)
{   //I use isnan() because I have tested others solutions with std::max(xxx,0) and in some cases crash
    //original cat02
    //x = ( 1.0978566 * r) - (0.277843 * g) + (0.179987 * b);
    //y = ( 0.455053 * r) + (0.473938 * g) + (0.0710096* b);
    //z = ( 0.000000 * r) - (0.000000 * g) + (1.000000 * b);
    float pl = plum;
    if(c16 == 1) {
        x = ( 0.99015849f * r) - (0.00838772f * g) + (0.018229217f * b); //Changjun Li
        y = ( 0.239565979f * r) + (0.758664642f * g) + (0.001770137f * b);
        z = ( 0.000000f * r) - (0.000000f * g) + (1.000000f * b);
    } else if(c16 == 16){//cat16
        x = ( 1.86206786f * r) - (1.01125463f * g) + (0.14918677f * b); //Cat16
        y = ( 0.38752654f * r) + (0.62144744f * g) + (-0.00897398f * b);
        z = ( -0.0158415f * r) - (0.03412294f * g) + (1.04996444f * b);
    }else if(c16 == 21){//cam16 PQ
        float lp = ( 1.86206786f * r) - (1.01125463f * g) + (0.14918677f * b); //Cat16
        float mp = ( 0.38752654f * r) + (0.62144744f * g) + (-0.00897398f * b);
        float sp = ( -0.0158415f * r) - (0.03412294f * g) + (1.04996444f * b);
        lp *= 0.01f;
        float tmp = pow_F(lp, Jzazbz_pi);
        if(std::isnan(tmp)) {//to avoid crash
            tmp = 0.f;
        }

        float prov = (Jzazbz_c1 - tmp) / ((Jzazbz_c3 * tmp) - Jzazbz_c2);
        x =  pl * pow_F(prov, Jzazbz_ni);
        if(std::isnan(x)) {//to avoid crash
            x = 0.f;
        }
        x *= 100.f;
        mp *= 0.01f;
        tmp = pow_F(mp, Jzazbz_pi);
        if(std::isnan(tmp)) {//to avoid crash
            tmp = 0.f;
        }
        prov = (Jzazbz_c1 - tmp) / ((Jzazbz_c3 * tmp) - Jzazbz_c2);
        y =  pl * pow_F(prov, Jzazbz_ni);
        if(std::isnan(y)) {
            y = 0.f;
        }
        y *= 100.f;
        sp *= 0.01f;
        tmp = pow_F(sp, Jzazbz_pi);
        if(std::isnan(tmp)) {//to avoid crash
            tmp = 0.f;
        }
        prov = (Jzazbz_c1 - tmp) / ((Jzazbz_c3 * tmp) - Jzazbz_c2);
        z =  pl * pow_F(prov, Jzazbz_ni);
        if(std::isnan(z)) {
            z = 0.;
        }
        z *= 100.f;
        
    }

}
#ifdef __SSE2__
void Ciecam02::cat02_to_xyzfloat ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b, int c16, vfloat plum )
{   //I use isnan() because I have tested others solutions with std::max(xxx,0) and in some cases crash
    vfloat plv = plum;
    //gamut correction M.H.Brill S.Susstrunk
    if(c16 == 1) {//cat02
        x = ( F2V (0.99015849f) * r) - (F2V (0.00838772f) * g) + (F2V (0.018229217f) * b); //Changjun Li
        y = ( F2V (0.239565979f) * r) + (F2V (0.758664642f) * g) + (F2V (0.001770137f) * b);
        z = b;
        } else  if(c16 == 16) {
        //cat16  
        x = ( F2V (1.86206786f) * r) - (F2V (1.01125463f) * g) + (F2V (0.14918677f) * b);
        y = ( F2V (0.38752654f) * r) + (F2V (0.621447744f) * g) - (F2V (0.00897398f) * b);
        z = -(F2V(0.0158415f) * r) - (F2V(0.03412294f) * g) + (F2V(1.04996444f) * b);
    }else if(c16 == 21){//cam16 PQ
            vfloat lp = ( F2V (1.86206786f) * r) - (F2V (1.01125463f) * g) + (F2V (0.14918677f) * b);
            vfloat mp = ( F2V (0.38752654f) * r) + (F2V (0.621447744f) * g) - (F2V (0.00897398f) * b);
            vfloat sp = -(F2V(0.0158415f) * r) - (F2V(0.03412294f) * g) + (F2V(1.04996444f) * b);
            float XX,YY,ZZ;
            vfloat Jzazbz_c1v = F2V(Jzazbz_c1);
            vfloat Jzazbz_c2v = F2V(Jzazbz_c2);
            vfloat Jzazbz_c3v = F2V(Jzazbz_c3);
            vfloat Jzazbz_piv = F2V(Jzazbz_pi);
            vfloat Jzazbz_niv = F2V(Jzazbz_ni);
            vfloat mulone = F2V(0.01f);
            vfloat mulhund = F2V(100.f);
            lp *= mulone;
            float pro;
            vfloat tmp = pow_F(lp, Jzazbz_piv);
            STVF(XX, tmp);
            if(std::isnan(XX)) {//to avoid crash
                tmp = F2V(0.f);;
            }
            vfloat prov = (Jzazbz_c1v - tmp) / ((Jzazbz_c3v * tmp) - Jzazbz_c2v);
            x =  plv * pow_F(prov, Jzazbz_niv);
            STVF(XX, x);
            if(std::isnan(XX)) {//to avoid crash
                x = F2V(0.f);;
            }
            x *= mulhund;
            mp *= mulone;
            tmp = pow_F(mp, Jzazbz_piv);
            STVF(YY, tmp);
            if(std::isnan(YY)) {//to avoid crash
                tmp = F2V(0.f);;
            }
            prov = (Jzazbz_c1v - tmp) / ((Jzazbz_c3v * tmp) - Jzazbz_c2v);
            y = plv * pow_F(prov, Jzazbz_niv);
            STVF(YY, y);
            if(std::isnan(YY)) {//to avoid crash
                y = F2V(0.f);;
            }
            y *= mulhund;
            sp *= mulone;
            tmp = pow_F(sp, Jzazbz_piv);
            STVF(ZZ, tmp);
            if(std::isnan(ZZ)) {//to avoid crash
                tmp = F2V(0.f);;
            }
            prov = (Jzazbz_c1v - tmp) / ((Jzazbz_c3v * tmp) - Jzazbz_c2v);
            STVF(pro, prov);
            z = plv * pow_F(prov, Jzazbz_niv);
            STVF(ZZ, z);
            if(std::isnan(ZZ)) {//to avoid crash
                z = F2V(0.f);;
            }
            z *= mulhund;
            
    }

}
#endif

void Ciecam02::hpe_to_xyzfloat ( float &x, float &y, float &z, float r, float g, float b, int c16)
{
        x = (1.910197f * r) - (1.112124f * g) + (0.201908f * b);
        y = (0.370950f * r) + (0.629054f * g) - (0.000008f * b);
        z = b;
}
#ifdef __SSE2__
void Ciecam02::hpe_to_xyzfloat ( vfloat &x, vfloat &y, vfloat &z, vfloat r, vfloat g, vfloat b, int c16)
{
        x = (F2V (1.910197f) * r) - (F2V (1.112124f) * g) + (F2V (0.201908f) * b);
        y = (F2V (0.370950f) * r) + (F2V (0.629054f) * g) - (F2V (0.000008f) * b);
        z = b;
}
#endif

void Ciecam02::cat02_to_hpefloat ( float &rh, float &gh, float &bh, float r, float g, float b, int c16)
{
// original cat02
//        rh = ( 0.7409792f * r) + (0.2180250f * g) + (0.0410058f * b);
//        gh = ( 0.2853532f * r) + (0.6242014f * g) + (0.0904454f * b);
//        bh = (-0.0096280f * r) - (0.0056980f * g) + (1.0153260f * b);
    if(c16 == 1) {//cat02
        rh = ( 0.550930835f * r) + (0.519435987f * g) - ( 0.070356303f * b);
        gh = ( 0.055954056f * r) + (0.89973132f * g) + (0.044315524f * b);
        bh = (0.0f * r) - (0.0f * g) + (1.0f * b);
    } else {//cat16
        rh = ( 1.f * r) + (0.f * g) + ( 0.f * b);
        gh = ( 0.f * r) + (1.f * g) + (0.f * b);
        bh = (0.0f * r) + (0.0f * g) + (1.0f * b);
    }

}

#ifdef __SSE2__
void Ciecam02::cat02_to_hpefloat ( vfloat &rh, vfloat &gh, vfloat &bh, vfloat r, vfloat g, vfloat b, int c16)
{
    if(c16 == 1) {
    //Changjun Li
        rh = ( F2V (0.550930835f) * r) + (F2V (0.519435987f) * g) - ( F2V (0.070356303f) * b);
        gh = ( F2V (0.055954056f) * r) + (F2V (0.89973132f) * g) + (F2V (0.044315524f) * b);
        bh = b;
    } else {//cat16
        rh = ( F2V (1.f) * r) + (F2V (0.f) * g) + ( F2V (0.f) * b);
        gh = ( F2V (0.f) * r) + (F2V (1.f) * g) + (F2V (0.f) * b);
        bh = b;
        
    }
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

void Ciecam02::initcam1float (float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                              float &cz, float &aw, float &wh, float &pfl, float &fl, float c, int c16, float plum)
{
    n = yb / yw;

    if (pilotd == 2.f) {
        d = d_factorfloat ( f, la );
    } else {
        d = pilotd;
    }

    fl = calculate_fl_from_la_ciecam02float ( la );
    nbb = ncb = 0.725f * pow_F ( 1.0f / n, 0.2f );
    cz = 1.48f + sqrt ( n );
    aw = achromatic_response_to_whitefloat ( xw, yw, zw, d, fl, nbb, c16, plum);
    wh = ( 4.0f / c ) * ( aw + 4.0f ) * pow_F ( fl, 0.25f );
    pfl = pow_F ( fl, 0.25f );
}

void Ciecam02::initcam2float (float yb, float pilotd, float f, float la, float xw, float yw, float zw, float &n, float &d, float &nbb, float &ncb,
                              float &cz, float &aw, float &fl, int c16, float plum)
{
    n = yb / yw;

    if (pilotd == 2.f) {
        d = d_factorfloat ( f, la );
    } else {
        d = pilotd;
    }

//   d = d_factorfloat( f, la );
    fl = calculate_fl_from_la_ciecam02float ( la );
    nbb = ncb = 0.725f * pow_F ( 1.0f / n, 0.2f );
    cz = 1.48f + sqrt ( n );
    aw = achromatic_response_to_whitefloat ( xw, yw, zw, d, fl, nbb, c16, plum);
}


void Ciecam02::xyz2jzczhz ( double &Jz, double &az, double &bz, double x, double y, double z, double pl, double &Lp, double &Mp, double &Sp, bool zcam)
{   //from various web 
    double Xp, Yp, Zp, L, M, S, Iz;
    double peakLum = 1. / pl;
    //I change 10000 for peaklum function of la (absolute luminance)- default 10000
    Xp = Jzazbz_b * x - ((Jzazbz_b - 1.) * z);
    Yp = Jzazbz_g * y - ((Jzazbz_g - 1.) * x);
    Zp = z;

    L = 0.41478972 * Xp + 0.579999 * Yp + 0.0146480 * Zp;
    M = -0.2015100 * Xp + 1.120649 * Yp + 0.0531008 * Zp;
    S = -0.0166008 * Xp + 0.264800 * Yp + 0.6684799 * Zp;
 
 //I use isnan() because I have tested others solutions with std::max(xxx,0) and in some cases crash 
 //   Lp = pow((Jzazbz_c1 + Jzazbz_c2 * pow(std::max((L * peakLum), 0.), Jzazbz_n)) / (1. + Jzazbz_c3 * pow((L * peakLum), Jzazbz_n)), Jzazbz_p);
 //   Mp = pow((Jzazbz_c1 + Jzazbz_c2 * pow(std::max((M * peakLum),0.), Jzazbz_n)) / (1. + Jzazbz_c3 * pow((M * peakLum), Jzazbz_n)), Jzazbz_p);
 //   Sp = pow((Jzazbz_c1 + Jzazbz_c2 * pow(std::max((S * peakLum), 0.), Jzazbz_n)) / (1. + Jzazbz_c3 * pow((S * peakLum), Jzazbz_n)), Jzazbz_p);
    double temp = pow(L * peakLum, Jzazbz_n);
    if(std::isnan(temp)) {//to avoid crash
        temp = 0.;
    }
    Lp = pow((Jzazbz_c1 + Jzazbz_c2 * temp) / (1. + Jzazbz_c3 * temp), Jzazbz_p);
    if(std::isnan(Lp)) {//to avoid crash
        Lp = 0.;
    }
    
    temp = pow(M * peakLum, Jzazbz_n);
    if(std::isnan(temp)) {//to avoid crash
        temp = 0.;
    }
    Mp = pow((Jzazbz_c1 + Jzazbz_c2 * temp) / (1. + Jzazbz_c3 * temp), Jzazbz_p);
    if(std::isnan(Mp)) {//to avoid crash
        Mp = 0.;
    }

    temp = pow(S * peakLum, Jzazbz_n);
    if(std::isnan(temp)) {//to avoid crash
        temp = 0.;
    }
    Sp = pow((Jzazbz_c1 + Jzazbz_c2 * temp) / (1. + Jzazbz_c3 * temp), Jzazbz_p);
    if(std::isnan(Sp)) {//to avoid crash
        Sp = 0.;
    }

    Iz = 0.5 * Lp + 0.5 * Mp;
    az = 3.524000 * Lp - 4.066708 * Mp + 0.542708 * Sp;
    bz = 0.199076 * Lp + 1.096799 * Mp - 1.295875 * Sp;
    if(!zcam) {
        Jz = (((1. + Jzazbz_d) * Iz) / (1. + Jzazbz_d * Iz)) - Jzazbz_d0;
      //  Jz = std::max((((1. + Jzazbz_d) * Iz) / (1. + Jzazbz_d * Iz)) - Jzazbz_d0, 0.);
    } else {
    //or if we use ZCAM Jz = Mp  - Jzazbz_d0
        Jz = Mp  - Jzazbz_d0;
    }
}


void Ciecam02::jzczhzxyz (double &x, double &y, double &z, double jz, double az, double bz, double pl, double &L, double &M, double &S, bool zcam)
{ //from various web
  //I use isnan() because I have tested others solutions with std::max(xxx,0) and in some cases crash 

    double Xp, Yp, Zp, Lp, Mp, Sp, Iz, tmp;

    if(!zcam) {
      //  Iz = std::max((jz + Jzazbz_d0) / (1. + Jzazbz_d - Jzazbz_d * (jz + Jzazbz_d0)), 0.);
        Iz = (jz + Jzazbz_d0) / (1. + Jzazbz_d - Jzazbz_d * (jz + Jzazbz_d0));
    } else {
    //or if we use ZCAM Iz = Jz  + Jzazbz_d0
        Iz = jz  + Jzazbz_d0;  
    }

    Lp = Iz + 0.138605043271539 * az + 0.0580473161561189 * bz;
    Mp = Iz - 0.138605043271539 * az - 0.0580473161561189 * bz;
    Sp = Iz - 0.0960192420263189 * az - 0.811891896056039 * bz;
    //I change optionally 10000 for pl function of la(absolute luminance) default 10000
   
    tmp = pow(Lp, Jzazbz_pi);
    if(std::isnan(tmp)) {//to avoid crash
        tmp = 0.;
    }
    L = pl * pow((Jzazbz_c1 - tmp) / ((Jzazbz_c3 * tmp) - Jzazbz_c2), Jzazbz_ni);
    if(std::isnan(L)) {//to avoid crash
        L = 0.;
    }

    tmp = pow(Mp, Jzazbz_pi);
    if(std::isnan(tmp)) {//to avoid crash
        tmp = 0.;
    }
    M = pl * pow((Jzazbz_c1 - tmp) / ((Jzazbz_c3 * tmp) - Jzazbz_c2), Jzazbz_ni);
    if(std::isnan(M)) {//to avoid crash
        M = 0.;
    }

    tmp = pow(Sp, Jzazbz_pi);
    if(std::isnan(tmp)) {//to avoid crash
        tmp = 0.;
    }
    S = pl * pow((Jzazbz_c1 - tmp) / ((Jzazbz_c3 * tmp) - Jzazbz_c2), Jzazbz_ni);
    if(std::isnan(S)) {//to avoid crash
        S = 0.;
    }

    Xp = 1.9242264357876067 * L - 1.0047923125953657 * M + 0.0376514040306180 * S;
    Yp = 0.3503167620949991 * L + 0.7264811939316552 * M - 0.0653844229480850 * S;
    Zp = -0.0909828109828475 * L - 0.3127282905230739 * M + 1.5227665613052603 * S;

    x = (Xp + (Jzazbz_b - 1.) * Zp) / Jzazbz_b;

    if(std::isnan(x)) {//to avoid crash
        x = 0.;
    }
    y = (Yp + (Jzazbz_g - 1.) * x) / Jzazbz_g;
    if(std::isnan(y)) {
        y = 0.;
    }
    z = Zp;
    if(std::isnan(z)) {
        z = 0.;
    }
}

void Ciecam02::xyz2jchqms_ciecam02float ( float &J, float &C, float &h, float &Q, float &M, float &s, float aw, float fl, float wh,
        float x, float y, float z, float xw, float yw, float zw,
        float c, float nc, float pow1, float nbb, float ncb, float pfl, float cz, float d, int c16, float plum)

{
    float r, g, b;
    float rw, gw, bw;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    float a, ca, cb;
    float e, t;
    float myh;
    xyz_to_cat02float ( r, g, b, x, y, z, c16, plum);
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw, c16, plum);
    rc = r * (((yw * d) / rw) + (1.f - d));
    gc = g * (((yw * d) / gw) + (1.f - d));
    bc = b * (((yw * d) / bw) + (1.f - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, c16);

    if(c16 == 1) {//cat02
        rp = MAXR (rp, 0.0f);
        gp = MAXR (gp, 0.0f);
        bp = MAXR (bp, 0.0f);
    } else {//cat16
        rp = MAXR (rc, 0.0f);
        gp = MAXR (gc, 0.0f);
        bp = MAXR (bc, 0.0f);
    }

    rpa = nonlinear_adaptationfloat ( rp, fl );
    gpa = nonlinear_adaptationfloat ( gp, fl );
    bpa = nonlinear_adaptationfloat ( bp, fl );

    ca = rpa - ((12.0f * gpa) - bpa) / 11.0f;
    cb = (0.11111111f) * (rpa + gpa - (2.0f * bpa));

    myh = xatan2f ( cb, ca );

    if ( myh < 0.0f ) {
        myh += (2.f * rtengine::RT_PI_F);
    }

    a = ((2.0f * rpa) + gpa + (0.05f * bpa) - 0.305f) * nbb;

    a = MAXR (a, 0.0f); //gamut correction M.H.Brill S.Susstrunk

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
        vfloat c, vfloat nc, vfloat pow1, vfloat nbb, vfloat ncb, vfloat pfl, vfloat cz, vfloat d, int c16, vfloat plum)

{
    vfloat r, g, b;
    vfloat rw, gw, bw;
    vfloat rc, gc, bc;
    vfloat rp, gp, bp;
    vfloat rpa, gpa, bpa;
    vfloat a, ca, cb;
    vfloat e, t;

    xyz_to_cat02float ( r, g, b, x, y, z, c16, plum);
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw, c16, plum);
    vfloat onev = F2V (1.f);
    rc = r * (((yw * d) / rw) + (onev - d));
    gc = g * (((yw * d) / gw) + (onev - d));
    bc = b * (((yw * d) / bw) + (onev - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, c16);

    //gamut correction M.H.Brill S.Susstrunk
    if(c16 == 1) {//cat02
        rp = vmaxf (rp, ZEROV);
        gp = vmaxf (gp, ZEROV);
        bp = vmaxf (bp, ZEROV);
    } else {//cat16
        rp = vmaxf (rc, ZEROV);
        gp = vmaxf (gc, ZEROV);
        bp = vmaxf (bc, ZEROV);
    }

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
    a = vmaxf (a, ZEROV);  //gamut correction M.H.Brill S.Susstrunk

    J = pow_F ( a / aw, c * cz * F2V (0.5f));

    e = ((F2V (961.53846f)) * nc * ncb) * (xcosf ( myh + F2V (2.0f) ) + F2V (3.8f));
    t = (e * vsqrtf ( (ca * ca) + (cb * cb) )) / (rpa + gpa + (F2V (1.05f) * bpa));

    C = pow_F ( t, F2V (0.9f) ) * J * pow1;

    Q = wh * J;
    J *= J * F2V (100.0f);
    M = C * pfl;
    Q = vmaxf (Q, F2V (0.0001f)); // avoid division by zero
    s = F2V (100.0f) * vsqrtf ( M / Q );
    h = (myh * F2V (180.f)) / F2V (rtengine::RT_PI);
}
#endif

void Ciecam02::xyz2jch_ciecam02float ( float &J, float &C, float &h, float aw, float fl,
                                       float x, float y, float z, float xw, float yw, float zw,
                                       float c, float nc, float pow1, float nbb, float ncb, float cz, float d, int c16, float plum)

{
    float r, g, b;
    float rw, gw, bw;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    float a, ca, cb;
    float e, t;
    float myh;
    xyz_to_cat02float ( r, g, b, x, y, z, c16, plum);
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw, c16, plum);
    rc = r * (((yw * d) / rw) + (1.f - d));
    gc = g * (((yw * d) / gw) + (1.f - d));
    bc = b * (((yw * d) / bw) + (1.f - d));

    cat02_to_hpefloat ( rp, gp, bp, rc, gc, bc, c16);

    if(c16 == 1) {//cat02
        rp = MAXR (rp, 0.0f);
        gp = MAXR (gp, 0.0f);
        bp = MAXR (bp, 0.0f);
    } else {//cat16
        rp = MAXR (rc, 0.0f);
        gp = MAXR (gc, 0.0f);
        bp = MAXR (bc, 0.0f);
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
        myh += (2.f * rtengine::RT_PI_F);
    }

    a = ((2.0f * rpa) + gpa + (0.05f * bpa) - 0.305f) * nbb;

    a = MAXR (a, 0.0f); //gamut correction M.H.Brill S.Susstrunk

    J = pow_F ( a / aw, c * cz * 0.5f);

    e = ((961.53846f) * nc * ncb) * (xcosf ( myh + 2.0f ) + 3.8f);
    t = (e * sqrtf ( (ca * ca) + (cb * cb) )) / (rpa + gpa + (1.05f * bpa));

    C = pow_F ( t, 0.9f ) * J * pow1;

    J *= J * 100.0f;
    h = (myh * 180.f) / (float)rtengine::RT_PI;
}

void Ciecam02::jch2xyz_ciecam02float ( float &x, float &y, float &z, float J, float C, float h,
                                       float xw, float yw, float zw,
                                       float c, float nc, float pow1, float nbb, float ncb, float fl, float cz, float d, float aw, int c16, float plum)
{
    float r, g, b;
    float rc, gc, bc;
    float rp, gp, bp;
    float rpa, gpa, bpa;
    float rw, gw, bw;
    float a, ca, cb;
    float e, t;
    xyz_to_cat02float(rw, gw, bw, xw, yw, zw, c16, plum);
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

    if(c16 == 1) {//cat02
        hpe_to_xyzfloat(x, y, z, rp, gp, bp, c16);
        xyz_to_cat02float(rc, gc, bc, x, y, z, c16, plum);

        r = rc / (((yw * d) / rw) + (1.0f - d));
        g = gc / (((yw * d) / gw) + (1.0f - d));
        b = bc / (((yw * d) / bw) + (1.0f - d));
    } else {//cat16
        r = rp / (((yw * d) / rw) + (1.0f - d));
        g = gp / (((yw * d) / gw) + (1.0f - d));
        b = bp / (((yw * d) / bw) + (1.0f - d));
    }

    cat02_to_xyzfloat(x, y, z, r, g, b, c16, plum);
}

#ifdef __SSE2__
void Ciecam02::jch2xyz_ciecam02float ( vfloat &x, vfloat &y, vfloat &z, vfloat J, vfloat C, vfloat h,
                                       vfloat xw, vfloat yw, vfloat zw,
                                       vfloat nc, vfloat pow1, vfloat nbb, vfloat ncb, vfloat fl, vfloat d, vfloat aw, vfloat reccmcz, int c16, vfloat plum)
{
    vfloat r, g, b;
    vfloat rc, gc, bc;
    vfloat rp, gp, bp;
    vfloat rpa, gpa, bpa;
    vfloat rw, gw, bw;
    vfloat a, ca, cb;
    vfloat e, t;
    xyz_to_cat02float ( rw, gw, bw, xw, yw, zw, c16, plum);
    e = ((F2V (961.53846f)) * nc * ncb) * (xcosf ( ((h * F2V (rtengine::RT_PI)) / F2V (180.0f)) + F2V (2.0f) ) + F2V (3.8f));
    a = pow_F ( J / F2V (100.0f), reccmcz ) * aw;
    t = pow_F ( F2V (10.f) * C / (vsqrtf ( J ) * pow1), F2V (1.1111111f) );

    calculate_abfloat ( ca, cb, h, e, t, nbb, a );
    Aab_to_rgbfloat ( rpa, gpa, bpa, a, ca, cb, nbb );

    rp = inverse_nonlinear_adaptationfloat ( rpa, fl );
    gp = inverse_nonlinear_adaptationfloat ( gpa, fl );
    bp = inverse_nonlinear_adaptationfloat ( bpa, fl );

    if(c16 == 1) {//cat02
        hpe_to_xyzfloat ( x, y, z, rp, gp, bp, c16);
        xyz_to_cat02float ( rc, gc, bc, x, y, z, c16, plum );

        r = rc / (((yw * d) / rw) + (F2V (1.0f) - d));
        g = gc / (((yw * d) / gw) + (F2V (1.0f) - d));
        b = bc / (((yw * d) / bw) + (F2V (1.0f) - d));
    } else {//cat16
        r = rp / (((yw * d) / rw) + (F2V (1.0f) - d));
        g = gp / (((yw * d) / gw) + (F2V (1.0f) - d));
        b = bp / (((yw * d) / bw) + (F2V (1.0f) - d));
    }

    cat02_to_xyzfloat ( x, y, z, r, g, b, c16, plum );
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
    c = vminf ( c, F2V (399.99f));
    return (F2V (100.0f) / fl) * pow_F ( (F2V (27.13f) * c) / (F2V (400.0f) - c), F2V (2.38095238f) );
}
#endif
}
