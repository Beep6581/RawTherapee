/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
 *  Small adaptation to Rawtherapee Locallab October 2019
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

/* film grain emulation.
 * Ported from darktable (src/iop/grain.c). Original copyright/license follows
 */
/*
    This file is part of darktable,
    copyright (c) 2010-2012 Henrik Andersson.
    adaptation to Rawtherapee 2021 Jacques Desmis jdesmis@gmail.com

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "imagefloat.h"
#include "improcfun.h"
#include "rt_math.h"


namespace rtengine {

namespace {

constexpr float GRAIN_LIGHTNESS_STRENGTH_SCALE = 0.15f;
constexpr float GRAIN_SCALE_FACTOR = 213.2f;

constexpr int GRAIN_LUT_SIZE = 128;
constexpr float GRAIN_LUT_DELTA_MAX = 2.0f;
constexpr float GRAIN_LUT_DELTA_MIN = 0.0001f;
constexpr float GRAIN_LUT_PAPER_GAMMA = 1.0f;


const int grad3[12][3] = { { 1, 1, 0 },
                           { -1, 1, 0 },
                           { 1, -1, 0 },
                           { -1, -1, 0 },
                           { 1, 0, 1 },
                           { -1, 0, 1 },
                           { 1, 0, -1 },
                           { -1, 0, -1 },
                           { 0, 1, 1 },
                           { 0, -1, 1 },
                           { 0, 1, -1 },
                           { 0, -1, -1 } };

const int permutation[]
    = { 151, 160, 137, 91,  90,  15,  131, 13,  201, 95,  96,  53,  194, 233, 7,   225, 140, 36,  103, 30,
        69,  142, 8,   99,  37,  240, 21,  10,  23,  190, 6,   148, 247, 120, 234, 75,  0,   26,  197, 62,
        94,  252, 219, 203, 117, 35,  11,  32,  57,  177, 33,  88,  237, 149, 56,  87,  174, 20,  125, 136,
        171, 168, 68,  175, 74,  165, 71,  134, 139, 48,  27,  166, 77,  146, 158, 231, 83,  111, 229, 122,
        60,  211, 133, 230, 220, 105, 92,  41,  55,  46,  245, 40,  244, 102, 143, 54,  65,  25,  63,  161,
        1,   216, 80,  73,  209, 76,  132, 187, 208, 89,  18,  169, 200, 196, 135, 130, 116, 188, 159, 86,
        164, 100, 109, 198, 173, 186, 3,   64,  52,  217, 226, 250, 124, 123, 5,   202, 38,  147, 118, 126,
        255, 82,  85,  212, 207, 206, 59,  227, 47,  16,  58,  17,  182, 189, 28,  42,  223, 183, 170, 213,
        119, 248, 152, 2,   44,  154, 163, 70,  221, 153, 101, 155, 167, 43,  172, 9,   129, 22,  39,  253,
        19,  98,  108, 110, 79,  113, 224, 232, 178, 185, 112, 104, 218, 246, 97,  228, 251, 34,  242, 193,
        238, 210, 144, 12,  191, 179, 162, 241, 81,  51,  145, 235, 249, 14,  239, 107, 49,  192, 214, 31,
        181, 199, 106, 157, 184, 84,  204, 176, 115, 121, 50,  45,  127, 4,   150, 254, 138, 236, 205, 93,
        222, 114, 67,  29,  24,  72,  243, 141, 128, 195, 78,  66,  215, 61,  156, 180 };


class GrainEvaluator {
public:
    GrainEvaluator(int offset_x, int offset_y, int full_width, int full_height, double scale, float divgr, int call, int fww, int fhh):
        ox(offset_x),
        oy(offset_y),
        fw(full_width),
        fh(full_height),
        scale(scale)
    {
        simplex_noise_init();
        constexpr float mb = 100.f;// * divgr;
        evaluate_grain_lut(mb, divgr);
    }
    
    void operator()(int isogr, int strengr, int scalegr, float divgr, Imagefloat *lab, bool multithread, int call, int fww, int fhh)
    {
        const double strength = (strengr / 100.0);
        const double octaves = 3.;
        const double wd = std::min(fw, fh);
        const double wdf = std::min(fww, fhh);
        
        const double zoom = (1.0 + 8 * (double(isogr) / GRAIN_SCALE_FACTOR) / 100.0) / 800.0;
        const double s = std::max(scale / 3.0, 1.0) / (double(std::max(scalegr, 1)) / 100.0);
        const int W = lab->getWidth();
        const int H = lab->getHeight();
        float **lab_L = lab->g.ptrs;
        double wddf = wd;
        if (call == 1 || call == 3) {
            wddf = wdf;
        }


#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int j = 0; j < H; ++j) {
            double wy = oy + j;
            double y = wy / wddf;
            for (int i = 0; i < W; ++i) {
                double wx = ox + i;
                double x = wx / wddf;
                double noise = simplex_2d_noise(x, y, octaves, zoom) / s;
                lab_L[j][i] += lut_lookup(noise * strength * GRAIN_LIGHTNESS_STRENGTH_SCALE, lab_L[j][i] / 32768.f);
            }
        }
    }

private:
    void simplex_noise_init()
    {
        for(int i = 0; i < 512; i++) perm[i] = permutation[i & 255];
    }

    double dot(const int *g, double x, double y, double z)
    {
        return g[0] * x + g[1] * y + g[2] * z;
    }

    float FASTFLOOR(float x)
    {
        return (x > 0 ? (int)(x) : (int)(x)-1);
    }

    double simplex_noise(double xin, double yin, double zin)
    {
        double n0, n1, n2, n3; // Noise contributions from the four corners
        // Skew the input space to determine which simplex cell we're in
        const double F3 = 1.0 / 3.0;
        const double s = (xin + yin + zin) * F3; // Very nice and simple skew factor for 3D
        const int i = FASTFLOOR(xin + s);
        const int j = FASTFLOOR(yin + s);
        const int k = FASTFLOOR(zin + s);
        const double G3 = 1.0 / 6.0; // Very nice and simple unskew factor, too
        const double t = (i + j + k) * G3;
        const double X0 = i - t; // Unskew the cell origin back to (x,y,z) space
        const double Y0 = j - t;
        const double Z0 = k - t;
        const double x0 = xin - X0; // The x,y,z distances from the cell origin
        const double y0 = yin - Y0;
        const double z0 = zin - Z0;
        // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
        // Determine which simplex we are in.
        int i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
        int i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
        if(x0 >= y0)
        {
            if(y0 >= z0)
            {
                i1 = 1; // X Y Z order
                j1 = 0;
                k1 = 0;
                i2 = 1;
                j2 = 1;
                k2 = 0;
            }
            else if(x0 >= z0)
            {
                i1 = 1; // X Z Y order
                j1 = 0;
                k1 = 0;
                i2 = 1;
                j2 = 0;
                k2 = 1;
            }
            else
            {
                i1 = 0; // Z X Y order
                j1 = 0;
                k1 = 1;
                i2 = 1;
                j2 = 0;
                k2 = 1;
            }
        }
        else // x0<y0
        {
            if(y0 < z0)
            {
                i1 = 0; // Z Y X order
                j1 = 0;
                k1 = 1;
                i2 = 0;
                j2 = 1;
                k2 = 1;
            }
            else if(x0 < z0)
            {
                i1 = 0; // Y Z X order
                j1 = 1;
                k1 = 0;
                i2 = 0;
                j2 = 1;
                k2 = 1;
            }
            else
            {
                i1 = 0; // Y X Z order
                j1 = 1;
                k1 = 0;
                i2 = 1;
                j2 = 1;
                k2 = 0;
            }
        }
        //  A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
        //  a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
        //  a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
        //  c = 1/6.
        const double x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
        const double y1 = y0 - j1 + G3;
        const double z1 = z0 - k1 + G3;
        const double x2 = x0 - i2 + 2.0 * G3; // Offsets for third corner in (x,y,z) coords
        const double y2 = y0 - j2 + 2.0 * G3;
        const double z2 = z0 - k2 + 2.0 * G3;
        const double x3 = x0 - 1.0 + 3.0 * G3; // Offsets for last corner in (x,y,z) coords
        const double y3 = y0 - 1.0 + 3.0 * G3;
        const double z3 = z0 - 1.0 + 3.0 * G3;
        // Work out the hashed gradient indices of the four simplex corners
        const int ii = i & 255;
        const int jj = j & 255;
        const int kk = k & 255;
        const int gi0 = perm[ii + perm[jj + perm[kk]]] % 12;
        const int gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]] % 12;
        const int gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]] % 12;
        const int gi3 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]] % 12;
        // Calculate the contribution from the four corners
        double t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
        if(t0 < 0)
            n0 = 0.0;
        else
        {
            t0 *= t0;
            n0 = t0 * t0 * dot(grad3[gi0], x0, y0, z0);
        }
        double t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
        if(t1 < 0)
            n1 = 0.0;
        else
        {
            t1 *= t1;
            n1 = t1 * t1 * dot(grad3[gi1], x1, y1, z1);
        }
        double t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
        if(t2 < 0)
            n2 = 0.0;
        else
        {
            t2 *= t2;
            n2 = t2 * t2 * dot(grad3[gi2], x2, y2, z2);
        }
        double t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
        if(t3 < 0)
            n3 = 0.0;
        else
        {
            t3 *= t3;
            n3 = t3 * t3 * dot(grad3[gi3], x3, y3, z3);
        }
        // Add contributions from each corner to get the final noise value.
        // The result is scaled to stay just inside [-1,1]
        return 32.0 * (n0 + n1 + n2 + n3);
    }

    double simplex_2d_noise(double x, double y, uint32_t octaves, double z)
    {
        double total = 0;

        // parametrization of octaves to match power spectrum of real grain scans
        static double f[] = {0.4910, 0.9441, 1.7280};
        static double a[] = {0.2340, 0.7850, 1.2150};

        for(uint32_t o = 0; o < octaves; o++)
        {
            total += (simplex_noise(x * f[o] / z, y * f[o] / z, o) * a[o]);
        }
        return total;
    }

    float paper_resp(float exposure, float mb, float gp, float divgr)
    {   
        float dived = 1.f;
        if(divgr > 1.8f) {
            dived = 1.f + (divgr - 1.8f);
        }
        const float delta = dived * GRAIN_LUT_DELTA_MAX * expf((mb / 100.0f) * logf(GRAIN_LUT_DELTA_MIN / dived));
        const float density = (1.0f + 2.0f * delta) / (1.0f + expf( (4.0f * gp * (0.5f - exposure)) / (1.0f + 2.0f * delta) )) - delta;
        return density;
    }

    float paper_resp_inverse(float density, float mb, float gp, float divgr)
    {
        float dived = 1.f;
        if(divgr > 1.8f) {
            dived = 1.f + (divgr - 1.8f);
        }
        const float delta =  dived * GRAIN_LUT_DELTA_MAX * expf((mb / 100.0f) * logf(GRAIN_LUT_DELTA_MIN / dived));
        const float exposure = -logf((1.0f + 2.0f * delta) / (density + delta) - 1.0f) * (1.0f + 2.0f * delta) / (4.0f * gp) + 0.5f;
        return exposure;
    }

    void evaluate_grain_lut(const float mb, float divgr)
    {
        for(int i = 0; i < GRAIN_LUT_SIZE; i++)
        {
            for(int j = 0; j < GRAIN_LUT_SIZE; j++)
            {
                float gu = (float)i / (GRAIN_LUT_SIZE - 1) - 0.5;
                float l = (float)j / (GRAIN_LUT_SIZE - 1);
                float divg = divgr; //1.f
                grain_lut[j * GRAIN_LUT_SIZE + i] = 32768.f * (paper_resp(gu + paper_resp_inverse(l, mb, divg * GRAIN_LUT_PAPER_GAMMA, divgr), mb, divg * GRAIN_LUT_PAPER_GAMMA, divgr) - l);
            }
        }
    }

    float lut_lookup(const float x, const float y)
    {
        const float _x = LIM((x + 0.5f) * (GRAIN_LUT_SIZE - 1), 0.f, float(GRAIN_LUT_SIZE - 1));
        const float _y = LIM(y * (GRAIN_LUT_SIZE - 1), 0.f, float(GRAIN_LUT_SIZE - 1));

        const int _x0 = _x < GRAIN_LUT_SIZE - 2 ? _x : GRAIN_LUT_SIZE - 2;
        const int _y0 = _y < GRAIN_LUT_SIZE - 2 ? _y : GRAIN_LUT_SIZE - 2;

        const int _x1 = _x0 + 1;
        const int _y1 = _y0 + 1;

        const float x_diff = _x - _x0;
        const float y_diff = _y - _y0;

        const float l00 = grain_lut[_y0 * GRAIN_LUT_SIZE + _x0];
        const float l01 = grain_lut[_y0 * GRAIN_LUT_SIZE + _x1];
        const float l10 = grain_lut[_y1 * GRAIN_LUT_SIZE + _x0];
        const float l11 = grain_lut[_y1 * GRAIN_LUT_SIZE + _x1];

        const float xy0 = (1.0 - y_diff) * l00 + l10 * y_diff;
        const float xy1 = (1.0 - y_diff) * l01 + l11 * y_diff;
        return xy0 * (1.0f - x_diff) + xy1 * x_diff;
    }

    int ox;
    int oy;
    int fw;
    int fh;
    double scale;
    int perm[512];
    float grain_lut[GRAIN_LUT_SIZE * GRAIN_LUT_SIZE];
};

} // namespace


void ImProcFunctions::filmGrain(Imagefloat *rgb, int isogr, int strengr, int scalegr, float divgr, int bfw, int bfh, int call, int fw, int fh)
{
    if (settings->verbose) {
        printf("iso=%i strength=%i scale=%i gamma=%f\n", isogr, strengr, scalegr, divgr);
    }

    GrainEvaluator ge(0, 0, bfw, bfh, scale, divgr, call, fw, fh);
    ge(isogr, strengr, scalegr, divgr, rgb, multiThread, call, fw, fh);
}

} // namespace rtengine
