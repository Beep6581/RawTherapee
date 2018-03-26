/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

/**
 * Adapted from PhotoFlow, Copyright (C) 2014 Ferrero Andrea
 * also distributed under the GPL V3+
 */

#include "gamutwarning.h"
#include <iostream>

namespace rtengine {

GamutWarning::GamutWarning(cmsHPROFILE iprof, cmsHPROFILE gamutprof, RenderingIntent intent, bool gamutbpc):
    lab2ref(nullptr),
    lab2softproof(nullptr),
    softproof2ref(nullptr)
{
    if (cmsIsMatrixShaper(gamutprof) && !cmsIsCLUT(gamutprof, intent, LCMS_USED_AS_OUTPUT)) {
        cmsHPROFILE aces = ICCStore::getInstance()->getProfile("ACES");
        if (aces) {
            lab2ref = cmsCreateTransform(iprof, TYPE_Lab_FLT, aces, TYPE_RGB_FLT, INTENT_ABSOLUTE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
            lab2softproof = cmsCreateTransform(iprof, TYPE_Lab_FLT, gamutprof, TYPE_RGB_FLT, INTENT_ABSOLUTE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
            softproof2ref = cmsCreateTransform(gamutprof, TYPE_RGB_FLT, aces, TYPE_RGB_FLT, INTENT_ABSOLUTE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | (gamutbpc ? cmsFLAGS_BLACKPOINTCOMPENSATION : 0));
        }
    } else {
        lab2ref = nullptr;
        lab2softproof = cmsCreateTransform(iprof, TYPE_Lab_FLT, gamutprof, TYPE_RGB_FLT, INTENT_ABSOLUTE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE);
        softproof2ref = cmsCreateTransform(gamutprof, TYPE_RGB_FLT, iprof, TYPE_Lab_FLT, INTENT_ABSOLUTE_COLORIMETRIC, cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE | (gamutbpc ? cmsFLAGS_BLACKPOINTCOMPENSATION : 0));
    }

    if (!softproof2ref) {
        if (lab2softproof) {
            cmsDeleteTransform(lab2softproof);
            lab2softproof = nullptr;
        }
    } else if (!lab2softproof) {
        if (softproof2ref) {
            cmsDeleteTransform(softproof2ref);
            softproof2ref = nullptr;
        }
    }
}


GamutWarning::~GamutWarning()
{
    if (softproof2ref) {
        cmsDeleteTransform(softproof2ref);
    }
    if (lab2ref) {
        cmsDeleteTransform(lab2ref);
    }
    if (lab2softproof) {
        cmsDeleteTransform(lab2softproof);
    }
}


void GamutWarning::markLine(Image8 *image, int y, float *srcbuf, float *buf1, float *buf2)
{
    if (softproof2ref) {
        const int width = image->getWidth();
        
        float delta_max = lab2ref ? 0.0001f : 4.9999f;
        cmsDoTransform(lab2softproof, srcbuf, buf2, width);
        cmsDoTransform(softproof2ref, buf2, buf1, width);
        
        float *proofdata = buf1;
        float *refdata = srcbuf;
        
        if (lab2ref) {
            cmsDoTransform(lab2ref, srcbuf, buf2, width);
            refdata = buf2;

            int iy = 0;
            for (int j = 0; j < width; ++j) {
                float delta = max(std::abs(proofdata[iy] - refdata[iy]), std::abs(proofdata[iy+1] - refdata[iy+1]), std::abs(proofdata[iy+2] - refdata[iy+2]));
                iy += 3;

                if (delta > delta_max) {
                    mark(image, y, j);
                }
            }
        } else {
            int iy = 0;
            for (int j = 0; j < width; ++j) {
                cmsCIELab lab1 = { proofdata[iy], proofdata[iy+1], proofdata[iy+2] };
                cmsCIELab lab2 = { refdata[iy], refdata[iy+1], refdata[iy+2] };
                iy += 3;
                float delta = cmsDeltaE(&lab1, &lab2);
                if (delta > delta_max) {
                    mark(image, y, j);
                }
            }
        }
    }
}


inline void GamutWarning::mark(Image8 *image, int y, int x)
{
    image->r(y, x) = 0;
    image->g(y, x) = 255;
    image->b(y, x) = 255;
}


} // namespace rtengine
