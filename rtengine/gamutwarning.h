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

#pragma once

#include "iccstore.h"
#include "noncopyable.h"
#include "image8.h"

namespace rtengine {

class GamutWarning: public NonCopyable {
public:
    GamutWarning(cmsHPROFILE iprof, cmsHPROFILE gamutprof, bool bpc);
    ~GamutWarning();
    void markLine(Image8 *image, int y, float *srcbuf, float *buf1, float *buf2);
    
private:
    void mark(Image8 *image, int i, int j);
    
    cmsHTRANSFORM lab2ref;
    cmsHTRANSFORM lab2softproof;
    cmsHTRANSFORM softproof2ref;
};

} // namespace rtengine
