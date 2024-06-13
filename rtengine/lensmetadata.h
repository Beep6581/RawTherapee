/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2024 Rawtherapee developers
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

#include <memory>

#include "lcp.h"
#include "metadata.h"
#include "procparams.h"
#include "rtengine.h"

namespace rtengine
{

/* MetadataLensCorrection is an abstract class for various lens correction based on raw file metadata
 this metadata is vendor dependent */
class MetadataLensCorrection : public LensCorrection,
                               public NonCopyable
{
public:
    virtual void initCorrections(int width, int height, const procparams::CoarseTransformParams &coarse, int rawRotationDeg) = 0;
};

/* MetadataLensCorrectionFinder tries to find and return MetadataLensCorrection for the provided metadata */
class MetadataLensCorrectionFinder
{
public:
    static std::unique_ptr<MetadataLensCorrection> findCorrection(const FramesMetaData *meta);
};

} // namespace rtengine
