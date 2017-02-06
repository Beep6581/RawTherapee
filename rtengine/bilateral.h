/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2017 RawTherapee development team
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
#pragma once

#include <memory>

#include "noncopyable.h"
#include "LUT.h"

namespace rtengine
{

class Bilateral final :
    public NonCopyable
{
public:
    Bilateral(
        float** source,
        float** destination,
        float** buffer,
        int width,
        int height,
        double sigma,
        double sensitivity
    );
    ~Bilateral();

    void compute();

private:
    class Implementation;

    std::unique_ptr<Implementation> implementation;
};

}
