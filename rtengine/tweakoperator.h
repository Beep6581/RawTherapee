/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

namespace rtengine
{

namespace procparams
{

class ProcParams;

}

/** This class can let objects alter the collected values of the ProcParams for a specific
 *  purpose, e.g. displaying a preview image at a specific point in the pipeline or with
 *  disabled tools. Before starting the rendering, the engine will call the TweakOperator
 *  (if set) to modify the ProcParams. The untweaked one will still exist as a backup, and
 *  can be sent back if necessary.  */
class TweakOperator
{
public:
    virtual ~TweakOperator() {}

    /** Callback that will alter the ProcParams before the image is computed. */
    virtual void tweakParams(procparams::ProcParams& pparams) = 0;
};

}
