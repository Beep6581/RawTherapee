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
#ifndef _PROCESSINGJOB_
#define _PROCESSINGJOB_

#include "rtengine.h"

namespace rtengine
{

class ProcessingJobImpl : public ProcessingJob
{

public:
    Glib::ustring fname;
    bool isRaw;
    InitialImage* initialImage;
    procparams::ProcParams pparams;

    ProcessingJobImpl (const Glib::ustring& fn, bool iR, const procparams::ProcParams& pp)
        : fname(fn), isRaw(iR), initialImage(NULL), pparams(pp) {}

    ProcessingJobImpl (InitialImage* iImage, const procparams::ProcParams& pp)
        : fname(""), isRaw(true), initialImage(iImage), pparams(pp)
    {
        iImage->increaseRef();
    }

    ~ProcessingJobImpl ()
    {
        if (initialImage) {
            initialImage->decreaseRef();
        }
    }
};

}

#endif
