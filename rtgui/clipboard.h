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
#ifndef _CLIPBOARD_
#define _CLIPBOARD_

#include <vector>
#include <rtengine.h>
#include <mydiagonalcurve.h>

class Clipboard {
    
    bool _hasIPTC;
    rtengine::MetadataList iptc;
    bool _hasProcParams;
    rtengine::procparams::ProcParams procParams;
    DiagonalCurveType hasCurveDataType;
    std::vector<double> curve;


    public:
        void                                               setIPTC (const rtengine::MetadataList& iptcc) { iptc = iptcc; _hasIPTC = true;}
        const rtengine::MetadataList&                      getIPTC () { return iptc;  }
        bool                                               hasIPTC () { return _hasIPTC; }
        
        void                                               setProcParams (const rtengine::procparams::ProcParams& pparams) { procParams = pparams; _hasProcParams = true; }
        const rtengine::procparams::ProcParams&            getProcParams ()                                                { return procParams; }
        bool                                               hasProcParams () { return _hasProcParams; }

        void                                               setCurveData (std::vector<double>& p, DiagonalCurveType type ) { curve = p;  hasCurveDataType = type; return; }
        const std::vector<double> &                        getCurveData ()                       { return curve; }
        DiagonalCurveType                                  hasCurveData ()                       { return hasCurveDataType; }


};

extern Clipboard clipboard;

#endif
