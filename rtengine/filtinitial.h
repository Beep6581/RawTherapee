/*
 * filtinitial.h
 *
 *  Created on: Sep 22, 2010
 *      Author: gabor
 */

#ifndef FILTINITIAL_H_
#define FILTINITIAL_H_

#include "filter.h"
#include "imagesource.h"

namespace rtengine {

class RawInitialFilterDescriptor : public FilterDescriptor {

    public:
        RawInitialFilterDescriptor ();
        void createAndAddToList (Filter* tail) const;
};

class StdInitialFilterDescriptor : public FilterDescriptor {

    public:
        StdInitialFilterDescriptor ();
        void createAndAddToList (Filter* tail) const;
};

class InitialFilter : public Filter {

        ImageSource* imgsrc;

    public:
        InitialFilter (ImageSource* imgs);
        Dim  getFullImageSize ();
        void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTINITIAL_H_ */
