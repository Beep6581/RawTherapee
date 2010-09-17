#ifndef _IMSOURCEFILTER_H_
#define _IMSOURCEFILTER_H_

#include <filter.h>

namespace rtengine {

class ImageSourceFilter : public Filter {

public:

	virtual void getFullImageSize (int& w, int& h) = 0;
};

}

#ifndef
