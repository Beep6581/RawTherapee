/*
 * rtcommon.h
 *
 *  Created on: Jan 19, 2011
 *      Author: gabor
 */

#ifndef RTCOMMON_H_
#define RTCOMMON_H_

#include <glibmm.h>
#include <cairomm/cairomm.h>
#include <vector>

namespace rtengine {

typedef Glib::ustring String;

typedef std::vector<float> 	FloatList;
typedef std::vector<int> 	IntList;
typedef std::vector<String> StringList;
typedef Cairo::RefPtr<Cairo::ImageSurface> DisplayImage;

}

#endif /* RTCOMMON_H_ */
