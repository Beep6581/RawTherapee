/*
 * filterfactory.h
 *
 *  Created on: Aug 18, 2010
 *      Author: gabor
 */

#ifndef FILTERFACTORY_H_
#define FILTERFACTORY_H_

#include <map>
#include "filter.h"

namespace rtengine {

class FilterFactory {

	std::map<String, FilterDescriptor*> filterDescriptors;

public:
	FilterFactory();
	void registerFilterDescriptor (FilterDescriptor* descr);
	FilterDescriptor* getFilterDescriptor (const String& name);
};

extern FilterFactory* filterFactory;
}
#endif /* FILTERFACTORY_H_ */
