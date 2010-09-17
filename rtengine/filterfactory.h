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

	std::map<std::string, FilterDescriptor*> filterDescriptors;

public:
	FilterFactory();
	void registerFilterDescriptor (FilterDescriptor* descr);
	void createFilterAddToList (const std::string& name, Filter* tail); // creates a new filter instance and adds it to the linked list pointed by "tail"
};

extern FilterFactory filterFactory;
}
#endif /* FILTERFACTORY_H_ */
