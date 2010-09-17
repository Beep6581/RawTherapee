#ifndef _FILTER_H_
#define _FILTER_H_

#include <set>
#include <string>
#include "imageview.h"
#include "multiimage.h"
#include "buffer.h"

namespace rtengine {

class FilterDescriptor {

protected:
	std::string 	  		name;
	std::set<ProcEvent> 	myEvents;
	bool 		   	  		forceOutCache;
	MultiImage::ColorSpace 	inputColorSpace, outputColorSpace;

public:
	FilterDescriptor (const std::string name, MultiImage::ColorSpace ics, MultiImage::ColorSpace ocs, bool forceCache = false);
	void addTriggerEvent (ProcEvent ev);

	std::string 			getName () const { return name; }
	MultiImage::ColorSpace  getInputColorSpace () const { return inputColorSpace; }
	MultiImage::ColorSpace  getOutputColorSpace () const { return outputColorSpace; }
	bool					forceOutputCache ();
	bool					myTriggerEvent (ProcEvent ev) const;
	virtual void    		createAndAddToList (Filter* tail) const = 0; // creates a new filter instance and adds it to the linked list pointed by "tail"
};

class FilterChain;
class Filter {

	friend class FilterChain;

private:

	Filter* next;
	Filter* prev;
	Filter* parent;
	FilterChain* myFilterChain;

	bool 		hasOutputCache;
	bool 		forceOutputCache;
	MultiImage* outputCache;

	ImageView   sourceImageView;
	ImageView   targetImageView;

	ProgressListener* plistener;

	void setupCache ();
	void setProcParams (ProcParams* pparams);

protected:

	bool multiThread;
	ProcParams* procParams;
	const FilterDescriptor* descriptor;

	const ImageView&    getSourceImageView () { return sourceImageView; }
	const ImageView&    getTargetImageView () { return targetImageView; }
	const Filter*	    getParentFilter    () { return parent; }
	FilterChain*		getFilterChain     () { return myFilterChain; }
    ProgressListener*   getProgressListener() { return plistener; }
	double              getScale           ();

public:

	void addNext (Filter* f);

public:

	Filter (FilterDescriptor* descr);
	virtual ~Filter ();

	// Return the size of the required buffer. One common (shared among all the filters) buffer is maintained and passed
	// to "process".
	virtual void getReqiredBufferSize (int& w, int& h);
	// return full image size if this filter was the last one
	virtual void getFullImageSize (int& w, int& h);
	// return the coordinates (xv,yv) corresponding to the "source" side given coordinates (x,y) corresponding to the target side
	virtual void reverseTransPoint (int x, int y, int& xv, int& yv);
	// returns true if set "events" contains at least one event that invalidates the result of this filter
	virtual bool isTriggerEvent (const std::set<ProcEvent>& events);
	// The main procedure of the filter, it must be overridden
	virtual void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int> buffer) = 0;
    // Calculate the image view (rectangle) provided by this filter when the given image view is requested from by the next filter
    // Most filters return "requestedImView" (this is the default implementation).
    // Some filters like demosaicing methods can calculate full image only, thus
    // they return the image view corresponding to the full image.
    virtual ImageView calculateTargetImageView (const ImageView& requestedImView);
    // Calculate the image view (rectangle) required by this filter from the previous filter to obtain the given image view
    // If there is no transform in geometry, it returns the target image view calculated above (default implementation).
    // Rotation / resize etc filters need to override it.
    virtual ImageView calculateSourceImageView (const ImageView& requestedImView);

    const FilterDescriptor* getDescriptor () const { return descriptor; }

};

}

#endif
