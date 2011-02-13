#ifndef _FILTER_H_
#define _FILTER_H_

#include <set>
#include <string>
#include "imageview.h"
#include "multiimage.h"
#include "buffer.h"
#include "procevents.h"
#include "rtengine.h"
#include "imageview.h"

namespace rtengine {

class Filter;
class FilterDescriptor {

protected:
	String		 	  		name;
	std::set<ProcEvent> 	myEvents;
	bool 		   	  		forceOutCache;
	MultiImage::ColorSpace 	inputColorSpace, outputColorSpace;
    bool                    applyOnRawImage;
    bool                    applyOnStdImage;
    bool                    applyOnThumbnail;

public:
	FilterDescriptor (const String& name, MultiImage::ColorSpace ics, MultiImage::ColorSpace ocs, bool forceCache = false);
	void addTriggerEvent (ProcEvent ev);

	virtual void			getDefaultParameters (ProcParams& defProcParams) const {}

	String		 			getName () const { return name; }
	MultiImage::ColorSpace  getInputColorSpace () const { return inputColorSpace; }
	MultiImage::ColorSpace  getOutputColorSpace () const { return outputColorSpace; }
	bool					forceOutputCache () const { return forceOutCache; }
	bool					myTriggerEvent (ProcEvent ev) const;
    bool                    isAppliedOnRawImage ()  const { return applyOnRawImage; }
    bool                    isAppliedOnStdImage ()  const { return applyOnStdImage; }
    bool                    isAppliedOnThumbnail () const { return applyOnThumbnail; }
	virtual void    		createAndAddToList (Filter* tail) const = 0; // creates a new filter instance and adds it to the linked list pointed by "tail"
};

class FilterChain;
class Filter {

	friend class FilterChain;

private:

	Filter* next;
	Filter* prev;
    Filter* shortCutPrev;
	Filter* parent;
	FilterChain* myFilterChain;

	bool        valid;
	bool 		hasOutputCache;
	bool 		forceOutputCache;
	MultiImage* outputCache;

	ImageView   sourceImageView;
	ImageView   targetImageView;
    ImageView   scaledSourceImageView;
    ImageView   scaledTargetImageView;

	ProgressListener* plistener;

	void setupCache ();
	void setProcParams (ProcParams* pparams);

protected:

	bool multiThread;
	ProcParams* procParams;
	const FilterDescriptor* descriptor;

	const ImageView&    getSourceImageView () const { return sourceImageView; }
	const ImageView&    getTargetImageView () const { return targetImageView; }
    Filter*             getPreviousFilter  () { return prev; }
	FilterChain*		getFilterChain     () { return myFilterChain; }
    ProgressListener*   getProgressListener() { return plistener; }
    const ImageView&    getScaledSourceImageView () { return scaledSourceImageView; }
    const ImageView&    getScaledTargetImageView () { return scaledTargetImageView; }

public:
    Filter*             getParentFilter    () { return parent; }

	void addNext (Filter* f);

public:

	Filter (FilterDescriptor* descr);
	virtual ~Filter ();

	// Return the size of the required buffer. One common (shared among all the filters) buffer is maintained and passed
	// to "process".
	virtual Dim getReqiredBufferSize ();
	// return full image size if this filter was the last one
	virtual Dim getFullImageSize ();

    // return the scale factor of the image at the target side of the filter
	virtual double getScale ();
	// return scale factor by assuming the given skip
	virtual void reverseTransPoint (int x, int y, int& xv, int& yv);
    // returns true if set "events" contains at least one event that invalidates the result of this filter
	virtual bool isTriggerEvent (const std::set<ProcEvent>& events);
	// The main procedure of the filter, it must be overridden
	virtual void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) = 0;
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

    // these two functions support scale pre-calculation before imageview setup (to support gui)
    // filters have to override these when:
    // a) they override the requested skip (e.g. demosaicing filter that can work only on the whole image (skip=1) whatever skip is requested by the next filter
    // b) the filter changes the scale of the result image (e.g. the resize filter when it is not applied on thumbnails)
    // returns target skip parameter assuming the next filter requires "nextInSkip" skip
    virtual int getTargetSkip (int nextInSkip);
    // returns target scale assuming the given skip
    virtual double getTargetScale (int skip);
};

}

#endif
