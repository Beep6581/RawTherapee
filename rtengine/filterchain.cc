#include "filterchain.h"
#include "imageview.h"
#include "settings.h"
#include "filterfactory.h"
#include "filtinitial.h"
#include "iccstore.h"
#include "curves.h"
#include <lcms2.h>
#include <iostream>
#include "macros.h"

namespace rtengine {

FilterChain::FilterChain (ImProcListener* listener, ImageSource* imgSource, ProcParams* params, bool multiThread, bool oneShot)
	: listener(listener), imgSource(imgSource), procParams(params), multiThread(multiThread), invalidated(true), oneShot(oneShot) {

    last = first = new InitialFilter (imgSource);

	if (params->getBoolean ("FilterOrder", "Custom") && !params->getStringList ("FilterOrder", "FilterList").empty())
		filterOrder = params->getStringList ("FilterOrder", "FilterList");
	else
		filterOrder = Settings::settings->filterList;
	setupChain (NULL);
	
	imgSource->increaseRef ();
}

FilterChain::FilterChain (ImProcListener* listener, FilterChain* previous)
	: listener(listener), imgSource(previous->imgSource),
	  procParams(previous->procParams), multiThread(previous->multiThread), invalidated(true) {

    last = first = new InitialFilter (imgSource);

    filterOrder = previous->filterOrder;
	setupChain (previous);
	
	imgSource->increaseRef ();
}

void FilterChain::setupChain (FilterChain* previous) {

    // First filter already there, it is the initial filter. We go on from the second one.
	Filter* parent = previous ? previous->first->next : NULL;
	for (int i=0; i<filterOrder.size(); i++)
		if (!oneShot && filterOrder[i] == "--Cache--")
			last->forceOutputCache = true;
		else {
		    FilterDescriptor* fDescr = filterFactory->getFilterDescriptor (filterOrder[i]);
		    if (fDescr && (
		            (fDescr->isAppliedOnThumbnail() && imgSource->isThumbnail()) ||
		            (fDescr->isAppliedOnRawImage() && imgSource->isRaw()) ||
		            (fDescr->isAppliedOnStdImage() && !imgSource->isRaw()))) {
                fDescr->createAndAddToList (last);
                while (last->next) {
                    last = last->next;
                    if (parent) {
                        parent = parent->next;
                        last->parent = parent;
                    }
                    last->multiThread = multiThread;
                    last->myFilterChain = this;
                }
		    }
		}
	for (Filter* curr = first; curr; curr = curr->next) {
	    curr->setProcParams (procParams);
	    curr->myFilterChain = this;
	}
}

void FilterChain::setNextChain (FilterChain* other) {

    if (!other)
        for (Filter* curr = first; curr; curr = curr->next)
            curr->parent = NULL;
    else
        for (Filter *curr = first, *pcurr = other->first; curr; curr = curr->next, pcurr = pcurr->next)
            curr->parent = pcurr;
}

FilterChain::~FilterChain () {

	Filter* curr = first;
	while (curr) {
		Filter* tmp = curr;
		curr = curr->next;
		delete tmp;
	}
	imgSource->decreaseRef ();
}

void FilterChain::invalidate () {

    invalidated = true;
}

void FilterChain::setupProcessing (const std::set<ProcEvent>& events, bool useShortCut) {

	// tell the listener the full size of the image with the current settings and let it decide the portion to refresh
	Dim fullSize = getFullImageSize ();
	ImageView reqView;
	if (listener)
	    reqView = listener->getViewToProcess (fullSize);
	else {
        reqView.w = fullSize.width;
        reqView.h = fullSize.height;
	}
	// walk through the list and find first filter that needs to be refreshed because of the changed view
	Filter* curr = last;
	firstToUpdate = last;
	ImageView view = reqView;
	while (curr) {
	    ImageView tiv = curr->calculateTargetImageView (view);
        if (curr->targetImageView != tiv)
            firstToUpdate = curr;
        curr->targetImageView = tiv;
        view  = curr->calculateSourceImageView (view);
        curr->sourceImageView = view;
        curr = curr->prev;
	}

	// calculate and set source and target image sizes of the filters
	for (curr = first; curr; curr = curr->next) {
        ImageView scSourceIV = curr->sourceImageView.getScaled (curr==first ? imgSource->getScale () : curr->prev->getScale ());
        ImageView scTargetIV = curr->targetImageView.getScaled (curr->getScale ());
        if (scSourceIV != curr->scaledSourceImageView || scTargetIV != curr->scaledTargetImageView) {
            curr->scaledSourceImageView = scSourceIV;
            curr->scaledTargetImageView = scTargetIV;
            curr->valid = false;
        }
	}

	// set mandatory cache points
	first->hasOutputCache = false;
	for (curr = first->next; curr != last; curr = curr->next)
		if (curr->forceOutputCache || curr->sourceImageView != curr->targetImageView
		        || curr->scaledSourceImageView != curr->scaledTargetImageView
		        || (curr->next && (curr->targetImageView != curr->next->sourceImageView || curr->scaledTargetImageView != curr->next->scaledSourceImageView)))
			curr->hasOutputCache = true;
		else
			curr->hasOutputCache = false;
	last->hasOutputCache = true;

	// walk further to find first filter that needs to be recalculated
	for (curr = firstToUpdate; curr; curr = curr->prev)
		if (invalidated || !curr->valid || curr->isTriggerEvent (events))
			firstToUpdate = curr;

	// walk further to find closest cache
	for (curr = firstToUpdate; curr; curr = curr->prev)
		if (!curr->prev || curr->prev->hasOutputCache) {
			firstToUpdate = curr;
			break;
		}

	// set every filter from firstToUpdate to be invalid
	for (curr = firstToUpdate; curr; curr = curr->next)
	    curr->valid = false;

	// detect shortcut possibilities in the filter group
	if (useShortCut) {
	    // clear shortcut pointers
	    for (curr = first; curr; curr = curr->next)
	        curr->shortCutPrev = NULL;
	    // find possible shortcut filter
	    for (curr = last; curr && curr != firstToUpdate; curr = curr->prev) {
	        Filter* p = curr->prev->parent;
	        while (p) {
                double skip = curr->getScale() / p->getScale();
	            if (p->hasOutputCache && fabs(skip-round(skip))<1e-12 && curr->sourceImageView.isPartOf (p->targetImageView)) {
	                curr->shortCutPrev = p;
	                firstToUpdate = curr;
	                break;
	            }
	            p = p->parent;
	        }
	        if (curr->shortCutPrev)
	            break;
	    }
	}

	// set up caches
	for (curr = first; curr; curr = curr->next)
		curr->setupCache ();
}

void FilterChain::process (const std::set<ProcEvent>& events, Buffer<float>* buffer, MultiImage* worker) {

    for (Filter* curr = firstToUpdate; curr; curr = curr->next) {

	    if (Settings::settings->verbose)
	        std::cout << "Applying filter " << String2StdString(curr->getDescriptor()->getName()) << " "
                << curr->sourceImageView << ", " << curr->targetImageView << " in progress...";
	    std::flush (std::cout);

	    MultiImage* sourceImage = NULL;
		MultiImage* targetImage = NULL;

		// set up source image for the filter
		if (curr != first) {
		    if (curr->shortCutPrev) {
		        sourceImage = worker;
                worker->setDimensions (curr->scaledSourceImageView.w, curr->scaledSourceImageView.h);
                int skip = (int)round (curr->getScale() / curr->shortCutPrev->getScale());
                worker->copyFrom (curr->shortCutPrev->outputCache, curr->scaledSourceImageView.x - curr->shortCutPrev->scaledTargetImageView.x * skip, curr->scaledSourceImageView.y - curr->shortCutPrev->scaledTargetImageView.y * skip, skip);
		    }
		    else if (curr->prev->hasOutputCache && curr->prev->targetImageView == curr->sourceImageView && curr->prev->scaledTargetImageView == curr->scaledSourceImageView)
				sourceImage = curr->prev->outputCache;
			else {
				sourceImage = worker;
				worker->setDimensions (curr->scaledSourceImageView.w, curr->scaledSourceImageView.h);
				if (curr->prev->hasOutputCache && curr->prev->targetImageView != curr->sourceImageView || curr->prev->scaledTargetImageView != curr->scaledSourceImageView) {
				    // There is a cache, but image views do not fit. Assume that sourceImageView is the part of prev->targetImageView.
                    int skip = (int)round (curr->getScale() / curr->prev->getScale());
				    worker->copyFrom (curr->prev->outputCache, curr->scaledSourceImageView.x - curr->prev->scaledTargetImageView.x * skip, curr->scaledSourceImageView.y - curr->prev->scaledTargetImageView.y * skip, skip);
				}
			}
			sourceImage->convertTo (curr->descriptor->getInputColorSpace());
		}

		// set up target image for the filter
		if (curr->hasOutputCache)
			targetImage = curr->outputCache;
		else {
			worker->setDimensions (curr->scaledTargetImageView.w, curr->scaledTargetImageView.h);
			targetImage = worker;
		}
		targetImage->switchTo (curr->descriptor->getOutputColorSpace());

		// apply filter
		curr->process (events, sourceImage, targetImage, buffer);
		curr->valid = true;

		if (Settings::settings->verbose)
            std::cout << "ready." << std::endl;
	}
	invalidated = false;
}

Dim FilterChain::getReqiredBufferSize () {

	Dim bufferSize;
	for (Filter* curr = first; curr; curr = curr->next)
		bufferSize.setMax (curr->getReqiredBufferSize ());

	return bufferSize;
}

Dim FilterChain::getReqiredWorkerSize () {

    Dim workerSize;
    for (Filter* curr = firstToUpdate; curr; curr = curr->next)
        if (!curr->hasOutputCache || (curr->prev && (curr->prev->targetImageView != curr->sourceImageView || curr->prev->scaledTargetImageView != curr->scaledSourceImageView)))
            workerSize.setMax (curr->scaledSourceImageView.getSize ());

    return workerSize;
}

Dim FilterChain::getFullImageSize () {

	if (!last)
		return Dim();
	else
		return last->getFullImageSize ();
}

// this method is called by the gui to determine the image view to request from rtengine such that
// the size of the result is given by a gui window
double FilterChain::getScale (int skip) {

    if (!last)
        return 1.0;
    else
        return last->getTargetScale (skip);
}

ImageView FilterChain::getLastImageView () {

    if (last)
        return last->targetImageView;
    else
        return ImageView();
}

double FilterChain::getLastScale () {

    return last ? last->getScale() : 1.0;
}

FinalImage16* FilterChain::getFinalImage () {

	String outputProfile  = procParams->getString ("ColorManagement", "OutputProfile");
	String workingProfile = procParams->getString ("ColorManagement", "WorkingProfile");

    last->outputCache->convertTo (MultiImage::RGB, true, workingProfile);

    // calculate crop rectangle
	int cx = 0, cy = 0, cw = last->outputCache->width, ch = last->outputCache->height;
	bool crenabled = procParams->getBoolean ("Crop", "Enabled");
	if (crenabled) {
		cx = procParams->getInteger ("Crop", "RectX");
		cy = procParams->getInteger ("Crop", "RectY");
		cw = procParams->getInteger ("Crop", "RectW");
		ch = procParams->getInteger ("Crop", "RectH");
	}

    // adjust to valid values
    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw > last->outputCache->width)  cw = last->outputCache->width - cx;
    if (cy+ch > last->outputCache->height) ch = last->outputCache->height - cy;

    // obtain cropped image
    Image* final = new Image (cw, ch);
    unsigned char* idata = final->getData ();
    int pitch = final->getScanLineSize ();

	bool ready = false;
    if (outputProfile!="") {
        // custom output icm profile specified, call lcms
        cmsHPROFILE oprof = iccStore->getProfile (outputProfile);
        cmsHPROFILE iprof = iccStore->getProfile (workingProfile);
        if (oprof && iprof) {
        	// TODO: do not create cmstransform every time, store it and re-create it only when it has been changed
        	cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16, oprof, TYPE_RGB_16, Settings::settings->colorimetricIntent, cmsFLAGS_NOCACHE);
			#pragma omp parallel for if (multiThread)
			for (int i=0; i<ch; i++) {
				FIRGB16* pixel  = (FIRGB16*)(idata + (ch-i-1)*pitch);
				for (int j=0; j<cw; j++) {
					pixel[j].red   = CLIP (65535.0 * last->outputCache->r[i+cy][j+cx]);
					pixel[j].green = CLIP (65535.0 * last->outputCache->r[i+cy][j+cx]);
					pixel[j].blue  = CLIP (65535.0 * last->outputCache->r[i+cy][j+cx]);
				}
				cmsDoTransform (hTransform, pixel, pixel, cw);
			}
        	cmsDeleteTransform(hTransform);
        	ready = true;
			ProfileContent outProfileContent = iccStore->getContent (outputProfile);
			final->setEmbeddedICCProfile (outProfileContent.length, outProfileContent.data);
        }
    }
    if (!ready && workingProfile != "sRGB") {
        // convert it to sRGB
        Matrix33 m;
        m.multiply (iccStore->workingSpaceMatrix ("sRGB"));
        m.multiply (iccStore->workingSpaceInverseMatrix (workingProfile));
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<ch; i++) {
			FIRGB16* pixel  = (FIRGB16*)(idata + (ch-i-1)*pitch);
			float r, g, b;
            for (int j=0; j<cw; j++) {
                m.transform (last->outputCache->r[i+cy][j+cx], last->outputCache->g[i+cy][j+cx], last->outputCache->b[i+cy][j+cx], r, g, b);
				pixel[j].red   = CLIP (65535.0 * CurveFactory::gamma2(r));
				pixel[j].green = CLIP (65535.0 * CurveFactory::gamma2(g));
				pixel[j].blue  = CLIP (65535.0 * CurveFactory::gamma2(b));
			}
		}
    }
	else {
		#pragma omp parallel for if (multiThread)
		for (int i=0; i<ch; i++) {
			FIRGB16* pixel  = (FIRGB16*)(idata + (ch-i-1)*pitch);
			for (int j=0; j<cw; j++) {
				pixel[j].red   = CLIP (65535.0 * CurveFactory::gamma2(last->outputCache->r[i+cy][j+cx]));
				pixel[j].green = CLIP (65535.0 * CurveFactory::gamma2(last->outputCache->g[i+cy][j+cx]));
				pixel[j].blue  = CLIP (65535.0 * CurveFactory::gamma2(last->outputCache->b[i+cy][j+cx]));
			}
		}
	}
	// set metadata
	Exiv2::ExifData ed;
	Exiv2::IptcData id;
	Exiv2::XmpData xd;
	getMetadataToSave (ed, id, xd);
	final->setMetadata (ed, id, xd);
	
    return final;
}

DisplayImage FilterChain::getDisplayImage () {

	String workingProfile = procParams->getString ("ColorManagement", "WorkingProfile");

    last->outputCache->convertTo(MultiImage::RGB, true, workingProfile);

#ifdef QTBUILD
    DisplayImage final (last->outputCache->width, last->outputCache->height, QImage::Format_RGB888);
    unsigned char* idata = final.bits ();
    int pitch = final.bytesPerLine ();
#else
    DisplayImage final = Cairo::ImageSurface::create (Cairo::FORMAT_RGB24, last->outputCache->width, last->outputCache->height);
    unsigned char* idata = final->get_data ();
    int pitch = final->get_stride ();
#endif

    if (Settings::settings->monitorProfile != "") {
        // custom output icm profile specified, call lcms
        cmsHPROFILE oprof = iccStore->getProfile (Settings::settings->monitorProfile);
        cmsHPROFILE iprof = iccStore->getProfile (workingProfile);
        if (oprof && iprof) {
			// TODO: do not create cmstransform every time, store it and re-create it only when it has been changed
			cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_8, oprof, TYPE_RGB_8, Settings::settings->colorimetricIntent, cmsFLAGS_NOCACHE);
			// TODO: do not convert to 8 bit before color space transformation for higher accuracy
			#pragma omp parallel for if (multiThread)
			for (int i=0; i<last->outputCache->height; i++) {
				int ix = i*pitch;
				for (int j=0; j<last->outputCache->width; j++) {
					idata[ix++] = CLIP (255.0 * last->outputCache->r[i][j]);
					idata[ix++] = CLIP (255.0 * last->outputCache->r[i][j]);
					idata[ix++] = CLIP (255.0 * last->outputCache->r[i][j]);
				}
				cmsDoTransform (hTransform, idata+ix, idata+ix, last->outputCache->width);
			}
			cmsDeleteTransform(hTransform);
			return final;
        }
    }
    
    if (workingProfile == "sRGB") {
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<last->outputCache->height; i++) {
            int ix = i * pitch;
            for (int j=0; j<last->outputCache->width; j++) {
                idata[ix++] = CLIPTO (CurveFactory::gamma2 (last->outputCache->r[i][j]) * 255.0, 0, 255);
                idata[ix++] = CLIPTO (CurveFactory::gamma2 (last->outputCache->g[i][j]) * 255.0, 0, 255);
                idata[ix++] = CLIPTO (CurveFactory::gamma2 (last->outputCache->b[i][j]) * 255.0, 0, 255);
            }
        }
    }
    else {
        Matrix33 m;
        m.multiply (iccStore->workingSpaceMatrix ("sRGB"));
        m.multiply (iccStore->workingSpaceInverseMatrix (workingProfile));
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<last->outputCache->height; i++) {
            int ix = i * pitch;
			float r, g, b;
            for (int j=0; j<last->outputCache->width; j++) {
                m.transform (last->outputCache->r[i][j], last->outputCache->g[i][j], last->outputCache->b[i][j], r, g, b);
                idata[ix++] = CLIPTO (CurveFactory::gamma2 (r) * 255.0, 0, 255);
                idata[ix++] = CLIPTO (CurveFactory::gamma2 (g) * 255.0, 0, 255);
                idata[ix++] = CLIPTO (CurveFactory::gamma2 (b) * 255.0, 0, 255);
            }
        }
    }
    return final;
}

void FilterChain::getMetadataToSave (Exiv2::ExifData& ed, Exiv2::IptcData& id, Exiv2::XmpData& xd) {

	if (imgSource->getMetaData()) {
		ed.clear (); 
		id.clear ();
		xd.clear ();
		
		// Assemble exif data. All tags belonging to image and subimage sections are deleted. (The tags describing the image representation
		// will be added by freeimage)
		for (Exiv2::ExifData::const_iterator it = imgSource->getMetaData()->getExifData().begin(); it!= imgSource->getMetaData()->getExifData().end (); it++)
			if (it->key () == "Exif.Image.ExifTag" || (it->groupName() != "Image" && it->groupName() != "SubImage1" && it->groupName() != "SubImage2"))
				ed.add (*it);

		// Add essential data to Image IFD
		Exiv2::ExifData::const_iterator it;
		if ((it=imgSource->getMetaData()->getExifData().findKey (Exiv2::ExifKey ("Exif.Image.Make"))) != imgSource->getMetaData()->getExifData().end())
			ed.add (*it);
		if ((it=imgSource->getMetaData()->getExifData().findKey (Exiv2::ExifKey ("Exif.Image.Model"))) != imgSource->getMetaData()->getExifData().end())
			ed.add (*it);
		if ((it=imgSource->getMetaData()->getExifData().findKey (Exiv2::ExifKey ("Exif.Image.DateTime"))) != imgSource->getMetaData()->getExifData().end())
			ed.add (*it);
		if ((it=imgSource->getMetaData()->getExifData().findKey (Exiv2::ExifKey ("Exif.Image.DateTimeOriginal"))) != imgSource->getMetaData()->getExifData().end())
			ed.add (*it);
	
		Exiv2::AsciiValue soft ("RawTherapee");
		ed.add (Exiv2::Exifdatum (Exiv2::ExifKey ("Exif.Image.Software"), &soft));

		// Apply changes requested by the GUI
		for (int i=0; i<procParams->exif.size (); i++) {
			if (procParams->exif[i].value == "#delete") {
				Exiv2::ExifData::iterator it = ed.findKey (Exiv2::ExifKey (String2StdString(procParams->exif[i].field)));
				if (it != ed.end())
					ed.erase (it);
			}
			else {
				Exiv2::AsciiValue nValue (String2StdString(procParams->exif[i].value));
				ed.add (Exiv2::Exifdatum (Exiv2::ExifKey (String2StdString(procParams->exif[i].field)), &nValue));
			}
		}
		
		// Assemble IPTC data
        for (int i=0; i<procParams->iptc.size (); i++)
			for (int j=0; j<procParams->iptc[i].values.size(); j++) {
				Exiv2::StringValue value (String2StdString(procParams->iptc[i].values[j]));
				id.add (Exiv2::IptcKey (String2StdString(procParams->iptc[i].field)), &value);
			}

		// Assemble XMP data
		procParams->addToXmp (xd);
	}
}


}
