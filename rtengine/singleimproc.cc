/*
 * singleimproc.cc
 *
 *  Created on: Oct 6, 2010
 *      Author: gabor
 */

#include "rtengine.h"
#include "processingjob.h"
#include "filterchain.h"
#include "buffer.h"
#include "multiimage.h"

namespace rtengine {

IImage16* SingleImageProcessor::process (ProcessingJob* pJob, ProgressListener* pListener, FinalImageListener* fiListener, int& errorCode) {

    if (fiListener) {
        ProcessingJob* currentJob = pJob;

        while (currentJob) {
            int errorCode;
            IImage16* img = SingleImageProcessor::process (currentJob, pListener, NULL, errorCode);
            if (errorCode && pListener)
                pListener->error ("Can not load input image.");
            currentJob = fiListener->imageReady (img);
        }
    }
    else {

        if (pListener) {
            pListener->setProgressStr ("Processing...");
            pListener->setProgress (0.0);
        }

        ProcessingJobImpl* job = (ProcessingJobImpl*)pJob;

        // load image, if not loaded so far
        InitialImage* ii = job->initialImage;
        if (!ii) {
            ii = InitialImage::load (job->fname, job->isRaw, &errorCode);
            if (errorCode) {
                ii->decreaseRef ();
                delete job;
                return NULL;
            }
        }

        // set up filter chain for processing
        FilterChain* fChain = new FilterChain (NULL, ii->getImageSource(), &job->pparams, true);
        std::set<ProcEvent> ev;
        ev.insert (EvAll);
        fChain->setupProcessing (ev, true);

        // create buffer, if necessary
        Buffer<int>* buffer = NULL;
        Dim bSize = fChain->getReqiredBufferSize ();
        if (bSize.nonZero())
            buffer = new Buffer<int> (bSize.width, bSize.height);

        // create worker image
        MultiImage* worker = NULL;
        Dim wSize = fChain->getReqiredWorkerSize ();
        if (wSize.nonZero())
            worker = new MultiImage (wSize.width, wSize.height);

        // perform processing
        fChain->process (ev, buffer, worker);

        return fChain->getFinalImage ();
    }
}
}
