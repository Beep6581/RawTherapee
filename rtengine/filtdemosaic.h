/*
 * filtdemosaic.h
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#ifndef FILTDEMOSAIC_H_
#define FILTDEMOSAIC_H_

#include "filter.h"

//
// D e m o s a i c   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

namespace rtengine {

class DemosaicFilterDescriptor : public FilterDescriptor {

	public:
        DemosaicFilterDescriptor ();
		void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern DemosaicFilterDescriptor demosaicFilterDescriptor;

class DemosaicFilter : public Filter {

        int border;

        void hphd_demosaic (MultiImage* si, MultiImage* ti, Buffer<float>* hpmap);
        void hphd_green (MultiImage* si, MultiImage* ti, Buffer<float>* hpmap);
        void hphd_horizontal (MultiImage* si, Buffer<float>* hpmap, int row_from, int row_to);
        void hphd_vertical (MultiImage* si, Buffer<float>* hpmap, int col_from, int col_to);
        void interpolate_rb_bilinear (MultiImage* si, MultiImage* ti);

        void correction_YIQ_LQ  (MultiImage* im, int times);
        void correction_YIQ_LQ_ (MultiImage* im, int row_from, int row_to);
        void convert_row_to_YIQ (float* r, float* g, float* b, float* Y, float* I, float* Q, int W);
        void convert_row_to_RGB (float* r, float* g, float* b, float* Y, float* I, float* Q, int W);

	public:
        DemosaicFilter ();

    	void      process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
        ImageView calculateTargetImageView (const ImageView& requestedImView);
        ImageView calculateSourceImageView (const ImageView& requestedImView);
        Dim       getFullImageSize ();
        Dim       getReqiredBufferSize ();
        void      reverseTransPoint (int x, int y, int& xv, int& yv);
        int       getTargetSkip (int nextInSkip);
};

}
#endif /* FILTDEMOSAIC_H_ */
