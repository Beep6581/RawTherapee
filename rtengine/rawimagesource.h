/*
 * rawimagesource.h
 *
 *  Created on: Aug 22, 2010
 *      Author: gabor
 */

#ifndef RAWIMAGESOURCE_H_
#define RAWIMAGESOURCE_H_

#include "imagesource.h"
#include "rawimage.h"
#include "imagedata.h"
#include "colortemp.h"
#include "multiimage.h"
#include "imageview.h"

namespace rtengine {

class RawImageSource : public ImageSource {

    private:
	 	int border;
        String fileName;
        RawImage* img;
        ImageData* idata;
        ColorTemp autoWB;
        bool autoWBComputed;
        cmsHPROFILE embProfile;

	public:

        RawImageSource ();
        virtual ~RawImageSource ();

		// inherited from InitialImage
		String 		 		 getFileName () { return fileName; }
		cmsHPROFILE 		 getEmbeddedProfile ();
		const ImageMetaData* getMetaData () { return idata; }

      // inherited from ImageSource
        int 		load (const String& fileName, ProgressListener* listener = NULL);

        ColorTemp   getCamWB    ();
        ColorTemp   getAutoWB   ();
        ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue);
        void        getAEHistogram (unsigned int* histogram, int& histcompr);
        Matrix33	getCamToRGBMatrix ();
        Matrix33	getRGBToCamMatrix ();

    	Dim  getFullImageSize ();
    	void getImage (const ImageView& targetImageView, MultiImage* targetImage);

        bool isRaw ()       { return true; }
        bool isThumbnail () { return false; }
        double getScale ()  { return 1.0; }
};

}

#endif /* RAWIMAGESOURCE_H_ */
