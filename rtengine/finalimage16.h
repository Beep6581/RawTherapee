/*
 * finalimage16.h
 *
 *  Created on: Feb 5, 2011
 *      Author: gabor
 */

#ifndef FINALIMAGE16_H_
#define FINALIMAGE16_H_

#include "rtcommon.h"

namespace rtengine {

class FinalImage16 {

	public:
		enum ErrorCodes {NoError=0, InvalidImage=1, UnknownFileExtension=2, LoadFailed=3, SaveFailed=4};
		enum JPEGSubSampling {JPEGSubSampling_411=1, JPEGSubSampling_420=2, JPEGSubSampling_422=3, JPEGSubSampling_444=4};
		enum PNGCompression {PNGDefault=1, PNGZBestSpeed=2, PNGZDefaultCompression=3, PNGZBestCompression=4, PNGZNoCompression=5};
		enum TIFFCompression {TIFFNoCompression=1, TIFFLZWCompression=2, TIFFDeflateCompression=3};

		virtual int getWidth () = 0;
		virtual int getHeight () = 0;
		virtual int getScanLineSize () = 0;
        virtual unsigned char* getData () = 0;

		virtual int saveAsPNG  (const String& fname, PNGCompression compr = PNGZDefaultCompression, bool bps16=true) = 0;
		virtual int saveAsJPEG (const String& fname, int quality = -1, JPEGSubSampling ss = JPEGSubSampling_420) = 0;
		virtual int saveAsTIFF (const String& fname, TIFFCompression compr = TIFFLZWCompression, bool bps16=true) = 0;
};
}
#endif /* FINALIMAGE16_H_ */
