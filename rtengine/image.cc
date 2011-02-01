#include "image.h"

namespace rtengine {

Image::Image (int width, int height) {

	bitmap = FreeImage_AllocateT (FIT_RGB16, width, height);
}

Image::~Image () {

	if (bitmap)
		FreeImage_Unload (bitmap);
}

Image::Image (FIBITMAP* img) : bitmap (img) {
}

void Image::getEmbeddedICCProfile (int& length, unsigned char*& pdata) {
	
	if (bitmap) {
		FIICCPROFILE* profile = FreeImage_GetICCProfile (bitmap);
		if (profile->data) {
			pdata = (unsigned char*)profile->data;
			length = profile->size;
			return;
		}
	}
	pdata = NULL;
	length = 0;
}

void Image::setEmbeddedICCProfile (int length, unsigned char* pdata) {
	
	if (bitmap)
		FreeImage_CreateICCProfile (bitmap, pdata, length);
}

int Image::getWidth () {
	
	if (bitmap)
		return FreeImage_GetWidth (bitmap);
	else
		return 0;
}

int Image::getScanLineSize () {
	
	if (bitmap)
		return FreeImage_GetPitch (bitmap);
	else
		return 0;
}

int Image::getHeight () {

	if (bitmap)
		return FreeImage_GetHeight (bitmap);
	else
		return 0;
}
        
unsigned char* Image::getData () {
	
	if (bitmap)
		return FreeImage_GetBits (bitmap);
	else
		return NULL;
}

Image* Image::load (const String& fname, bool fast) {

	// check if format is supported and load it
	FREE_IMAGE_FORMAT fif = FreeImage_GetFileType (fname.c_str(), 0);
	if (fif == FIF_UNKNOWN)
		fif = FreeImage_GetFIFFromFilename (fname.c_str());

	FIBITMAP* bitmap = NULL;

	if (fif == FIF_PNG || fif == FIF_TIFF)
		bitmap = FreeImage_Load (fif, fname.c_str(), 0);
	else if (fif == FIF_JPEG)
		bitmap = FreeImage_Load (FIF_JPEG, fname.c_str(), fast ? JPEG_FAST : JPEG_ACCURATE);
		
	if (!bitmap)
		return NULL;

	// convert the image to RGB16 representation
	unsigned width  = FreeImage_GetWidth (bitmap);
	unsigned height = FreeImage_GetHeight (bitmap);
	unsigned pitch  = FreeImage_GetPitch (bitmap);
	unsigned char* bits = FreeImage_GetBits (bitmap);
	
	FIBITMAP* bmp16 = NULL;
	// RGB16, just as we like it
	if (FreeImage_GetImageType(bitmap) == FIT_RGB16)
		bmp16 = bitmap;
	// RGBA16, remove alpha channel
	else if (FreeImage_GetImageType(bitmap) == FIT_RGBA16) {
		bmp16 = FreeImage_AllocateT (FIT_RGB16, width, height);
		unsigned pitch16 = FreeImage_GetPitch (bmp16);
		unsigned char* bits16 = FreeImage_GetBits (bmp16);
		for (int y = 0; y < height; y++) {
			FIRGBA16 *pixela = (FIRGBA16*)bits;
			FIRGB16  *pixel  = (FIRGB16*)bits16;
			for (int x = 0; x < width; x++) {
				pixel[x].red   = pixela[x].red;
				pixel[x].green = pixela[x].green;
				pixel[x].blue  = pixela[x].blue;
			}
			bits16 += pitch16;
			bits += pitch;
		}
		FIICCPROFILE* profile = FreeImage_GetICCProfile (bitmap);
		if (profile->data)
			FreeImage_CreateICCProfile (bmp16, profile->data, profile->size);
	}
	// RGB8, convert it to RGB16
	else if (FreeImage_GetImageType(bitmap) == FIT_BITMAP && FreeImage_GetColorType(bitmap) == FIC_RGB && FreeImage_GetBPP(bitmap) == 24) {
		bmp16 = FreeImage_AllocateT (FIT_RGB16, width, height);
		unsigned pitch16 = FreeImage_GetPitch (bmp16);
		unsigned char* bits16 = FreeImage_GetBits (bmp16);
		for (int y = 0; y < height; y++) {
			RGBTRIPLE *pixel8 = (RGBTRIPLE*)bits;
			FIRGB16  *pixel16  = (FIRGB16*)bits16;
			for (int x = 0; x < width; x++) {
				pixel16[x].red   = pixel8[x].rgbtRed << 8;
				pixel16[x].green = pixel8[x].rgbtGreen << 8;
				pixel16[x].blue  = pixel8[x].rgbtBlue << 8;
			}
			bits16 += pitch16;
			bits += pitch;
		}
		FIICCPROFILE* profile = FreeImage_GetICCProfile (bitmap);
		if (profile->data)
			FreeImage_CreateICCProfile (bmp16, profile->data, profile->size);
	}
	// RGBA8, convert it to RGB16
	else if (FreeImage_GetImageType(bitmap) == FIT_BITMAP && FreeImage_GetColorType(bitmap) == FIC_RGBALPHA && FreeImage_GetBPP(bitmap) == 32) {
		bmp16 = FreeImage_AllocateT (FIT_RGB16, width, height);
		unsigned pitch16 = FreeImage_GetPitch (bmp16);
		unsigned char* bits16 = FreeImage_GetBits (bmp16);
		for (int y = 0; y < height; y++) {
			RGBQUAD *pixel8 = (RGBQUAD*)bits;
			FIRGB16 *pixel16 = (FIRGB16*)bits16;
			for (int x = 0; x < width; x++) {
				pixel16[x].red   = pixel8[x].rgbRed << 8;
				pixel16[x].green = pixel8[x].rgbGreen << 8;
				pixel16[x].blue  = pixel8[x].rgbBlue << 8;
			}
			bits16 += pitch16;
			bits += pitch;
		}
		FIICCPROFILE* profile = FreeImage_GetICCProfile (bitmap);
		if (profile->data)
			FreeImage_CreateICCProfile (bmp16, profile->data, profile->size);
	}
	
	// if none, return NULL
	if (bmp16 == NULL) {
		FreeImage_Unload (bitmap);
		return NULL;
	}
	// if conversion was needed, delete the original
	else if (bmp16 != bitmap) {
		FreeImage_Unload (bitmap);
		bitmap = bmp16;
	}
	
	// return
	return new Image (bitmap);
}

int Image::save (const String& fname) {
	
	FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename (fname.c_str());

	if (fif == FIF_PNG)
		return saveAsPNG (fname);
	else if (fif == FIF_JPEG)
		return saveAsJPEG (fname);
	else if (fif == FIF_TIFF)
		return saveAsTIFF (fname);
	else
		return Image::UnknownFileExtension;
}

FIBITMAP* Image::convertTo24bpp () {
	
	FIBITMAP* bmp8 = FreeImage_Allocate (getWidth(), getHeight(), 24);

	unsigned width   = FreeImage_GetWidth (bitmap);
	unsigned height  = FreeImage_GetHeight (bitmap);
	unsigned pitch16 = FreeImage_GetPitch (bitmap);
	unsigned pitch8  = FreeImage_GetPitch (bmp8);

	unsigned char *bits16 = FreeImage_GetBits (bitmap);
	unsigned char *bits8  = FreeImage_GetBits (bmp8);
	
	for (int y = 0; y < height; y++) {
		FIRGB16 *pixel16  = (FIRGB16*)bits16;
		RGBTRIPLE *pixel8 = (RGBTRIPLE*)bits8;
		for (int x = 0; x < width; x++) {
			pixel8[x].rgbtRed   = pixel16[x].red >> 8;
			pixel8[x].rgbtGreen = pixel16[x].green >> 8;
			pixel8[x].rgbtBlue  = pixel16[x].blue >> 8;
		}
		bits16 += pitch16;
		bits8 += pitch8;
	}
	
	FIICCPROFILE* profile = FreeImage_GetICCProfile (bitmap);
	if (profile->data)
		FreeImage_CreateICCProfile (bmp8, profile->data, profile->size);
	
	return bmp8;
}
		
int Image::saveAsPNG  (const String& fname, PNGCompression compr, bool bps16) {
	
	if (bitmap) {
		FIBITMAP* bmp;
		if (bps16)
			bmp = bitmap;
		else
			bmp = convertTo24bpp ();
		
		int flags = 0;
		if (compr == PNGDefault)
			flags |= PNG_DEFAULT;
		else if (compr == PNGZBestSpeed)
			flags |= PNG_Z_BEST_SPEED;
		else if (compr == PNGZDefaultCompression)
			flags |= PNG_Z_DEFAULT_COMPRESSION;
		else if (compr == PNGZBestCompression)
			flags |= PNG_Z_BEST_COMPRESSION;
		else if (compr == PNGZNoCompression)
			flags |= PNG_Z_NO_COMPRESSION;
		
		bool success = FreeImage_Save (FIF_PNG, bmp, fname.c_str(), flags);

		if (bmp != bitmap)
			FreeImage_Unload (bmp);

		if (success) {
			writeMetadata (fname);
			return Image::NoError;
		}
		else 
			return Image::SaveFailed;
	}
	else
		return Image::InvalidImage;
}

int Image::saveAsJPEG (const String& fname, int quality, JPEGSubSampling ss) {
	
	if (bitmap) {
		FIBITMAP* bmp8 = convertTo24bpp ();
		
		int flags = 0;
		if (quality == -1)
			flags |= JPEG_DEFAULT;
		else
			flags += quality;
		if (ss == Image::JPEGSubSampling_411)
			flags |= JPEG_SUBSAMPLING_411;
		else if (ss == Image::JPEGSubSampling_420)
			flags |= JPEG_SUBSAMPLING_420;
		else if (ss == Image::JPEGSubSampling_422)
			flags |= JPEG_SUBSAMPLING_422;
		else if (ss == Image::JPEGSubSampling_444)
			flags |= JPEG_SUBSAMPLING_444;
		
		bool success = FreeImage_Save (FIF_JPEG, bmp8, fname.c_str(), flags);
		FreeImage_Unload (bmp8);
		
		if (success) {
			writeMetadata (fname);
			return Image::NoError;
		}
		else 
			return Image::SaveFailed;
	}
	else
		return Image::InvalidImage;
}

int Image::saveAsTIFF (const String& fname, TIFFCompression compr, bool bps16) {

	if (bitmap) {
		FIBITMAP* bmp;
		if (bps16)
			bmp = bitmap;
		else
			bmp = convertTo24bpp ();
		
		int flags = 0;
		if (compr == TIFFNoCompression)
			flags |= TIFF_NONE;
		else if (compr == TIFFLZWCompression)
			flags |= TIFF_LZW;
		else if (compr == TIFFDeflateCompression)
			flags |= TIFF_DEFLATE;
		
		bool success = FreeImage_Save (FIF_TIFF, bmp, fname.c_str(), flags);

		if (bmp != bitmap)
			FreeImage_Unload (bmp);

		if (success) {
			writeMetadata (fname);
			return Image::NoError;
		}
		else 
			return Image::SaveFailed;
	}
	else
		return Image::InvalidImage;
}

void Image::setMetadata (const Exiv2::ExifData& ed, const Exiv2::IptcData& id, const Exiv2::XmpData& xd) {
	
	exifData = ed;
	iptcData = id;
	xmpData  = xd;
}

void Image::writeMetadata (const String& fname) {

	try {
		// open image
		Exiv2::Image::AutoPtr image = Exiv2::ImageFactory::open (fname);
		
		// read original metadata
		image->readMetadata();
		Exiv2::ExifData& origExifData = image->exifData ();
		Exiv2::IptcData& origIptcData = image->iptcData ();
		Exiv2::XmpData&  origXmpData  = image->xmpData ();
		
		// update exif metadata
		for (Exiv2::ExifData::const_iterator i = exifData.begin(); i!= exifData.end (); i++)
			origExifData.add (*i);
		
		// update iptc metadata
		for (Exiv2::IptcData::const_iterator i = iptcData.begin(); i!= iptcData.end (); i++)
			origIptcData.add (*i);
		
		// update xmp metadata
		for (Exiv2::XmpData::const_iterator i = xmpData.begin(); i!= xmpData.end (); i++)
			origXmpData.add (*i);
		
		// write back metadata
		origExifData.sortByTag ();
		image->writeMetadata ();
	}
	catch (const Exiv2::Error& e) {
	}
}

}
