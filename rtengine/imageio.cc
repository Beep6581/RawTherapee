/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <png.h>
#include <glib/gstdio.h>
#include <imageio.h>
#include <tiff.h>
#include <tiffio.h>
#include <stdio.h>
#include <string.h>
extern "C" {
#include "jdatasrc.c"
#include <jpeglib.h>
#include <iccjpeg.h>
}
#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif
#include <iptcpairs.h>
#include <libiptcdata/iptc-jpeg.h>

using namespace rtengine;
using namespace rtengine::procparams;

Glib::Mutex& IImage::getMutex () {}
cmsHPROFILE IImage::getProfile () {}
int IImage::getWidth ()  {}
int IImage::getHeight () {}
int IImage::getBitsPerPixel () {}
int IImage::saveToFile (Glib::ustring fname) {}
int IImage::saveAsPNG  (Glib::ustring fname, int compression, int bps) { }
int IImage::saveAsJPEG (Glib::ustring fname, int quality) {}
int IImage::saveAsTIFF (Glib::ustring fname, int bps) {}
void IImage::setSaveProgressListener (ProgressListener* pl) {} 
void IImage::free () {}

Glib::ustring ImageIO::errorMsg[6] = {"Success", "Cannot read file.", "Invalid header.","Error while reading header.","File reading error", "Image format not supported."};

void ImageIO::setMetadata (const rtexif::TagDirectory* eroot, const std::vector<ExifPair>& exif, const std::vector<IPTCPair>& iptcc) {

    // store exif info
    exifChange.resize (exif.size());
    for (int i=0; i<exif.size(); i++) {
        exifChange[i].first  = exif[i].field;
        exifChange[i].second = exif[i].value;
    }
    delete exifRoot;

    exifRoot = NULL;
    if (eroot)
        exifRoot = ((rtexif::TagDirectory*)eroot)->clone (NULL);

    if (iptc) 
        iptc_data_free (iptc);
    iptc = NULL;        
    
    // build iptc structures for libiptcdata
    if (iptcc.size()==0)
        return;
        
    iptc = iptc_data_new ();
    for (int i=0; i<iptcc.size(); i++) {
        if (iptcc[i].field == "Keywords" && iptcc[i].values.size()>0) {
            for (int j=0; j<iptcc[i].values.size(); j++) {
                IptcDataSet * ds = iptc_dataset_new ();
                iptc_dataset_set_tag (ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS);
                std::string loc = Glib::locale_from_utf8(iptcc[i].values[j]);
                iptc_dataset_set_data (ds, (unsigned char*)loc.c_str(), MIN(64,loc.size()), IPTC_DONT_VALIDATE);
                iptc_data_add_dataset (iptc, ds);
                iptc_dataset_unref (ds);
            }
            continue;
        }
        else if (iptcc[i].field == "SupplementalCategories" && iptcc[i].values.size()>0) {
            for (int j=0; j<iptcc[i].values.size(); j++) {
                IptcDataSet * ds = iptc_dataset_new ();
                iptc_dataset_set_tag (ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY);
                std::string loc = Glib::locale_from_utf8(iptcc[i].values[j]);
                iptc_dataset_set_data (ds, (unsigned char*)loc.c_str(), MIN(32,loc.size()), IPTC_DONT_VALIDATE);
                iptc_data_add_dataset (iptc, ds);
                iptc_dataset_unref (ds);
            }
            continue;
        }
        for (int j=0; j<16; j++)
            if (iptcc[i].field == strTags[j].field && iptcc[i].values.size()>0) {
                IptcDataSet * ds = iptc_dataset_new ();
                iptc_dataset_set_tag (ds, IPTC_RECORD_APP_2, strTags[j].tag);
                std::string loc = Glib::locale_from_utf8(iptcc[i].values[0]);
                iptc_dataset_set_data (ds, (unsigned char*)loc.c_str(), MIN(strTags[j].size,loc.size()), IPTC_DONT_VALIDATE);
                iptc_data_add_dataset (iptc, ds);
                iptc_dataset_unref (ds);
            }
    }
    iptc_data_sort (iptc);
}

void ImageIO::setOutputProfile  (char* pdata, int plen) {

    delete [] profileData;
    if (pdata) {
        profileData = new char [plen];
        memcpy (profileData, pdata, plen);
    }
    else
        profileData = NULL;
    profileLength = plen;
}

ImageIO::~ImageIO () {

    if (embProfile)
        cmsCloseProfile(embProfile);
    delete loadedProfileData;
    delete exifRoot;
    delete profileData;
}

void png_read_data(png_struct_def  *png_ptr, unsigned char *data, size_t length);
void png_write_data(png_struct_def *png_ptr, unsigned char *data, size_t length);
void png_flush(png_struct_def *png_ptr);

int ImageIO::loadPNG  (Glib::ustring fname) {

    FILE *file = g_fopen (fname.c_str(),"rb");
    if (!file) 
      return IMIO_CANNOTREADFILE;

    if (pl) {
      pl->setProgressStr ("Loading PNG file...");
      pl->setProgress (0.0);
    }

	//reading PNG header
	unsigned char header[8];
	fread (header, 1, 8, file);
	if (!png_check_sig (header, 8)) {
		fclose(file);
		return IMIO_HEADERERROR;
	}
	//initializing main structures
	png_structp png = png_create_read_struct (PNG_LIBPNG_VER_STRING, 0, 0, 0);
	if (!png) {
		fclose (file);
		return IMIO_HEADERERROR;
	}
	png_infop info = png_create_info_struct (png);
	png_infop end_info = png_create_info_struct (png);
	if (!end_info || !info) {
		png_destroy_read_struct (&png, &info, &end_info);
		fclose (file);
		return IMIO_HEADERERROR;
	}

	if (setjmp (png_jmpbuf(png))) {
		png_destroy_read_struct (&png, &info, &end_info);
		fclose (file);
		return IMIO_READERROR;
    }

	//set up png read
    png_set_read_fn (png, file, png_read_data);
	png_set_sig_bytes (png,8);

	png_read_info(png,info);

    embProfile = NULL;

	//retrieving image information
	unsigned long width,height;
	int bit_depth,color_type,interlace_type,compression_type,filter_method;
	png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

	//converting to 32bpp format
	if (color_type==PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
	
	if (color_type==PNG_COLOR_TYPE_GRAY || color_type==PNG_COLOR_TYPE_GRAY_ALPHA)
          png_set_gray_to_rgb(png);

	if (png_get_valid(png,info,PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);

	if (interlace_type!=PNG_INTERLACE_NONE) {
    	  png_destroy_read_struct (&png, &info, &end_info);
		  fclose (file);
          return IMIO_VARIANTNOTSUPPORTED;
    }

    if (color_type & PNG_COLOR_MASK_ALPHA)
        png_set_strip_alpha(png);

	//setting gamma
	double gamma;
	if (png_get_gAMA(png,info,&gamma))
		png_set_gamma(png, 2.0, gamma);
	else
		png_set_gamma(png,2.0, 0.45455);

    int bps = getBPS ();


//	if (bps==8 && bit_depth==16) png_set_strip_16(png);

	//updating png info struct
	png_read_update_info(png,info);
	png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

        if (color_type & PNG_COLOR_MASK_ALPHA)
          png_set_strip_alpha(png);

	png_read_update_info(png,info);
	png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,&interlace_type,
		&compression_type, &filter_method);

    allocate (width, height);

    int rowlen = width*3*bit_depth/8;
    unsigned char *row = new unsigned char [rowlen];

	for (unsigned int i=0;i<height;i++) {

  	    png_read_row (png, (png_byte*)row, NULL);
  	    if (bit_depth==16) {  // convert scanline to host byte order
  	        unsigned short* srow = (unsigned short*)row;
  	        for (int j=0; j<width*3; j++)
  	            srow[j] = ntohs (srow[j]);
  	    }
        setScanline (i, row, bit_depth);

//        if (bps==16 && bit_depth==8)
//            setScanline (i, row, 8);
//        else
//            setScanline (i, row, bps);

        if (pl && !(i%100))
            pl->setProgress ((double)(i+1)/height);
    }

	png_read_end (png, 0);
	png_destroy_read_struct (&png, &info, &end_info);
	
	delete [] row;
	fclose(file);
    if (pl) {
        pl->setProgressStr ("Ready.");
        pl->setProgress (1.0);
    }
    return IMIO_SUCCESS;
}

extern jmp_buf jpeg_jmp_buf;

int ImageIO::loadJPEG (Glib::ustring fname) {

	FILE *file=g_fopen(fname.c_str(),"rb");
	if (!file) 
        return IMIO_CANNOTREADFILE;

    jpeg_decompress_struct cinfo;
    jpeg_error_mgr jerr;
    cinfo.err = my_jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    my_jpeg_stdio_src (&cinfo,file);

    if (pl) {
      pl->setProgressStr ("Loading JPEG file...");
      pl->setProgress (0.0);

    }

    setup_read_icc_profile (&cinfo);

    if (!setjmp(jpeg_jmp_buf)) {
    	jpeg_stdio_src(&cinfo,file);
    	jpeg_read_header(&cinfo, TRUE);

        unsigned int proflen;    	
        delete loadedProfileData;
        loadedProfileData = NULL;
        bool hasprofile = read_icc_profile (&cinfo, (JOCTET**)&loadedProfileData, (unsigned int*)&loadedProfileLength);
    	if (hasprofile) 
    	    embProfile = cmsOpenProfileFromMem (loadedProfileData, loadedProfileLength);
    	else 
            embProfile = NULL;
    	
	    jpeg_start_decompress(&cinfo);

        int width = cinfo.output_width;
        int height = cinfo.output_height;

        allocate (width, height);

	    unsigned char *row=new unsigned char[width*3];
	    while (cinfo.output_scanline < height) {
		    if (jpeg_read_scanlines(&cinfo,&row,1) < 1) {
  	          jpeg_finish_decompress(&cinfo);
 	          jpeg_destroy_decompress(&cinfo);
 	          delete [] row;
              return IMIO_READERROR;
            }
            setScanline (cinfo.output_scanline-1, row, 8);

            if (pl && !(cinfo.output_scanline%100))
              pl->setProgress ((double)(cinfo.output_scanline)/cinfo.output_height);
	    }
	    delete [] row;

	    jpeg_finish_decompress(&cinfo);
	    jpeg_destroy_decompress(&cinfo);
	    fclose(file);
        if (pl) {
          pl->setProgressStr ("Ready.");
          pl->setProgress (1.0);
        }
        return IMIO_SUCCESS;
    }
    else {
        jpeg_destroy_decompress(&cinfo);
        return IMIO_READERROR;
    }
}

int ImageIO::loadTIFF (Glib::ustring fname) {

#ifdef WIN32
    wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
    TIFF* in = TIFFOpenW (wfilename, "r");
    g_free (wfilename);
#else
    TIFF* in = TIFFOpen(fname.c_str(), "r");
#endif
	if (in == NULL) 
          return IMIO_CANNOTREADFILE;

    if (pl) {
        pl->setProgressStr ("Loading TIFF file...");
        pl->setProgress (0.0);
    }
    
    int width, height;
	TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &width);
	TIFFGetField(in, TIFFTAG_IMAGELENGTH, &height);

    uint16 bitspersample, samplesperpixel;
	TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
    uint16 photometric;
	if (!TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photometric) ||
	    photometric != PHOTOMETRIC_RGB || samplesperpixel < 3) {
        TIFFClose(in);
		return IMIO_VARIANTNOTSUPPORTED;
	}

    uint16 config;
	TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);
	if (config != PLANARCONFIG_CONTIG) {
        TIFFClose(in);
		return IMIO_VARIANTNOTSUPPORTED;
	}

    char* profdata;
	delete loadedProfileData;
    loadedProfileData = NULL;
   	if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &loadedProfileLength, &profdata)) {
   	    embProfile = cmsOpenProfileFromMem (profdata, loadedProfileLength);
        loadedProfileData = new char [loadedProfileLength];
        memcpy (loadedProfileData, profdata, loadedProfileLength);
    }
   	else 
        embProfile = NULL;
        

    allocate (width, height);

    unsigned char* linebuffer = new unsigned char[TIFFScanlineSize(in)];
    for (int row = 0; row < height; row++) {
        if (TIFFReadScanline(in, linebuffer, row, 0) <0) {
          TIFFClose(in);
          delete [] linebuffer;
          return IMIO_READERROR;
        }
        if (samplesperpixel>3) 
            for (int i=0; i<width; i++) 
                memcpy (linebuffer+i*3*bitspersample/8, linebuffer+i*samplesperpixel*bitspersample/8, 3*bitspersample/8);
        setScanline (row, linebuffer, bitspersample);
              
        if (pl && !(row%100))
            pl->setProgress ((double)(row+1)/height);
    }
    TIFFClose(in);
    delete [] linebuffer;

    if (pl) {
      pl->setProgressStr ("Ready.");
      pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}


int ImageIO::savePNG  (Glib::ustring fname, int compression, int bps) {

    FILE* file=g_fopen(fname.c_str (),"wb");

    if (!file) 
      return IMIO_CANNOTREADFILE;

    if (pl) {
      pl->setProgressStr ("Saving PNG file...");
      pl->setProgress (0.0);
    }

	png_structp png = png_create_write_struct (PNG_LIBPNG_VER_STRING,0,0,0);
	if (!png) {
		fclose (file);
		return IMIO_HEADERERROR;
	}
	png_infop info = png_create_info_struct(png);
	if (!info) {
		png_destroy_write_struct (&png,0);
		fclose (file);
		return IMIO_HEADERERROR;
        }

	if (setjmp(png_jmpbuf(png))) {
		png_destroy_write_struct (&png,&info);
		fclose(file);
		return IMIO_READERROR;
    }

	png_set_write_fn (png, file, png_write_data, png_flush);	

	png_set_compression_level(png,compression);

    int width = getW ();
    int height = getH ();
    if (bps<0)
        bps = getBPS ();

	png_set_IHDR(png, info, width, height, bps, PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_BASE);


    int rowlen = width*3*bps/8;
    unsigned char *row = new unsigned char [rowlen];

	png_write_info(png,info);
	for (unsigned int i=0;i<height;i++) {
        getScanline (i, row, bps);
        if (bps==16) {
            // convert to network byte order
            for (int j=0; j<width*6; j+=2) {
                unsigned char tmp = row[j];
                row[j] = row[j+1];
                row[j+1] = tmp;
            }
        }
        png_write_row (png, (png_byte*)row);
        
        if (pl && !(i%100))
            pl->setProgress ((double)(i+1)/height);
    }

	png_write_end(png,info);
	png_destroy_write_struct(&png,&info);

    delete [] row;
	fclose (file);

    if (pl) {
        pl->setProgressStr ("Ready.");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}


int ImageIO::saveJPEG (Glib::ustring fname, int quality) {

	jpeg_compress_struct cinfo;
	jpeg_error_mgr jerr;
	
	cinfo.err = jpeg_std_error (&jerr);
	jpeg_create_compress (&cinfo);

	FILE *file = g_fopen (fname.c_str (), "wb");

	if (!file)
          return IMIO_CANNOTREADFILE;

    if (pl) {
        pl->setProgressStr ("Saving JPEG file...");
        pl->setProgress (0.0);
    }

	jpeg_stdio_dest (&cinfo, file);

	int width = getW ();
    int height = getH ();

    cinfo.image_width  = width;
	cinfo.image_height = height;
	cinfo.in_color_space = JCS_RGB;
	cinfo.input_components = 3;
	jpeg_set_defaults (&cinfo);
    cinfo.write_JFIF_header = FALSE;
    

	if (quality>=0 && quality<=100) 
	    jpeg_set_quality (&cinfo, quality, true);

	jpeg_start_compress(&cinfo, TRUE);

    // buffer for exif and iptc markers
    unsigned char buffer[165535];
    unsigned int size;
    // assemble and write exif marker
   if (exifRoot) {
        int size = rtexif::ExifManager::createJPEGMarker (exifRoot, exifChange, cinfo.image_width, cinfo.image_height, buffer);
        if (size>0 && size<65530)
            jpeg_write_marker(&cinfo, JPEG_APP0+1, buffer, size);
    }
    // assemble and write iptc marker
    if (iptc) {
        unsigned char* iptcdata;
        bool error = false;
        if (iptc_data_save (iptc, &iptcdata, &size)) {
            if (iptcdata)
                iptc_data_free_buf (iptc, iptcdata);
            error = true;
        }
        int bytes = 0;
        if (!error && (bytes = iptc_jpeg_ps3_save_iptc (NULL, 0, iptcdata, size, buffer, 65532)) < 0) {
            if (iptcdata)
                iptc_data_free_buf (iptc, iptcdata);
            error = true;
        }
        if (!error)
            jpeg_write_marker(&cinfo, JPEG_APP0+13, buffer, bytes);
    }
    // write icc profile to the output
    if (profileData)
        write_icc_profile (&cinfo, (JOCTET*)profileData, profileLength);

    // write image data
    int rowlen = width*3;
    unsigned char *row = new unsigned char [rowlen];

	while (cinfo.next_scanline < cinfo.image_height) {
        
        getScanline (cinfo.next_scanline, row, 8);
        
		if (jpeg_write_scanlines (&cinfo, &row, 1) < 1) {
            jpeg_finish_compress (&cinfo);
	        jpeg_destroy_compress (&cinfo);
	        fclose (file);
            return IMIO_READERROR;
        }

        if (pl && !(cinfo.next_scanline%100))
            pl->setProgress ((double)(cinfo.next_scanline)/cinfo.image_height);
	}

	jpeg_finish_compress (&cinfo);
	jpeg_destroy_compress (&cinfo);

    delete [] row;
	fclose (file);

    if (pl) {
        pl->setProgressStr ("Ready.");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

int ImageIO::saveTIFF (Glib::ustring fname, int bps) {

	int width = getW ();
    int height = getH ();
    
    if (bps<0)
        bps = getBPS ();

    int lineWidth = width*3*bps/8;
    unsigned char* linebuffer = new unsigned char[lineWidth];

    if (exifRoot) {
        FILE *file = g_fopen (fname.c_str (), "wb");

        if (!file)
            return IMIO_CANNOTREADFILE;           
            
        if (pl) {
            pl->setProgressStr ("Saving TIFF file...");
            pl->setProgress (0.0);
        }
        
        // buffer for the exif and iptc
        unsigned char buffer[165535];
        unsigned char* iptcdata = NULL;
        unsigned int iptclen = 0;
        if (iptc && iptc_data_save (iptc, &iptcdata, &iptclen) && iptcdata) {
            iptc_data_free_buf (iptc, iptcdata);
            iptcdata = NULL;
        }
        int size = rtexif::ExifManager::createTIFFHeader (exifRoot, exifChange, width, height, bps, profileData, profileLength, (char*)iptcdata, iptclen, buffer);
        if (iptcdata) 
            iptc_data_free_buf (iptc, iptcdata);
        if (size>0 && size<165530)
            fwrite (buffer, size, 1, file);

        bool needsReverse = bps==16 && exifRoot->getOrder()==rtexif::MOTOROLA;
        
        for (int i=0; i<height; i++) {
            getScanline (i, linebuffer, bps);
            if (needsReverse)
                for (int i=0; i<lineWidth; i+=2) {
                    char c = linebuffer[i];
                    linebuffer[i] = linebuffer[i+1];
                    linebuffer[i+1] = c;
                }
            fwrite (linebuffer, lineWidth, 1, file);
            if (pl && !(i%100))
                pl->setProgress ((double)(i+1)/height);
        }
        
        fclose (file);
    }
    else {
        #ifdef WIN32
        wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
        TIFF* out = TIFFOpenW (wfilename, "w");
        g_free (wfilename);
        #else
        TIFF* out = TIFFOpen(fname.c_str(), "w");
        #endif
        if (!out) 
            return IMIO_CANNOTREADFILE;

        if (pl) {
            pl->setProgressStr ("Saving TIFF file...");
            pl->setProgress (0.0);
        }
    
        TIFFSetField (out, TIFFTAG_IMAGEWIDTH, width);
        TIFFSetField (out, TIFFTAG_IMAGELENGTH, height);
        TIFFSetField (out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField (out, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField (out, TIFFTAG_BITSPERSAMPLE, bps);
        TIFFSetField (out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
        TIFFSetField (out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        TIFFSetField (out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);

        if (profileData) 
            TIFFSetField (out, TIFFTAG_ICCPROFILE, profileLength, profileData);   

        for (int row = 0; row < height; row++) {
            getScanline (row, linebuffer, bps);
        
            if (TIFFWriteScanline (out, linebuffer, row, 0) < 0) {
                TIFFClose (out);
                delete [] linebuffer;
                return IMIO_READERROR;
            }
            if (pl && !(row%100))
                pl->setProgress ((double)(row+1)/height);
        }
        TIFFClose (out);
    }

    delete [] linebuffer;
    if (pl) {
        pl->setProgressStr ("Ready.");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

// PNG read and write routines:

void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length) {
   png_size_t check;

   /* fread() returns 0 on error, so it is OK to store this in a png_size_t
    * instead of an int, which is what fread() actually returns.
    */
   check = (png_size_t)fread(data, (png_size_t)1, length, (FILE *)png_ptr->io_ptr);

   if (check != length)
   {
      png_error(png_ptr, "Read Error");
   }
}

void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length) {
   png_uint_32 check;

   check = fwrite(data, 1, length, (FILE *)(png_ptr->io_ptr));
   if (check != length)
   {
      png_error(png_ptr, "Write Error");
   }
}

void png_flush(png_structp png_ptr) {
   FILE *io_ptr;
   io_ptr = (FILE *)CVT_PTR((png_ptr->io_ptr));
   if (io_ptr != NULL)
      fflush(io_ptr);
}

int ImageIO::load (Glib::ustring fname) {

  int lastdot = fname.find_last_of ('.');

  if (!fname.casefold().compare (lastdot, 4, ".png"))
    return loadPNG (fname);
  else if (!fname.casefold().compare (lastdot, 4, ".jpg"))
    return loadJPEG (fname);
  else if (!fname.casefold().compare (lastdot, 4, ".tif"))
    return loadTIFF (fname);
}

int ImageIO::save (Glib::ustring fname) {

  int lastdot = fname.find_last_of ('.');

  if (!fname.casefold().compare (lastdot, 4, ".png"))
    return savePNG (fname);
  else if (!fname.casefold().compare (lastdot, 4, ".jpg"))
    return saveJPEG (fname);
  else if (!fname.casefold().compare (lastdot, 4, ".tif"))
    return saveTIFF (fname);
}

