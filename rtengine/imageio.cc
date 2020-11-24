/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2010 Oliver Duis <www.oliverduis.de>
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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <png.h>
#include <glib/gstdio.h>
#include <tiff.h>
#include <tiffio.h>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <libiptcdata/iptc-jpeg.h>
#include "rt_math.h"
#include "procparams.h"
#include "utils.h"
#include "../rtgui/options.h"
#include "../rtgui/version.h"
#include "../rtexif/rtexif.h"

#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

#include "imageio.h"
#include "iptcpairs.h"
#include "iccjpeg.h"
#include "color.h"

#include "jpeg.h"

using namespace std;
using namespace rtengine;
using namespace rtengine::procparams;

namespace
{

// Opens a file for binary writing and request exclusive lock (cases were you need "wb" mode plus locking)
FILE* g_fopen_withBinaryAndLock(const Glib::ustring& fname)
{

#ifdef WIN32

    // Use native function to disallow sharing, i.e. lock the file for exclusive access.
    // This is important to e.g. prevent Windows Explorer from crashing RT due to concurrently scanning an image file.
    std::unique_ptr<wchar_t, GFreeFunc> wfname (reinterpret_cast<wchar_t*>(g_utf8_to_utf16 (fname.c_str (), -1, NULL, NULL, NULL)), g_free);

    HANDLE hFile = CreateFileW ( wfname.get (), GENERIC_READ | GENERIC_WRITE, 0 /* no sharing allowed */, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    FILE* f = nullptr;

    if (hFile != INVALID_HANDLE_VALUE) {
        f = _fdopen (_open_osfhandle ((intptr_t)hFile, 0), "wb");
    }

#else

    FILE* f = ::g_fopen (fname.c_str (), "wb");

#endif

    return f;
}

}

Glib::ustring ImageIO::errorMsg[6] = {"Success", "Cannot read file.", "Invalid header.", "Error while reading header.", "File reading error", "Image format not supported."};

// For only copying the raw input data
void ImageIO::setMetadata (const rtexif::TagDirectory* eroot)
{
    if (exifRoot != nullptr) {
        delete exifRoot;
        exifRoot = nullptr;
    }

    if (eroot) {
        rtexif::TagDirectory* td = eroot->clone (nullptr);

        // make IPTC and XMP pass through
        td->keepTag(0x83bb);  // IPTC
        td->keepTag(0x02bc);  // XMP

        exifRoot = td;
    }
}

// For merging with RT specific data
void ImageIO::setMetadata (const rtexif::TagDirectory* eroot, const rtengine::procparams::ExifPairs& exif, const rtengine::procparams::IPTCPairs& iptcc)
{

    // store exif info
    exifChange->clear();
    *exifChange = exif;

    if (exifRoot != nullptr) {
        delete exifRoot;
        exifRoot = nullptr;
    }

    if (eroot) {
        exifRoot = eroot->clone (nullptr);
    }

    if (iptc != nullptr) {
        iptc_data_free (iptc);
        iptc = nullptr;
    }

    // build iptc structures for libiptcdata
    if (iptcc.empty()) {
        return;
    }

    iptc = iptc_data_new ();

    const unsigned char utf8Esc[] = {0x1B, '%', 'G'};
    IptcDataSet * ds = iptc_dataset_new ();
    iptc_dataset_set_tag (ds, IPTC_RECORD_OBJECT_ENV, IPTC_TAG_CHARACTER_SET);
    iptc_dataset_set_data (ds, utf8Esc, 3, IPTC_DONT_VALIDATE);
    iptc_data_add_dataset (iptc, ds);
    iptc_dataset_unref (ds);

    for (rtengine::procparams::IPTCPairs::const_iterator i = iptcc.begin(); i != iptcc.end(); ++i) {
        if (i->first == "Keywords" && !(i->second.empty())) {
            for (unsigned int j = 0; j < i->second.size(); j++) {
                IptcDataSet * ds = iptc_dataset_new ();
                iptc_dataset_set_tag (ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS);
                iptc_dataset_set_data (ds, (const unsigned char*)i->second.at(j).c_str(), min(static_cast<size_t>(64), i->second.at(j).bytes()), IPTC_DONT_VALIDATE);
                iptc_data_add_dataset (iptc, ds);
                iptc_dataset_unref (ds);
            }

            continue;
        } else if (i->first == "SupplementalCategories" && !(i->second.empty())) {
            for (unsigned int j = 0; j < i->second.size(); j++) {
                IptcDataSet * ds = iptc_dataset_new ();
                iptc_dataset_set_tag (ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY);
                iptc_dataset_set_data (ds, (const unsigned char*)i->second.at(j).c_str(), min(static_cast<size_t>(32), i->second.at(j).bytes()), IPTC_DONT_VALIDATE);
                iptc_data_add_dataset (iptc, ds);
                iptc_dataset_unref (ds);
            }

            continue;
        }

        for (int j = 0; j < 16; j++)
            if (i->first == strTags[j].field && !(i->second.empty())) {
                IptcDataSet * ds = iptc_dataset_new ();
                iptc_dataset_set_tag (ds, IPTC_RECORD_APP_2, strTags[j].tag);
                iptc_dataset_set_data (ds, (const unsigned char*)i->second.at(0).c_str(), min(strTags[j].size, i->second.at(0).bytes()), IPTC_DONT_VALIDATE);
                iptc_data_add_dataset (iptc, ds);
                iptc_dataset_unref (ds);
            }
    }

    iptc_data_sort (iptc);
}

void ImageIO::setOutputProfile  (const char* pdata, int plen)
{

    delete [] profileData;

    if (pdata) {
        profileData = new char [plen];
        memcpy (profileData, pdata, plen);
    } else {
        profileData = nullptr;
    }

    profileLength = plen;
}

ImageIO::ImageIO() :
    pl(nullptr),
    embProfile(nullptr),
    profileData(nullptr),
    profileLength(0),
    loadedProfileData(nullptr),
    loadedProfileDataJpg(false),
    loadedProfileLength(0),
    exifChange(new procparams::ExifPairs),
    iptc(nullptr),
    exifRoot(nullptr),
    sampleFormat(IIOSF_UNKNOWN),
    sampleArrangement(IIOSA_UNKNOWN)
{
}

ImageIO::~ImageIO ()
{

    if (embProfile) {
        cmsCloseProfile(embProfile);
    }

    deleteLoadedProfileData();
    delete exifRoot;
    delete [] profileData;
}

void png_read_data(png_struct_def  *png_ptr, unsigned char *data, size_t length);
void png_write_data(png_struct_def *png_ptr, unsigned char *data, size_t length);
void png_flush(png_struct_def *png_ptr);

int ImageIO::getPNGSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement)
{
    FILE *file = g_fopen (fname.c_str (), "rb");

    if (!file) {
        return IMIO_CANNOTREADFILE;
    }

    //reading PNG header
    unsigned char header[8];

    if (fread (header, 1, 8, file) != 8 || png_sig_cmp (header, 0, 8)) {
        fclose(file);
        return IMIO_HEADERERROR;
    }

    //initializing main structures
    png_structp png = png_create_read_struct (PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

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
    png_set_sig_bytes (png, 8);

    png_read_info(png, info);

    //retrieving image information
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;
    png_get_IHDR(png, info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    png_destroy_read_struct (&png, &info, &end_info);
    fclose (file);

    if (interlace_type != PNG_INTERLACE_NONE) {
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (bit_depth == 8) {
        sArrangement = IIOSA_CHUNKY;
        sFormat = IIOSF_UNSIGNED_CHAR;
        return IMIO_SUCCESS;
    } else if (bit_depth == 16) {
        sArrangement = IIOSA_CHUNKY;
        sFormat = IIOSF_UNSIGNED_SHORT;
        return IMIO_SUCCESS;
    } else {
        sArrangement = IIOSA_UNKNOWN;
        sFormat = IIOSF_UNKNOWN;
        return IMIO_VARIANTNOTSUPPORTED;
    }
}

int ImageIO::loadPNG  (const Glib::ustring &fname)
{

    FILE *file = g_fopen (fname.c_str (), "rb");

    if (!file) {
        return IMIO_CANNOTREADFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_LOADPNG");
        pl->setProgress (0.0);
    }

    //reading PNG header
    unsigned char header[8];

    if (fread (header, 1, 8, file) != 8 || png_sig_cmp (header, 0, 8)) {
        fclose(file);
        return IMIO_HEADERERROR;
    }

    //initializing main structures
    png_structp png = png_create_read_struct (PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

    if (!png) {
        fclose (file);
        return IMIO_HEADERERROR;
    }

    // silence the warning about "invalid" sRGB profiles -- see #4260
#if defined(PNG_SKIP_sRGB_CHECK_PROFILE) && defined(PNG_SET_OPTION_SUPPORTED)
    png_set_option(png, PNG_SKIP_sRGB_CHECK_PROFILE, PNG_OPTION_ON);
#endif

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
    png_set_sig_bytes (png, 8);

    png_read_info(png, info);

    embProfile = nullptr;

    //retrieving image information
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;
    png_get_IHDR(png, info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    if (color_type == PNG_COLOR_TYPE_PALETTE || interlace_type != PNG_INTERLACE_NONE )  {
        // we don't support interlaced png or png with palette
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        printf("%s uses an unsupported feature: <palette-indexed colors|interlacing>. Skipping.\n", fname.data());
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(png);
    }

    if (png_get_valid(png, info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(png);
    }

    if (color_type & PNG_COLOR_MASK_ALPHA) {
        png_set_strip_alpha(png);
    }

    // reading the embedded ICC profile if any
    if (png_get_valid(png, info, PNG_INFO_iCCP)) {
        png_charp name;
        int compression_type;
#if PNG_LIBPNG_VER < 10500
        png_charp profdata;
#else
        png_bytep profdata;
#endif
        png_uint_32 proflen;
        png_get_iCCP(png, info, &name, &compression_type, &profdata, &proflen);
        embProfile = cmsOpenProfileFromMem(profdata, proflen);
        loadedProfileData = new char[proflen];
        memcpy(loadedProfileData, profdata, proflen);
    }

    //setting gamma
    double gamma;

    if (png_get_gAMA(png, info, &gamma)) {
        png_set_gamma(png, 1.0 / gamma, gamma);    // use gamma from metadata
    } else {
        png_set_gamma(png, 2.2, 1.0 / 2.2);    // no gamma in metadata, suppose gamma 2.2
    }

//  if (bps==8 && bit_depth==16) png_set_strip_16(png);

    //updating png info struct
    png_read_update_info(png, info);
    png_get_IHDR(png, info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    allocate (width, height);

    int rowlen = width * 3 * bit_depth / 8;
    unsigned char *row = new unsigned char [rowlen];

    // set a new jump point to avoid memory leak
    if (setjmp (png_jmpbuf(png))) {
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        delete [] row;
        return IMIO_READERROR;
    }

    for (unsigned int i = 0; i < height; i++) {

        png_read_row (png, (png_byte*)row, nullptr);

        if (bit_depth == 16) { // convert scanline to host byte order
            unsigned short* srow = (unsigned short*)row;

            for (unsigned int j = 0; j < width * 3; j++) {
                srow[j] = ntohs (srow[j]);
            }
        }

        setScanline (i, row, bit_depth);

        if (pl && !(i % 100)) {
            pl->setProgress ((double)(i + 1) / height);
        }
    }

    png_read_end (png, nullptr);
    png_destroy_read_struct (&png, &info, &end_info);

    delete [] row;
    fclose(file);

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

typedef struct  {
    struct jpeg_error_mgr pub;  /* "public" fields */
    jmp_buf setjmp_buffer;  /* for return to caller */
} my_error_mgr;

void my_error_exit (j_common_ptr cinfo)
{
    /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
    my_error_mgr *myerr = (my_error_mgr*) cinfo->err;
    /* Always display the message. */
    /* We could postpone this until after returning, if we chose. */
    (*cinfo->err->output_message) (cinfo);

    /* Return control to the setjmp point */
#if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)
    __builtin_longjmp(myerr->setjmp_buffer, 1);
#else
    longjmp(myerr->setjmp_buffer, 1);
#endif
}


int ImageIO::loadJPEGFromMemory (const char* buffer, int bufsize)
{
    jpeg_decompress_struct cinfo;
    jpeg_create_decompress(&cinfo);
    jpeg_memory_src (&cinfo, (const JOCTET*)buffer, bufsize);

    /* We use our private extension JPEG error handler.
       Note that this struct must live as long as the main JPEG parameter
       struct, to avoid dangling-pointer problems.
    */
    my_error_mgr jerr;
    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
#if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)

    if (__builtin_setjmp(jerr.setjmp_buffer)) {
#else

    if (setjmp(jerr.setjmp_buffer)) {
#endif
        /* If we get here, the JPEG code has signaled an error.
           We need to clean up the JPEG object and return.
        */
        jpeg_destroy_decompress(&cinfo);
        return IMIO_READERROR;
    }


    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_LOADJPEG");
        pl->setProgress (0.0);

    }

    setup_read_icc_profile (&cinfo);

    jpeg_read_header(&cinfo, TRUE);

    deleteLoadedProfileData();
    loadedProfileDataJpg = true;
    bool hasprofile = read_icc_profile (&cinfo, (JOCTET**)&loadedProfileData, (unsigned int*)&loadedProfileLength);

    if (hasprofile) {
        embProfile = cmsOpenProfileFromMem (loadedProfileData, loadedProfileLength);
    } else {
        embProfile = nullptr;
    }

    jpeg_start_decompress(&cinfo);

    unsigned int width = cinfo.output_width;
    unsigned int height = cinfo.output_height;

    allocate (width, height);

    unsigned char *row = new unsigned char[width * 3];

    while (cinfo.output_scanline < height) {
        if (jpeg_read_scanlines(&cinfo, &row, 1) < 1) {
            jpeg_finish_decompress(&cinfo);
            jpeg_destroy_decompress(&cinfo);
            delete [] row;
            return IMIO_READERROR;
        }

        setScanline (cinfo.output_scanline - 1, row, 8, cinfo.num_components);

        if (pl && !(cinfo.output_scanline % 100)) {
            pl->setProgress ((double)(cinfo.output_scanline) / cinfo.output_height);
        }
    }

    delete [] row;

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

int ImageIO::loadJPEG (const Glib::ustring &fname)
{
    FILE *file = g_fopen(fname.c_str (), "rb");

    if (!file) {
        return IMIO_CANNOTREADFILE;
    }

    jpeg_decompress_struct cinfo;
    jpeg_error_mgr jerr;
    cinfo.err = my_jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);

    my_jpeg_stdio_src (&cinfo, file);

#if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)
    if ( __builtin_setjmp((reinterpret_cast<rt_jpeg_error_mgr*>(cinfo.src))->error_jmp_buf) == 0 ) {
#else
    if ( setjmp((reinterpret_cast<rt_jpeg_error_mgr*>(cinfo.src))->error_jmp_buf) == 0 ) {
#endif
        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_LOADJPEG");
            pl->setProgress (0.0);
        }

        setup_read_icc_profile (&cinfo);

        //jpeg_stdio_src(&cinfo,file);
        jpeg_read_header(&cinfo, TRUE);

        //if JPEG is CMYK, then abort reading
        if (cinfo.jpeg_color_space == JCS_CMYK || cinfo.jpeg_color_space == JCS_YCCK) {
            jpeg_destroy_decompress(&cinfo);
            return IMIO_READERROR;
        }

        cinfo.out_color_space = JCS_RGB;

        deleteLoadedProfileData();
        loadedProfileDataJpg = true;
        bool hasprofile = read_icc_profile (&cinfo, (JOCTET**)&loadedProfileData, (unsigned int*)&loadedProfileLength);

        if (hasprofile) {
            embProfile = cmsOpenProfileFromMem (loadedProfileData, loadedProfileLength);
        } else {
            embProfile = nullptr;
        }

        jpeg_start_decompress(&cinfo);

        unsigned int width = cinfo.output_width;
        unsigned int height = cinfo.output_height;

        allocate (width, height);

        unsigned char *row = new unsigned char[width * 3];

        while (cinfo.output_scanline < height) {
            if (jpeg_read_scanlines(&cinfo, &row, 1) < 1) {
                jpeg_finish_decompress(&cinfo);
                jpeg_destroy_decompress(&cinfo);
                delete [] row;
                return IMIO_READERROR;
            }

            setScanline (cinfo.output_scanline - 1, row, 8);

            if (pl && !(cinfo.output_scanline % 100)) {
                pl->setProgress ((double)(cinfo.output_scanline) / cinfo.output_height);
            }
        }

        delete [] row;

        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);
        fclose(file);

        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_READY");
            pl->setProgress (1.0);
        }

        return IMIO_SUCCESS;
    } else {
        jpeg_destroy_decompress(&cinfo);
        return IMIO_READERROR;
    }
}

int ImageIO::getTIFFSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement)
{
#ifdef WIN32
    wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
    TIFF* in = TIFFOpenW (wfilename, "r");
    g_free (wfilename);
#else
    TIFF* in = TIFFOpen(fname.c_str(), "r");
#endif

    if (in == nullptr) {
        return IMIO_CANNOTREADFILE;
    }

    uint16 bitspersample = 0, samplesperpixel = 0, sampleformat = 0;
    int hasTag = TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    hasTag &= TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);

    if (!hasTag) {
        // These are needed
        TIFFClose(in);
        sFormat = IIOSF_UNKNOWN;
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (!TIFFGetField(in, TIFFTAG_SAMPLEFORMAT, &sampleformat))
        /*
         * WARNING: This is a dirty hack!
         * We assume that files which doesn't contain the TIFFTAG_SAMPLEFORMAT tag
         * (which is the case with uncompressed TIFFs produced by RT!) are RGB files,
         * but that may be not true.   --- Hombre
         */
    {
        sampleformat = SAMPLEFORMAT_UINT;
    }

    uint16 config;
    TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);

    if (config == PLANARCONFIG_CONTIG) {
        sArrangement = IIOSA_CHUNKY;
    } else {
        sFormat = IIOSF_UNKNOWN;
        sArrangement = IIOSA_UNKNOWN;
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    uint16 photometric;

    if (!TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photometric)) {
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    uint16 compression;

    if (photometric == PHOTOMETRIC_LOGLUV)
        if (!TIFFGetField(in, TIFFTAG_COMPRESSION, &compression)) {
            compression = COMPRESSION_NONE;
        }

    TIFFClose(in);

    if (photometric == PHOTOMETRIC_RGB || photometric == PHOTOMETRIC_MINISBLACK) {
        if ((samplesperpixel == 1 || samplesperpixel == 3 || samplesperpixel == 4) && sampleformat == SAMPLEFORMAT_UINT) {
            if (bitspersample == 8) {
                sFormat = IIOSF_UNSIGNED_CHAR;
                return IMIO_SUCCESS;
            }

            if (bitspersample == 16) {
                sFormat = IIOSF_UNSIGNED_SHORT;
                return IMIO_SUCCESS;
            }
        } else if ((samplesperpixel == 3 || samplesperpixel == 4) && sampleformat == SAMPLEFORMAT_IEEEFP) {
            if (bitspersample==16) {
                sFormat = IIOSF_FLOAT16;
                return IMIO_SUCCESS;
            }
            if (bitspersample == 24) {
                sFormat = IIOSF_FLOAT24;
                return IMIO_SUCCESS;
            }
            if (bitspersample == 32) {
                sFormat = IIOSF_FLOAT32;
                return IMIO_SUCCESS;
            }
        }
    } else if ((samplesperpixel == 3 || samplesperpixel == 4) && photometric == PHOTOMETRIC_LOGLUV) {
        if (compression == COMPRESSION_SGILOG24) {
            sFormat = IIOSF_LOGLUV24;
            return IMIO_SUCCESS;
        } else if (compression == COMPRESSION_SGILOG) {
            sFormat = IIOSF_LOGLUV32;
            return IMIO_SUCCESS;
        }
    }

    return IMIO_VARIANTNOTSUPPORTED;
}

int ImageIO::loadTIFF (const Glib::ustring &fname)
{

    static MyMutex thumbMutex;
    MyMutex::MyLock lock(thumbMutex);

    if(!options.serializeTiffRead) {
        lock.release();
    }

#ifdef WIN32
    wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
    TIFF* in = TIFFOpenW (wfilename, "r");
    g_free (wfilename);
#else
    TIFF* in = TIFFOpen(fname.c_str(), "r");
#endif

    if (in == nullptr) {
        return IMIO_CANNOTREADFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_LOADTIFF");
        pl->setProgress (0.0);
    }

    int width, height;
    TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(in, TIFFTAG_IMAGELENGTH, &height);

    uint16 bitspersample, samplesperpixel;
    int hasTag = TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    hasTag &= TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);

    if (!hasTag) {
        // These are needed
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    uint16 config;
    TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);

    if (config != PLANARCONFIG_CONTIG) {
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (sampleFormat & (IIOSF_LOGLUV24 | IIOSF_LOGLUV32)) {
        TIFFSetField(in, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
    }

    /*
     * We could use the min/max values set in TIFFTAG_SMINSAMPLEVALUE and
     * TIFFTAG_SMAXSAMPLEVALUE, but for now, we normalize the image to the
     * effective minimum and maximum values
     */
    if (settings->verbose) {
        printf("Information of \"%s\":\n", fname.c_str());
        uint16 tiffDefaultScale, tiffBaselineExposure, tiffLinearResponseLimit;
        if (TIFFGetField(in, TIFFTAG_DEFAULTSCALE, &tiffDefaultScale)) {
            printf("   DefaultScale: %d\n", tiffDefaultScale);
        }
        else
            printf("   No DefaultScale value!\n");
        if (TIFFGetField(in, TIFFTAG_BASELINEEXPOSURE, &tiffBaselineExposure)) {
            printf("   BaselineExposure: %d\n", tiffBaselineExposure);
        }
        else
            printf("   No BaselineExposure value!\n");
        if (TIFFGetField(in, TIFFTAG_LINEARRESPONSELIMIT, &tiffLinearResponseLimit)) {
            printf("   LinearResponseLimit: %d\n", tiffLinearResponseLimit);
        }
        else
            printf("   No LinearResponseLimit value!\n");

        uint16 tiffMinValue, tiffMaxValue;
        if (TIFFGetField(in, TIFFTAG_SMINSAMPLEVALUE, &tiffMinValue)) {
            printf("   MinValue: %d\n", tiffMinValue);
        }
        else
            printf("   No minimum value!\n");
        if (TIFFGetField(in, TIFFTAG_SMAXSAMPLEVALUE, &tiffMaxValue)) {
            printf("   MaxValue: %d\n\n", tiffMaxValue);
        }
        else
            printf("   No maximum value!\n\n");
        printf("   Those values are not taken into account, the image data are normalized to a [0;1] range\n\n");
    }

    char* profdata;
    deleteLoadedProfileData();
    loadedProfileDataJpg = false;

    if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &loadedProfileLength, &profdata)) {
        embProfile = cmsOpenProfileFromMem (profdata, loadedProfileLength);
        loadedProfileData = new char [loadedProfileLength];
        memcpy (loadedProfileData, profdata, loadedProfileLength);
    } else {
        embProfile = nullptr;
    }

    allocate (width, height);

    unsigned char* linebuffer = new unsigned char[TIFFScanlineSize(in) * (samplesperpixel == 1 ? 3 : 1)];

    for (int row = 0; row < height; row++) {
        if (TIFFReadScanline(in, linebuffer, row, 0) < 0) {
            TIFFClose(in);
            delete [] linebuffer;
            return IMIO_READERROR;
        }

        if (samplesperpixel > 3) {
            for (int i = 0; i < width; i++) {
                memcpy (linebuffer + i * 3 * bitspersample / 8, linebuffer + i * samplesperpixel * bitspersample / 8, 3 * bitspersample / 8);
            }
        }
        else if (samplesperpixel == 1) {
            const size_t bytes = bitspersample / 8;
            for (int i = width - 1; i >= 0; --i) {
                const unsigned char* const src = linebuffer + i * bytes;
                unsigned char* const dest = linebuffer + i * 3 * bytes;
                memcpy(dest + 2 * bytes, src, bytes);
                memcpy(dest + 1 * bytes, src, bytes);
                memcpy(dest + 0 * bytes, src, bytes);
            }
        }

        setScanline (row, linebuffer, bitspersample);

        if (pl && !(row % 100)) {
            pl->setProgress ((double)(row + 1) / height);
        }
    }

    TIFFClose(in);
    delete [] linebuffer;

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

int ImageIO::loadPPMFromMemory(const char* buffer, int width, int height, bool swap, int bps)
{
    allocate (width, height);

    int line_length(width * 3 * (bps / 8));

    if ( swap && bps > 8 ) {
        char swapped[line_length];

        for ( int row = 0; row < height; ++row ) {
            ::rtengine::swab(((const char*)buffer) + (row * line_length), swapped, line_length);
            setScanline(row, (unsigned char*)&swapped[0], bps);
        }
    } else {
        for ( int row = 0; row < height; ++row ) {
            setScanline(row, ((const unsigned char*)buffer) + (row * line_length), bps);
        }
    }

    return IMIO_SUCCESS;
}


namespace {

// Taken from Darktable -- src/imageio/format/png.c
//
/* Write EXIF data to PNG file.
 * Code copied from DigiKam's libs/dimg/loaders/pngloader.cpp.
 * The EXIF embedding is defined by ImageMagicK.
 * It is documented in the ExifTool page:
 * http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/PNG.html
 *
 * ..and in turn copied from ufraw. thanks to udi and colleagues
 * for making useful code much more readable and discoverable ;)
 */

void PNGwriteRawProfile(png_struct *ping, png_info *ping_info, const char *profile_type, guint8 *profile_data, png_uint_32 length)
{
    png_textp text;
    long i;
    guint8 *sp;
    png_charp dp;
    png_uint_32 allocated_length, description_length;

    const guint8 hex[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
    text = static_cast<png_textp>(png_malloc(ping, sizeof(png_text)));
    description_length = strlen(profile_type);
    allocated_length = length * 2 + (length >> 5) + 20 + description_length;

    text[0].text = static_cast<png_charp>(png_malloc(ping, allocated_length));
    text[0].key = static_cast<png_charp>(png_malloc(ping, 80));
    text[0].key[0] = '\0';

    g_strlcat(text[0].key, "Raw profile type ", 80);
    g_strlcat(text[0].key, profile_type, 80);

    sp = profile_data;
    dp = text[0].text;
    *dp++ = '\n';

    g_strlcpy(dp, profile_type, allocated_length);

    dp += description_length;
    *dp++ = '\n';
    *dp = '\0';

    g_snprintf(dp, allocated_length - strlen(text[0].text), "%8lu ", static_cast<unsigned long int>(length));

    dp += 8;

    for(i = 0; i < long(length); i++)
    {
        if(i % 36 == 0) *dp++ = '\n';

        *(dp++) = hex[((*sp >> 4) & 0x0f)];
        *(dp++) = hex[((*sp++) & 0x0f)];
    }

    *dp++ = '\n';
    *dp = '\0';
    text[0].text_length = (dp - text[0].text);
    text[0].compression = -1;

    if(text[0].text_length <= allocated_length) png_set_text(ping, ping_info, text, 1);

    png_free(ping, text[0].text);
    png_free(ping, text[0].key);
    png_free(ping, text);
}

} // namespace

int ImageIO::savePNG  (const Glib::ustring &fname, int bps) const
{
    if (getWidth() < 1 || getHeight() < 1) {
        return IMIO_HEADERERROR;
    }

    FILE* const file = g_fopen_withBinaryAndLock (fname);

    if (!file) {
        return IMIO_CANNOTWRITEFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_SAVEPNG");
        pl->setProgress (0.0);
    }

    png_structp png = png_create_write_struct (PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

    if (!png) {
        fclose (file);
        return IMIO_HEADERERROR;
    }

    // silence the warning about "invalid" sRGB profiles -- see #4260
#if defined(PNG_SKIP_sRGB_CHECK_PROFILE) && defined(PNG_SET_OPTION_SUPPORTED)
    png_set_option(png, PNG_SKIP_sRGB_CHECK_PROFILE, PNG_OPTION_ON);
#endif
    
    png_infop info = png_create_info_struct(png);

    if (!info) {
        png_destroy_write_struct (&png, nullptr);
        fclose (file);
        return IMIO_HEADERERROR;
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct (&png, &info);
        fclose(file);
        return IMIO_CANNOTWRITEFILE;
    }

    png_set_write_fn (png, file, png_write_data, png_flush);

    png_set_filter(png, 0, PNG_FILTER_PAETH);
    png_set_compression_level(png, 6);
    png_set_compression_strategy(png, 3);

    int width = getWidth ();
    int height = getHeight ();

    if (bps < 0) {
        bps = getBPS ();
    }
    if (bps > 16) {
        bps = 16;
    }

    png_set_IHDR(png, info, width, height, bps, PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_BASE);

    if (profileData) {
#if PNG_LIBPNG_VER < 10500
        png_charp profdata = reinterpret_cast<png_charp>(profileData);
#else
        png_bytep profdata = reinterpret_cast<png_bytep>(profileData);
#endif
        png_set_iCCP(png, info, const_cast<png_charp>("icc"), 0, profdata, profileLength);
    }

    {
        // buffer for the exif and iptc
        unsigned int bufferSize;
        unsigned char* buffer = nullptr; // buffer will be allocated in createTIFFHeader
        unsigned char* iptcdata = nullptr;
        unsigned int iptclen = 0;

        if (iptc && iptc_data_save (iptc, &iptcdata, &iptclen) && iptcdata) {
            iptc_data_free_buf (iptc, iptcdata);
            iptcdata = nullptr;
        }

        int size = rtexif::ExifManager::createPNGMarker(exifRoot, *exifChange, width, height, bps, (char*)iptcdata, iptclen, buffer, bufferSize);

        if (iptcdata) {
            iptc_data_free_buf (iptc, iptcdata);
        }
        if (buffer && size) {
            PNGwriteRawProfile(png, info, "exif", buffer, size);
            delete[] buffer;
        }
    }


    int rowlen = width * 3 * bps / 8;
    unsigned char *row = new unsigned char [rowlen];

    png_write_info(png, info);

    for (int i = 0; i < height; i++) {
        getScanline (i, row, bps);

        if (bps == 16) {
            // convert to network byte order
#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
            for (int j = 0; j < width * 6; j += 2) {
                unsigned char tmp = row[j];
                row[j] = row[j + 1];
                row[j + 1] = tmp;
            }

#endif
        }

        png_write_row (png, (png_byte*)row);

        if (pl && !(i % 100)) {
            pl->setProgress ((double)(i + 1) / height);
        }
    }

    png_write_end(png, info);
    png_destroy_write_struct(&png, &info);

    delete [] row;
    fclose (file);

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}



// Quality 0..100, subsampling: 1=low quality, 2=medium, 3=high
int ImageIO::saveJPEG (const Glib::ustring &fname, int quality, int subSamp) const
{
    if (getWidth() < 1 || getHeight() < 1) {
        return IMIO_HEADERERROR;
    }

    FILE* const file = g_fopen_withBinaryAndLock (fname);

    if (!file) {
        return IMIO_CANNOTWRITEFILE;
    }

    jpeg_compress_struct cinfo;
    /* We use our private extension JPEG error handler.
       Note that this struct must live as long as the main JPEG parameter
       struct, to avoid dangling-pointer problems.
    */
    my_error_mgr jerr;
    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
#if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)

    if (__builtin_setjmp(jerr.setjmp_buffer)) {
#else

    if (setjmp(jerr.setjmp_buffer)) {
#endif
        /* If we get here, the JPEG code has signaled an error.
           We need to clean up the JPEG object, close the file, remove the already saved part of the file and return.
        */
        jpeg_destroy_compress(&cinfo);
        fclose(file);
        g_remove (fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }

    jpeg_create_compress (&cinfo);



    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_SAVEJPEG");
        pl->setProgress (0.0);
    }

    jpeg_stdio_dest (&cinfo, file);

    int width = getWidth ();
    int height = getHeight ();

    cinfo.image_width  = width;
    cinfo.image_height = height;
    cinfo.in_color_space = JCS_RGB;
    cinfo.input_components = 3;
    jpeg_set_defaults (&cinfo);
    cinfo.write_JFIF_header = FALSE;

    // compute optimal Huffman coding tables for the image. Bit slower to generate, but size of result image is a bit less (default was FALSE)
    cinfo.optimize_coding = TRUE;

    // Since math coprocessors are common these days, FLOAT should be a bit more accurate AND fast (default is ISLOW)
    // (machine dependency is not really an issue, since we all run on x86 and having exactly the same file is not a requirement)
    cinfo.dct_method = JDCT_FLOAT;

    if (quality >= 0 && quality <= 100) {
        jpeg_set_quality (&cinfo, quality, true);
    }

    cinfo.comp_info[1].h_samp_factor = cinfo.comp_info[1].v_samp_factor = 1;
    cinfo.comp_info[2].h_samp_factor = cinfo.comp_info[2].v_samp_factor = 1;

    if (subSamp == 1) {
        // Best compression, default of the JPEG library:  2x2, 1x1, 1x1 (4:2:0)
        cinfo.comp_info[0].h_samp_factor = cinfo.comp_info[0].v_samp_factor = 2;
    } else if (subSamp == 2) {
        // Widely used normal ratio 2x1, 1x1, 1x1 (4:2:2)
        cinfo.comp_info[0].h_samp_factor = 2;
        cinfo.comp_info[0].v_samp_factor = 1;
    } else if (subSamp == 3) {
        // Best quality 1x1 1x1 1x1 (4:4:4)
        cinfo.comp_info[0].h_samp_factor = cinfo.comp_info[0].v_samp_factor = 1;
    }

    jpeg_start_compress(&cinfo, TRUE);

    // buffer for exif and iptc markers
    unsigned char* buffer = new unsigned char[165535]; //FIXME: no buffer size check so it can be overflowed in createJPEGMarker() for large tags, and then software will crash
    unsigned int size;

    // assemble and write exif marker
    if (exifRoot) {
        int size = rtexif::ExifManager::createJPEGMarker (exifRoot, *exifChange, cinfo.image_width, cinfo.image_height, buffer);

        if (size > 0 && size < 65530) {
            jpeg_write_marker(&cinfo, JPEG_APP0 + 1, buffer, size);
        }
    }

    // assemble and write iptc marker
    if (iptc) {
        unsigned char* iptcdata;
        bool error = false;

        if (iptc_data_save (iptc, &iptcdata, &size)) {
            if (iptcdata) {
                iptc_data_free_buf (iptc, iptcdata);
            }

            error = true;
        }

        int bytes = 0;

        if (!error && (bytes = iptc_jpeg_ps3_save_iptc (nullptr, 0, iptcdata, size, buffer, 65532)) < 0) {
            error = true;
        }

        if (iptcdata) {
            iptc_data_free_buf (iptc, iptcdata);
        }

        if (!error) {
            jpeg_write_marker(&cinfo, JPEG_APP0 + 13, buffer, bytes);
        }
    }

    delete [] buffer;

    // write icc profile to the output
    if (profileData) {
        write_icc_profile (&cinfo, (JOCTET*)profileData, profileLength);
    }

    // write image data
    int rowlen = width * 3;
    unsigned char *row = new unsigned char [rowlen];

    /* To avoid memory leaks we establish a new setjmp return context for my_error_exit to use. */
#if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)

    if (__builtin_setjmp(jerr.setjmp_buffer)) {
#else

    if (setjmp(jerr.setjmp_buffer)) {
#endif
        /* If we get here, the JPEG code has signaled an error.
           We need to clean up the JPEG object, close the file, remove the already saved part of the file and return.
        */
        delete [] row;
        jpeg_destroy_compress(&cinfo);
        fclose(file);
        g_remove (fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }

    while (cinfo.next_scanline < cinfo.image_height) {

        getScanline (cinfo.next_scanline, row, 8);

        if (jpeg_write_scanlines (&cinfo, &row, 1) < 1) {
            jpeg_destroy_compress (&cinfo);
            delete [] row;
            fclose (file);
            g_remove (fname.c_str());
            return IMIO_CANNOTWRITEFILE;
        }

        if (pl && !(cinfo.next_scanline % 100)) {
            pl->setProgress ((double)(cinfo.next_scanline) / cinfo.image_height);
        }
    }

    jpeg_finish_compress (&cinfo);
    jpeg_destroy_compress (&cinfo);

    delete [] row;

    fclose (file);

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

int ImageIO::saveTIFF (const Glib::ustring &fname, int bps, bool isFloat, bool uncompressed) const
{
    if (getWidth() < 1 || getHeight() < 1) {
        return IMIO_HEADERERROR;
    }

    bool writeOk = true;
    int width = getWidth ();
    int height = getHeight ();

    if (bps < 0) {
        bps = getBPS ();
    }

    int lineWidth = width * 3 * bps / 8;
    unsigned char* linebuffer = new unsigned char[lineWidth];

    // little hack to get libTiff to use proper byte order (see TIFFClienOpen()):
    const char *mode = !exifRoot ? "w" : (exifRoot->getOrder() == rtexif::INTEL ? "wl" : "wb");
#ifdef WIN32
    FILE *file = g_fopen_withBinaryAndLock (fname);
    int fileno = _fileno(file);
    int osfileno = _get_osfhandle(fileno);
    TIFF* out = TIFFFdOpen (osfileno, fname.c_str(), mode);
#else
    TIFF* out = TIFFOpen(fname.c_str(), mode);
    int fileno = TIFFFileno (out);
#endif

    if (!out) {
        delete [] linebuffer;
        return IMIO_CANNOTWRITEFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_SAVETIFF");
        pl->setProgress (0.0);
    }

    bool applyExifPatch = false;

    if (exifRoot) {
        rtexif::TagDirectory* cl = (const_cast<rtexif::TagDirectory*> (exifRoot))->clone (nullptr);

        // ------------------ remove some unknown top level tags which produce warnings when opening a tiff (might be useless) -----------------

        rtexif::Tag *removeTag = cl->getTag (0x9003);

        if (removeTag) {
            removeTag->setKeep (false);
        }

        removeTag = cl->getTag (0x9211);

        if (removeTag) {
            removeTag->setKeep (false);
        }

        // ------------------ Apply list of change -----------------

        for (auto currExifChange : *exifChange) {
            cl->applyChange (currExifChange.first, currExifChange.second);
        }

        rtexif::Tag *tag = cl->getTag (TIFFTAG_EXIFIFD);

        if (tag && tag->isDirectory()) {
            rtexif::TagDirectory *exif = tag->getDirectory();

            if (exif)   {
                int exif_size = exif->calculateSize();
                unsigned char *buffer = new unsigned char[exif_size + 8];
                // TIFFOpen writes out the header and sets file pointer at position 8

                exif->write (8, buffer);

                write (fileno, buffer + 8, exif_size);

                delete [] buffer;
                // let libtiff know that scanlines or any other following stuff should go
                // at a different offset:
                TIFFSetWriteOffset (out, exif_size + 8);
                TIFFSetField (out, TIFFTAG_EXIFIFD, 8);
                applyExifPatch = true;
            }
        }

        //TODO Even though we are saving EXIF IFD - MakerNote still comes out screwed.

        if ((tag = cl->getTag (TIFFTAG_MODEL)) != nullptr) {
            TIFFSetField (out, TIFFTAG_MODEL, tag->getValue());
        }

        if ((tag = cl->getTag (TIFFTAG_MAKE)) != nullptr) {
            TIFFSetField (out, TIFFTAG_MAKE, tag->getValue());
        }

        if ((tag = cl->getTag (TIFFTAG_DATETIME)) != nullptr) {
            TIFFSetField (out, TIFFTAG_DATETIME, tag->getValue());
        }

        if ((tag = cl->getTag (TIFFTAG_ARTIST)) != nullptr) {
            TIFFSetField (out, TIFFTAG_ARTIST, tag->getValue());
        }

        if ((tag = cl->getTag (TIFFTAG_COPYRIGHT)) != nullptr) {
            TIFFSetField (out, TIFFTAG_COPYRIGHT, tag->getValue());
        }

        delete cl;
    }

    unsigned char* iptcdata = nullptr;
    unsigned int iptclen = 0;

    if (iptc && iptc_data_save (iptc, &iptcdata, &iptclen)) {
        if (iptcdata) {
            iptc_data_free_buf (iptc, iptcdata);
            iptcdata = nullptr;
        }
    }

#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
        bool needsReverse = exifRoot && exifRoot->getOrder() == rtexif::MOTOROLA;
#else
        bool needsReverse = exifRoot && exifRoot->getOrder() == rtexif::INTEL;
#endif
    if (iptcdata) {
        rtexif::Tag iptcTag(nullptr, rtexif::lookupAttrib (rtexif::ifdAttribs, "IPTCData"));
        iptcTag.initLongArray((char*)iptcdata, iptclen);
        if (needsReverse) {
            unsigned char *ptr = iptcTag.getValue();
            for (int a = 0; a < iptcTag.getCount(); ++a) {
                unsigned char cc;
                cc = ptr[3];
                ptr[3] = ptr[0];
                ptr[0] = cc;
                cc = ptr[2];
                ptr[2] = ptr[1];
                ptr[1] = cc;
                ptr += 4;
            }
        }
        TIFFSetField (out, TIFFTAG_RICHTIFFIPTC, iptcTag.getCount(), (long*)iptcTag.getValue());
        iptc_data_free_buf (iptc, iptcdata);
    }

    TIFFSetField (out, TIFFTAG_SOFTWARE, "RawTherapee " RTVERSION);
    TIFFSetField (out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField (out, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField (out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField (out, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField (out, TIFFTAG_ROWSPERSTRIP, height);
    TIFFSetField (out, TIFFTAG_BITSPERSAMPLE, bps);
    TIFFSetField (out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField (out, TIFFTAG_COMPRESSION, uncompressed ? COMPRESSION_NONE : COMPRESSION_ADOBE_DEFLATE);
    TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, (bps == 16 || bps == 32) && isFloat ? SAMPLEFORMAT_IEEEFP : SAMPLEFORMAT_UINT);

    [out]()
    {
        const std::vector<rtexif::Tag*> default_tags = rtexif::ExifManager::getDefaultTIFFTags(nullptr);

        TIFFSetField (out, TIFFTAG_XRESOLUTION, default_tags[2]->toDouble());
        TIFFSetField (out, TIFFTAG_YRESOLUTION, default_tags[3]->toDouble());
        TIFFSetField (out, TIFFTAG_RESOLUTIONUNIT, default_tags[4]->toInt());

        for (auto default_tag : default_tags) {
            delete default_tag;
        }
    }();

    if (!uncompressed) {
        TIFFSetField (out, TIFFTAG_PREDICTOR, (bps == 16 || bps == 32) && isFloat ? PREDICTOR_FLOATINGPOINT : PREDICTOR_HORIZONTAL);
    }
    if (profileData) {
        TIFFSetField (out, TIFFTAG_ICCPROFILE, profileLength, profileData);
    }

    for (int row = 0; row < height; row++) {
        getScanline (row, linebuffer, bps, isFloat);

        if (bps == 16) {
            if(needsReverse && !uncompressed && isFloat) {
                for(int i = 0; i < lineWidth; i += 2) {
                    char temp = linebuffer[i];
                    linebuffer[i] = linebuffer[i + 1];
                    linebuffer[i + 1] = temp;
                }
            }
        } else if (bps == 32) {
            if(needsReverse && !uncompressed) {
                for(int i = 0; i < lineWidth; i += 4) {
                    char temp = linebuffer[i];
                    linebuffer[i] = linebuffer[i + 3];
                    linebuffer[i + 3] = temp;
                    temp = linebuffer[i + 1];
                    linebuffer[i + 1] = linebuffer[i + 2];
                    linebuffer[i + 2] = temp;
                }
            }
        }

        if (TIFFWriteScanline (out, linebuffer, row, 0) < 0) {
            TIFFClose (out);
            delete [] linebuffer;
            return IMIO_CANNOTWRITEFILE;
        }

        if (pl && !(row % 100)) {
            pl->setProgress ((double)(row + 1) / height);
        }
    }

    if (TIFFFlush(out) != 1) {
        writeOk = false;
    }

    /************************************************************************************************************
     *
     * Hombre: This is a dirty hack to update the Exif tag data type to 0x0004 so that Windows can understand it.
     *         libtiff will set this data type to 0x000d and doesn't provide any mechanism to update it before
     *         dumping to the file.
     *
     */
    if (applyExifPatch) {
        unsigned char b[10];
        uint16 tagCount = 0;
        lseek(fileno, 4, SEEK_SET);
        read(fileno, b, 4);
        uint32 ifd0Offset = rtexif::sget4(b, exifRoot->getOrder());
        lseek(fileno, ifd0Offset, SEEK_SET);
        read(fileno, b, 2);
        tagCount = rtexif::sget2(b, exifRoot->getOrder());
        for (size_t i = 0; i < tagCount ; ++i) {
            uint16 tagID = 0;
            read(fileno, b, 2);
            tagID = rtexif::sget2(b, exifRoot->getOrder());
            if (tagID == 0x8769) {
                rtexif::sset2(4, b, exifRoot->getOrder());
                write(fileno, b, 2);
                break;
            } else {
                read(fileno, b, 10);
            }
        }
    }
    /************************************************************************************************************/


    TIFFClose (out);
#ifdef WIN32
    fclose (file);
#endif

    delete [] linebuffer;

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    if(writeOk) {
        return IMIO_SUCCESS;
    } else {
        g_remove (fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }
}

// PNG read and write routines:

void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
    png_size_t check;

    /* fread() returns 0 on error, so it is OK to store this in a png_size_t
     * instead of an int, which is what fread() actually returns.
     */
    check = (png_size_t)fread(data, (png_size_t)1, length, (FILE *)png_get_io_ptr(png_ptr));

    if (check != length) {
        png_error(png_ptr, "Read Error");
    }
}

void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
    png_uint_32 check;

    check = fwrite(data, 1, length, (FILE *)png_get_io_ptr(png_ptr));

    if (check != length) {
        png_error(png_ptr, "Write Error");
    }
}

void png_flush(png_structp png_ptr)
{
    FILE *io_ptr;
    io_ptr = (FILE *)(png_get_io_ptr(png_ptr));

    if (io_ptr != nullptr) {
        fflush(io_ptr);
    }
}

int ImageIO::load (const Glib::ustring &fname)
{

    if (hasPngExtension(fname)) {
        return loadPNG (fname);
    } else if (hasJpegExtension(fname)) {
        return loadJPEG (fname);
    } else if (hasTiffExtension(fname)) {
        return loadTIFF (fname);
    } else {
        return IMIO_FILETYPENOTSUPPORTED;
    }
}

int ImageIO::save (const Glib::ustring &fname) const
{
    if (hasPngExtension(fname)) {
        return savePNG (fname);
    } else if (hasJpegExtension(fname)) {
        return saveJPEG (fname);
    } else if (hasTiffExtension(fname)) {
        return saveTIFF (fname);
    } else {
        return IMIO_FILETYPENOTSUPPORTED;
    }
}

void ImageIO::setProgressListener (ProgressListener* l)
{
    pl = l;
}

void ImageIO::setSampleFormat(IIOSampleFormat sFormat)
{
    sampleFormat = sFormat;
}

IIOSampleFormat ImageIO::getSampleFormat() const
{
    return sampleFormat;
}

void ImageIO::setSampleArrangement(IIOSampleArrangement sArrangement)
{
    sampleArrangement = sArrangement;
}

IIOSampleArrangement ImageIO::getSampleArrangement() const
{
    return sampleArrangement;
}

cmsHPROFILE ImageIO::getEmbeddedProfile () const
{
    return embProfile;
}

void ImageIO::getEmbeddedProfileData (int& length, unsigned char*& pdata) const
{
    length = loadedProfileLength;
    pdata = (unsigned char*)loadedProfileData;
}

MyMutex& ImageIO::mutex ()
{
    return imutex;
}

void ImageIO::deleteLoadedProfileData( )
{
    if(loadedProfileData) {
        if(loadedProfileDataJpg) {
            free(loadedProfileData);
        } else {
            delete[] loadedProfileData;
        }
    }

    loadedProfileData = nullptr;
}
