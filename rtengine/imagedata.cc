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
#include <imagedata.h>
#include <iptcpairs.h>
#include <glib/gstdio.h>
#ifdef RAWZOR_SUPPORT
#include <rwz_sdk.h>
#endif

using namespace rtengine;

extern "C" IptcData *iptc_data_new_from_jpeg_file (FILE* infile);

ImageMetaData* ImageMetaData::fromFile (const Glib::ustring& fname, RawMetaDataLocation* rml) {

    return new ImageData (fname, rml);
}

ImageData::ImageData (Glib::ustring fname, RawMetaDataLocation* ri) {

    int dotpos = fname.find_last_of ('.');
    root = NULL;
    iptc = NULL;
#ifdef RAWZOR_SUPPORT
    // RAWZOR support begin
    if (dotpos<fname.size()-3 && !fname.casefold().compare (dotpos, 4, ".rwz") && ri && (ri->exifBase>=0 || ri->ciffBase>=0)) {
        FILE* f = g_fopen (fname.c_str (), "rb");
        if (f) {
        	fseek (f, 0, SEEK_END);
        	int rzwSize = ftell (f);
            char* rzwData = new char [rzwSize];
        	fseek (f, 0, SEEK_SET);
	        fread (rzwData, 1, rzwSize, f);
            fclose(f);
            int rawSize;
            if (!m_rwz_check (rzwData, rzwSize, &rawSize)) {
                char* rawData = new char [rawSize];
                if (!m_rwz_get_meta_only (rzwData, rzwSize, rawData, rawSize)) {
                    std::string tfname;
                    int fd = Glib::file_open_tmp (tfname, "");
                    FILE* tf = fdopen (fd, "w+b");
                    fwrite (rawData, 1, rawSize, tf);
                    if (ri->exifBase>=0)
                        root = rtexif::ExifManager::parse (tf, ri->exifBase);
                    else if (ri->ciffBase>=0)
                        root = rtexif::ExifManager::parseCIFF (tf, ri->ciffBase, ri->ciffLength);
                    fclose (tf);
                    ::g_remove (tfname.c_str());
                    extractInfo ();
                }
                delete [] rawData;
                
            }
            delete [] rzwData;
        }
    }    
    // RAWZOR support end
    else 
#endif
    if (ri && (ri->exifBase>=0 || ri->ciffBase>=0)) {
        FILE* f = g_fopen (fname.c_str(), "rb");
        if (f) {
            if (ri->exifBase>=0) {
                root = rtexif::ExifManager::parse (f, ri->exifBase);
                if (root) {
                    rtexif::Tag* t = root->getTag (0x83BB);
                    if (t) 
                        iptc = iptc_data_new_from_data ((unsigned char*)t->getValue (), (unsigned)t->getValueSize ());
                }
            }
            else if (ri->ciffBase>=0)
                root = rtexif::ExifManager::parseCIFF (f, ri->ciffBase, ri->ciffLength);
            fclose (f);
            extractInfo ();
        }
    }    
    else if (dotpos<fname.size()-3 && !fname.casefold().compare (dotpos, 4, ".jpg")) {
        FILE* f = g_fopen (fname.c_str (), "rb");
        if (f) {
            root = rtexif::ExifManager::parseJPEG (f);
            extractInfo ();
            fclose (f);
            FILE* ff = g_fopen (fname.c_str (), "rb");
            iptc = iptc_data_new_from_jpeg_file (ff);
            fclose (ff);
        }
    }    
    else if ((dotpos<fname.size()-3 && !fname.casefold().compare (dotpos, 4, ".tif")) || (dotpos<fname.size()-4 && !fname.casefold().compare (dotpos, 5, ".tiff"))) {
        FILE* f = g_fopen (fname.c_str (), "rb");
        if (f) {
            root = rtexif::ExifManager::parseTIFF (f);
            fclose (f);
            extractInfo ();
            if (root) {
                rtexif::Tag* t = root->getTag (0x83BB);
                if (t) 
                    iptc = iptc_data_new_from_data ((unsigned char*)t->getValue (), (unsigned)t->getValueSize ());
            }
        }
    }    
    else {
        root = new rtexif::TagDirectory ();
        shutter = 0;
        aperture = 0;
        iso_speed = 0;
        lens = "Unknown";
        make = "Unknown";
        model = "Unknown";
        focal_len = 0;
        memset (&time, 0, sizeof(time));
    }
}

void ImageData::extractInfo () {

  if (!root) 
    return;

  char buffer[256];

  make = "";
  model = "";
  shutter = 0;
  aperture = 0;
  focal_len = 0;
  iso_speed = 0;
  memset (&time, 0, sizeof(time));
  
  if (root->getTag ("Make"))
    make = root->getTag ("Make")->valueToString ();  
  if (root->getTag ("Model"))
    model = root->getTag ("Model")->valueToString ();  

  rtexif::TagDirectory* exif = NULL;
  if (root->getTag ("Exif"))
    exif = root->getTag ("Exif")->getDirectory ();
  
  if (exif) {
  
    // standard exif tags
    if (exif->getTag ("ShutterSpeedValue"))
        shutter = exif->getTag ("ShutterSpeedValue")->toDouble ();
    if (exif->getTag ("ExposureTime"))
        shutter = exif->getTag ("ExposureTime")->toDouble ();
    if (exif->getTag ("ApertureValue"))
        aperture = exif->getTag ("ApertureValue")->toDouble ();
    if (exif->getTag ("FNumber"))
        aperture = exif->getTag ("FNumber")->toDouble ();
    if (exif->getTag ("FocalLength"))
        focal_len = exif->getTag ("FocalLength")->toDouble ();
    if (exif->getTag ("ISOSpeedRatings"))
        iso_speed = exif->getTag ("ISOSpeedRatings")->toDouble ();
    if (exif->getTag ("DateTimeOriginal")) {
        if (sscanf ((const char*)exif->getTag("DateTimeOriginal")->getValue(), "%d:%d:%d %d:%d:%d", &time.tm_year, &time.tm_mon, &time.tm_mday, &time.tm_hour, &time.tm_min, &time.tm_sec) == 6) {
            time.tm_year -= 1900;
            time.tm_mon -= 1;
        }
    }
    // guess lens...
    lens = "Unknown";

    if (exif->getTag ("MakerNote")) {
        rtexif::TagDirectory* mnote = exif->getTag ("MakerNote")->getDirectory();
        if (mnote && !make.compare (0, 5, "NIKON")) {
            bool lensOk = false;
            if (mnote->getTag ("LensData")) {
                std::string ldata = mnote->getTag ("LensData")->valueToString ();
                int pos;
                if (ldata.size()>10 && (pos=ldata.find ("Lens = "))!=Glib::ustring::npos) {
                    lens = ldata.substr (pos + 7);
                    if (lens.compare (0, 7, "Unknown"))
                        lensOk = true;
                }
            }
            if (!lensOk && mnote->getTag ("Lens")) {
                std::string ldata = mnote->getTag ("Lens")->valueToString ();
                int i=0, j=0; 
                double n[4];
                for (int m=0; m<4; m++) {
                    while (i<ldata.size() && ldata[i]!='/') i++;
                    int nom = atoi(ldata.substr(j, i).c_str());
                    j = i+1; i++;
                    while (i<ldata.size() && ldata[i]!=',') i++;
                    int den = atoi(ldata.substr(j, i).c_str());
                    j = i+2; i+=2;
                    n[m] = (double) nom/den;
                }
                std::ostringstream str;
                if (n[0]==n[1])
                    str << "Unknown " << n[0] << "mm F/" << n[2];
                else if (n[2]==n[3])
                    str << "Unknown " << n[0] << "-" << n[1] << "mm F/" << n[2];
                else
                    str << "Unknown " << n[0] << "-" << n[1] << "mm F/" << n[2] << "-" << n[3];
                lens = str.str();
            }
        }
        else if (mnote && !make.compare (0, 5, "Canon")) {
            bool lensOk = false;
            if (mnote->getTag ("LensType")) {
                std::string ldata = mnote->getTag ("LensType")->valueToString ();
                if (ldata.size()>1) {
                    lens = ldata;
                    lensOk = true;
                }
            }
            if (!lensOk && mnote->getTag ("CanonCameraSettings")) {
                std::string ccs = mnote->getTag ("CanonCameraSettings")->valueToString ();
                int i = ccs.find ("LongFocal = ");
                double a = 0;
                if (i!=ccs.npos) {
                    i += 12;
                    int j = i;
                    while (j!=ccs.npos && ccs[j]!='\n' && ccs[j]!=' ') j++;
                    a = atof (ccs.substr (i, j-i).c_str());
                }
                i = ccs.find ("ShortFocal = ");
                double b = 0;
                if (i!=ccs.npos) {
                    i += 13;
                    int j = i;
                    while (j!=ccs.npos && ccs[j]!='\n' && ccs[j]!=' ') j++;
                    b = atof (ccs.substr (i, j-i).c_str());
                }
                if (a>0 && b>0) {
                    std::ostringstream str;
                    if (a==b)
                        str << "Unknown " << a << "mm";
                    else 
                        str << "Unknown " << b << "-" << a << "mm";
                    lens = str.str();
                }                
            }
        }
        else if (mnote && !make.compare (0, 6, "PENTAX")) {
            if (mnote->getTag ("LensType")) 
                lens = mnote->getTag ("LensType")->valueToString ();
        }
        else if (mnote && (!make.compare (0, 4, "SONY") || !make.compare (0, 6, "KONICA"))) {
            if (mnote->getTag ("LensID")) 
                lens = mnote->getTag ("LensID")->valueToString ();
        }
        else if (mnote && !make.compare (0, 7, "OLYMPUS")) {
            if (mnote->getTag ("Equipment"))  {
                rtexif::TagDirectory* eq = mnote->getTag ("Equipment")->getDirectory ();
                if (eq->getTag ("LensType"))
                    lens = eq->getTag ("LensType")->valueToString ();
            }
        }
    }
  }
}

ImageData::~ImageData () {

    delete root;
    if (iptc) 
        iptc_data_free (iptc);
}

const std::vector<procparams::IPTCPair> ImageData::getIPTCData () const {

    std::vector<procparams::IPTCPair> iptcc;
    if (!iptc)
        return iptcc;
   
    unsigned char buffer[2100];
    for (int i=0; i<16; i++) {
        IptcDataSet* ds = iptc_data_get_next_dataset (iptc, NULL, IPTC_RECORD_APP_2, strTags[i].tag);
        if (ds) {
            iptc_dataset_get_data (ds, buffer, 2100);
            procparams::IPTCPair ic;
            ic.field = strTags[i].field;
            try {
                ic.values.push_back (Glib::locale_to_utf8((char*)buffer));
            }
            catch (const Glib::ConvertError& e) {
                ic.values.push_back (Glib::convert_with_fallback((char*)buffer, "UTF8", "LATIN1","?"));
            }
            iptcc.push_back (ic);
            iptc_dataset_unref (ds);
        }
    }
    IptcDataSet* ds = NULL;
    procparams::IPTCPair ickw;
    ickw.field = "Keywords";
    while (ds=iptc_data_get_next_dataset (iptc, ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS)) {
        iptc_dataset_get_data (ds, buffer, 2100);
        try {
            ickw.values.push_back (Glib::locale_to_utf8((char*)buffer));
        }
        catch (const Glib::ConvertError& e) {
            ickw.values.push_back (Glib::convert_with_fallback((char*)buffer, "UTF8", "LATIN1","?"));
        }
    }
    iptcc.push_back (ickw);
    ds = NULL;
    procparams::IPTCPair icsc;
    icsc.field = "SupplementalCategories";
    while (ds=iptc_data_get_next_dataset (iptc, ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY)) {
        iptc_dataset_get_data (ds, buffer, 2100);
        try {
            icsc.values.push_back (Glib::locale_to_utf8((char*)buffer));
        }
        catch (const Glib::ConvertError& e) {
            icsc.values.push_back (Glib::convert_with_fallback((char*)buffer, "UTF8", "LATIN1","?"));
        }
        iptc_dataset_unref (ds);
    }
    iptcc.push_back (icsc);
}

//------inherited functions--------------//


std::string ImageMetaData::apertureToString (double aperture) {

    char buffer[256];
    sprintf (buffer, "%0.1f", aperture);
    return buffer;
}

std::string ImageMetaData::shutterToString (double shutter) {

    char buffer[256];
    if (shutter > 0.0 && shutter < 0.9)
        sprintf (buffer, "1/%0.0f", 1.0 / shutter);
    else
        sprintf (buffer, "%0.1f", shutter);
    return buffer;
}

double ImageMetaData::shutterFromString (std::string s) {

    int i = s.find_first_of ('/');
    if (i==std::string::npos)
        return atof (s.c_str());
    else 
        return atof (s.substr(0,i).c_str()) / atof (s.substr(i+1).c_str());
}

double ImageMetaData::apertureFromString (std::string s) {

    return atof (s.c_str());
}

extern "C" {

#include <libiptcdata/iptc-data.h>
#include <libiptcdata/iptc-jpeg.h>

struct _IptcDataPrivate
{
	unsigned int ref_count;

	IptcLog *log;
	IptcMem *mem;
};

IptcData *
iptc_data_new_from_jpeg_file (FILE *infile)
{
	IptcData *d;
	unsigned char * buf;
	int buf_len = 256*256;
	int len, offset;
        unsigned int iptc_len;

	if (!infile)
		return NULL;

	d = iptc_data_new ();
	if (!d) 
		return NULL;

	buf = (unsigned char*)iptc_mem_alloc (d->priv->mem, buf_len);
	if (!buf) {
		iptc_data_unref (d);
		return NULL;
	}

	len = iptc_jpeg_read_ps3 (infile, buf, buf_len);

	if (len <= 0) {
		goto failure;
	}

	offset = iptc_jpeg_ps3_find_iptc (buf, len, &iptc_len);
	if (offset <= 0) {
		goto failure;
	}

	iptc_data_load (d, buf + offset, iptc_len);

	iptc_mem_free (d->priv->mem, buf);
	return d;

failure:
	iptc_mem_free (d->priv->mem, buf);
	iptc_data_unref (d);
	return NULL;
}

}
