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
#include "imagedata.h"
#include "iptcpairs.h"
#include <glib/gstdio.h>
#include "safegtk.h"

#ifndef GLIBMM_EXCEPTIONS_ENABLED
#include <memory>
#endif

using namespace rtengine;

extern "C" IptcData *iptc_data_new_from_jpeg_file (FILE* infile);

ImageMetaData* ImageMetaData::fromFile (const Glib::ustring& fname, RawMetaDataLocation* rml) {

    return new ImageData (fname, rml);
}

ImageData::ImageData (Glib::ustring fname, RawMetaDataLocation* ri) {

    size_t dotpos = fname.find_last_of ('.');
    root = NULL;
    iptc = NULL;

    if (ri && (ri->exifBase>=0 || ri->ciffBase>=0)) {
        FILE* f = safe_g_fopen (fname, "rb");
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
    else if ((dotpos<fname.size()-3 && !fname.casefold().compare (dotpos, 4, ".jpg")) || (dotpos<fname.size()-4 && !fname.casefold().compare (dotpos, 5, ".jpeg"))) {
        FILE* f = safe_g_fopen (fname, "rb");
        if (f) {
            root = rtexif::ExifManager::parseJPEG (f);
            extractInfo ();
            fclose (f);
            FILE* ff = safe_g_fopen (fname, "rb");
            iptc = iptc_data_new_from_jpeg_file (ff);
            fclose (ff);
        }
    }
    else if ((dotpos<fname.size()-3 && !fname.casefold().compare (dotpos, 4, ".tif")) || (dotpos<fname.size()-4 && !fname.casefold().compare (dotpos, 5, ".tiff"))) {
        FILE* f = safe_g_fopen (fname, "rb");
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
        orientation = "Unknown";
        expcomp = 0;
        focal_len = 0;
        memset (&time, 0, sizeof(time));
    }
}

void ImageData::extractInfo () {

  if (!root)
    return;

  make = "";
  model = "";
  serial = "";
  orientation = "";
  expcomp = 0;
  shutter = 0;
  aperture = 0;
  focal_len = focal_len35mm = 0;
  focus_dist = 0;
  iso_speed = 0;
  memset (&time, 0, sizeof(time));
  timeStamp = 0;

  if (root->getTag ("Make")){
     make = root->getTag ("Make")->valueToString ();
     // same dcraw treatment
     static const char *corp[] =
       { "Canon", "NIKON", "EPSON", "KODAK", "Kodak", "OLYMPUS", "PENTAX", "RICOH",
         "MINOLTA", "Minolta", "Konica", "CASIO", "Sinar", "Phase One",
         "SAMSUNG", "Mamiya", "MOTOROLA", "Leaf" };
     for (size_t i=0; i < (sizeof(corp)/sizeof(*corp)); i++)
       if ( make.find( corp[i] ) != std::string::npos ){		/* Simplify company names */
   	     make = corp[i];
   	     break;
       }
       make.erase( make.find_last_not_of(' ')+1 );
  }

  if (root->getTag ("Model"))
     model = root->getTag ("Model")->valueToString ();
  if (!(model.size()==0)) {
     std::size_t i=0;
     if (  make.find("KODAK") != std::string::npos ){
   	  if( (i = model.find(" DIGITAL CAMERA")) !=  std::string::npos ||
   	      (i = model.find(" Digital Camera")) !=  std::string::npos ||
   	      (i = model.find("FILE VERSION")) )
        model.resize( i );
     }

     model.erase( model.find_last_not_of(' ')+1 );

     //if( (i=model.find( make )) != std::string::npos )
     if( !strncasecmp (model.c_str(), make.c_str(), make.size()) )
     	if( model.at( make.size() )==' ')
     	   model.erase(0,make.size()+1);
     if( model.find( "Digital Camera ") != std::string::npos )
     	model.erase(0,15);
  }
  else {
      model = "Unknown";
  }
  if (root->getTag ("Orientation")){
     orientation = root->getTag ("Orientation")->valueToString ();
  }

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
    if (exif->getTag ("ExposureBiasValue"))
        expcomp = exif->getTag ("ExposureBiasValue")->toDouble ();
    if (exif->getTag ("FocalLength"))
        focal_len = exif->getTag ("FocalLength")->toDouble ();
    if (exif->getTag ("FocalLengthIn35mmFilm"))
        focal_len35mm = exif->getTag ("FocalLengthIn35mmFilm")->toDouble ();

    // Focus distance from EXIF or XMP. MakerNote ones are scattered and partly encrypted
    int num=-3, denom=-3;

    // First try, offical EXIF. Set by Adobe on some DNGs
    rtexif::Tag* pDst=exif->getTag("SubjectDistance");
    if (pDst) {
        int num, denom;
        pDst->toRational(num,denom);
    } else {
        // Second try, XMP data
        char sXMPVal[64];
        if (root->getXMPTagValue("aux:ApproximateFocusDistance",sXMPVal)) { sscanf(sXMPVal,"%d/%d",&num,&denom); }
    }

    if (num!=-3) {
        if ((denom==1 && num>=10000) || num<0 || denom<0)
            focus_dist=10000;  // infinity
        else if (denom>0) {
            focus_dist=(float)num/denom;
        }
    }

    if (exif->getTag ("ISOSpeedRatings"))
        iso_speed = exif->getTag ("ISOSpeedRatings")->toDouble ();
    if (exif->getTag ("DateTimeOriginal")) {
        if (sscanf ((const char*)exif->getTag("DateTimeOriginal")->getValue(), "%d:%d:%d %d:%d:%d", &time.tm_year, &time.tm_mon, &time.tm_mday, &time.tm_hour, &time.tm_min, &time.tm_sec) == 6) {
            time.tm_year -= 1900;
            time.tm_mon -= 1;
            time.tm_isdst = -1;
            timeStamp = mktime(&time);
        }
    }
    rtexif::Tag *snTag = exif->findTag ("SerialNumber");
    if(!snTag)
    	snTag = exif->findTag ("InternalSerialNumber");
    if ( snTag )
        serial = snTag->valueToString();
    // guess lens...
    lens = "Unknown";

    // Sometimes (e.g. DNG) EXIF already contains lens data

	if(!make.compare (0, 8, "FUJIFILM")) {
		if(exif->getTag ("LensModel")) {
			lens = exif->getTag ("LensModel")->valueToString ();
		}
	} else if (root->findTag("MakerNote")) {
        rtexif::TagDirectory* mnote = root->findTag("MakerNote")->getDirectory();
        if (mnote && !make.compare (0, 5, "NIKON")) {
            // ISO at max value supported, check manufacturer specific
            if (iso_speed == 65535 || iso_speed == 0) {
                rtexif::Tag* isoTag = mnote->getTagP("ISOInfo/ISO");
                if (isoTag)
                    iso_speed = isoTag->toInt();
            }
            bool lensOk = false;
            if (mnote->getTag ("LensData")) {
                std::string ldata = mnote->getTag ("LensData")->valueToString ();
                int pos;
                if (ldata.size()>10 && (pos=ldata.find ("Lens = "))!=Glib::ustring::npos) {
                    lens = ldata.substr (pos + 7);
                    if (lens.compare (0, 7, "Unknown"))
                        lensOk = true;
					else {
						int pos = lens.find("$FL$");		// is there a placeholder for focallength?
						if(pos != Glib::ustring::npos) {				// then fill in focallength
							lens = lens.replace(pos,4,exif->getTag ("FocalLength")->valueToString ());
							if(mnote->getTag ("LensType")) {
								std::string ltype = mnote->getTag ("LensType")->valueToString ();
								if(ltype.find("MF = Yes")!=Glib::ustring::npos)		// check, whether it's a MF lens, should be always
									lens = lens.replace(0,7,"MF");
								lensOk = true;
							}
						}
					}
                }
            }
            if (!lensOk && mnote->getTag ("Lens")) {
                std::string ldata = mnote->getTag ("Lens")->valueToString ();
                size_t i=0, j=0;
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
                // Look whether it's MF or AF
				if(mnote->getTag ("LensType")) {
					std::string ltype = mnote->getTag ("LensType")->valueToString ();
					if(ltype.find("MF = Yes")!=Glib::ustring::npos)		// check, whether it's a MF lens
						lens = lens.replace(0,7,"MF");					// replace 'Unknwon' with 'MF'
					else
						lens = lens.replace(0,7,"AF");					// replace 'Unknwon' with 'AF'
				}
            }
        }
        else if (mnote && !make.compare (0, 5, "Canon")) {
            // ISO at max value supported, check manufacturer specific
            if (iso_speed == 65535 || iso_speed == 0) {
                rtexif::Tag* baseIsoTag = mnote->getTagP("CanonShotInfo/BaseISO");
                if (baseIsoTag)
                    iso_speed = baseIsoTag->toInt();
            }
            int found=false;
            // canon EXIF have a string for lens model
            rtexif::Tag *lt = mnote->getTag("LensType");
            if ( lt ) {
                std::string ldata = lt->valueToString ();
                if (ldata.size()>1) {
                    found=true;
                    lens = "Canon " + ldata;
                }
            }
            if( !found || lens.substr(lens.find(' ')).length() < 7 ){
                lt = mnote->findTag("LensID");
                if ( lt ) {
                    std::string ldata = lt->valueToString ();
                    if (ldata.size()>1) {
                        lens = ldata;
                    }
                }
            }
        }
        else if (mnote && (!make.compare (0, 6, "PENTAX") || (!make.compare (0, 5, "RICOH") && !model.compare (0, 6, "PENTAX")))) {
            if (mnote->getTag ("LensType"))
                lens = mnote->getTag ("LensType")->valueToString ();

            // Try to get the FocalLength from the LensInfo structure, where length below 10mm will be correctly set
            rtexif::Tag* flt=mnote->getTagP ("LensInfo/FocalLength");
            if (flt)
                focal_len = flt->toDouble ();
            else if ((flt = mnote->getTagP ("FocalLength"))) {
                rtexif::Tag* flt = mnote->getTag ("FocalLength");
                focal_len = flt->toDouble ();
            }

            if (mnote->getTag ("FocalLengthIn35mmFilm"))
                focal_len35mm = mnote->getTag ("FocalLengthIn35mmFilm")->toDouble ();
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
    } else if (exif->getTag ("DNGLensInfo")) {
        lens = exif->getTag ("DNGLensInfo")->valueToString ();
    } else if (exif->getTag ("LensModel")) {
        lens = exif->getTag ("LensModel")->valueToString ();
    } else if (exif->getTag ("LensInfo")) {
        lens = exif->getTag ("LensInfo")->valueToString ();
    }
  }
}

ImageData::~ImageData () {

    delete root;
    if (iptc)
        iptc_data_free (iptc);
}

const procparams::IPTCPairs ImageData::getIPTCData () const {

    procparams::IPTCPairs iptcc;
    if (!iptc)
        return iptcc;

    unsigned char buffer[2100];
    for (int i=0; i<16; i++) {
        IptcDataSet* ds = iptc_data_get_next_dataset (iptc, NULL, IPTC_RECORD_APP_2, strTags[i].tag);
        if (ds) {
            iptc_dataset_get_data (ds, buffer, 2100);
            std::vector<Glib::ustring> icValues;
            icValues.push_back (safe_locale_to_utf8((char*)buffer));

            iptcc[strTags[i].field] = icValues;
            iptc_dataset_unref (ds);
        }
    }
    IptcDataSet* ds = NULL;
    std::vector<Glib::ustring> keywords;
    while ((ds=iptc_data_get_next_dataset (iptc, ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS))) {
        iptc_dataset_get_data (ds, buffer, 2100);
        keywords.push_back (safe_locale_to_utf8((char*)buffer));
    }
    iptcc["Keywords"] = keywords;
    ds = NULL;
    std::vector<Glib::ustring> suppCategories;
    while ((ds=iptc_data_get_next_dataset (iptc, ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY))) {
        iptc_dataset_get_data (ds, buffer, 2100);
        suppCategories.push_back (safe_locale_to_utf8((char*)buffer));
        iptc_dataset_unref (ds);
    }
    iptcc["SupplementalCategories"] = suppCategories;
    return iptcc;
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

std::string ImageMetaData::expcompToString (double expcomp, bool maskZeroexpcomp) {

    char buffer[256];
    if (maskZeroexpcomp==true){
        if (expcomp!=0.0){
    	    sprintf (buffer, "%0.2f", expcomp);
    	    return buffer;
        }
        else
    	    return "";
    }
    else{
    	sprintf (buffer, "%0.2f", expcomp);
    	return buffer;
    }
}

double ImageMetaData::shutterFromString (std::string s) {

    size_t i = s.find_first_of ('/');
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
