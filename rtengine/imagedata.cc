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

using namespace rtengine;

ImageMetaData* ImageMetaData::fromFile (const String& fname) {

    return new ImageData (fname);
}

ImageData::ImageData (const String& fname) 
	: expTime (0,1), fNumber (0), iso (0), lens ("Unknown"), defRot (0),
	  make ("Unknown"), model ("Unknown"), focalLen (0) {

	memset (&time, 0, sizeof(time));

	try {
		Exiv2::Image::AutoPtr image = Exiv2::ImageFactory::open (fname);
		image->readMetadata();
		exifData = image->exifData ();
		iptcData = image->iptcData ();
		xmpData = image->xmpData ();
		extractInfo ();
/*		std::cout << "Make: " << make << std::endl;
		std::cout << "Model: " << model << std::endl;
		std::cout << "Focal: " << focalLen << std::endl;
		std::cout << "ISO: " << iso << std::endl;
		std::cout << "ExpTime: " << exposureTimeToString(expTime) << std::endl;
		std::cout << "FNumber: " << fNumberToString(fNumber) << std::endl;
		std::cout << "Rotation: " << defRot << std::endl;
		std::cout << "Lens: " << lens << std::endl;
		std::cout << "Time: " << asctime (&time) << std::endl;*/
	}
	catch (const Exiv2::Error& e) {
	}
}

void ImageData::extractInfo () {

	Exiv2::ExifData::const_iterator i;

	std::stringstream ss;
	if ((i=Exiv2::make (exifData)) != exifData.end())
		i->write (ss, &exifData);
	make  = ss.str();
	
	ss.str(std::string());
	if ((i=Exiv2::model (exifData)) != exifData.end())
		i->write (ss, &exifData);
	model  = ss.str();

	ss.str(std::string());
	if ((i=Exiv2::lensName (exifData)) != exifData.end())
		i->write (ss, &exifData);
	lens  = ss.str();

	ss.str(std::string());
	if ((i=Exiv2::isoSpeed (exifData)) != exifData.end())
		i->write (ss, &exifData);
	iso  = atoi (ss.str().c_str());

	if ((i=Exiv2::focalLength (exifData)) != exifData.end())
		focalLen = i->toFloat ();
		
	if ((i=Exiv2::exposureTime (exifData)) != exifData.end())
		expTime = i->toRational ();
		
	if ((i=Exiv2::fNumber (exifData)) != exifData.end())
		fNumber = i->toFloat ();
	
	int ori = 0;
	if ((i=Exiv2::orientation (exifData)) != exifData.end())
		ori = i->toLong ();
	
	if (ori==6)
		defRot = 270;
	else if (ori==8)
		defRot = 90;
	else if (ori==3)
		defRot = 180;
	else 
		defRot = 0;
		
	if ((i=exifData.findKey (Exiv2::ExifKey ("Exif.Image.DateTimeOriginal"))) != exifData.end()) {
		std::istringstream ss (i->toString());
		char c;
		ss >> time.tm_year >> c;
		ss >> time.tm_mon >> c;
		ss >> time.tm_mday;
		ss >> time.tm_hour >> c;
		ss >> time.tm_min >> c;
		ss >> time.tm_sec;
		time.tm_year -= 1900;
		time.tm_mon -= 1;
	}
}

std::string ImageMetaData::fNumberToString (float fNumber) {

	std::stringstream ss;
	ss << std::setprecision(2) << fNumber;
	return ss.str ();
}

std::string ImageMetaData::exposureTimeToString (Exiv2::Rational expTime) {

	std::stringstream ss;
	if (expTime.first > 1 && expTime.second >= expTime.first) {
		expTime.first = 1;
		expTime.second = expTime.second / expTime.first;
	}
	else if (expTime.first > 1 && expTime.second < expTime.first) {
		expTime.second = 1;
		expTime.second = expTime.first / expTime.second;
	}
	if (expTime.second == 1)
		ss << expTime.first;
	else
		ss << expTime.first << '/' << expTime.second;
	return ss.str ();
}

Exiv2::Rational ImageMetaData::exposureTimeFromString (const std::string& s) {

    int i = s.find_first_of ('/');
    if (i==std::string::npos)
        return Exiv2::Rational (atof (s.c_str()), 1);
    else 
        return Exiv2::Rational (atof (s.substr(0,i).c_str()), atof (s.substr(i+1).c_str()));
}

float ImageMetaData::fNumberFromString (const std::string& s) {

    return atof (s.c_str());
}

std::string ImageData::getIptcKey (const std::string& rtIptcKey) {

	if (rtIptcKey == "Caption") 
		return String ("Iptc.Application2.Caption");
	else if (rtIptcKey == "CaptionWriter")
		return String ("Iptc.Application2.Writer");
	else if (rtIptcKey == "Headline") 
		return String ("Iptc.Application2.Headline");
	else if (rtIptcKey == "Instructions")
		return String ("Iptc.Application2.SpecialInstructions");
	else if (rtIptcKey == "Keywords") 
		return String ("Iptc.Application2.Keywords");
	else if (rtIptcKey == "Category")
		return String ("Iptc.Application2.Category");
	else if (rtIptcKey == "SupplementalCategories") 
		return String ("Iptc.Application2.SuppCategory");
	else if (rtIptcKey == "Author")
		return String ("Iptc.Application2.Byline");
	else if (rtIptcKey == "AuthorsPosition")
		return String ("Iptc.Application2.BylineTitle");
	else if (rtIptcKey == "Credit")
		return String ("Iptc.Application2.Credit");
	else if (rtIptcKey == "Source")
		return String ("Iptc.Application2.Source");
	else if (rtIptcKey == "Copyright")
		return String ("Iptc.Application2.Copyright");
	else if (rtIptcKey == "City")
		return String ("Iptc.Application2.City");
	else if (rtIptcKey == "Province")
		return String ("Iptc.Application2.ProvinceState");
	else if (rtIptcKey == "Country")
		return String ("Iptc.Application2.CountryName");
	else if (rtIptcKey == "Title")
		return String ("Iptc.Application2.ObjectName");
	else if (rtIptcKey == "DateCreated")
		return String ("Iptc.Application2.DateCreated");
	else if (rtIptcKey == "TransReference")
		return String ("Iptc.Application2.TransmissionReference");
}
