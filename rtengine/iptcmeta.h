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
#ifndef __IPTCMETA_H__
#define __IPTCMETA_H__


#include <map>
#include <glibmm.h>
#include <exiv2/exiv2.hpp>

namespace rtengine {

extern const char *kIPTCArtworkRights;
extern const char *kIPTCArtworkCreator;
extern const char *kIPTCArtworkDate;
extern const char *kIPTCArtworkSource;
extern const char *kIPTCArtworkNumber;
extern const char *kIPTCArtworkTitle;
extern const char *kIPTCAuthorPos;
extern const char *kIPTCCategory;
extern const char *kIPTCCity;
extern const char *kIPTCCountry;
extern const char *kIPTCCountryCode;
extern const char *kIPTCCreator;
extern const char *kIPTCCreatorAdrCity;
extern const char *kIPTCCreatorAdrCtry;
extern const char *kIPTCCreatorExtadr;
extern const char *kIPTCCreatorPcode;
extern const char *kIPTCCreatorRegion;
extern const char *kIPTCCreatorEmail;
extern const char *kIPTCCreatorTel;
extern const char *kIPTCCreatorUrl;
extern const char *kIPTCCredit;
extern const char *kIPTCCVterm;
extern const char *kIPTCDate;
extern const char *kIPTCDescription;
extern const char *kIPTCEvent;
extern const char *kIPTCGenre;
extern const char *kIPTCGUID;
extern const char *kIPTCGUIDSupplier;
extern const char *kIPTCHeadline;
extern const char *kIPTCImageCreator;
extern const char *kIPTCImageSupplier;
extern const char *kIPTCInstruct;
extern const char *kIPTCKeywords;
extern const char *kIPTCLicensor;
extern const char *kIPTCLocation;
extern const char *kIPTCLocationCity;
extern const char *kIPTCLocationCode;
extern const char *kIPTCLocationCtry;
extern const char *kIPTCLocationState;
extern const char *kIPTCLocationSubloc;
extern const char *kIPTCLocationRegion;
extern const char *kIPTCLocCreateCity;
extern const char *kIPTCLocCreateCode;
extern const char *kIPTCLocCreateCtry;
extern const char *kIPTCLocCreateState;
extern const char *kIPTCLocCreateSubloc;
extern const char *kIPTCLocCreateRegion;
extern const char *kIPTCMaxHeight;
extern const char *kIPTCMaxWidth;
extern const char *kIPTCMinorDisclosure;
extern const char *kIPTCModelAge;
extern const char *kIPTCModelInfo;
extern const char *kIPTCModelReleaseID;
extern const char *kIPTCModelReleaseSt;
extern const char *kIPTCOrganization;
extern const char *kIPTCPerson;
extern const char *kIPTCPlusVersion;
extern const char *kIPTCReference;
extern const char *kIPTCPropertyRelID;
extern const char *kIPTCPropertyRelSt;
extern const char *kIPTCRegistryID;
extern const char *kIPTCRegistryOrgID;
extern const char *kIPTCRights;
extern const char *kIPTCRightsOwner;
extern const char *kIPTCScene;
extern const char *kIPTCSource;
extern const char *kIPTCState;
extern const char *kIPTCSubjCode;
extern const char *kIPTCSuppCateg;
extern const char *kIPTCTitle;
extern const char *kIPTCUrgency;
extern const char *kIPTCUsageTerms;
extern const char *kIPTCWriter;


typedef std::map< std::string, Glib::ustring > IPTCPairList_t;

/** This class defines standard CORE and extended metadata for photos*/
class IPTCMeta
{
public:
	Exiv2::TypeId arrType; // { xaNone, xaAlt, xaBag, xaSeq }
	bool isStructArray;
	std::string key;
	Glib::ustring guiName;
	Glib::ustring description;

	IPTCMeta():arrType(Exiv2::xmpText){}
	IPTCMeta( const std::string &vkey, const Glib::ustring &vName, const Glib::ustring &vdesc, Exiv2::TypeId arr=Exiv2::xmpText,bool isArr= false ): key(vkey), guiName(vName), description(vdesc), arrType(arr),isStructArray(isArr){}

	/** Extract a Xmp key formatted for exiv2; index is meaningful only for array types */
	std::string getXmpKey( int index=0 ) const;

	typedef std::map< std::string, IPTCMeta > IPTCtagsList_t;

	static IPTCtagsList_t IPTCtags;  // Mapping between key and information
	static IPTCPairList_t IPTCScene; // Scene codes
	static IPTCPairList_t IPTCSubject; // Subject codes
	static IPTCPairList_t IPTCGenre; // Genre
	static IPTCPairList_t IPTCReleaseStatus; // Release status values
	static IPTCPairList_t IPTCWorldRegion;
	static IPTCPairList_t IPTCISO3166;

	static void initIPTCMeta();
};


} // namespace


#endif
