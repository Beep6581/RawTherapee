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
#include <iptcmeta.h>

namespace rtengine{


IPTCMeta::IPTCtagsList_t IPTCMeta::IPTCtags;
IPTCPairList_t IPTCMeta::IPTCScene;
IPTCPairList_t IPTCMeta::IPTCSubject;
IPTCPairList_t IPTCMeta::IPTCGenre;
IPTCPairList_t IPTCMeta::IPTCReleaseStatus;
IPTCPairList_t IPTCMeta::IPTCCopyrightStatus;
IPTCPairList_t IPTCMeta::IPTCWorldRegion;
IPTCPairList_t IPTCMeta::IPTCISO3166;
IPTCPairList_t IPTCMeta::IPTCBoolean;

const char *kIPTCArtworkRights  = "Iptc4xmpExt:ArtworkOrObjectDetails/Iptc4xmpExt:AOCopyrightNotice";
const char *kIPTCArtworkCreator = "Iptc4xmpExt:ArtworkOrObjectDetails/Iptc4xmpExt:AOCreator";
const char *kIPTCArtworkDate    = "Iptc4xmpExt:ArtworkOrObjectDetails/Iptc4xmpExt:AODateCreated";
const char *kIPTCArtworkSource  = "Iptc4xmpExt:ArtworkOrObjectDetails/Iptc4xmpExt:AOSource";
const char *kIPTCArtworkNumber  = "Iptc4xmpExt:ArtworkOrObjectDetails/Iptc4xmpExt:AOSourceInvNo";
const char *kIPTCArtworkTitle   = "Iptc4xmpExt:ArtworkOrObjectDetails/Iptc4xmpExt:AOTitle";
const char *kIPTCAuthorPos  = "photoshop:AuthorsPosition";
const char *kIPTCCategory   = "photoshop:Category";
const char *kIPTCCity       = "photoshop:City";
const char *kIPTCCountry    = "photoshop:Country";
const char *kIPTCCountryCode= "Iptc4xmpCore:CountryCode";
const char *kIPTCCreator    = "dc:creator";
const char *kIPTCCreatorAdrCity = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiAdrCity";
const char *kIPTCCreatorAdrCtry = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiAdrCtry";
const char *kIPTCCreatorEmail   = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiEmailWork";
const char *kIPTCCreatorExtadr  = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiAdrExtadr";
const char *kIPTCCreatorPcode   = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiAdrPcode";
const char *kIPTCCreatorRegion  = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiAdrRegion";
const char *kIPTCCreatorTel		= "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiTelWork";
const char *kIPTCCreatorUrl     = "Iptc4xmpCore:CreatorContactInfo/Iptc4xmpCore:CiUrlWork";
const char *kIPTCCredit     = "photoshop:Credit";
const char *kIPTCCVterm     = "Iptc4xmpExt:CVterm";
const char *kIPTCDate       = "photoshop:DateCreated";
const char *kIPTCDescription= "dc:description";
const char *kIPTCEvent      = "Iptc4xmpExt:Event";
const char *kIPTCGenre      = "Iptc4xmpCore:IntellectualGenre";
const char *kIPTCGUID       = "Iptc4xmpExt:DigImageGUID";
const char *kIPTCGUIDSupplier   = "plus:ImageSupplierImageID";
const char *kIPTCHeadline   = "photoshop:Headline";
const char *kIPTCImageCreator   = "plus:ImageCreator/plus:ImageCreatorName";
const char *kIPTCImageSupplier  = "plus:ImageSupplier/plus:ImageSupplierName";
const char *kIPTCInstruct   = "photoshop:Instructions";
const char *kIPTCKeywords   = "dc:subject";
const char *kIPTCLicensor   = "plus:Licensor/plus:LicensorName";
const char *kIPTCLocation   = "Iptc4xmpCore:Location";
const char *kIPTCLocationCity   = "Iptc4xmpExt:LocationShown/Iptc4xmpExt:City";
const char *kIPTCLocationCode   = "Iptc4xmpExt:LocationShown/Iptc4xmpExt:CV-code";
const char *kIPTCLocationCtry   = "Iptc4xmpExt:LocationShown/Iptc4xmpExt:CountryName";
const char *kIPTCLocationState  = "Iptc4xmpExt:LocationShown/Iptc4xmpExt:ProvinceState";
const char *kIPTCLocationSubloc = "Iptc4xmpExt:LocationShown/Iptc4xmpExt:Sublocation";
const char *kIPTCLocationRegion = "Iptc4xmpExt:LocationShown/Iptc4xmpExt:WorldRegion";
const char *kIPTCLocCreateCity  = "Iptc4xmpExt:LocationCreated/Iptc4xmpExt:City";
const char *kIPTCLocCreateCode  = "Iptc4xmpExt:LocationCreated/Iptc4xmpExt:CV-code";
const char *kIPTCLocCreateCtry  = "Iptc4xmpExt:LocationCreated/Iptc4xmpExt:CountryName";
const char *kIPTCLocCreateState = "Iptc4xmpExt:LocationCreated/Iptc4xmpExt:ProvinceState";
const char *kIPTCLocCreateSubloc= "Iptc4xmpExt:LocationCreated/Iptc4xmpExt:Sublocation";
const char *kIPTCLocCreateRegion= "Iptc4xmpExt:LocationCreated/Iptc4xmpExt:WorldRegion";
const char *kIPTCMinorDisclosure= "plus:MinorModelAgeDisclosure";
const char *kIPTCMaxHeight  = "Iptc4xmpExt:MaxAvailHeight";
const char *kIPTCMaxWidth   = "Iptc4xmpExt:MaxAvailWidth";
const char *kIPTCModelAge   = "Iptc4xmpExt:ModelAge";
const char *kIPTCModelInfo  = "Iptc4xmpExt:AddlModelInfo";
const char *kIPTCModelReleaseID = "plus:ModelReleaseID";
const char *kIPTCModelReleaseSt = "plus:ModelReleaseStatus";
const char *kIPTCOrganization   = "Iptc4xmpExt:OrganisationInImageCode";
const char *kIPTCPerson     = "Iptc4xmpExt:PersonInImage";
const char *kIPTCPlusVersion    = "plus:Version";
const char *kIPTCPropertyRelID  = "plus:PropertyReleaseID";
const char *kIPTCPropertyRelSt  = "plus:PropertyReleaseStatus";
const char *kIPTCReference  = "photoshop:TransmissionReference";
const char *kIPTCRegistryID = "Iptc4xmpExt:RegistryEntryDetails/Iptc4xmpExt:RegItemId";
const char *kIPTCRegistryOrgID  = "Iptc4xmpExt:RegistryEntryDetails/Iptc4xmpExt:RegOrgId";
const char *kIPTCRights     = "dc:rights";
const char *kIPTCRightsCertif   = "xmpRights:Certificate";
const char *kIPTCRightsMarked   = "xmpRights:Marked";
const char *kIPTCRightsOwner    = "plus:CopyrightOwner/plus:CopyrightOwnerName";
const char *kIPTCRightsStatement= "xmpRights:WebStatement";
const char *kIPTCRightsStatus   = "plus:CopyrightStatus";
const char *kIPTCScene      = "Iptc4xmpCore:Scene";
const char *kIPTCSource     = "photoshop:Source";
const char *kIPTCState      = "photoshop:State";
const char *kIPTCSubjCode   = "Iptc4xmpCore:SubjectCode";
const char *kIPTCSuppCateg  = "photoshop:SupplementalCategories";
const char *kIPTCTitle      = "dc:title";
const char *kIPTCUrgency    = "photoshop:Urgency";
const char *kIPTCUsageTerms = "xmpRights:UsageTerms";
const char *kIPTCWriter     = "photoshop:CaptionWriter";

std::string IPTCMeta::getXmpKey( int index ) const
{
	std::string xKey(key);
	std::string::size_type i = xKey.find_first_of(':');
	if( i != std::string::npos )
		xKey.replace(i,1,".");
	xKey = std::string("Xmp.")+ xKey;

	if( (arrType == Exiv2::xmpBag || arrType == Exiv2::xmpSeq ) && index>0){
		xKey = xKey + Glib::ustring::compose( "[%1]",index );
	}
	return xKey;
}

void IPTCMeta::initIPTCMeta()
{
	IPTCtags[kIPTCHeadline] = IPTCMeta(kIPTCHeadline, "IPTC_HEADLINE", "IPTC_HEADLINE_HINT");
	IPTCtags[kIPTCCity]     = IPTCMeta(kIPTCCity,     "IPTC_CITY", "IPTC_CITY_HINT");
	IPTCtags[kIPTCCountry]  = IPTCMeta(kIPTCCountry,  "IPTC_COUNTRY", "IPTC_COUNTRY_HINT");
    IPTCtags[kIPTCState]    = IPTCMeta(kIPTCState,    "IPTC_STATE", "IPTC_STATE_HINT");
    IPTCtags[kIPTCAuthorPos]= IPTCMeta(kIPTCAuthorPos,"IPTC_AUTHOR_POS","IPTC_AUTHOR_POS_HINT");

    IPTCtags[kIPTCWriter]   = IPTCMeta(kIPTCWriter, "IPTC_WRITER","IPTC_WRITER_HINT");
    IPTCtags[kIPTCInstruct] = IPTCMeta(kIPTCInstruct, "IPTC_INSTRUCTIONS", "IPTC_INSTRUCTIONS_HINT");
    IPTCtags[kIPTCReference]= IPTCMeta(kIPTCReference,"IPTC_JOBID","IPTC_JOBID_HINT");
    IPTCtags[kIPTCCredit]   = IPTCMeta(kIPTCCredit, "IPTC_CREDIT" ,"IPTC_CREDIT_HINT");
    IPTCtags[kIPTCSource]   = IPTCMeta(kIPTCSource, "IPTC_SOURCE", "IPTC_SOURCE_HINT");
    IPTCtags[kIPTCUrgency]  = IPTCMeta(kIPTCUrgency, "IPTC_URGENCY", "IPTC_URGENCY_HINT");
    IPTCtags[kIPTCCategory] = IPTCMeta(kIPTCCategory,"IPTC_CATEGORY", "IPTC_CATEGORY_HINT");
    IPTCtags[kIPTCSuppCateg]= IPTCMeta(kIPTCSuppCateg,"IPTC_SUPPLCATEGORIES","IPTC_SUPPLCATEGORIES_HINT");

    IPTCtags[kIPTCLocation]     = IPTCMeta(kIPTCLocation, "IPTC_SUBLOCATION", "IPTC_SUBLOCATION_HINT");
    IPTCtags[kIPTCCountryCode]  = IPTCMeta(kIPTCCountryCode, "IPTC_COUNTRYCODE","IPTC_COUNTRYCODE_HINT");
    IPTCtags[kIPTCGenre]        = IPTCMeta(kIPTCGenre,"IPTC_GENRE" ,"IPTC_GENRE_HINT");
    IPTCtags[kIPTCScene]        = IPTCMeta(kIPTCScene, "IPTC_SCENE", "IPTC_SCENE_HINT", Exiv2::xmpBag );
    IPTCtags[kIPTCSubjCode]     = IPTCMeta(kIPTCSubjCode,"IPTC_SUBJECT","IPTC_SUBJECT_HINT", Exiv2::xmpBag );

    IPTCtags[kIPTCKeywords]     = IPTCMeta(kIPTCKeywords,"IPTC_KEYWORDS","IPTC_KEYWORDS_HINT", Exiv2::xmpBag );
    IPTCtags[kIPTCCreator]      = IPTCMeta(kIPTCCreator,"IPTC_AUTHOR","IPTC_AUTHOR_HINT", Exiv2::langAlt );
    IPTCtags[kIPTCTitle]        = IPTCMeta(kIPTCTitle, "IPTC_TITLE","IPTC_TITLE_HINT", Exiv2::langAlt );
    IPTCtags[kIPTCRights]       = IPTCMeta(kIPTCRights,"IPTC_RIGHTS", "IPTC_RIGHTS_HINT", Exiv2::langAlt );
    IPTCtags[kIPTCDescription]  = IPTCMeta(kIPTCDescription, "IPTC_DESCRIPTION", "IPTC_DESCRIPTION_HINT", Exiv2::langAlt );

    IPTCtags[kIPTCUsageTerms]   = IPTCMeta(kIPTCUsageTerms,"IPTC_USAGETERMS","IPTC_USAGETERMS_HINT",Exiv2::langAlt);

    IPTCtags[kIPTCCreatorAdrCity] = IPTCMeta(kIPTCCreatorAdrCity,"IPTC_CREATORCITY","IPTC_CREATORCITY_HINT" );
    IPTCtags[kIPTCCreatorAdrCtry] = IPTCMeta(kIPTCCreatorAdrCtry,"IPTC_CREATORCOUNTRY","IPTC_CREATORCOUNTRY_HINT");
    IPTCtags[kIPTCCreatorExtadr]  = IPTCMeta(kIPTCCreatorExtadr,"IPTC_CREATORADDRESS","IPTC_CREATORADDRESS_HINT");
    IPTCtags[kIPTCCreatorPcode]   = IPTCMeta(kIPTCCreatorPcode,"IPTC_CREATORPCODE","IPTC_CREATORPCODE_HINT");
    IPTCtags[kIPTCCreatorRegion]  = IPTCMeta(kIPTCCreatorRegion,"IPTC_CREATORSTATE","IPTC_CREATORSTATE_HINT");
    IPTCtags[kIPTCCreatorEmail]   = IPTCMeta(kIPTCCreatorEmail,"IPTC_CREATOREMAIL","IPTC_CREATOREMAIL_HINT");
    IPTCtags[kIPTCCreatorTel]     = IPTCMeta(kIPTCCreatorTel,"IPTC_CREATORPHONE","IPTC_CREATORPHONE_HINT");
    IPTCtags[kIPTCCreatorUrl]     = IPTCMeta(kIPTCCreatorUrl,"IPTC_CREATORURL","IPTC_CREATORURL_HINT");

    IPTCtags[kIPTCModelInfo]      = IPTCMeta(kIPTCModelInfo,"IPTC_MODELINFO","IPTC_MODELINFO_HINT");
    IPTCtags[kIPTCModelAge]       = IPTCMeta(kIPTCModelAge,"IPTC_MODELAGE","IPTC_MODELAGE_HINT",Exiv2::xmpBag);
    IPTCtags[kIPTCPerson]         = IPTCMeta(kIPTCPerson,"IPTC_PERSON","IPTC_PERSON_HINT",Exiv2::xmpBag);
    IPTCtags[kIPTCOrganization]   = IPTCMeta(kIPTCOrganization,"IPTC_ORGANIZATION","IPTC_ORGANIZATION_HINT");
    IPTCtags[kIPTCCVterm]         = IPTCMeta(kIPTCCVterm,"IPTC_VOCABULARY","IPTC_VOCABULARY_HINT",Exiv2::xmpBag );
    IPTCtags[kIPTCEvent]          = IPTCMeta(kIPTCEvent,"IPTC_EVENT","IPTC_EVENT_HINT",Exiv2::langAlt);

    IPTCtags[kIPTCLocationCity]   = IPTCMeta(kIPTCLocationCity,"IPTC_LOCATIONCITY","IPTC_LOCATIONCITY_HINT");
    IPTCtags[kIPTCLocationCode]   = IPTCMeta(kIPTCLocationCode,"IPTC_LOCATIONCOUNTRYCODE","IPTC_LOCATIONCOUNTRYCODE_HINT");
    IPTCtags[kIPTCLocationCtry]   = IPTCMeta(kIPTCLocationCtry,"IPTC_LOCATIONCOUNTRY","IPTC_LOCATIONCOUNTRY_HINT");
    IPTCtags[kIPTCLocationState]  = IPTCMeta(kIPTCLocationState,"IPTC_LOCATIONSTATE","IPTC_LOCATIONSTATE_HINT");
    IPTCtags[kIPTCLocationSubloc] = IPTCMeta(kIPTCLocationSubloc,"IPTC_LOCATIONSUBLOC","IPTC_LOCATIONSUBLOC_HINT");
    IPTCtags[kIPTCLocationRegion] = IPTCMeta(kIPTCLocationRegion,"IPTC_LOCATIONREGION","IPTC_LOCATIONREGION_HINT");

    IPTCtags[kIPTCLocCreateCity]  = IPTCMeta(kIPTCLocCreateCity,"IPTC_CREATEDCITY","IPTC_CREATEDCITY_HINT");
    IPTCtags[kIPTCLocCreateCode]  = IPTCMeta(kIPTCLocCreateCode,"IPTC_CREATEDCOUNTRYCODE","IPTC_CREATEDCOUNTRYCODE_HINT");
    IPTCtags[kIPTCLocCreateCtry]  = IPTCMeta(kIPTCLocCreateCtry,"IPTC_CREATEDCOUNTRY","IPTC_CREATEDCOUNTRY_HINT");
    IPTCtags[kIPTCLocCreateState] = IPTCMeta(kIPTCLocCreateState,"IPTC_CREATEDSTATE","IPTC_CREATEDSTATE_HINT");
    IPTCtags[kIPTCLocCreateSubloc]= IPTCMeta(kIPTCLocCreateSubloc,"IPTC_CREATEDSUBLOC","IPTC_CREATEDSUBLOC_HINT");
    IPTCtags[kIPTCLocCreateRegion]= IPTCMeta(kIPTCLocCreateRegion,"IPTC_CREATEDREGION","IPTC_CREATEDREGION_HINT");

    IPTCtags[kIPTCArtworkRights]  = IPTCMeta(kIPTCArtworkRights,"IPTC_ARTWORKRIGHTS","IPTC_ARTWORKRIGHTS_HINT");
    IPTCtags[kIPTCArtworkCreator] = IPTCMeta(kIPTCArtworkCreator,"IPTC_ARTWORKCREATOR","IPTC_ARTWORKCREATOR_HINT");
    IPTCtags[kIPTCArtworkDate]    = IPTCMeta(kIPTCArtworkDate,"IPTC_ARTWORKDATE","IPTC_ARTWORKDATE_HINT");
    IPTCtags[kIPTCArtworkSource]  = IPTCMeta(kIPTCArtworkSource,"IPTC_ARTWORKSOURCE","IPTC_ARTWORKSOURCE_HINT");
    IPTCtags[kIPTCArtworkNumber]  = IPTCMeta(kIPTCArtworkNumber,"IPTC_ARTWORKNUMBER","IPTC_ARTWORKNUMBER_HINT");
    IPTCtags[kIPTCArtworkTitle]   = IPTCMeta(kIPTCArtworkTitle,"IPTC_ARTWORKTITLE","IPTC_ARTWORKTITLE_HINT");

    IPTCtags[kIPTCRightsOwner]    = IPTCMeta(kIPTCRightsOwner,"IPTC_RIGHTSOWNER","IPTC_RIGHTSOWNER_HINT",Exiv2::xmpBag);
    IPTCtags[kIPTCRightsStatus]   = IPTCMeta(kIPTCRightsStatus,"IPTC_RIGHTSSTATUS","IPTC_RIGHTSSTATUS_HINT" );
    IPTCtags[kIPTCRightsMarked]   = IPTCMeta(kIPTCRightsMarked,"IPTC_RIGHTSMARKED","IPTC_RIGHTSMARKED_HINT");
    IPTCtags[kIPTCRightsCertif]   = IPTCMeta(kIPTCRightsCertif,"IPTC_RIGHTSCERTIFICATE","IPTC_RIGHTSCERTIFICATE_HINT");
    IPTCtags[kIPTCRightsStatement]= IPTCMeta(kIPTCRightsStatement,"IPTC_RIGHTSSTATEMENT","IPTC_RIGHTSSTATEMENT_HINT");
    IPTCtags[kIPTCImageCreator]   = IPTCMeta(kIPTCImageCreator,"IPTC_IMAGECREATOR","IPTC_IMAGECREATOR_HINT");
    IPTCtags[kIPTCLicensor]       = IPTCMeta(kIPTCLicensor,"IPTC_LICENSOR","IPTC_LICENSOR_HINT");
    IPTCtags[kIPTCImageSupplier]  = IPTCMeta(kIPTCImageSupplier,"IPTC_IMAGESUPPLIER","IPTC_IMAGESUPPLIER_HINT");
    IPTCtags[kIPTCGUIDSupplier]   = IPTCMeta(kIPTCGUIDSupplier,"IPTC_IMAGESUPPLIERGUID","IPTC_IMAGESUPPLIERGUID_HINT");
    IPTCtags[kIPTCMinorDisclosure]= IPTCMeta(kIPTCMinorDisclosure,"IPTC_MINORDISCLOSURE","IPTC_MINORDISCLOSURE_HINT");
    IPTCtags[kIPTCModelReleaseID] = IPTCMeta(kIPTCModelReleaseID,"IPTC_MODELRELID","IPTC_MODELRELID_HINT", Exiv2::xmpBag);
    IPTCtags[kIPTCModelReleaseSt] = IPTCMeta(kIPTCModelReleaseSt,"IPTC_MODELRELSTATUS","IPTC_MODELRELSTATUS_HINT");
    IPTCtags[kIPTCPropertyRelID]  = IPTCMeta(kIPTCPropertyRelID,"IPTC_PROPERTYRELID","IPTC_PROPERTYRELID_HINT", Exiv2::xmpBag);
    IPTCtags[kIPTCPropertyRelSt]  = IPTCMeta(kIPTCPropertyRelSt,"IPTC_PROPERTYRELSTATUS","IPTC_PROPERTYRELSTATUS_HINT");
    IPTCtags[kIPTCRegistryID]     = IPTCMeta(kIPTCRegistryID,"IPTC_REGISTRYID","IPTC_REGISTRYID_HINT");
    IPTCtags[kIPTCRegistryOrgID]  = IPTCMeta(kIPTCRegistryOrgID,"IPTC_REGISTRYORG","IPTC_REGISTRYORG_HINT");

    // readonly
    IPTCtags[kIPTCDate]           = IPTCMeta(kIPTCDate,"IPTC_DATE" ,"IPTC_DATE_HINT");
    IPTCtags[kIPTCDate].readOnly  = true;
    IPTCtags[kIPTCGUID]           = IPTCMeta(kIPTCGUID,"IPTC_IMAGEGUID","IPTC_IMAGEGUID_HINT");
    IPTCtags[kIPTCGUID].readOnly  = true;
    IPTCtags[kIPTCMaxHeight]      = IPTCMeta(kIPTCMaxHeight,"IPTC_IMAGEMAXHEIGHT","IPTC_IMAGEMAXHEIGHT_HINT");
    IPTCtags[kIPTCMaxHeight].readOnly = true;
    IPTCtags[kIPTCMaxWidth]       = IPTCMeta(kIPTCMaxWidth,"IPTC_IMAGEMAXWIDTH","IPTC_IMAGEMAXWIDTH_HINT");
    IPTCtags[kIPTCMaxWidth].readOnly = true;
    // PLUS versione is 1.2.0
    IPTCtags[kIPTCPlusVersion]    = IPTCMeta(kIPTCPlusVersion,"IPTC_PLUSVERSION","IPTC_PLUSVERSION_HINT");
    IPTCtags[kIPTCPlusVersion].readOnly = true;

    /* Scene codes: http://cv.iptc.org/newscodes/scene/ */
    IPTCScene["010100"] = "headshot:a head only view of a person (or animal/s) or persons as in a montage";
    IPTCScene["010200"] = "half-length: a torso and head view of a person or persons";
    IPTCScene["010300"] = "full-length: a view from head to toe of a person or persons";
    IPTCScene["010400"] = "profile: a view of a person from the side";
    IPTCScene["010500"] = "rear view: a view of a person or persons from the rear";
    IPTCScene["010600"] = "single: a view of only one person, object or animal";
    IPTCScene["010700"] = "couple: a view of two people who are in a personal relationship, for example engaged, married or in a romantic partnership";
    IPTCScene["010800"] = "two: a view of two people";
    IPTCScene["010900"] = "group: a view of more than two people";
    IPTCScene["011000"] = "general view: an overall view of the subject and its surrounds";
    IPTCScene["011100"] = "panoramic view: a panoramic or wide angle view of a subject and its surrounds";
    IPTCScene["011200"] = "aerial view: a view taken from above";
    IPTCScene["011300"] = "under-water: a photo taken under water";
    IPTCScene["011400"] = "night scene: a photo taken during darkness";
    IPTCScene["011500"] = "satellite: a photo taken from a satellite in orbit";
    IPTCScene["011600"] = "exterior view: a photo that shows the exterior of a building or other object";
    IPTCScene["011700"] = "interior view: a scene or view of the interior of a building or other object";
    IPTCScene["011800"] = "close-up: a view of, or part of a person/object taken at close range in order to emphasize detail or accentuate mood. Macro photography";
    IPTCScene["011900"] = "action: a subject in motion such as children jumping, horse running";
    IPTCScene["012000"] = "performing: subject or subjects on a stage performing to an audience";
    IPTCScene["012100"] = "posing: subject or subjects posing such as a 'victory' pose or other stance that symbolizes leadership";
    IPTCScene["012200"] = "symbolic: a posed picture symbolizing an event - two rings for marriage";
    IPTCScene["012300"] = "off-beat: an attractive, perhaps fun picture of everyday events - dog with sunglasses, people cooling off in the fountain";
    IPTCScene["012400"] = "movie scene: photos taken during the shooting of a movie or TV production";

    /* Subject codes: http://cv.iptc.org/newscodes/subjectcode/    */
    IPTCSubject["01000000"] = "arts, culture and entertainment: Matters pertaining to the advancement and refinement of the human mind, of interests, skills, tastes and emotions";
    IPTCSubject["01001000"] = "archaeology: Probing the past through ruins and artefacts";
    IPTCSubject["01002000"] = "architecture: Designing of buildings, monuments and the spaces around them";
    IPTCSubject["01003000"] = "bullfighting: Classical contest pitting man against the bull";
    IPTCSubject["01004000"] = "festive event (including carnival): Parades, parties, celebrations and the like not necessarily tied to a fixed occasion or date";
    IPTCSubject["01005000"] = "cinema: Stories related to cinema as art and entertainment";
    IPTCSubject["01005001"] = "film festival: Stories about national and international motion pictures festivals, selections, festival juries, nominations, awards etc.";
    IPTCSubject["01006000"] = "dance: The expression of emotion or message through movement";
    IPTCSubject["01007000"] = "fashion: The design of clothing and accessories";
    IPTCSubject["01007001"] = "jewelry: Accessories to clothing";
    IPTCSubject["01008000"] = "language: The means by which people communicate with each other";
    IPTCSubject["01009000"] = "library and museum: Edifices used to house collections of books, music, art, or objects from the past and present for public use and display";
    IPTCSubject["01010000"] = "literature: The use of pamphlets, books or other printed matter to convey ideas, stories or other messages for the public";
    IPTCSubject["01010001"] = "fiction: Structured stories that are usually not based on fact but are the creation of the authors imagination";
    IPTCSubject["01010002"] = "poetry: The art, structure, forms of poetic expression";
    IPTCSubject["01011000"] = "music: Expressing emotion or message through instruments or voice using different sounds, tones, harmonies and the like";
    IPTCSubject["01011001"] = "classical music: Music that follows classic structures of rhythm and harmony";
    IPTCSubject["01011002"] = "folk music: Music that developed from folk cultures, often based on story-telling";
    IPTCSubject["01011003"] = "jazz music: A music of diverse harmonics, often improvised";
    IPTCSubject["01011004"] = "popular music: The latest fad in music, generally aimed at the younger generation";
    IPTCSubject["01011005"] = "country music: Similar to folk but is unique to the United States and is less about story telling than about loves sought and lost";
    IPTCSubject["01011006"] = "rock and roll music: Popular dance music developed in the 1950s";
    IPTCSubject["01012000"] = "painting: Using the mediums of oils, watercolour, pastel, pencils, chalk, crayon etc on various grounds to express emotion or message";
    IPTCSubject["01013000"] = "photography:	Mechanical means of creating images of objects by use of light and light sensitive materials with chemicals or by digitals means";
    IPTCSubject["01014000"] = "radio: Stories related to radio as art and entertainment";
    IPTCSubject["01015000"] = "sculpture: Representation of forms in clays, stone, woods, metals or other materials";
    IPTCSubject["01015001"] = "plastic art:	Forms of hand created art including installations";
    IPTCSubject["01016000"] = "television: Stories related to television as art and entertainment";
    IPTCSubject["01017000"] = "theatre: Telling of a story or idea through dialogue, music and physical expression in a space or building designed for it";
    IPTCSubject["01017001"] = "music theatre: Opera, operetta, music revues etc";
    IPTCSubject["01018000"] = "monument and heritage site: Areas containing commemorative objects for historical people or events";
    IPTCSubject["01019000"] = "customs and tradition: A particular way of behaving, or observances that have developed over time by a group of people";
    IPTCSubject["01020000"] = "arts (general): The collective expression of message or emotion through music, literature, painting, theatre or other means";
    IPTCSubject["01021000"] = "entertainment (general):	The collective use of television, radio, theatre, music and the like for the amusement of people";
    IPTCSubject["01022000"] = "culture (general): The ideas, customs, arts, skills of a particular group";
    IPTCSubject["01022001"] = "cultural development: The history of the development of art and culture such as the rise of cave paintings, pre-Colombian art, Chinese paper-making, anything non-political";
    IPTCSubject["01023000"] = "nightclub: A commercial establishment providing music, or other entertainment along with food and drink to selected clientele";
    IPTCSubject["01024000"] = "cartoon:	Still images such as editorial cartoons and comic strips";
    IPTCSubject["01025000"] = "animation: Animation, including full-length and short cinema, artists and merchandising of goods featuring animation characters.";
    IPTCSubject["01026000"] = "mass media: Television, radio, magazines, newspapers etc";
    IPTCSubject["01026001"] = "periodicals:	Written material that is usually published weekly, bi-weekly, monthly or annually for a long time";
    IPTCSubject["01026002"] = "news media: Television, wire services, radio that collect facts about incidents, developing and presenting them to audiences as a whole story";
    IPTCSubject["01026003"] = "newspapers: Daily or weekly publications that present the day to day history of the world, as well as features, comics etc";
    IPTCSubject["01026004"] = "reviews:	A critical look at someone else's work, whether film, theatre or writing";
    // ... others

    /* Genre codes http://cv.iptc.org/newscodes/genre/ */
    IPTCGenre["Actuality"] = "The object contains the recording of the event.";
    IPTCGenre["Advice"] = "The object contains advice, typically letters and answers about personal problems, that are publishable.";
    IPTCGenre["Almanac"] = "List of data, including birthdays of famous people and items of historical significance, for a given day";
    IPTCGenre["Analysis"] = "The object contains data and conclusions drawn by a journalist who has researched the story in depth.";
    IPTCGenre["Anniversary"] = "Stories about the anniversary of some important event that took place in recent history, usually bringing a short review of the event itself.";
    IPTCGenre["Archive material"] = "The object contains material distributed previously that has been selected from the originator's archives.";
    IPTCGenre["Background"] = "The object provides some scene setting and explanation for the event being reported.";
    IPTCGenre["Current"] = "The object content is about events taking place at the time of the report.";
    IPTCGenre["Curtain Raiser"] = "The object contains information about the staging and outcome of an immediately upcoming event.";
    IPTCGenre["Daybook"] = "Items filed on a regular basis that are lists of upcoming events with time and place, designed to inform others of events for planning purposes.";
    IPTCGenre["Exclusive"] = "Information content, in any form, that is unique to a specific information provider.";
    IPTCGenre["Feature"] = "The object content is about a particular event or individual that may not be significant to the current breaking news.";
    // ...other

    // Release status
    IPTCReleaseStatus["None"] = "No release is available";
    IPTCReleaseStatus["Not Applicable"] = "There are no items requiring a release";
    IPTCReleaseStatus["Unlimited Property Releases"] = "Releases are available for all";
    IPTCReleaseStatus["Incomplete Model Releases"] = "There are releases only for some";

    IPTCCopyrightStatus["Unknown"]="Unknown status";
    IPTCCopyrightStatus["Protected"] = "Resource is protected by copyright";
    IPTCCopyrightStatus["Public Domain"] ="Public Domain resource";

    // World regions/continent
    IPTCWorldRegion["r001"] = "World: The whole world";
    IPTCWorldRegion["r002"] = "Africa: 	The continent of Africa";
    IPTCWorldRegion["r005"] = "South America: The southern region of the Americas";
    IPTCWorldRegion["r009"] = "Oceania: The continent of Oceania";
    IPTCWorldRegion["r021"] = "North America: The northern part of the Americas";
    IPTCWorldRegion["r142"] = "Asia: The continent of Asia";
    IPTCWorldRegion["r150"] = "Europe: The continent of Europe";
    IPTCWorldRegion["r901"] = "Antarctica: The continent of the Antarctica";

    IPTCBoolean["True"]="";
    IPTCBoolean["False"]="";

    IPTCISO3166["ABW"] = "Aruba";
    IPTCISO3166["AFG"] = "Afghanistan";
    IPTCISO3166["AGO"] = "Angola";
    IPTCISO3166["AIA"] = "Anguilla";
    IPTCISO3166["ALA"] = "Aland Islands"; //TODO
    IPTCISO3166["ALB"] = "Albania";
    IPTCISO3166["AND"] = "Andorra";
    IPTCISO3166["ARE"] = "United Arab Emirates";
    IPTCISO3166["ARG"] = "Argentina";
    IPTCISO3166["ARM"] = "Armenia";
    IPTCISO3166["ASM"] = "American Samoa";
    IPTCISO3166["ATA"] = "Antarctica";
    IPTCISO3166["ATF"] = "French Southern Territories";
    IPTCISO3166["ATG"] = "Antigua and Barbuda";
    IPTCISO3166["AUS"] = "Australia";
    IPTCISO3166["AUT"] = "Austria";
    IPTCISO3166["AZE"] = "Azerbaijan";
    IPTCISO3166["BDI"] = "Burundi";
    IPTCISO3166["BEL"] = "Belgium";
    IPTCISO3166["BEN"] = "Benin";
    IPTCISO3166["BES"] = "Bonaire, Sint Eustatius and Saba";
    IPTCISO3166["BFA"] = "Burkina Faso";
    IPTCISO3166["BGD"] = "Bangladesh";
    IPTCISO3166["BGR"] = "Bulgaria";
    IPTCISO3166["BHR"] = "Bahrain";
    IPTCISO3166["BHS"] = "Bahamas";
    IPTCISO3166["BIH"] = "Bosnia and Herzegovina";
    IPTCISO3166["BLM"] = "Saint Barthelemy"; //TODO
    IPTCISO3166["BLR"] = "Belarus";
    IPTCISO3166["BLZ"] = "Belize";
    IPTCISO3166["BMU"] = "Bermuda";
    IPTCISO3166["BOL"] = "Bolivia, Plurinational State of";
    IPTCISO3166["BRA"] = "Brazil";
    IPTCISO3166["BRB"] = "Barbados";
    IPTCISO3166["BRN"] = "Brunei Darussalam";
    IPTCISO3166["BTN"] = "Bhutan";
    IPTCISO3166["BVT"] = "Bouvet Island";
    IPTCISO3166["BWA"] = "Botswana";
    IPTCISO3166["CAF"] = "Central African Republic";
    IPTCISO3166["CAN"] = "Canada";
    IPTCISO3166["CCK"] = "Cocos (Keeling) Islands";
    IPTCISO3166["CHE"] = "Switzerland";
    IPTCISO3166["CHL"] = "Chile";
    IPTCISO3166["CHN"] = "China";
    IPTCISO3166["CIV"] = "Cote d'Ivoire"; //TODO
    IPTCISO3166["CMR"] = "Cameroon";
    IPTCISO3166["COD"] = "Congo, the Democratic Republic of the";
    IPTCISO3166["COG"] = "Congo";
    IPTCISO3166["COK"] = "Cook Islands";
    IPTCISO3166["COL"] = "Colombia";
    IPTCISO3166["COM"] = "Comoros";
    IPTCISO3166["CPV"] = "Cape Verde";
    IPTCISO3166["CRI"] = "Costa Rica";
    IPTCISO3166["CUB"] = "Cuba";
    IPTCISO3166["CUW"] = "Curacao"; //TODO
    IPTCISO3166["CXR"] = "Christmas Island";
    IPTCISO3166["CYM"] = "Cayman Islands";
    IPTCISO3166["CYP"] = "Cyprus";
    IPTCISO3166["CZE"] = "Czech Republic";
    IPTCISO3166["DEU"] = "Germany";
    IPTCISO3166["DJI"] = "Djibouti";
    IPTCISO3166["DMA"] = "Dominica";
    IPTCISO3166["DNK"] = "Denmark";
    IPTCISO3166["DOM"] = "Dominican Republic";
    IPTCISO3166["DZA"] = "Algeria";
    IPTCISO3166["ECU"] = "Ecuador";
    IPTCISO3166["EGY"] = "Egypt";
    IPTCISO3166["ERI"] = "Eritrea";
    IPTCISO3166["ESH"] = "Western Sahara";
    IPTCISO3166["ESP"] = "Spain";
    IPTCISO3166["EST"] = "Estonia";
    IPTCISO3166["ETH"] = "Ethiopia";
    IPTCISO3166["FIN"] = "Finland";
    IPTCISO3166["FJI"] = "Fiji";
    IPTCISO3166["FLK"] = "Falkland Islands (Malvinas)";
    IPTCISO3166["FRA"] = "France";
    IPTCISO3166["FRO"] = "Faer Oer Islands"; //TODO
    IPTCISO3166["FSM"] = "Micronesia, Federated States of";
    IPTCISO3166["GAB"] = "Gabon";
    IPTCISO3166["GBR"] = "United Kingdom";
    IPTCISO3166["GEO"] = "Georgia";
    IPTCISO3166["GGY"] = "Guernsey";
    IPTCISO3166["GHA"] = "Ghana";
    IPTCISO3166["GIB"] = "Gibraltar";
    IPTCISO3166["GIN"] = "Guinea";
    IPTCISO3166["GLP"] = "Guadeloupe";
    IPTCISO3166["GMB"] = "Gambia";
    IPTCISO3166["GNB"] = "Guinea-Bissau";
    IPTCISO3166["GNQ"] = "Equatorial Guinea";
    IPTCISO3166["GRC"] = "Greece";
    IPTCISO3166["GRD"] = "Grenada";
    IPTCISO3166["GRL"] = "Greenland";
    IPTCISO3166["GTM"] = "Guatemala";
    IPTCISO3166["GUF"] = "French Guiana";
    IPTCISO3166["GUM"] = "Guam";
    IPTCISO3166["GUY"] = "Guyana";
    IPTCISO3166["HKG"] = "Hong Kong";
    IPTCISO3166["HMD"] = "Heard Island and McDonald Islands";
    IPTCISO3166["HND"] = "Honduras";
    IPTCISO3166["HRV"] = "Croatia";
    IPTCISO3166["HTI"] = "Haiti";
    IPTCISO3166["HUN"] = "Hungary";
    IPTCISO3166["IDN"] = "Indonesia";
    IPTCISO3166["IMN"] = "Isle of Man";
    IPTCISO3166["IND"] = "India";
    IPTCISO3166["IOT"] = "British Indian Ocean Territory";
    IPTCISO3166["IRL"] = "Ireland";
    IPTCISO3166["IRN"] = "Iran, Islamic Republic of";
    IPTCISO3166["IRQ"] = "Iraq";
    IPTCISO3166["ISL"] = "Iceland";
    IPTCISO3166["ISR"] = "Israel";
    IPTCISO3166["ITA"] = "Italy";
    IPTCISO3166["JAM"] = "Jamaica";
    IPTCISO3166["JEY"] = "Jersey";
    IPTCISO3166["JOR"] = "Jordan";
    IPTCISO3166["JPN"] = "Japan";
    IPTCISO3166["KAZ"] = "Kazakhstan";
    IPTCISO3166["KEN"] = "Kenya";
    IPTCISO3166["KGZ"] = "Kyrgyzstan";
    IPTCISO3166["KHM"] = "Cambodia";
    IPTCISO3166["KIR"] = "Kiribati";
    IPTCISO3166["KNA"] = "Saint Kitts and Nevis";
    IPTCISO3166["KOR"] = "Korea, Republic of";
    IPTCISO3166["KWT"] = "Kuwait";
    IPTCISO3166["LAO"] = "Lao People's Democratic Republic";
    IPTCISO3166["LBN"] = "Lebanon";
    IPTCISO3166["LBR"] = "Liberia";
    IPTCISO3166["LBY"] = "Libyan Arab Jamahiriya";
    IPTCISO3166["LCA"] = "Saint Lucia";
    IPTCISO3166["LIE"] = "Liechtenstein";
    IPTCISO3166["LKA"] = "Sri Lanka";
    IPTCISO3166["LSO"] = "Lesotho";
    IPTCISO3166["LTU"] = "Lithuania";
    IPTCISO3166["LUX"] = "Luxembourg";
    IPTCISO3166["LVA"] = "Latvia";
    IPTCISO3166["MAC"] = "Macao";
    IPTCISO3166["MAF"] = "Saint Martin (French part)";
    IPTCISO3166["MAR"] = "Morocco";
    IPTCISO3166["MCO"] = "Monaco";
    IPTCISO3166["MDA"] = "Moldova, Republic of";
    IPTCISO3166["MDG"] = "Madagascar";
    IPTCISO3166["MDV"] = "Maldives";
    IPTCISO3166["MEX"] = "Mexico";
    IPTCISO3166["MHL"] = "Marshall Islands";
    IPTCISO3166["MKD"] = "Macedonia, the former Yugoslav Republic of";
    IPTCISO3166["MLI"] = "Mali";
    IPTCISO3166["MLT"] = "Malta";
    IPTCISO3166["MMR"] = "Myanmar";
    IPTCISO3166["MNE"] = "Montenegro";
    IPTCISO3166["MNG"] = "Mongolia";
    IPTCISO3166["MNP"] = "Northern Mariana Islands";
    IPTCISO3166["MOZ"] = "Mozambique";
    IPTCISO3166["MRT"] = "Mauritania";
    IPTCISO3166["MSR"] = "Montserrat";
    IPTCISO3166["MTQ"] = "Martinique";
    IPTCISO3166["MUS"] = "Mauritius";
    IPTCISO3166["MWI"] = "Malawi";
    IPTCISO3166["MYS"] = "Malaysia";
    IPTCISO3166["MYT"] = "Mayotte";
    IPTCISO3166["NAM"] = "Namibia";
    IPTCISO3166["NCL"] = "New Caledonia";
    IPTCISO3166["NER"] = "Niger";
    IPTCISO3166["NFK"] = "Norfolk Island";
    IPTCISO3166["NGA"] = "Nigeria";
    IPTCISO3166["NIC"] = "Nicaragua";
    IPTCISO3166["NIU"] = "Niue";
    IPTCISO3166["NLD"] = "Netherlands";
    IPTCISO3166["NOR"] = "Norway";
    IPTCISO3166["NPL"] = "Nepal";
    IPTCISO3166["NRU"] = "Nauru";
    IPTCISO3166["NZL"] = "New Zealand";
    IPTCISO3166["OMN"] = "Oman";
    IPTCISO3166["PAK"] = "Pakistan";
    IPTCISO3166["PAN"] = "Panama";
    IPTCISO3166["PCN"] = "Pitcairn";
    IPTCISO3166["PER"] = "Peru";
    IPTCISO3166["PHL"] = "Philippines";
    IPTCISO3166["PLW"] = "Palau";
    IPTCISO3166["PNG"] = "Papua New Guinea";
    IPTCISO3166["POL"] = "Poland";
    IPTCISO3166["PRI"] = "Puerto Rico";
    IPTCISO3166["PRK"] = "Korea, Democratic People's Republic of";
    IPTCISO3166["PRT"] = "Portugal";
    IPTCISO3166["PRY"] = "Paraguay";
    IPTCISO3166["PSE"] = "Palestinian Territory, Occupied";
    IPTCISO3166["PYF"] = "French Polynesia";
    IPTCISO3166["QAT"] = "Qatar";
    IPTCISO3166["REU"] = "Reunion"; //TODO
    IPTCISO3166["ROU"] = "Romania";
    IPTCISO3166["RUS"] = "Russian Federation";
    IPTCISO3166["RWA"] = "Rwanda";
    IPTCISO3166["SAU"] = "Saudi Arabia";
    IPTCISO3166["SDN"] = "Sudan";
    IPTCISO3166["SEN"] = "Senegal";
    IPTCISO3166["SGP"] = "Singapore";
    IPTCISO3166["SGS"] = "South Georgia and the South Sandwich Islands";
    IPTCISO3166["SHN"] = "Saint Helena, Ascension and Tristan da Cunha";
    IPTCISO3166["SJM"] = "Svalbard and Jan Mayen";
    IPTCISO3166["SLB"] = "Solomon Islands";
    IPTCISO3166["SLE"] = "Sierra Leone";
    IPTCISO3166["SLV"] = "El Salvador";
    IPTCISO3166["SMR"] = "San Marino";
    IPTCISO3166["SOM"] = "Somalia";
    IPTCISO3166["SPM"] = "Saint Pierre and Miquelon";
    IPTCISO3166["SRB"] = "Serbia";
    IPTCISO3166["SSD"] = "South Sudan";
    IPTCISO3166["STP"] = "Sao Tome and Principe";
    IPTCISO3166["SUR"] = "Suriname";
    IPTCISO3166["SVK"] = "Slovakia";
    IPTCISO3166["SVN"] = "Slovenia";
    IPTCISO3166["SWE"] = "Sweden";
    IPTCISO3166["SWZ"] = "Swaziland";
    IPTCISO3166["SXM"] = "Sint Maarten (Dutch part)";
    IPTCISO3166["SYC"] = "Seychelles";
    IPTCISO3166["SYR"] = "Syrian Arab Republic";
    IPTCISO3166["TCA"] = "Turks and Caicos Islands";
    IPTCISO3166["TCD"] = "Chad";
    IPTCISO3166["TGO"] = "Togo";
    IPTCISO3166["THA"] = "Thailand";
    IPTCISO3166["TJK"] = "Tajikistan";
    IPTCISO3166["TKL"] = "Tokelau";
    IPTCISO3166["TKM"] = "Turkmenistan";
    IPTCISO3166["TLS"] = "Timor-Leste";
    IPTCISO3166["TON"] = "Tonga";
    IPTCISO3166["TTO"] = "Trinidad and Tobago";
    IPTCISO3166["TUN"] = "Tunisia";
    IPTCISO3166["TUR"] = "Turkey";
    IPTCISO3166["TUV"] = "Tuvalu";
    IPTCISO3166["TWN"] = "Taiwan, Province of China";
    IPTCISO3166["TZA"] = "Tanzania, United Republic of";
    IPTCISO3166["UGA"] = "Uganda";
    IPTCISO3166["UKR"] = "Ukraine";
    IPTCISO3166["UMI"] = "United States Minor Outlying Islands";
    IPTCISO3166["URY"] = "Uruguay";
    IPTCISO3166["USA"] = "United States";
    IPTCISO3166["UZB"] = "Uzbekistan";
    IPTCISO3166["VAT"] = "Holy See (Vatican City State)";
    IPTCISO3166["VCT"] = "Saint Vincent and the Grenadines";
    IPTCISO3166["VEN"] = "Venezuela, Bolivarian Republic of";
    IPTCISO3166["VGB"] = "Virgin Islands, British";
    IPTCISO3166["VIR"] = "Virgin Islands, U.S.";
    IPTCISO3166["VNM"] = "Viet Nam";
    IPTCISO3166["VUT"] = "Vanuatu";
    IPTCISO3166["WLF"] = "Wallis and Futuna";
    IPTCISO3166["WSM"] = "Samoa";
    IPTCISO3166["YEM"] = "Yemen";
    IPTCISO3166["ZAF"] = "South Africa";
    IPTCISO3166["ZMB"] = "Zambia";
    IPTCISO3166["ZWE"] = "Zimbabwe";
}

}
