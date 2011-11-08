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
IPTCPairList_t IPTCMeta::IPTCWorldRegion;
IPTCPairList_t IPTCMeta::IPTCISO3166;

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
const char *kIPTCRightsOwner    = "plus:CopyrightOwner/plus:CopyrightOwnerName";
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

	IPTCtags[kIPTCHeadline] = IPTCMeta(kIPTCHeadline, "Headline", "A brief synopsis of the caption.");
	IPTCtags[kIPTCCity]     = IPTCMeta(kIPTCCity,     "City", "Name of the city of the location shown in the image. IIM 2:90");
	IPTCtags[kIPTCCountry]  = IPTCMeta(kIPTCCountry,  "Country", "Full name of the country of the location shown in the image. IIM 2:101");
    IPTCtags[kIPTCState]    = IPTCMeta(kIPTCState,    "Province/State", "Name of the subregion of a country of the location shown in the image. IIM 2:95");
    IPTCtags[kIPTCAuthorPos]= IPTCMeta(kIPTCAuthorPos,"Creator's job title","The job title of the photographer IIM 2:85 By line title");

    IPTCtags[kIPTCWriter]   = IPTCMeta(kIPTCWriter, "Description Writer","The name of the person involved in writing, editing or correcting the description of the image IIM 2:122");
    IPTCtags[kIPTCInstruct] = IPTCMeta(kIPTCInstruct, "Instructions", "information about embargoes, or other restrictions not covered by the Rights Usage IIM 2:40 ");
    IPTCtags[kIPTCReference]= IPTCMeta(kIPTCReference,"Job Id","A number or identifier needed for workflow control or tracking IIM 2:103");
    IPTCtags[kIPTCCredit]   = IPTCMeta(kIPTCCredit, "Credit Line" ,"The credit to person(s) and/or organisation(s) required by the supplier of the image to be used when published IIM 2:110");
    IPTCtags[kIPTCSource]   = IPTCMeta(kIPTCSource, "Source", "Identifies the original owner of the copyright for the intellectual content of the image. This could be an agency, a member of an agency or an individual. Source could be different from Creator and from the entities in the CopyrightNotice. IIM 2:115");
    IPTCtags[kIPTCUrgency]  = IPTCMeta(kIPTCUrgency, "Urgency", "IIM  2:09 Urgency (deprecated)");
    IPTCtags[kIPTCCategory] = IPTCMeta(kIPTCCategory,"Category", "IIM 2:15 Category (deprecated)");
    IPTCtags[kIPTCSuppCateg]= IPTCMeta(kIPTCSuppCateg,"Supplemental Categories"," IIM 2:20 Supplemental Category (deprecated)");

    IPTCtags[kIPTCLocation]     = IPTCMeta(kIPTCLocation, "Sublocation", "Exact name of the sublocation shown in the image. This sublocation name could either be the name of a sublocation to a city or the name of a well known location or (natural) monument outside a city. IIM 2:92");
    IPTCtags[kIPTCCountryCode]  = IPTCMeta(kIPTCCountryCode, "Country Code","The code should be taken from ISO 3166 two or three letter code. The full name of a country should go to the 'Country' element IIM 2:100");
    IPTCtags[kIPTCGenre]        = IPTCMeta(kIPTCGenre,"Intellectual Genre" ,"Describe the nature of the image in terms of its intellectual or journalistic characteristics, such as daybook, or feature (examples at: http://www.newscodes.org/) IIM 2:04 ");
    IPTCtags[kIPTCScene]        = IPTCMeta(kIPTCScene, "Scene code", "Describes the scene of a photo content. Specifies one ore more terms from the IPTC 'Scene-NewsCodes' (http://www.newscodes.org/). Each Scene is represented as a string of 6 digits in an unordered list.", Exiv2::xmpBag );
    IPTCtags[kIPTCSubjCode]     = IPTCMeta(kIPTCSubjCode,"Subject Code","Specifies one or more Subjects from the IPTC 'Subject-NewsCodes' (http://www.newscodes.org/) taxonomy to categorise the image. Each Subject is represented as a string of 8 digits in an unordered list.", Exiv2::xmpBag );

    IPTCtags[kIPTCKeywords]     = IPTCMeta(kIPTCKeywords,"Keywords","Keywords to express the subject of the image. Any number of keywords, terms or phrases used to express the subject matter in the image. IIM 2:25", Exiv2::xmpBag );
    IPTCtags[kIPTCCreator]      = IPTCMeta(kIPTCCreator,"Author","Contains the name of the photographer, but in cases where the photographer should not be identified the name of a company or organisation may be appropriate. IIM 2:80 By line", Exiv2::langAlt );
    IPTCtags[kIPTCTitle]        = IPTCMeta(kIPTCTitle, "Title","A shorthand reference for the digital image; this may be the file name. IIM 2:05 Object Name", Exiv2::langAlt );
    IPTCtags[kIPTCRights]       = IPTCMeta(kIPTCRights,"Copyright Notice", "Contains any necessary copyright notice for claiming the intellectual property for this photograph and should identify the current owner of the copyright for the photograph IIM 2:116 Copyright Notice", Exiv2::langAlt );
    IPTCtags[kIPTCDescription]  = IPTCMeta(kIPTCDescription, "Description", "Enter a 'caption' describing the who, what, and why of what is happening in this image, this might include names of people, and/or their role in the action that is taking place within the image. IIM 2:120 Caption/Abstract", Exiv2::langAlt );

    IPTCtags[kIPTCUsageTerms]   = IPTCMeta(kIPTCUsageTerms,"Usage terms","Instructions on how this image can legally be used",Exiv2::langAlt);

    IPTCtags[kIPTCCreatorAdrCity] = IPTCMeta(kIPTCCreatorAdrCity,"Creator city","The city for the address of the person that created this image" );
    IPTCtags[kIPTCCreatorAdrCtry] = IPTCMeta(kIPTCCreatorAdrCtry,"Creator country","The country name for the address of the person that created this image");
    IPTCtags[kIPTCCreatorExtadr]  = IPTCMeta(kIPTCCreatorExtadr,"Creator address","The address for the person that created this image");
    IPTCtags[kIPTCCreatorPcode]   = IPTCMeta(kIPTCCreatorPcode,"Creator postal code","The postal code for the address of the person that created this image");
    IPTCtags[kIPTCCreatorRegion]  = IPTCMeta(kIPTCCreatorRegion,"Creator State/Province","The state or province for the address of the person that created this image");
    IPTCtags[kIPTCCreatorEmail]   = IPTCMeta(kIPTCCreatorEmail,"Creator Email(s)","The work Email address(es) for the person that created this image");
    IPTCtags[kIPTCCreatorTel]     = IPTCMeta(kIPTCCreatorTel,"Creator Phone num.","The work Phone number(s) for the person that created this image, using the international format");
    IPTCtags[kIPTCCreatorUrl]     = IPTCMeta(kIPTCCreatorUrl,"Creator Web URL","The work Web URL(s) for the person that created this image, such as http://www.domain.com/");

    IPTCtags[kIPTCModelInfo]      = IPTCMeta(kIPTCModelInfo,"Model information","Information like ethnicity or other details about the model(s) in this image");
    IPTCtags[kIPTCModelAge]       = IPTCMeta(kIPTCModelAge,"Model(s) Age","The age of the human model(s) at the time this image was made",Exiv2::xmpBag);
    IPTCtags[kIPTCPerson]         = IPTCMeta(kIPTCPerson,"Person shown","Name of a person shown in the image",Exiv2::xmpBag);
    IPTCtags[kIPTCOrganization]   = IPTCMeta(kIPTCOrganization,"Code of Organization","Enter an identifier for the controlled vocabulary, then a colon, and finally the code from the vocabulary assigned to the organisation shown in this image (e.g. nasdaq:companyA)");
    IPTCtags[kIPTCCVterm]         = IPTCMeta(kIPTCCVterm,"Controlled Vocabulary Term","Term to describe the content of the image by a value from a Controlled Vocabulary: an identifier for the controlled vocabulary, then a colon, and finally the code from the vocabulary assigned to the term",Exiv2::xmpBag );
    IPTCtags[kIPTCEvent]          = IPTCMeta(kIPTCEvent,"Event","the name or description of the event where this image was taken",Exiv2::langAlt);

    IPTCtags[kIPTCLocationCity]   = IPTCMeta(kIPTCLocationCity,"Location City","Details about the City which is shown in this image");
    IPTCtags[kIPTCLocationCode]   = IPTCMeta(kIPTCLocationCode,"Location Country code","Enter the 2 or 3 letter ISO 3166 Country Code of the Country");
    IPTCtags[kIPTCLocationCtry]   = IPTCMeta(kIPTCLocationCtry,"Location Country","Details about the Country which is shown in this image");
    IPTCtags[kIPTCLocationState]  = IPTCMeta(kIPTCLocationState,"Location State/Province","Details about the State/Province which is shown in this image");
    IPTCtags[kIPTCLocationSubloc] = IPTCMeta(kIPTCLocationSubloc,"Location Sublocation","Name of a sublocation. This sublocation name could either be the name of a sublocation to a city or the name of a well known location or (natural) monument outside a city. In the sense of a sublocation to a city this element is at the fifth level of a top-down geographical hierarchy.");
    IPTCtags[kIPTCLocationRegion] = IPTCMeta(kIPTCLocationRegion,"Location World region","The name of a world region of a location. This element is at the first (top I) level of a top-down geographical hierarchy.");

    IPTCtags[kIPTCLocCreateCity]  = IPTCMeta(kIPTCLocCreateCity,"Created City","The location the content of the item was created: details about a location where this image was created");
    IPTCtags[kIPTCLocCreateCode]  = IPTCMeta(kIPTCLocCreateCode,"Created Country code","The location the content of the item was created: details about a location where this image was created");
    IPTCtags[kIPTCLocCreateCtry]  = IPTCMeta(kIPTCLocCreateCtry,"Created Country","The location the content of the item was created: details about a location where this image was created");
    IPTCtags[kIPTCLocCreateState] = IPTCMeta(kIPTCLocCreateState,"Created State/Province","The location the content of the item was created: details about a location where this image was created");
    IPTCtags[kIPTCLocCreateSubloc]= IPTCMeta(kIPTCLocCreateSubloc,"Created Sublocation","The location the content of the item was created: details about a location where this image was created");
    IPTCtags[kIPTCLocCreateRegion]= IPTCMeta(kIPTCLocCreateRegion,"Created World region","The location the content of the item was created: details about a location where this image was created");

    IPTCtags[kIPTCArtworkRights]  = IPTCMeta(kIPTCArtworkRights,"Artwork copyrights","Contains any necessary copyright notice for claiming the intellectual property for artwork or an object in the image and should identify the current owner of the copyright of this work with associated intellectual property rights");
    IPTCtags[kIPTCArtworkCreator] = IPTCMeta(kIPTCArtworkCreator,"Artwork creator","Contains the name of the artist who has created artwork or an object in the image. In cases where the artist could or should not be identified the name of a company or organisation may be appropriate.");
    IPTCtags[kIPTCArtworkDate]    = IPTCMeta(kIPTCArtworkDate,"Artwork creation date","Designates the date and optionally the time the artwork or object in the image was created. This relates to artwork or objects with associated intellectual property rights");
    IPTCtags[kIPTCArtworkSource]  = IPTCMeta(kIPTCArtworkSource,"Artwork source","The organisation or body holding and registering the artwork or object in the image for inventory purposes.");
    IPTCtags[kIPTCArtworkNumber]  = IPTCMeta(kIPTCArtworkNumber,"Source Inventory Number","The inventory number issued by the organisation or body holding and registering the artwork or object in the image.");
    IPTCtags[kIPTCArtworkTitle]   = IPTCMeta(kIPTCArtworkTitle,"Artwork title","The verbal and human readable name of the artwork or object in this image");

    IPTCtags[kIPTCRightsOwner]    = IPTCMeta(kIPTCRightsOwner,"Copyright Owner","The owner or owners of the copyright in the licensed image",Exiv2::xmpSeq);
    IPTCtags[kIPTCImageCreator]   = IPTCMeta(kIPTCImageCreator,"Image Creator","Creator or creators of the image");
    IPTCtags[kIPTCLicensor]       = IPTCMeta(kIPTCLicensor,"Licensor","A person or company that should be contacted to obtain a licence for using the item or who has licensed the item");
    IPTCtags[kIPTCImageSupplier]  = IPTCMeta(kIPTCImageSupplier,"Image supplier","The identifier for the most recent supplier of this image - note that this might not be the creator or owner of the image");
    IPTCtags[kIPTCGUIDSupplier]   = IPTCMeta(kIPTCGUIDSupplier,"Image supplier GUID","Optional identifier assigned by the Image Supplier to the image");
    IPTCtags[kIPTCMinorDisclosure]= IPTCMeta(kIPTCMinorDisclosure,"Minor Model Age Disclosure","The age of the youngest model pictured in this image, at the time that this image was made.");
    IPTCtags[kIPTCModelReleaseID] = IPTCMeta(kIPTCModelReleaseID,"Model Release ID","Identifier associated with each Model Release", Exiv2::xmpBag);
    IPTCtags[kIPTCModelReleaseSt] = IPTCMeta(kIPTCModelReleaseSt,"Model Release Status","The availability and scope of model releases authorising usage of the likenesses of persons appearing in the photograph.");
    IPTCtags[kIPTCPropertyRelID]  = IPTCMeta(kIPTCPropertyRelID,"Property Release Id","Optional identifier associated with each Property Release", Exiv2::xmpBag);
    IPTCtags[kIPTCPropertyRelSt]  = IPTCMeta(kIPTCPropertyRelSt,"Property Release Status","The availability and scope of property releases authorising usage of the properties appearing in the photograph.");
    IPTCtags[kIPTCRegistryID]     = IPTCMeta(kIPTCRegistryID,"Registry Image Id","A unique identifier created by a registry and applied by the creator of the digital image. This value shall not be changed after being applied.");
    IPTCtags[kIPTCRegistryOrgID]  = IPTCMeta(kIPTCRegistryOrgID,"Registry Organisation Id","An identifier for the registry which issued the corresponding Registry Image Id");

    // readonly
    IPTCtags[kIPTCDate]           = IPTCMeta(kIPTCDate,     "Date Created" ,"The Date the image was taken IIM 2:55");
    IPTCtags[kIPTCGUID]           = IPTCMeta(kIPTCGUID,"Image GUID","Globally unique identifier for this digital image");
    IPTCtags[kIPTCMaxHeight]      = IPTCMeta(kIPTCMaxHeight,"Max Avail Height","The maximum available height in pixels of the original photo from which this photo has been derived by downsizing.");
    IPTCtags[kIPTCMaxWidth]       = IPTCMeta(kIPTCMaxWidth,"Max Avail Width","The maximum available width in pixels of the original photo from which this photo has been derived by downsizing.");
    // PLUS versione is 1.2.0
    IPTCtags[kIPTCPlusVersion]    = IPTCMeta(kIPTCPlusVersion,"PLUS Version","The version number of the PLUS standards in place at the time of the transaction");

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

    // World regions/continent
    IPTCWorldRegion["r001"] = "World: The whole world";
    IPTCWorldRegion["r002"] = "Africa: 	The continent of Africa";
    IPTCWorldRegion["r005"] = "South America: The southern region of the Americas";
    IPTCWorldRegion["r009"] = "Oceania: The continent of Oceania";
    IPTCWorldRegion["r021"] = "North America: The northern part of the Americas";
    IPTCWorldRegion["r142"] = "Asia: The continent of Asia";
    IPTCWorldRegion["r150"] = "Europe: The continent of Europe";
    IPTCWorldRegion["r901"] = "Antarctica: The continent of the Antarctica";

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
