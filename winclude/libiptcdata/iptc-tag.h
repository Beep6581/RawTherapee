/* iptc-tag.h
 *
 * Copyright © 2001 David Moore <dcm@acm.org>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details. 
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __IPTC_TAG_H__
#define __IPTC_TAG_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


typedef enum {
	IPTC_RECORD_OBJECT_ENV		= 1,
	IPTC_RECORD_APP_2		= 2,
	IPTC_RECORD_APP_3		= 3,
	IPTC_RECORD_APP_4		= 4,
	IPTC_RECORD_APP_5		= 5,
	IPTC_RECORD_APP_6		= 6,
	IPTC_RECORD_PREOBJ_DATA		= 7,
	IPTC_RECORD_OBJ_DATA		= 8,
	IPTC_RECORD_POSTOBJ_DATA	= 9
} IptcRecord;

typedef enum {
	IPTC_TAG_MODEL_VERSION		= 0,	/* Begin record 1 tags */
	IPTC_TAG_DESTINATION		= 5,
	IPTC_TAG_FILE_FORMAT		= 20,
	IPTC_TAG_FILE_VERSION		= 22,
	IPTC_TAG_SERVICE_ID		= 30,
	IPTC_TAG_ENVELOPE_NUM		= 40,
	IPTC_TAG_PRODUCT_ID		= 50,
	IPTC_TAG_ENVELOPE_PRIORITY	= 60,
	IPTC_TAG_DATE_SENT		= 70,
	IPTC_TAG_TIME_SENT		= 80,
	IPTC_TAG_CHARACTER_SET		= 90,
	IPTC_TAG_UNO			= 100,
	IPTC_TAG_ARM_ID			= 120,
	IPTC_TAG_ARM_VERSION		= 122,	/* End record 1 tags */
	IPTC_TAG_RECORD_VERSION		= 0,	/* Begin record 2 tags */
	IPTC_TAG_OBJECT_TYPE		= 3,
	IPTC_TAG_OBJECT_ATTRIBUTE	= 4,
	IPTC_TAG_OBJECT_NAME		= 5,
	IPTC_TAG_EDIT_STATUS		= 7,
	IPTC_TAG_EDITORIAL_UPDATE	= 8,
	IPTC_TAG_URGENCY		= 10,
	IPTC_TAG_SUBJECT_REFERENCE	= 12,
	IPTC_TAG_CATEGORY		= 15,
	IPTC_TAG_SUPPL_CATEGORY		= 20,
	IPTC_TAG_FIXTURE_ID		= 22,
	IPTC_TAG_KEYWORDS		= 25,
	IPTC_TAG_CONTENT_LOC_CODE	= 26,
	IPTC_TAG_CONTENT_LOC_NAME	= 27,
	IPTC_TAG_RELEASE_DATE		= 30,
	IPTC_TAG_RELEASE_TIME		= 35,
	IPTC_TAG_EXPIRATION_DATE	= 37,
	IPTC_TAG_EXPIRATION_TIME	= 38,
	IPTC_TAG_SPECIAL_INSTRUCTIONS	= 40,
	IPTC_TAG_ACTION_ADVISED		= 42,
	IPTC_TAG_REFERENCE_SERVICE	= 45,
	IPTC_TAG_REFERENCE_DATE		= 47,
	IPTC_TAG_REFERENCE_NUMBER	= 50,
	IPTC_TAG_DATE_CREATED		= 55,
	IPTC_TAG_TIME_CREATED		= 60,
	IPTC_TAG_DIGITAL_CREATION_DATE	= 62,
	IPTC_TAG_DIGITAL_CREATION_TIME	= 63,
	IPTC_TAG_ORIGINATING_PROGRAM	= 65,
	IPTC_TAG_PROGRAM_VERSION	= 70,
	IPTC_TAG_OBJECT_CYCLE		= 75,
	IPTC_TAG_BYLINE			= 80,
	IPTC_TAG_BYLINE_TITLE		= 85,
	IPTC_TAG_CITY			= 90,
	IPTC_TAG_SUBLOCATION		= 92,
	IPTC_TAG_STATE			= 95,
	IPTC_TAG_COUNTRY_CODE		= 100,
	IPTC_TAG_COUNTRY_NAME		= 101,
	IPTC_TAG_ORIG_TRANS_REF		= 103,
	IPTC_TAG_HEADLINE		= 105,
	IPTC_TAG_CREDIT			= 110,
	IPTC_TAG_SOURCE			= 115,
	IPTC_TAG_COPYRIGHT_NOTICE	= 116,
	IPTC_TAG_PICASA_UNKNOWN		= 117,
	IPTC_TAG_CONTACT		= 118,
	IPTC_TAG_CAPTION		= 120,
	IPTC_TAG_WRITER_EDITOR		= 122,
	IPTC_TAG_RASTERIZED_CAPTION	= 125,
	IPTC_TAG_IMAGE_TYPE		= 130,
	IPTC_TAG_IMAGE_ORIENTATION	= 131,
	IPTC_TAG_LANGUAGE_ID		= 135,
	IPTC_TAG_AUDIO_TYPE		= 150,
	IPTC_TAG_AUDIO_SAMPLING_RATE	= 151,
	IPTC_TAG_AUDIO_SAMPLING_RES	= 152,
	IPTC_TAG_AUDIO_DURATION		= 153,
	IPTC_TAG_AUDIO_OUTCUE		= 154,
	IPTC_TAG_PREVIEW_FORMAT		= 200,
	IPTC_TAG_PREVIEW_FORMAT_VER	= 201,
	IPTC_TAG_PREVIEW_DATA		= 202,	/* End record 2 tags */
	IPTC_TAG_SIZE_MODE		= 10,	/* Begin record 7 tags */
	IPTC_TAG_MAX_SUBFILE_SIZE	= 20,
	IPTC_TAG_SIZE_ANNOUNCED		= 90,
	IPTC_TAG_MAX_OBJECT_SIZE	= 95,	/* End record 7 tags */
	IPTC_TAG_SUBFILE		= 10,	/* Record 8 tags */
	IPTC_TAG_CONFIRMED_DATA_SIZE	= 10	/* Record 9 tags */
} IptcTag;

typedef enum {
	IPTC_OPTIONAL = 0,
	IPTC_MANDATORY = 1
} IptcMandatory;

typedef enum {
	IPTC_NOT_REPEATABLE = 0,
	IPTC_REPEATABLE = 1
} IptcRepeatable;

typedef enum {
	IPTC_FORMAT_UNKNOWN,
	IPTC_FORMAT_BINARY,
	IPTC_FORMAT_BYTE,
	IPTC_FORMAT_SHORT,
	IPTC_FORMAT_LONG,
	IPTC_FORMAT_STRING,
	IPTC_FORMAT_NUMERIC_STRING,
	IPTC_FORMAT_DATE,
	IPTC_FORMAT_TIME
} IptcFormat;

typedef struct _IptcTagInfo IptcTagInfo;

struct _IptcTagInfo {
	IptcRecord	record;
	IptcTag		tag;
	const char     *name;
	const char     *title;
	const char     *description;
	IptcFormat	format;
	IptcMandatory	mandatory;
	IptcRepeatable	repeatable;
	unsigned int	minbytes;
	unsigned int	maxbytes;
};

const char     *iptc_tag_get_name        (IptcRecord record, IptcTag tag);
char           *iptc_tag_get_title       (IptcRecord record, IptcTag tag);
char           *iptc_tag_get_description (IptcRecord record, IptcTag tag);
const IptcTagInfo *iptc_tag_get_info     (IptcRecord record, IptcTag tag);
char           *iptc_format_get_name	 (IptcFormat format);

int iptc_tag_find_by_name (const char * name, IptcRecord * record, IptcTag * tag);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_TAG_H__ */
