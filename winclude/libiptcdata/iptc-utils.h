/* iptc-utils.h
 *
 * Copyright © 2005 David Moore <dcm@acm.org>
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

#ifndef __IPTC_UTILS_H__
#define __IPTC_UTILS_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <libiptcdata/_stdint.h>

typedef enum {
	IPTC_BYTE_ORDER_MOTOROLA,
	IPTC_BYTE_ORDER_INTEL
} IptcByteOrder;


/* If these definitions don't work for you, please let us fix the 
 * macro generating _stdint.h */
	
typedef char		IptcByte;          /* 1 byte  */
typedef uint16_t	IptcShort;         /* 2 bytes */
typedef uint32_t	IptcLong;          /* 4 bytes */
typedef int32_t		IptcSLong;         /* 4 bytes */


IptcShort     iptc_get_short     (const unsigned char *b, IptcByteOrder order);
IptcLong      iptc_get_long      (const unsigned char *b, IptcByteOrder order);
IptcSLong     iptc_get_slong     (const unsigned char *b, IptcByteOrder order);

void iptc_set_short     (unsigned char *b, IptcByteOrder order,
			 IptcShort value);
void iptc_set_long      (unsigned char *b, IptcByteOrder order,
			 IptcLong value);
void iptc_set_slong     (unsigned char *b, IptcByteOrder order,
			 IptcSLong value);

#undef  MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

/* For compatibility with older versions */
#define IPTC_TAG_SUBSEC_TIME IPTC_TAG_SUB_SEC_TIME

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_UTILS_H__ */
