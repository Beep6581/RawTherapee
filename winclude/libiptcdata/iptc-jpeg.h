/* iptc-jpeg.h
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

#ifndef __IPTC_JPEG_H__
#define __IPTC_JPEG_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <stdio.h>
#include <libiptcdata/iptc-data.h>

int iptc_jpeg_read_ps3 (FILE * infile, unsigned char * buf, unsigned int size);
int iptc_jpeg_ps3_find_iptc (const unsigned char * ps3,
		unsigned int ps3_size, unsigned int * iptc_len);

int iptc_jpeg_ps3_save_iptc (const unsigned char * ps3, unsigned int ps3_size,
		const unsigned char * iptc, unsigned int iptc_size,
		unsigned char * buf, unsigned int size);
int iptc_jpeg_save_with_ps3 (FILE * infile, FILE * outfile,
		const unsigned char * ps3, unsigned int ps3_size);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_JPEG_H__ */
