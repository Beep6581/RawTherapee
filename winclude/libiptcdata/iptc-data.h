/* iptc-data.h
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

#ifndef __IPTC_DATA_H__
#define __IPTC_DATA_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _IptcData        IptcData;
typedef struct _IptcDataPrivate IptcDataPrivate;

#include <libiptcdata/iptc-tag.h>
#include <libiptcdata/iptc-dataset.h>
#include <libiptcdata/iptc-mem.h>
#include <libiptcdata/iptc-log.h>
#include <stdio.h>

typedef enum {
	IPTC_ENCODING_UNKNOWN = 0,
	IPTC_ENCODING_UNSPECIFIED = 1,
	IPTC_ENCODING_UTF8 = 2
} IptcEncoding;

/* The version of the spec implemented by this library */
#define IPTC_IIM_VERSION	4

struct _IptcData
{
        IptcDataSet **datasets;
        unsigned int count;

	IptcDataPrivate *priv;
};

/* Lifecycle */
IptcData    *iptc_data_new     (void);
IptcData    *iptc_data_new_mem (IptcMem *mem);
IptcData    *iptc_data_new_from_jpeg (const char *path);
IptcData    *iptc_data_new_from_jpeg_file (FILE* infile);
IptcData    *iptc_data_new_from_data (const unsigned char *buf,
				   unsigned int size);
void         iptc_data_ref     (IptcData *data);
void         iptc_data_unref   (IptcData *data);
void         iptc_data_free    (IptcData *data);

int          iptc_data_load (IptcData *data, const unsigned char *buf, 
			       unsigned int size);
int          iptc_data_save (IptcData *data, unsigned char **buf,
			       unsigned int *size);
void         iptc_data_free_buf (IptcData *data, unsigned char *buf);
			       
int          iptc_data_add_dataset     (IptcData *data, IptcDataSet *ds);
int          iptc_data_add_dataset_before (IptcData *data, IptcDataSet *ds,
						IptcDataSet *newds);
int          iptc_data_add_dataset_after (IptcData *data, IptcDataSet *ds,
						IptcDataSet *newds);
int          iptc_data_remove_dataset  (IptcData *data, IptcDataSet *ds);
IptcDataSet *iptc_data_get_dataset     (IptcData *data, IptcRecord record,
						IptcTag tag);
IptcDataSet *iptc_data_get_next_dataset (IptcData *data, IptcDataSet *ds,
						IptcRecord record, IptcTag tag);

typedef void (* IptcDataForeachDataSetFunc) (IptcDataSet *dataset,
		void *user_data);
void         iptc_data_foreach_dataset (IptcData *data,
					 IptcDataForeachDataSetFunc func,
					 void *user_data);
void         iptc_data_sort (IptcData *data);
IptcEncoding iptc_data_get_encoding (IptcData *data);
int          iptc_data_set_encoding_utf8 (IptcData *data);

int          iptc_data_set_version (IptcData *data, unsigned int version);

int iptc_data_add_dataset_with_value (IptcData *data, IptcRecord record,
		IptcTag tag, unsigned int value, IptcValidate validate);
int iptc_data_add_dataset_with_contents (IptcData *data, IptcRecord record,
		IptcTag tag, const unsigned char * buf,
		unsigned int size, IptcValidate validate);

void iptc_data_dump  (IptcData *data, unsigned int indent);
void iptc_data_log  (IptcData *data, IptcLog *log);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_DATA_H__ */
