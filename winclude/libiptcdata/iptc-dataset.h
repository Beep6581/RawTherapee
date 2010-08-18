/* iptc-dataset.h
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

#ifndef __IPTC_DATASET_H__
#define __IPTC_DATASET_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef struct _IptcDataSet        IptcDataSet;
typedef struct _IptcDataSetPrivate IptcDataSetPrivate;

typedef enum {
	IPTC_DONT_VALIDATE = 0,
	IPTC_VALIDATE      = 1
} IptcValidate;

#include <libiptcdata/iptc-data.h>
#include <libiptcdata/iptc-mem.h>

struct _IptcDataSet {
	IptcRecord record;
        IptcTag tag;
        const IptcTagInfo * info;

        unsigned char *data;
        unsigned int size;

	/* Data containing this dataset */
	IptcData *parent;

	IptcDataSetPrivate *priv;
};


/* Lifecycle */
IptcDataSet  *iptc_dataset_new     (void);
IptcDataSet  *iptc_dataset_new_mem (IptcMem * mem);
IptcDataSet  *iptc_dataset_copy    (IptcDataSet *dataset);
void        iptc_dataset_ref     (IptcDataSet *dataset);
void        iptc_dataset_unref   (IptcDataSet *dataset);
void        iptc_dataset_free  (IptcDataSet *dataset);

void        iptc_dataset_set_tag (IptcDataSet *dataset, IptcRecord record, IptcTag tag);
IptcFormat  iptc_dataset_get_format (IptcDataSet *dataset);

int         iptc_dataset_get_data (IptcDataSet *dataset, unsigned char * buf,
				unsigned int size);
unsigned int iptc_dataset_get_value (IptcDataSet *dataset);
int	    iptc_dataset_get_date (IptcDataSet *dataset, int *year, int *month, int *day);
int	    iptc_dataset_get_time (IptcDataSet *dataset, int *hour, int *min, int *sec,
				int *tz);

int         iptc_dataset_set_data (IptcDataSet *dataset, const unsigned char * buf,
				unsigned int size, IptcValidate validate);
int         iptc_dataset_set_value (IptcDataSet *dataset, unsigned int value,
				IptcValidate validate);
int	    iptc_dataset_set_date (IptcDataSet *dataset, int year, int month, int day,
				IptcValidate validate);
int	    iptc_dataset_set_time (IptcDataSet *dataset, int hour, int min, int sec,
				int tz, IptcValidate validate);


/* For your convenience */
const char *iptc_dataset_get_as_str (IptcDataSet *dataset, char *buf,
				  unsigned int size);

void        iptc_dataset_dump      (IptcDataSet *dataset, unsigned int indent);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_DATASET_H__ */
