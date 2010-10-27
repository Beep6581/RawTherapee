/* iptc-mem.h
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

#ifndef __IPTC_MEM_H__
#define __IPTC_MEM_H__

#include <libiptcdata/iptc-utils.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* Should work like calloc: Needs to return initialized memory. */
typedef void * (* IptcMemAllocFunc)   (IptcLong);

typedef void * (* IptcMemReallocFunc) (void *, IptcLong);
typedef void   (* IptcMemFreeFunc)    (void *);

typedef struct _IptcMem IptcMem;

IptcMem *iptc_mem_new   (IptcMemAllocFunc, IptcMemReallocFunc,
			 IptcMemFreeFunc);
void     iptc_mem_ref   (IptcMem *);
void     iptc_mem_unref (IptcMem *);

void *iptc_mem_alloc   (IptcMem *, IptcLong);
void *iptc_mem_realloc (IptcMem *, void *, IptcLong);
void  iptc_mem_free    (IptcMem *, void *);

/* For your convenience */
IptcMem *iptc_mem_new_default (void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_MEM_H__ */
