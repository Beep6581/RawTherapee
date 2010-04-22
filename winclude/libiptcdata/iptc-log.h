/* iptc-log.h
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

#ifndef __IPTC_LOG_H__
#define __IPTC_LOG_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <libiptcdata/iptc-mem.h>
#include <stdarg.h>

typedef struct _IptcLog        IptcLog;

IptcLog *iptc_log_new     (void);
IptcLog *iptc_log_new_mem (IptcMem *);
void     iptc_log_ref     (IptcLog *log);
void     iptc_log_unref   (IptcLog *log);
void     iptc_log_free    (IptcLog *log);

typedef enum {
	IPTC_LOG_CODE_NONE,
	IPTC_LOG_CODE_DEBUG,
	IPTC_LOG_CODE_NO_MEMORY,
	IPTC_LOG_CODE_CORRUPT_DATA
} IptcLogCode;
const char *iptc_log_code_get_title   (IptcLogCode); /* Title for dialog   */
const char *iptc_log_code_get_message (IptcLogCode); /* Message for dialog */

typedef void (* IptcLogFunc) (IptcLog *log, IptcLogCode, const char *domain,
			      const char *format, va_list args, void *data);

void     iptc_log_set_func (IptcLog *log, IptcLogFunc func, void *data);

void     iptc_log  (IptcLog *log, IptcLogCode, const char *domain,
		    const char *format, ...)
#ifdef __GNUC__
			__attribute__((__format__(printf,4,5)))
#endif
;

void     iptc_logv (IptcLog *log, IptcLogCode, const char *domain,
		    const char *format, va_list args);

/* For your convenience */
#define IPTC_LOG_NO_MEMORY(l,d,s) iptc_log (l, IPTC_LOG_CODE_NO_MEMORY, d, "Could not allocate %i byte(s).", s)

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __IPTC_LOG_H__ */
