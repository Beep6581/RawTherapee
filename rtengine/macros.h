/*
 * macros.h
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#ifndef MACROS_H_
#define MACROS_H_

#undef MAX
#undef MIN
#undef CLIP
#undef CLIPTO
#undef CLIPTOC

#define MAX(a,b) ((a)<(b)?(b):(a))
#define MIN(a,b) ((a)>(b)?(b):(a))
#define CLIP(a) ((a)>0?((a)<65535?(a):65535):0)
#define CLIPTO(a,b,c) ((a)>(b)?((a)<(c)?(a):(c)):(b))
#define CLIPTOC(a,b,c,d) ((a)>=(b)?((a)<=(c)?(a):(d=true,(c))):(d=true,(b)))

#endif /* MACROS_H_ */
