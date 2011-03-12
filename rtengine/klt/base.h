/*********************************************************************
 * base.h
 *********************************************************************/

#ifndef _BASE_H_
#define _BASE_H_

#ifndef uchar
#define uchar unsigned char
#endif

#ifndef schar
#define schar signed char
#endif

#ifndef uint
#define uint unsigned int
#endif

#ifndef ushort
#define ushort unsigned short
#endif

#ifndef ulong
#define ulong unsigned long
#endif

#ifndef max
#define max(a,b)	((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)	((a) < (b) ? (a) : (b))
#endif
#define max3(a,b,c)	((a) > (b) ? max((a),(c)) : max((b),(c)))
#define min3(a,b,c)	((a) < (b) ? min((a),(c)) : min((b),(c)))

#endif

