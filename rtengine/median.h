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
#include "rt_math.h"

// middle 4 of 6 elements,
#define MIDDLE4OF6(s0,s1,s2,s3,s4,s5,d0,d1,d2,d3,d4,d5,temp) \
{\
d1 = min(s1,s2);\
d2 = max(s1,s2);\
d0 = min(s0,d2);\
d2 = max(s0,d2);\
temp = min(d0,d1);\
d1 = max(d0,d1);\
d0 = temp;\
d4 = min(s4,s5);\
d5 = max(s4,s5);\
d3 = min(s3,d5);\
d5 = max(s3,d5);\
temp = min(d3,d4);\
d4 = max(d3,d4);\
d3 = max(d0,temp);\
d2 = min(d2,d5);\
}

// middle 4 of 6 elements, vectorized
#define VMIDDLE4OF6(s0,s1,s2,s3,s4,s5,d0,d1,d2,d3,d4,d5,temp) \
{\
d1 = vminf(s1,s2);\
d2 = vmaxf(s1,s2);\
d0 = vminf(s0,d2);\
d2 = vmaxf(s0,d2);\
temp = vminf(d0,d1);\
d1 = vmaxf(d0,d1);\
d0 = temp;\
d4 = vminf(s4,s5);\
d5 = vmaxf(s4,s5);\
d3 = vminf(s3,d5);\
d5 = vmaxf(s3,d5);\
temp = vminf(d3,d4);\
d4 = vmaxf(d3,d4);\
d3 = vmaxf(d0,temp);\
d2 = vminf(d2,d5);\
}


#define MEDIAN7(s0,s1,s2,s3,s4,s5,s6,t0,t1,t2,t3,t4,t5,t6,median) \
{\
t0 = min(s0,s5);\
t5 = max(s0,s5);\
t3 = max(t0,s3);\
t0 = min(t0,s3);\
t1 = min(s1,s6);\
t6 = max(s1,s6);\
t2 = min(s2,s4);\
t4 = max(s2,s4);\
t1 = max(t0,t1);\
median = min(t3,t5);\
t5 = max(t3,t5);\
t3 = median;\
median = min(t2,t6);\
t6 = max(t2,t6);\
t3 = max(median,t3);\
t3 = min(t3,t6);\
t4 = min(t4,t5);\
median = min(t1,t4);\
t4 = max(t1,t4);\
t3 = max(median,t3);\
median = min(t3,t4);\
}

#define VMEDIAN7(s0,s1,s2,s3,s4,s5,s6,t0,t1,t2,t3,t4,t5,t6,median) \
{\
t0 = vminf(s0,s5);\
t5 = vmaxf(s0,s5);\
t3 = vmaxf(t0,s3);\
t0 = vminf(t0,s3);\
t1 = vminf(s1,s6);\
t6 = vmaxf(s1,s6);\
t2 = vminf(s2,s4);\
t4 = vmaxf(s2,s4);\
t1 = vmaxf(t0,t1);\
median = vminf(t3,t5);\
t5 = vmaxf(t3,t5);\
t3 = median;\
median = vminf(t2,t6);\
t6 = vmaxf(t2,t6);\
t3 = vmaxf(median,t3);\
t3 = vminf(t3,t6);\
t4 = vminf(t4,t5);\
median = vminf(t1,t4);\
t4 = vmaxf(t1,t4);\
t3 = vmaxf(median,t3);\
median = vminf(t3,t4);\
}
