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
#define SORT3(a1,a2,a3,b1,b2,b3) \
  {  \
   if ((a1)<(a2)) {   \
     if ((a2)<(a3)) { \
       (b1) = (a1); (b2) = (a2); (b3) = (a3); \
     } \
     else if ((a1)<(a3)) { \
       (b1) = (a1); (b2) = (a3); (b3) = (a2); \
     } \
     else { \
       (b1) = (a3); (b2) = (a1); (b3) = (a2); \
     } \
   } \
   else { \
     if ((a3)<(a2)) { \
       (b1) = (a3); (b2) = (a2); (b3) = (a1); \
     } \
     else if ((a3)<(a1)) { \
       (b1) = (a2); (b2) = (a3); (b3) = (a1); \
     } \
     else { \
       (b1) = (a2); (b2) = (a1); (b3) = (a3); \
     } \
   } \
 }

#define MERGESORT(a1,a2,a3,b1,b2,b3,c1,c2,c3,c4,c5,c6) \
  {\
  if (a1<b1) { \
    c1 = a1; \
    if (a2<b1) {  \
      c2 = a2; \
      if (a3<b1) { \
        c3 = a3; c4 = b1; c5 = b2; c6 = b3;\
      }\
      else {\
        c3 = b1;\
        if (a3<b2) {\
          c4 = a3; c5 = b2; c6 = b3;\
        }\
        else {\
          c4 = b2;\
          if (a3<b3) {\
            c5 = a3; c6 = b3;\
          }\
          else {\
            c5 = b3; c6 = a3;\
          }\
        }\
      }\
    }\
    else {\
      c2 = b1;\
      if (a2<b2) {\
        c3 = a2;\
        if (a3<b2) {\
          c4 = a3; c5 = b2; c6 = b3;\
        }\
        else {\
          c4 = b2;\
          if (a3<b3) {\
            c5 = a3; c6 = b3;\
          }\
          else {\
            c5 = b3; c6 = a3;\
          }\
        }\
      }\
      else {\
        c3 = b2;\
        if (a2<b3) {\
          c4 = a2;\
          if (a3<b3) {\
            c5 = a3; c6 = b3;\
          }\
          else {\
            c5 = b3; c6 = a3;\
          }\
        }\
        else {\
          c4 = b3; c5 = a2; c6 = a3;\
        }\
      }\
    }\
  }\
  else {\
    c1 = b1;\
    if (a1<b2) {\
      c2 = a1;\
      if (a2<b2) {\
        c3 = a2;\
        if (a3<b2) {\
          c4 = a3; c5 = b2; c6 = b3;\
        }\
        else {\
          c4 = b2;\
          if (a3<b3) {\
            c5 = a3; c6 = b3;\
          }\
          else {\
            c5 = b3; c6 = a3;\
          }\
        }\
      }\
      else {\
        c3 = b2;\
        if (a2<b3) {\
          c4 = a2;\
          if (a3<b3) {\
            c5 = a3; c6 = b3;\
          }\
          else {\
            c5 = b3; c6 = a3;\
          }\
        }\
        else {\
          c4 = b3; c5 = a2; c6 = a3;\
        }\
      }\
    }\
    else {\
      c2 = b2;\
      if (a1<b3) {\
        c3 = a1;\
        if (a2<b3) {\
          c4 = a2;\
          if (a3<b3) {\
            c5 = a3; c6 = b3;\
          }\
          else {\
            c5 = b3; c6 = a3;\
          }\
        }\
        else {\
          c4 = b3; c5 = a2; c6 = a3;\
        }\
      }\
      else {\
        c3 = b3; c4 = a1; c5 = a2; c6 = a3;\
      }\
    }\
  }\
}

#define MEDIAN7(a1,a2,a3,b1,b2,b3,b4,median) \
  { \
  if (a1<b1) {\
    if (a2<b1) {\
      if (a3<b1) \
        median = b1; \
      else {\
        if (a3<b2) \
          median = a3;\
        else \
          median = b2;\
      }\
    }\
    else {\
      if (a2<b2) {\
        if (a3<b2) \
          median = a3;\
        else \
          median = b2;\
      }\
      else {\
        if (a2<b3) \
          median = a2;\
        else \
          median = b3;\
      }\
    }\
  }\
  else {\
    if (a1<b2) {\
      if (a2<b2) {\
        if (a3<b2) \
          median = a3;\
        else\
          median = b2;\
      }\
      else {\
        if (a2<b3)\
          median = a2;\
        else\
          median = b3;\
      }\
    }\
    else {\
      if (a1<b3) {\
        if (a2<b3)\
          median = a2;\
        else \
          median = b3; \
      }\
      else {\
        if (a1<b4) \
          median = a1;\
        else\
          median = b4;\
      } \
    }\
  }\
}

