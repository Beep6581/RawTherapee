/*
 * dim.h
 *
 *  Created on: Sep 30, 2010
 *      Author: gabor
 */

#ifndef DIM_H_
#define DIM_H_

#include <algorithm>

namespace rtengine {

class Dim {

    public:
        int width, height;

        Dim (int w=0, int h=0) : width(w), height(h) {}
        bool operator== (const Dim& other) { return width == other.width && height == other.height; }
        bool operator!= (const Dim& other) { return !(*this == other); }
        bool operator<  (const Dim& other) { return width < other.width && height < other.height; }
        void setMax (const Dim& other)     { width = std::max(width, other.width); height = std::max(height, other.height); }
        bool nonZero ()                    { return width > 0 && height > 0; }
};

}

#endif /* DIM_H_ */
