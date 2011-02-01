/*
 * buffer.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef BUFFER_H_
#define BUFFER_H_

namespace rtengine {

template <class T>
class Buffer {

        bool owner;

    public:
        int width;
        int height;
        T* data;
        T** rows;

        Buffer (int W, int H) : width(W), height(H), owner(true) {
            data = new T [width * height];
            rows = new T* [height];
            for (int i=0; i<height; i++)
                rows[i] = data + i*width;
		}
        Buffer (int W, int H, T* d, T** r)
            : width(W), height(H), data(d), rows(r), owner(false) {}
        ~Buffer () {
            if (owner) {
                delete [] data;
                delete [] rows;
            }
        }
};
}
#endif /* BUFFER_H_ */
