#include "cieimage.h"

#include <new>
#include <cstring>
namespace rtengine
{

CieImage::CieImage (int w, int h) : fromImage(false), W(w), H(h)
{
    J_p = new float*[H];
    Q_p = new float*[H];
    M_p = new float*[H];
    C_p = new float*[H];
    sh_p = new float*[H];
    h_p = new float*[H];

    // Initialize the pointers to zero
    for (unsigned int c = 0; c < 6; ++c) {
        data[c] = nullptr;
    }

    // Trying to allocate all in one block
    data[0] = new (std::nothrow) float [W * H * 6];

    if (data[0]) {
        float * index = data[0];

        for (int i = 0; i < H; i++) {
            J_p[i] = index + i * W;
        }

        index += W * H;

        for (int i = 0; i < H; i++) {
            Q_p[i] = index + i * W;
        }

        index += W * H;

        for (int i = 0; i < H; i++) {
            M_p[i] = index + i * W;
        }

        index += W * H;

        for (int i = 0; i < H; i++) {
            C_p[i] = index + i * W;
        }

        index += W * H;

        for (int i = 0; i < H; i++) {
            sh_p[i] = index + i * W;
        }

        index += W * H;

        //   for (int i=0; i<H; i++)
        //        ch_p[i] = index + i*W;
        //   index+=W*H;
        for (int i = 0; i < H; i++) {
            h_p[i] = index + i * W;
        }
    } else {
        // Allocating each plane separately
        for (unsigned int c = 0; c < 6; ++c) {
            data[c] = new float [W * H];
        }

        unsigned int c = 0;

        for (int i = 0; i < H; i++) {
            J_p[i] = data[c] + i * W;
        }

        ++c;

        for (int i = 0; i < H; i++) {
            Q_p[i] = data[c] + i * W;
        }

        ++c;

        for (int i = 0; i < H; i++) {
            M_p[i] = data[c] + i * W;
        }

        ++c;

        for (int i = 0; i < H; i++) {
            C_p[i] = data[c] + i * W;
        }

        ++c;

        for (int i = 0; i < H; i++) {
            sh_p[i] = data[c] + i * W;
        }

        ++c;

        for (int i = 0; i < H; i++) {
            h_p[i] = data[c] + i * W;
        }
    }
}

CieImage::~CieImage ()
{

    if (!fromImage) {
        delete [] J_p;
        delete [] Q_p;
        delete [] M_p;
        delete [] C_p;
        delete [] sh_p;
//      delete [] ch_p;
        delete [] h_p;

        for (unsigned int c = 0; c < 6; ++c)
            if (data[c]) {
                delete [] data[c];
            }
    }
}

void CieImage::CopyFrom(CieImage *Img)
{
    if (!data [1])
        // Only one allocated block
    {
        memcpy(data, Img->data, W * H * 6 * sizeof(float));
    } else

        // Separate allocation
        for (unsigned int c = 0; c < 6; ++c) {
            memcpy(data[c], Img->data[c], W * H * sizeof(float));
        }
}

}
