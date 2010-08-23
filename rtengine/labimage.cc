#include <labimage.h>
namespace rtengine {

LabImage::LabImage (int w, int h) : W(w), H(h), fromImage(false) {

    L = new unsigned short*[H];
    for (int i=0; i<H; i++)
        L[i] = new unsigned short[W];

    a = new short*[H];
    for (int i=0; i<H; i++)
        a[i] = new short[W];

    b = new short*[H];
    for (int i=0; i<H; i++)
        b[i] = new short[W];
}

LabImage::LabImage (Image16* im) {

    W = im->width;
    H = im->height;
    L = im->r;
    a = (short**) im->g;
    b = (short**) im->b;
    fromImage = true;
}

LabImage::~LabImage () {

    if (!fromImage) {
        for (int i=0; i<H; i++) {
            delete [] L[i];
            delete [] a[i];
            delete [] b[i];
        }
        delete [] L;
        delete [] a;
        delete [] b;
    }
}
}
