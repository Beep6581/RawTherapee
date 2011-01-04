#include <labimage.h>
namespace rtengine {

LabImage::LabImage (int w, int h) : fromImage(false), W(w), H(h) {

    L = new float*[H];
    for (int i=0; i<H; i++)
        L[i] = new float[W];

    a = new float*[H];
    for (int i=0; i<H; i++)
        a[i] = new float[W];

    b = new float*[H];
    for (int i=0; i<H; i++)
        b[i] = new float[W];
}

LabImage::LabImage (Image16* im) {

    W = im->width;
    H = im->height;
	for (int i=0; i<H; i++) 
		for (int j=0; j<W; j++) {
			L[i][j] = im->r[i][j];
			a[i][j] = im->g[i][j];
			b[i][j] = im->b[i][j];
		}
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
