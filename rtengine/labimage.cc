#include <labimage.h>
namespace rtengine {

LabImage::LabImage (int w, int h) : fromImage(false), W(w), H(h) {

    L = new float*[H];
    a = new float*[H];
    b = new float*[H];

    data = new float [W*H*3];
    float * index = data;
    for (int i=0; i<H; i++)
        L[i] = index + i*W;
    index+=W*H;
    for (int i=0; i<H; i++)
        a[i] = index + i*W;
    index+=W*H;

    for (int i=0; i<H; i++)
        b[i] = index + i*W;
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
        delete [] L;
        delete [] a;
        delete [] b;
        delete [] data;
    }
}
}
