#include "cieimage.h"
#include <memory.h>
namespace rtengine {

CieImage::CieImage (int w, int h) : fromImage(false), W(w), H(h) {
    J_p = new float*[H];
    Q_p = new float*[H];
    M_p = new float*[H];
    C_p = new float*[H];
 //   s_p = new float*[H];
    h_p = new float*[H];
	
    data = new float [W*H*5];
    float * index = data;
    for (int i=0; i<H; i++)
        J_p[i] = index + i*W;
    index+=W*H;
    for (int i=0; i<H; i++)
        Q_p[i] = index + i*W;
    index+=W*H;
   for (int i=0; i<H; i++)
        M_p[i] = index + i*W;
    index+=W*H;
    for (int i=0; i<H; i++)
        C_p[i] = index + i*W;
    index+=W*H;
/*    for (int i=0; i<H; i++)
        s_p[i] = index + i*W;
    index+=W*H;*/
    for (int i=0; i<H; i++)
        h_p[i] = index + i*W;
}

CieImage::~CieImage () {

    if (!fromImage) {
        delete [] J_p;       
		delete [] Q_p;
        delete [] M_p;
        delete [] C_p;
    //    delete [] s_p;
        delete [] h_p;
		
		
        delete [] data;
    }
}

void CieImage::CopyFrom(CieImage *Img){
	memcpy(data, Img->data, W*H*5*sizeof(float));
}

}
