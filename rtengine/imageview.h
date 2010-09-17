#ifndef _IMAGEVIEW_H_
#define _IMAGEVIEW_H_

namespace rtengine {

class ImageView {

public:
	int x, y, w, h, skip;

	ImageView (int x, int y, int w, int h, int skip=1) : x(x), y(y), w(w), h(h), skip(skip) {}

	bool operator== (const ImageView& other) { return x==other.x && y==other.y && w==other.w && h==other.h && skip==other.skip; }
	int getPixelWidth () { return w >> skip; }
	int getPixelHeight () { return h >> skip; }
};

}

#endif
