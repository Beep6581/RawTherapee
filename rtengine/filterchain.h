#ifndef _FILTERCHAIN_H
#define _FILTERCHAIN_H

#include <set>
#include <vector>
#include "improclistener.h"
#include "imagesource.h"
#include "procparams.h"

namespace rtengine {

using namespace procparams;
class FilterChain {

protected:

    bool multiThread;
	ImageSource* first;
	Filter* last;
	ImProcListener* listener;
	Filter* firstToUpdate;
	ProcParams* procParams;
	std::vector<Glib::ustring> filterOrder;

	void setupChain (FilterChain* previous);

public:

	FilterChain (ImProcListener* listener, ImageSource* imgSource, ProcParams* params, bool multiThread);
	FilterChain (ImProcListener* listener, FilterChain* previous);
	~FilterChain ();
	void getReqiredBufferSize (int& w, int& h);
	void getFullImageSize (int& w, int& h);

	void setupProcessing (const set<ProcEvent>& events, int fullW, int fullH, int& maxWorkerWidth, int& maxWorkerHeight);
	void process (Buffer<int>* buffer);

	ImProcListener* getListener () { return listener; }
	void invalidate ();

	ImageSource* getImageSource () { return first; }
};
}
#endif
