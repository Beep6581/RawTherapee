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
#include <batchqueue.h>
#include <glibmm.h>
#include <glib/gstdio.h>
#include <multilangmgr.h>
#include <filecatalog.h>
#include <batchqueuebuttonset.h>
#include <guiutils.h>
#include <safegtk.h>

using namespace rtengine;

BatchQueue::BatchQueue () : processing(NULL), listener(NULL)  {

    int p = 0;
    pmenu = new Gtk::Menu ();
    pmenu->attach (*(cancel = new Gtk::MenuItem (M("FILEBROWSER_POPUPCANCELJOB"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(head = new Gtk::MenuItem (M("FILEBROWSER_POPUPMOVEHEAD"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(tail = new Gtk::MenuItem (M("FILEBROWSER_POPUPMOVEEND"))), 0, 1, p, p+1); p++;
    pmenu->attach (*(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;
    pmenu->attach (*(selall = new Gtk::MenuItem (M("FILEBROWSER_POPUPSELECTALL"))), 0, 1, p, p+1); p++;
    pmenu->show_all ();

    cancel->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &BatchQueue::cancelItems), &selected));    
    head->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &BatchQueue::headItems), &selected));    
    tail->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &BatchQueue::tailItems), &selected));    
    selall->signal_activate().connect (sigc::mem_fun(*this, &BatchQueue::selectAll));
}

void BatchQueue::rightClicked (ThumbBrowserEntryBase* entry) {

    pmenu->popup (3, 0);
}

void BatchQueue::addEntries ( std::vector<BatchQueueEntry*> &entries, bool head)
{
	for( std::vector<BatchQueueEntry*>::iterator entry = entries.begin(); entry != entries.end();entry++ ){
		(*entry)->setParent (this);
		(*entry)->resize (options.thumbSize);
		Glib::ustring tempFile = getTempFilenameForParams( (*entry)->filename );

		// recovery save
		if( !(*entry)->params.save( tempFile ) )
		   (*entry)->savedParamsFile = tempFile;

		(*entry)->selected = false;
		if (!head)
			fd.push_back (*entry);
		else {
			std::vector<ThumbBrowserEntryBase*>::iterator pos;
			for (pos=fd.begin(); pos!=fd.end(); pos++)
				if (!(*pos)->processing) {
					fd.insert (pos, *entry);
					break;
				}
			if (pos==fd.end())
				fd.push_back (*entry);
		}
		if ((*entry)->thumbnail)
			(*entry)->thumbnail->imageEnqueued ();

		BatchQueueButtonSet* bqbs = new BatchQueueButtonSet (*entry);
		bqbs->setButtonListener (this);
		(*entry)->addButtonSet (bqbs);
	}
    saveBatchQueue( );
    arrangeFiles ();
    queue_draw ();
    notifyListener ();
}

bool BatchQueue::saveBatchQueue( )
{
    Glib::ustring savedQueueFile;
    savedQueueFile = options.rtdir+"/batch/queue";
    FILE *f = safe_g_fopen (savedQueueFile, "wt");

    if (f==NULL)
        return false;

    for (std::vector<ThumbBrowserEntryBase*>::iterator pos=fd.begin(); pos!=fd.end(); pos++){
    	BatchQueueEntry* bqe = reinterpret_cast<BatchQueueEntry*>(*pos);
    	fprintf(f,"%s;%s\n", bqe->filename.c_str(),bqe->savedParamsFile.c_str() );
    }
    fclose (f);
    return true;
}

bool BatchQueue::loadBatchQueue( )
{
    Glib::ustring savedQueueFile;
    savedQueueFile = options.rtdir+"/batch/queue";
    FILE *f = safe_g_fopen (savedQueueFile, "rt");

    if (f==NULL)
        return false;
    char *buffer = new char[1024];
    unsigned numLoaded=0;
    while (fgets (buffer, 1024, f)){
    	char *p = strchr(buffer,';' );
    	if( p ){
            char *le = buffer + strlen(buffer);
            while( --le > buffer && (*le == '\n' || *le == '\r') );
    	    Glib::ustring source(buffer, p-buffer );
    	    Glib::ustring paramsFile(p+1, (le +1)- (p+1) );

    	    rtengine::procparams::ProcParams pparams;
    	    if( pparams.load( paramsFile ) )
    	    	continue;

    	    ::Thumbnail *thumb = cacheMgr->getEntry( source );
    	    if( thumb ){
    	        rtengine::ProcessingJob* job = rtengine::ProcessingJob::create(source, thumb->getType() == FT_Raw, pparams);

				int prevh = options.maxThumbnailHeight;
				int prevw = prevh;
				guint8* prev = NULL;
				double tmpscale;
				rtengine::IImage8* img = thumb->processThumbImage(pparams, options.maxThumbnailHeight, tmpscale);
				if (img) {
					prevw = img->getWidth();
					prevh = img->getHeight();
					prev = new guint8[prevw * prevh * 3];
					memcpy(prev, img->getData(), prevw * prevh * 3);
					img->free();
				}
				BatchQueueEntry *entry = new BatchQueueEntry(job, pparams,source, prev, prevw, prevh, thumb);
				entry->setParent(this);
				entry->resize(options.thumbSize);
				entry->savedParamsFile = paramsFile;
				entry->selected = false;
				fd.push_back(entry);

				BatchQueueButtonSet* bqbs = new BatchQueueButtonSet(entry);
				bqbs->setButtonListener(this);
				entry->addButtonSet(bqbs);
				numLoaded++;
    	    }
    	}
    }
    delete [] buffer;
    fclose(f);
    arrangeFiles ();
    queue_draw ();
    return numLoaded > 0;
}

Glib::ustring BatchQueue::getTempFilenameForParams( const Glib::ustring filename )
{
    time_t rawtime;
    struct tm *timeinfo;
    char stringTimestamp [80];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (stringTimestamp,sizeof(stringTimestamp),"_%Y%m%d%H%M%S_",timeinfo);
    Glib::ustring savedParamPath;
    savedParamPath = options.rtdir+"/batch/";
    safe_g_mkdir_with_parents (savedParamPath, 0755);
    savedParamPath += Glib::path_get_basename (filename);
    savedParamPath += stringTimestamp;
    savedParamPath += paramFileExtension;
    return savedParamPath;
}

int deleteitem (void* data)
{
	safe_g_remove(  ((BatchQueueEntry*)data)->savedParamsFile );
    gdk_threads_enter ();
    delete (BatchQueueEntry*)data;
    gdk_threads_leave ();
    return 0;
}

void BatchQueue::cancelItems (std::vector<ThumbBrowserEntryBase*>* items) {

    for (int i=0; i<items->size(); i++) {
        BatchQueueEntry* entry = (BatchQueueEntry*)(*items)[i];
        if (entry->processing)
            continue;
        std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (fd.begin(), fd.end(), entry);
        if (pos!=fd.end()) {
            fd.erase (pos);
            rtengine::ProcessingJob::destroy (entry->job);
            if (entry->thumbnail)
                entry->thumbnail->imageRemovedFromQueue ();
            g_idle_add (deleteitem, entry);
        }
    }
    for (int i=0; i<fd.size(); i++) 
        fd[i]->selected = false;
    lastClicked = NULL;
    selected.clear ();

    saveBatchQueue( );

    redraw ();
    notifyListener ();
}

void BatchQueue::headItems (std::vector<ThumbBrowserEntryBase*>* items) {

    for (int i=items->size()-1; i>=0; i--) {
        BatchQueueEntry* entry = (BatchQueueEntry*)(*items)[i];
        if (entry->processing)
            continue;
        std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (fd.begin(), fd.end(), entry);
        if (pos!=fd.end() && pos!=fd.begin()) {
            fd.erase (pos);
            // find the first item that is not under processing
            for (pos=fd.begin(); pos!=fd.end(); pos++) 
                if (!(*pos)->processing) {
                    fd.insert (pos, entry);
                    break;
                }
        }
    }
    saveBatchQueue( );
    redraw ();
}

void BatchQueue::tailItems (std::vector<ThumbBrowserEntryBase*>* items) {

    for (int i=0; i<items->size(); i++) {
        BatchQueueEntry* entry = (BatchQueueEntry*)(*items)[i];
        if (entry->processing)
            continue;
        std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (fd.begin(), fd.end(), entry);
        if (pos!=fd.end()) {
            fd.erase (pos);
            fd.push_back (entry);
        }
    }
    saveBatchQueue( );
    redraw ();
}
   
void BatchQueue::selectAll () {

    lastClicked = NULL;
    selected.clear ();
    for (int i=0; i<fd.size(); i++) {
        if (fd[i]->processing)
            continue;
        fd[i]->selected = true;
        selected.push_back (fd[i]);
    }
    queue_draw ();
}
void BatchQueue::startProcessing () {

    if (!processing && fd.size()>0) {
        BatchQueueEntry* next = (BatchQueueEntry*)fd[0];
        // tag it as processing        
        next->processing = true;
        processing = next;
        // remove from selection
        if (processing->selected) {
            std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (selected.begin(), selected.end(), processing);
            if (pos!=selected.end())
                selected.erase (pos);
            processing->selected = false;
        }
        // remove button set
        next->removeButtonSet ();
        // start batch processing
        rtengine::startBatchProcessing (next->job, this);
        queue_draw ();
    }
}

rtengine::ProcessingJob* BatchQueue::imageReady (rtengine::IImage16* img) {

    gdk_threads_enter ();
    // save image img
    Glib::ustring fname;
    SaveFormat saveFormat;
    if (processing->outFileName=="") {   // auto file name
        Glib::ustring s = calcAutoFileNameBase (processing->filename);
        fname = autoCompleteFileName (s, options.saveFormat.format);
        saveFormat = options.saveFormat;
    }
    else {  // use the save-as filename with automatic completion for uniqueness
        fname = autoCompleteFileName (removeExtension(processing->outFileName), getExtension(processing->outFileName));
        saveFormat = processing->saveFormat;
    }
    //printf ("fname=%s, %s\n", fname.c_str(), removeExtension(fname).c_str());

    if (img && fname!="") {
        int err = 0;
        if (saveFormat.format=="tif")
            err = img->saveAsTIFF (fname, saveFormat.tiffBits,saveFormat.tiffUncompressed);
        else if (saveFormat.format=="png")
            err = img->saveAsPNG (fname, saveFormat.pngCompression, saveFormat.pngBits);
        else if (saveFormat.format=="jpg")
            err = img->saveAsJPEG (fname, saveFormat.jpegQuality);
        img->free ();

		if (err) throw "Unable to save output file";

        if (saveFormat.saveParams) {
			// We keep the extension to avoid overwriting the profile when we have
			// the same output filename with different extension
            //processing->params.save (removeExtension(fname) + paramFileExtension);
            processing->params.save (fname + paramFileExtension);
        }

        if (processing->thumbnail) {
            processing->thumbnail->imageDeveloped ();
            processing->thumbnail->imageRemovedFromQueue ();
            if (listener)
                listener->imageProcessingReady (processing->filename);
        }
    }
    // save temporary params file name: delete as last thing
    Glib::ustring processedParams = processing->savedParamsFile;
    
    // delete from the queue
    delete processing;
    processing = NULL;
    fd.erase (fd.begin());

    // return next job
    if (fd.size()==0) {
        if (listener)
            listener->queueEmpty ();
    }
    else if (listener && listener->canStartNext ()) {
        BatchQueueEntry* next = (BatchQueueEntry*)fd[0];
        // tag it as selected        
        next->processing = true;
        processing = next;
        // remove from selection
        if (processing->selected) {
            std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (selected.begin(), selected.end(), processing);
            if (pos!=selected.end())
                selected.erase (pos);
            processing->selected = false;
        }
        // remove button set
        next->removeButtonSet ();
    }
    if( saveBatchQueue( ) ){
       safe_g_remove( processedParams );
       // Delete all files in directory \batch when finished, just to be sure to remove zombies
       if( fd.size()==0 ){
    	    std::vector<Glib::ustring> names;
    	    Glib::ustring batchdir = options.rtdir+"/batch/";
    	    Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (batchdir);
    		safe_build_file_list (dir, names, batchdir);
    		for(std::vector<Glib::ustring>::iterator iter=names.begin(); iter != names.end();iter++ )
    			safe_g_remove( *iter );
       }
    }
    redraw ();
    notifyListener ();
    gdk_threads_leave ();
    return processing ? processing->job : NULL;
}

// Calculates automatic filename of processed batch entry, but just the base name
// example output: "c:\out\converted\dsc0121"
Glib::ustring BatchQueue::calcAutoFileNameBase (const Glib::ustring& origFileName) {

    std::vector<Glib::ustring> pa;
    std::vector<Glib::ustring> da;

    for (int i=0; i<origFileName.size(); i++) {
        while ((i<origFileName.size()) && (origFileName[i]=='\\' || origFileName[i]=='/'))
            i++;
        if (i>=origFileName.size())
            break;
        Glib::ustring tok = "";
        while ((i<origFileName.size()) && !(origFileName[i]=='\\' || origFileName[i]=='/')) 
            tok = tok + origFileName[i++];
        da.push_back (tok);
    }

    if (origFileName[0]=='/' || origFileName[0]=='\\')
        pa.push_back ("/" + da[0]);
    else
        pa.push_back (da[0]);
        
    for (int i=1; i<da.size(); i++)
        pa.push_back (pa[i-1] + "/" + da[i]);

//    for (int i=0; i<da.size(); i++)
//        printf ("da: %s\n", da[i].c_str());
//    for (int i=0; i<pa.size(); i++)
//        printf ("pa: %s\n", pa[i].c_str());

    // extracting filebase
    Glib::ustring filename;
    
    int extpos = origFileName.size()-1;
    for (; extpos>=0 && origFileName[extpos]!='.'; extpos--);
    for (int k=extpos-1; k>=0 && origFileName[k]!='/' && origFileName[k]!='\\'; k--)
        filename = origFileName[k] + filename;

//    printf ("%d, |%s|\n", extpos, filename.c_str());
 
    // constructing full output path
//    printf ("path=|%s|\n", options.savePath.c_str());

    Glib::ustring path="";
    if (options.saveUsePathTemplate) {
        int ix=0;
        while (options.savePathTemplate[ix]!=0) {
            if (options.savePathTemplate[ix]=='%') {
                ix++;
                if (options.savePathTemplate[ix]=='p') {
                    ix++;
                    int i = options.savePathTemplate[ix]-'0';
                    if (i<pa.size())
                        path = path + pa[pa.size()-i-1] + '/';
                    ix++;
                }
                else if (options.savePathTemplate[ix]=='d') {
                    ix++;
                    int i = options.savePathTemplate[ix]-'0';
                    if (i<da.size())
                        path = path + da[da.size()-i-1] + '/';
                    ix++;
                }
                else if (options.savePathTemplate[ix]=='f') {
                    path = path + filename;
                }
            }
            else
                path = path + options.savePathTemplate[ix];
            ix++;
        }
    }
    else
        path = Glib::build_filename (options.savePathFolder, filename);

    return path;
}

Glib::ustring BatchQueue::autoCompleteFileName (const Glib::ustring& fileName, const Glib::ustring& format) {

    // separate filename and the path to the destination directory
    Glib::ustring dstdir = Glib::path_get_dirname (fileName);
    Glib::ustring dstfname = Glib::path_get_basename (fileName);

    // create directory, if does not exist
    if (safe_g_mkdir_with_parents (dstdir, 0755) ) 
        return "";
	
	// In overwrite mode we TRY to delete the old file first.
	// if that's not possible (e.g. locked by viewer, R/O), we revert to the standard naming scheme
	bool inOverwriteMode=options.overwriteOutputFile;

    for (int tries=0; tries<100; tries++) {
        Glib::ustring fname;
        if (tries==0)
            fname = Glib::ustring::compose ("%1.%2", Glib::build_filename (dstdir,  dstfname), format);
        else
            fname = Glib::ustring::compose ("%1-%2.%3", Glib::build_filename (dstdir,  dstfname), tries, format);

		int fileExists=safe_file_test (fname, Glib::FILE_TEST_EXISTS); 
        
		if (inOverwriteMode && fileExists) {
			if (safe_g_remove(fname) == -1)
				inOverwriteMode = false;  // failed to delete- revert to old naming scheme
			else
				fileExists = false;  // deleted now
		}
		
		if (!fileExists) {
            return fname;
        }
    }
}

int bqredraw (void* p) {

    gdk_threads_enter ();
    ((BatchQueue*)p)->redraw();
    gdk_threads_leave ();
    return 0;
}

void BatchQueue::setProgress (double p) {

    if (processing)
        processing->progress = p;

    g_idle_add (bqredraw, this);
}

void BatchQueue::buttonPressed (LWButton* button, int actionCode, void* actionData) {
    
    std::vector<ThumbBrowserEntryBase*> bqe;
    bqe.push_back ((BatchQueueEntry*)actionData);

    if (actionCode==10)  // cancel
        cancelItems (&bqe);
    else if (actionCode==8)  // to head
        headItems (&bqe);
    else if (actionCode==9)  // to tail
        tailItems (&bqe);
}

struct NLParams {
    BatchQueueListener* listener;
    int qsize;
};

int bqnotifylistener (void* data) {

    gdk_threads_enter ();
    NLParams* params = (NLParams*)data;
    params->listener->queueSizeChanged (params->qsize);
    delete params;
    gdk_threads_leave ();
    return 0;
}

void BatchQueue::notifyListener () {

    if (listener) {
        NLParams* params = new NLParams;
        params->listener = listener;
        params->qsize = fd.size();
        g_idle_add (bqnotifylistener, params);
    }
}

void BatchQueue::redrawNeeded (LWButton* button) {
    
    queue_draw ();
}
