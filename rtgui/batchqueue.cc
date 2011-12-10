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
#include <processingjob.h>
#include <rtimage.h>
#include <cstring>

using namespace rtengine;

BatchQueue::BatchQueue () : processing(NULL), listener(NULL)  {

    int p = 0;
    pmenu = new Gtk::Menu ();
    pmenu->attach (*Gtk::manage(selall = new Gtk::MenuItem (M("FILEBROWSER_POPUPSELECTALL"))), 0, 1, p, p+1); p++;
    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;

    pmenu->attach (*Gtk::manage(head = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPMOVEHEAD"))), 0, 1, p, p+1); p++;
    head->set_image(*Gtk::manage(new RTImage ("toleftend.png")));

    pmenu->attach (*Gtk::manage(tail = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPMOVEEND"))), 0, 1, p, p+1); p++;
    tail->set_image(*Gtk::manage(new RTImage ("torightend.png")));

    pmenu->attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p+1); p++;

    pmenu->attach (*Gtk::manage(cancel = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPCANCELJOB"))), 0, 1, p, p+1); p++;
    cancel->set_image(*Gtk::manage(new RTImage ("gtk-close.png")));

    pmenu->show_all ();

    cancel->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &BatchQueue::cancelItems), &selected));    
    head->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &BatchQueue::headItems), &selected));    
    tail->signal_activate().connect (sigc::bind(sigc::mem_fun(*this, &BatchQueue::tailItems), &selected));    
    selall->signal_activate().connect (sigc::mem_fun(*this, &BatchQueue::selectAll));

    setArrangement (ThumbBrowserBase::TB_Vertical);
}

BatchQueue::~BatchQueue ()
{
	delete pmenu;
}

// Reduce the max size of a thumb, since thumb is processed synchronously on adding to queue
// leading to very long waiting when adding more images
int BatchQueue::calcMaxThumbnailHeight() {
    return MIN(options.maxThumbnailHeight, 200);
}

// Function for virtual override in thumbbrowser base
int BatchQueue::getMaxThumbnailHeight() const {
    return calcMaxThumbnailHeight();
}


void BatchQueue::rightClicked (ThumbBrowserEntryBase* entry) {

    pmenu->popup (3, this->eventTime);
}

void BatchQueue::addEntries ( std::vector<BatchQueueEntry*> &entries, bool head)
{
	{
		// TODO: Check for Linux
		#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
		#endif

	for( std::vector<BatchQueueEntry*>::iterator entry = entries.begin(); entry != entries.end();entry++ ){
		(*entry)->setParent (this);
		(*entry)->resize (MIN(options.thumbSize, getMaxThumbnailHeight()));  // batch queue might have smaller, restricted size

	    // Is not and already present snapshot, so we create a new one
		if( (*entry)->currentSnapshoId <0 ){
		    time_t rawtime;
		    struct tm *timeinfo;
		    char stringTimestamp [80];
		    time ( &rawtime );
		    timeinfo = localtime ( &rawtime );
		    strftime (stringTimestamp,sizeof(stringTimestamp),"Queued_%Y-%m-%d %H:%M:%S",timeinfo);

			int id = (*entry)->thumbnail->newSnapshot(stringTimestamp,(*entry)->params,true );
			(*entry)->currentSnapshoId = id;
		}else
			(*entry)->thumbnail->setQueued((*entry)->currentSnapshoId,true );


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

		BatchQueueButtonSet* bqbs = new BatchQueueButtonSet (*entry);
		bqbs->setButtonListener (this);
		(*entry)->addButtonSet (bqbs);
	}
    saveBatchQueue( );
	}

    redraw();
    notifyListener (false);
}

bool BatchQueue::saveBatchQueue( )
{
    Glib::ustring savedParamPath;
    savedParamPath = options.rtdir+"/batch/";
    safe_g_mkdir_with_parents (savedParamPath, 0755);

    Glib::ustring savedQueueFile;
    savedQueueFile = options.rtdir+"/batch/queue";
    FILE *f = safe_g_fopen (savedQueueFile, "wt");

    if (f==NULL)
        return false;

	// method is already running with entryLock, so no need to lock again
    for (std::vector<ThumbBrowserEntryBase*>::iterator pos=fd.begin(); pos!=fd.end(); pos++){
    	BatchQueueEntry* bqe = reinterpret_cast<BatchQueueEntry*>(*pos);
    	fprintf(f,"%s;%d\n", bqe->filename.c_str(),bqe->currentSnapshoId );
    }
    fclose (f);
    return true;
}

void BatchQueue::loadBatchQueue( )
{
    {
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
#endif

        Glib::ustring savedQueueFile;
        savedQueueFile = options.rtdir+"/batch/queue";
        FILE *f = safe_g_fopen (savedQueueFile, "rt");

        if (f!=NULL) {
            char *buffer = new char[1024];
            unsigned numLoaded=0;
            while (fgets (buffer, 1024, f)){
                char *p = strchr(buffer,';' );
                if( p ){
                    std::string _source(buffer, p-buffer );
                    int id=-1;
                    sscanf(p+1,"%d",&id );
                    Glib::ustring source(_source);

                    ::Thumbnail *thumb = cacheMgr->getEntry( source );
                    if( thumb ){
                    	SnapshotInfo si = thumb->getSnapshot( id );
                    	if( si.id != id ) // better checking for info returned
                    		continue;
                    	rtengine::ImageMetaData* md= new rtengine::ImageMetaData( *(thumb->getMetadata()) );
                    	rtengine::ProcessingJobImpl* job = (rtengine::ProcessingJobImpl*)rtengine::ProcessingJob::create(source, thumb->getType() == FT_Raw, si.params, md,options.outputMetaData);

                        int prevh = getMaxThumbnailHeight();
                        int prevw = prevh;
                        guint8* prev = NULL;
                        double tmpscale;
                        rtengine::IImage8* img = thumb->processThumbImage( si.params, prevh, tmpscale);
                        if (img) {
                            prevw = img->getWidth();
                            prevh = img->getHeight();
                            prev = new guint8[prevw * prevh * 3];
                            memcpy(prev, img->getData(), prevw * prevh * 3);
                            img->free();
                        }
                        BatchQueueEntry *entry = new BatchQueueEntry(job, si.params, source, prev, prevw, prevh, thumb);
                        entry->currentSnapshoId = id;
                        entry->setParent(this);
                        entry->resize(options.thumbSize);
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
        }
    }

    redraw();
    notifyListener(false);
}

int cancelItemUI (void* data)
{
    delete (BatchQueueEntry*)data;
    return 0;
}

void BatchQueue::cancelItems (std::vector<ThumbBrowserEntryBase*>* items) {
	{
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
#endif

        for (int i=0; i<items->size(); i++) {
            BatchQueueEntry* entry = (BatchQueueEntry*)(*items)[i];
            if (entry->processing)
                continue;
            std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (fd.begin(), fd.end(), entry);
            if (pos!=fd.end()) {
                fd.erase (pos);
                rtengine::ProcessingJob::destroy (entry->job);
                if (entry->thumbnail)
                    entry->thumbnail->setQueued( entry->currentSnapshoId,false );
                g_idle_add (cancelItemUI, entry);
            }
        }
        for (int i=0; i<fd.size(); i++) 
            fd[i]->selected = false;
        lastClicked = NULL;
        selected.clear ();

        saveBatchQueue( );
    }

    redraw ();
    notifyListener (false);
}

void BatchQueue::headItems (std::vector<ThumbBrowserEntryBase*>* items) {
	{
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
#endif
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
    }

    redraw ();
}

void BatchQueue::tailItems (std::vector<ThumbBrowserEntryBase*>* items) {
	{
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
#endif
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
    }

    redraw ();
}
   
void BatchQueue::selectAll () {
    {
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::ReaderLock l(entryRW);
#endif

        lastClicked = NULL;
        selected.clear ();
        for (int i=0; i<fd.size(); i++) {
            if (fd[i]->processing)
                continue;
            fd[i]->selected = true;
            selected.push_back (fd[i]);
        }
    }
    queue_draw ();
}

void BatchQueue::startProcessing () {
    if (!processing && fd.size()>0) {
        BatchQueueEntry* next;

        {
            // TODO: Check for Linux
	        #ifdef WIN32
	        Glib::RWLock::WriterLock l(entryRW);
	        #endif

            next = (BatchQueueEntry*)fd[0];
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
        }

        // start batch processing
        rtengine::startBatchProcessing (next->job, this );
        queue_draw ();
    }
}

rtengine::ProcessingJob* BatchQueue::imageReady (rtengine::IImage16* img) {

    // save image img
    Glib::ustring fname;
    SaveFormat saveFormat;
    if (processing->outFileName=="") {   // auto file name
        Glib::ustring s = calcAutoFileNameBase ( processing->thumbnail );
        saveFormat = options.saveFormatBatch;
        fname = autoCompleteFileName (s, saveFormat.format);
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

        if (processing->thumbnail) {
            processing->thumbnail->imageDeveloped ();
        	processing->thumbnail->setSaved(processing->currentSnapshoId,true, fname );
        }
    }
    // save temporary params file name: delete as last thing
    //Glib::ustring processedParams = processing->savedParamsFile;
    
    // delete from the queue
    delete processing; processing = NULL;
    bool queueEmptied=false;
	{
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
#endif

        fd.erase (fd.begin());

        // return next job
        if (fd.size()==0) {
            queueEmptied=true;
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
        if (saveBatchQueue( )) {
            //safe_g_remove( processedParams );
            // Delete all files in directory \batch when finished, just to be sure to remove zombies
            /*if( fd.size()==0 ){
                std::vector<Glib::ustring> names;
                Glib::ustring batchdir = options.rtdir+"/batch/";
                Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (batchdir);
                safe_build_file_list (dir, names, batchdir);
                for(std::vector<Glib::ustring>::iterator iter=names.begin(); iter != names.end();iter++ )
                    safe_g_remove( *iter );
            }*/
        }
    }

    redraw ();
    notifyListener (queueEmptied);

    return processing ? processing->job : NULL;
}

// Calculates automatic filename of processed batch entry, but just the base name
// example output: "c:\out\converted\dsc0121"
Glib::ustring BatchQueue::calcAutoFileNameBase ( ::Thumbnail *thumb ) {

    std::vector<Glib::ustring> pa;
    std::vector<Glib::ustring> da;
    Glib::ustring origFileName( thumb->getFileName() );

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
                        path = path + da[da.size()-i-1];
                }
                else if (options.savePathTemplate[ix]=='f') {
                    path = path + filename;
                }
                else if (options.savePathTemplate[ix]=='r') { // rank from pparams

                    int r = thumb->getRank();
                    char rank = r>=0 ? '0'+r:'x';
					/*rtengine::procparams::ProcParams pparams;
					if( pparams.load(origFileName + paramFileExtension)==0 ){
						if (!pparams.inTrash)
							rank = pparams.rank + '0';
						else
							rank = 'x';
					}
					else
						rank = '0'; // if param file not loaded (e.g. does not exist), default to rank=0
						*/

					path += rank;
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

int setProgressUI (void* p) {
    ((BatchQueue*)p)->redraw();
    return 0;
}

void BatchQueue::setProgress (double p) {

    if (processing)
        processing->progress = p;

    g_idle_add (setProgressUI, this);
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
    bool queueEmptied;
};

int bqnotifylistenerUI (void* data) {
    NLParams* params = (NLParams*)data;
    params->listener->queueSizeChanged (params->qsize, params->queueEmptied);
    delete params;
    return 0;
}

void BatchQueue::notifyListener (bool queueEmptied) {

    if (listener) {
        NLParams* params = new NLParams;
        params->listener = listener;
        params->qsize = fd.size();
        params->queueEmptied = queueEmptied;
        g_idle_add (bqnotifylistenerUI, params);
    }
}

void BatchQueue::redrawNeeded (LWButton* button) {
    
    queue_draw ();
}
