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
#include <glibmm.h>
#include <glib/gstdio.h>
#include <cstring>
#include "../rtengine/rt_math.h"

#include "batchqueue.h"
#include "multilangmgr.h"
#include "filecatalog.h"
#include "batchqueuebuttonset.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"

using namespace std;
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
    return std::min(options.maxThumbnailHeight, 200);
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
		(*entry)->resize (std::min(options.thumbSize, getMaxThumbnailHeight()));  // batch queue might have smaller, restricted size
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
	}

    redraw();
    notifyListener (false);
}

bool BatchQueue::saveBatchQueue( )
{
    Glib::ustring savedQueueFile;
    savedQueueFile = options.rtdir+"/batch/queue.csv";
    FILE *f = safe_g_fopen (savedQueueFile, "wt");

    if (f==NULL)
        return false;

    if (fd.size())
        // The column's header is mandatory (the first line will be skipped when loaded)
        fprintf(f,"input image full path|param file full path|output image full path|file format|jpeg quality|jpeg subsampling|"
                  "png bit depth|png compression|tiff bit depth|uncompressed tiff|save output params|force format options|<end of line>\n");

    // method is already running with entryLock, so no need to lock again
    for (std::vector<ThumbBrowserEntryBase*>::iterator pos=fd.begin(); pos!=fd.end(); pos++){
        BatchQueueEntry* bqe = reinterpret_cast<BatchQueueEntry*>(*pos);
        // Warning: for code's simplicity in loadBatchQueue, each field must end by the '|' character, safer than ';' or ',' since it can't be used in paths
        fprintf(f,"%s|%s|%s|%s|%d|%d|%d|%d|%d|%d|%d|%d|\n",
            bqe->filename.c_str(),bqe->savedParamsFile.c_str(), bqe->outFileName.c_str(), bqe->saveFormat.format.c_str(),
            bqe->saveFormat.jpegQuality, bqe->saveFormat.jpegSubSamp,
            bqe->saveFormat.pngBits, bqe->saveFormat.pngCompression,
            bqe->saveFormat.tiffBits, bqe->saveFormat.tiffUncompressed,
            bqe->saveFormat.saveParams, bqe->forceFormatOpts
        );
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
        savedQueueFile = options.rtdir+"/batch/queue.csv";
        FILE *f = safe_g_fopen (savedQueueFile, "rt");

        if (f!=NULL) {
            char *buffer = new char[1024];
            unsigned numLoaded=0;
            // skipping the first line
            bool firstLine=true;
            while (fgets (buffer, 1024, f)){

                if (firstLine) {
                    // skipping the column's title line
                    firstLine=false;
                    continue;
                }

                size_t pos;
                Glib::ustring source;
                Glib::ustring paramsFile;
                Glib::ustring outputFile;
                Glib::ustring saveFmt(options.saveFormat.format);
                int jpegQuality=options.saveFormat.jpegQuality, jpegSubSamp     =options.saveFormat.jpegSubSamp;
                int pngBits    =options.saveFormat.pngBits,     pngCompression  =options.saveFormat.pngCompression;
                int tiffBits   =options.saveFormat.tiffBits,    tiffUncompressed=options.saveFormat.tiffUncompressed;
                int saveParams =options.saveFormat.saveParams;
                int forceFormatOpts =options.forceFormatOpts;

                Glib::ustring currLine(buffer);
                int a = 0;
                if (currLine.rfind('\n') != Glib::ustring::npos) a++;
                if (currLine.rfind('\r') != Glib::ustring::npos) a++;
                if (a)
                    currLine = currLine.substr(0, currLine.length()-a);

                // Looking for the image's full path
                pos = currLine.find('|');
                if (pos != Glib::ustring::npos) {
                source = currLine.substr(0, pos);
                currLine = currLine.substr(pos+1);

                // Looking for the procparams' full path
                pos = currLine.find('|');
                if (pos != Glib::ustring::npos) {
                paramsFile = currLine.substr(0, pos);
                currLine = currLine.substr(pos+1);

                // Looking for the full output path; if empty, it'll use the template string
                pos = currLine.find('|');
                if (pos != Glib::ustring::npos) {
                    outputFile = currLine.substr(0, pos);
                    currLine = currLine.substr(pos+1);

                // No need to bother reading the last options, they will be ignored if outputFile is empty!
                if (!outputFile.empty()) {

                    // Looking for the saving format
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        saveFmt = currLine.substr(0, pos);
                        currLine = currLine.substr(pos+1);

                    // Looking for the jpeg quality
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        jpegQuality = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking for the jpeg subsampling
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        jpegSubSamp = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking for the png bit depth
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        pngBits = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking for the png compression
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        pngCompression = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking for the tiff bit depth
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        tiffBits = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking for the tiff uncompression
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        tiffUncompressed = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking out if we have to save the procparams
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        saveParams = atoi(currLine.substr(0, pos).c_str());
                        currLine = currLine.substr(pos+1);

                    // Looking out if we have to to use the format options
                    pos = currLine.find('|');
                    if (pos != Glib::ustring::npos) {
                        forceFormatOpts = atoi(currLine.substr(0, pos).c_str());
                        // currLine = currLine.substr(pos+1);

                }}}}}}}}}}}}}

                if( !source.empty() && !paramsFile.empty() ){
                    rtengine::procparams::ProcParams pparams;
                    if( pparams.load( paramsFile ) )
                        continue;

                    ::Thumbnail *thumb = cacheMgr->getEntry( source );
                    if( thumb ){
                        rtengine::ProcessingJob* job = rtengine::ProcessingJob::create(source, thumb->getType() == FT_Raw, pparams);

                        int prevh = getMaxThumbnailHeight();
                        int prevw = prevh;
                        guint8* prev = NULL;
                        double tmpscale;
                        rtengine::IImage8* img = thumb->processThumbImage(pparams, prevh, tmpscale);
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
                        entry->outFileName = outputFile;
                        if (!outputFile.empty()) {
                            entry->saveFormat.format = saveFmt;
                            entry->saveFormat.jpegQuality = jpegQuality;
                            entry->saveFormat.jpegSubSamp = jpegSubSamp;
                            entry->saveFormat.pngBits = pngBits;
                            entry->saveFormat.pngCompression = pngCompression;
                            entry->saveFormat.tiffBits = tiffBits;
                            entry->saveFormat.tiffUncompressed = tiffUncompressed!=0;
                            entry->saveFormat.saveParams = saveParams!=0;
                            entry->forceFormatOpts = forceFormatOpts!=0;
                        }
                        else
                            entry->forceFormatOpts = false;
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

int cancelItemUI (void* data)
{
    safe_g_remove( (static_cast<BatchQueueEntry*>(data))->savedParamsFile );
    delete static_cast<BatchQueueEntry*>(data);
    return 0;
}

void BatchQueue::cancelItems (std::vector<ThumbBrowserEntryBase*>* items) {
	{
        // TODO: Check for Linux
#ifdef WIN32
        Glib::RWLock::WriterLock l(entryRW);
#endif

	for (size_t i=0; i<items->size(); i++) {
            BatchQueueEntry* entry = (BatchQueueEntry*)(*items)[i];
            if (entry->processing)
                continue;
            std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (fd.begin(), fd.end(), entry);
            if (pos!=fd.end()) {
                fd.erase (pos);
                rtengine::ProcessingJob::destroy (entry->job);
                if (entry->thumbnail)
                    entry->thumbnail->imageRemovedFromQueue ();
                g_idle_add (cancelItemUI, entry);
            }
        }
        for (size_t i=0; i<fd.size(); i++)
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
	for (size_t i=0; i<items->size(); i++) {
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
	for (size_t i=0; i<fd.size(); i++) {
            if (fd[i]->processing)
                continue;
            fd[i]->selected = true;
            selected.push_back (fd[i]);
        }
    }
    queue_draw ();
}

void BatchQueue::startProcessing () {
    if (!processing && !fd.empty()) {
        BatchQueueEntry* next;

        {
            // TODO: Check for Linux
	        #ifdef WIN32
	        Glib::RWLock::WriterLock l(entryRW);
	        #endif

	    next = static_cast<BatchQueueEntry*>(fd[0]);
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
        rtengine::startBatchProcessing (next->job, this, options.tunnelMetaData);
        queue_draw ();
    }
}

rtengine::ProcessingJob* BatchQueue::imageReady (rtengine::IImage16* img) {

    // save image img
    Glib::ustring fname;
    SaveFormat saveFormat;
    if (processing->outFileName=="") {   // auto file name
        Glib::ustring s = calcAutoFileNameBase (processing->filename);
        saveFormat = options.saveFormatBatch;
        fname = autoCompleteFileName (s, saveFormat.format);
    }
    else {  // use the save-as filename with automatic completion for uniqueness
        if (processing->forceFormatOpts)
            saveFormat = processing->saveFormat;
        else
            saveFormat = options.saveFormatBatch;
        // The output filename's extension is forced to the current or selected output format,
        // despite what the user have set in the fielneame's field of the "Save as" dialgo box
        fname = autoCompleteFileName (removeExtension(processing->outFileName), saveFormat.format);
        //fname = autoCompleteFileName (removeExtension(processing->outFileName), getExtension(processing->outFileName));
    }
    //printf ("fname=%s, %s\n", fname.c_str(), removeExtension(fname).c_str());

    if (img && fname!="") {
        int err = 0;
        if (saveFormat.format=="tif")
            err = img->saveAsTIFF (fname, saveFormat.tiffBits,saveFormat.tiffUncompressed);
        else if (saveFormat.format=="png")
            err = img->saveAsPNG (fname, saveFormat.pngCompression, saveFormat.pngBits);
        else if (saveFormat.format=="jpg")
            err = img->saveAsJPEG (fname, saveFormat.jpegQuality, saveFormat.jpegSubSamp);
        img->free ();

		if (err) throw "Unable to save output file";

        if (saveFormat.saveParams) {
			// We keep the extension to avoid overwriting the profile when we have
			// the same output filename with different extension
            //processing->params.save (removeExtension(fname) + paramFileExtension);
            processing->params.save (fname + ".out" + paramFileExtension);
        }

        if (processing->thumbnail) {
            processing->thumbnail->imageDeveloped ();
            processing->thumbnail->imageRemovedFromQueue ();
        }
    }
    // save temporary params file name: delete as last thing
    Glib::ustring processedParams = processing->savedParamsFile;
    
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
        if (fd.empty()) {
            queueEmptied=true;
        }
        else if (listener && listener->canStartNext ()) {
	    BatchQueueEntry* next = static_cast<BatchQueueEntry*>(fd[0]);
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
            safe_g_remove( processedParams );
            // Delete all files in directory \batch when finished, just to be sure to remove zombies
            if( fd.empty() ){
                std::vector<Glib::ustring> names;
                Glib::ustring batchdir = options.rtdir+"/batch/";
                Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (batchdir);
                safe_build_file_list (dir, names, batchdir);
                for(std::vector<Glib::ustring>::iterator iter=names.begin(); iter != names.end();iter++ )
                    safe_g_remove( *iter );
            }
        }
    }

    redraw ();
    notifyListener (queueEmptied);

    return processing ? processing->job : NULL;
}

// Calculates automatic filename of processed batch entry, but just the base name
// example output: "c:\out\converted\dsc0121"
Glib::ustring BatchQueue::calcAutoFileNameBase (const Glib::ustring& origFileName) {

    std::vector<Glib::ustring> pa;
    std::vector<Glib::ustring> da;

    for (size_t i=0; i<origFileName.size(); i++) {
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
        
    for (size_t i=1; i<da.size(); i++)
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
                    char rank;
					rtengine::procparams::ProcParams pparams;
					if( pparams.load(origFileName + paramFileExtension)==0 ){
						if (!pparams.inTrash)
							rank = pparams.rank + '0';
						else
							rank = 'x';
					}
					else
						rank = '0'; // if param file not loaded (e.g. does not exist), default to rank=0
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
    Glib::ustring fname;

    // create directory, if does not exist
    if (safe_g_mkdir_with_parents (dstdir, 0755) ) 
        return "";
	
	// In overwrite mode we TRY to delete the old file first.
	// if that's not possible (e.g. locked by viewer, R/O), we revert to the standard naming scheme
	bool inOverwriteMode=options.overwriteOutputFile;

    for (int tries=0; tries<100; tries++) {
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
    return "";
}

int setProgressUI (void* p) {
  (static_cast<BatchQueue*>(p))->redraw();
    return 0;
}

void BatchQueue::setProgress (double p) {

    if (processing)
        processing->progress = p;

    g_idle_add (setProgressUI, this);
}

void BatchQueue::buttonPressed (LWButton* button, int actionCode, void* actionData) {
    
    std::vector<ThumbBrowserEntryBase*> bqe;
    bqe.push_back (static_cast<BatchQueueEntry*>(actionData));

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
    NLParams* params = static_cast<NLParams*>(data);
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
