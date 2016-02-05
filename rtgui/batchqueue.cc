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

#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>

#include "thumbnail.h"
#include "batchqueue.h"
#include "multilangmgr.h"
#include "filecatalog.h"
#include "batchqueuebuttonset.h"
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"

using namespace std;
using namespace rtengine;

BatchQueue::BatchQueue (FileCatalog* aFileCatalog) : processing(NULL), fileCatalog(aFileCatalog), sequence(0), listener(NULL)
{

    location = THLOC_BATCHQUEUE;

    int p = 0;

    pmenu.attach (*Gtk::manage(open = new Gtk::MenuItem (M("FILEBROWSER_POPUPOPENINEDITOR"))), 0, 1, p, p + 1);
    p++;
    pmenu.attach (*Gtk::manage(selall = new Gtk::MenuItem (M("FILEBROWSER_POPUPSELECTALL"))), 0, 1, p, p + 1);
    p++;
    pmenu.attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p + 1);
    p++;

    pmenu.attach (*Gtk::manage(head = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPMOVEHEAD"))), 0, 1, p, p + 1);
    p++;
    head->set_image(*Gtk::manage(new RTImage ("toleftend.png")));

    pmenu.attach (*Gtk::manage(tail = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPMOVEEND"))), 0, 1, p, p + 1);
    p++;
    tail->set_image(*Gtk::manage(new RTImage ("torightend.png")));

    pmenu.attach (*Gtk::manage(new Gtk::SeparatorMenuItem ()), 0, 1, p, p + 1);
    p++;

    pmenu.attach (*Gtk::manage(cancel = new Gtk::ImageMenuItem (M("FILEBROWSER_POPUPCANCELJOB"))), 0, 1, p, p + 1);
    p++;
    cancel->set_image(*Gtk::manage(new RTImage ("gtk-close.png")));

    pmenu.show_all ();

    // Accelerators
    pmaccelgroup = Gtk::AccelGroup::create ();
    pmenu.set_accel_group (pmaccelgroup);
    open->add_accelerator ("activate", pmaccelgroup, GDK_e, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    selall->add_accelerator ("activate", pmaccelgroup, GDK_a, Gdk::CONTROL_MASK, Gtk::ACCEL_VISIBLE);
    head->add_accelerator ("activate", pmaccelgroup, GDK_Home, (Gdk::ModifierType)0, Gtk::ACCEL_VISIBLE);
    tail->add_accelerator ("activate", pmaccelgroup, GDK_End, (Gdk::ModifierType)0, Gtk::ACCEL_VISIBLE);
    cancel->add_accelerator ("activate", pmaccelgroup, GDK_Delete, (Gdk::ModifierType)0, Gtk::ACCEL_VISIBLE);

    open->signal_activate().connect(sigc::mem_fun(*this, &BatchQueue::openLastSelectedItemInEditor));
    cancel->signal_activate().connect (std::bind (&BatchQueue::cancelItems, this, std::ref (selected)));
    head->signal_activate().connect (std::bind (&BatchQueue::headItems, this, std::ref (selected)));
    tail->signal_activate().connect (std::bind (&BatchQueue::tailItems, this, std::ref (selected)));
    selall->signal_activate().connect (sigc::mem_fun(*this, &BatchQueue::selectAll));

    setArrangement (ThumbBrowserBase::TB_Vertical);
}

BatchQueue::~BatchQueue ()
{
    MYWRITERLOCK(l, entryRW);

    // The listener merges parameters with old values, so delete afterwards
    for (size_t i = 0; i < fd.size(); i++) {
        delete fd.at(i);
    }

    fd.clear ();
}

void BatchQueue::resizeLoadedQueue()
{
    MYWRITERLOCK(l, entryRW);

    const auto height = getThumbnailHeight ();

    for (const auto entry : fd)
        entry->resize(height);
}

// Reduce the max size of a thumb, since thumb is processed synchronously on adding to queue
// leading to very long waiting when adding more images
int BatchQueue::calcMaxThumbnailHeight()
{
    return std::min(options.maxThumbnailHeight, 200);
}

// Function for virtual override in thumbbrowser base
int BatchQueue::getMaxThumbnailHeight() const
{
    return calcMaxThumbnailHeight();
}

void BatchQueue::saveThumbnailHeight (int height)
{
    options.thumbSizeQueue = height;
}

int BatchQueue::getThumbnailHeight ()
{
    // The user could have manually forced the option to a too big value
    return std::max(std::min(options.thumbSizeQueue, 200), 10);
}

void BatchQueue::rightClicked (ThumbBrowserEntryBase* entry)
{
    pmenu.popup (3, this->eventTime);
}

void BatchQueue::doubleClicked(ThumbBrowserEntryBase* entry)
{
    openItemInEditor(entry);
}

bool BatchQueue::keyPressed (GdkEventKey* event)
{
    bool ctrl  = event->state & GDK_CONTROL_MASK;

    if ((event->keyval == GDK_A || event->keyval == GDK_a) && ctrl) {
        selectAll ();
        return true;
    } else if ((event->keyval == GDK_E || event->keyval == GDK_e) && ctrl) {
        openLastSelectedItemInEditor();
        return true;
    } else if (event->keyval == GDK_Home) {
        headItems (selected);
        return true;
    } else if (event->keyval == GDK_End) {
        tailItems (selected);
        return true;
    } else if (event->keyval == GDK_Delete) {
        cancelItems (selected);
        return true;
    }

    return false;
}

void BatchQueue::addEntries (const std::vector<BatchQueueEntry*>& entries, bool head, bool save)
{
    {
        MYWRITERLOCK(l, entryRW);

        for (const auto entry : entries) {

            entry->setParent (this);

            // BatchQueueButtonSet have to be added before resizing to take them into account
            const auto bqbs = new BatchQueueButtonSet (entry);
            bqbs->setButtonListener (this);
            entry->addButtonSet (bqbs);

            // batch queue might have smaller, restricted size
            entry->resize (getThumbnailHeight());

            // recovery save
            const auto tempFile = getTempFilenameForParams (entry->filename);

            if (!entry->params.save (tempFile))
                entry->savedParamsFile = tempFile;

            entry->selected = false;

            // insert either at the end, or before the first non-processing entry
            auto pos = fd.end ();

            if (head)
                pos = std::find_if (fd.begin (), fd.end (), [] (const ThumbBrowserEntryBase* fdEntry) { return !fdEntry->processing; });

            fd.insert (pos, entry);

            if (entry->thumbnail)
                entry->thumbnail->imageEnqueued ();
        }
    }

    if (save)
        saveBatchQueue ();

    redraw ();
    notifyListener (false);
}

bool BatchQueue::saveBatchQueue ()
{
    const auto fileName = Glib::build_filename (options.rtdir, "batch", "queue.csv");

    std::ofstream file (fileName, std::ios::binary | std::ios::trunc);

    if (!file.is_open ())
        return false;

    {
        MYREADERLOCK(l, entryRW);

        if (fd.empty ())
            return true;

        // The column's header is mandatory (the first line will be skipped when loaded)
        file << "input image full path|param file full path|output image full path|file format|jpeg quality|jpeg subsampling|"
             << "png bit depth|png compression|tiff bit depth|uncompressed tiff|save output params|force format options|<end of line>"
             << std::endl;

        // method is already running with entryLock, so no need to lock again
        for (const auto fdEntry : fd) {

            const auto entry = static_cast<BatchQueueEntry*> (fdEntry);
            const auto& saveFormat = entry->saveFormat;

            // Warning: for code's simplicity in loadBatchQueue, each field must end by the '|' character, safer than ';' or ',' since it can't be used in paths
            file << entry->filename << '|' << entry->savedParamsFile << '|' << entry->outFileName << '|' << saveFormat.format << '|'
                 << saveFormat.jpegQuality << '|' << saveFormat.jpegSubSamp << '|'
                 << saveFormat.pngBits << '|' << saveFormat.pngCompression << '|'
                 << saveFormat.tiffBits << '|'  << saveFormat.tiffUncompressed << '|'
                 << saveFormat.saveParams << '|' << entry->forceFormatOpts << '|'
                 << std::endl;
        }
    }

    return true;
}

bool BatchQueue::loadBatchQueue ()
{
    const auto fileName = Glib::build_filename (options.rtdir, "batch", "queue.csv");

    std::ifstream file (fileName, std::ios::binary);

    if (file.is_open ()) {
        // Yes, it's better to get the lock for the whole file reading,
        // to update the list in one shot without any other concurrent access!
        MYWRITERLOCK(l, entryRW);

        std::string row, column;
        std::vector<std::string> values;

        // skipping the first row
        std::getline (file, row);

        while (std::getline (file, row)) {

            std::istringstream line (row);

            values.clear ();

            while (std::getline(line, column, '|')) {
                values.push_back (column);
            }

            auto value = values.begin ();

            const auto nextStringOr = [&] (const Glib::ustring& defaultValue) -> Glib::ustring
            {
                return value != values.end () ? Glib::ustring(*value++) : defaultValue;
            };
            const auto nextIntOr = [&] (int defaultValue) -> int
            {
                try {
                    return value != values.end () ? std::stoi(*value++) : defaultValue;
                }
                catch (std::exception&) {
                    return defaultValue;
                }
            };

            const auto source = nextStringOr (Glib::ustring ());
            const auto paramsFile = nextStringOr (Glib::ustring ());

            if (source.empty () || paramsFile.empty ())
                continue;

            const auto outputFile = nextStringOr (Glib::ustring ());
            const auto saveFmt = nextStringOr (options.saveFormat.format);
            const auto jpegQuality = nextIntOr (options.saveFormat.jpegQuality);
            const auto jpegSubSamp = nextIntOr (options.saveFormat.jpegSubSamp);
            const auto pngBits = nextIntOr (options.saveFormat.pngBits);
            const auto pngCompression = nextIntOr (options.saveFormat.pngCompression);
            const auto tiffBits = nextIntOr (options.saveFormat.tiffBits);
            const auto tiffUncompressed = nextIntOr (options.saveFormat.tiffUncompressed);
            const auto saveParams = nextIntOr (options.saveFormat.saveParams);
            const auto forceFormatOpts = nextIntOr (options.forceFormatOpts);

            rtengine::procparams::ProcParams pparams;

            if (pparams.load (paramsFile))
                continue;

            auto thumb = CacheManager::getInstance ()->getEntry (source);

            if (!thumb)
                continue;

            auto job = rtengine::ProcessingJob::create (source, thumb->getType () == FT_Raw, pparams);

            auto prevh = getMaxThumbnailHeight ();
            auto prevw = prevh;
            thumb->getThumbnailSize (prevw, prevh, &pparams);

            auto entry = new BatchQueueEntry (job, pparams, source, prevw, prevh, thumb);
            thumb->decreaseRef ();  // Removing the refCount acquired by cacheMgr->getEntry
            entry->setParent (this);

            // BatchQueueButtonSet have to be added before resizing to take them into account
            auto bqbs = new BatchQueueButtonSet (entry);
            bqbs->setButtonListener (this);
            entry->addButtonSet (bqbs);

            entry->savedParamsFile = paramsFile;
            entry->selected = false;
            entry->outFileName = outputFile;

            if (!outputFile.empty ()) {
                auto& saveFormat = entry->saveFormat;
                saveFormat.format = saveFmt;
                saveFormat.jpegQuality = jpegQuality;
                saveFormat.jpegSubSamp = jpegSubSamp;
                saveFormat.pngBits = pngBits;
                saveFormat.pngCompression = pngCompression;
                saveFormat.tiffBits = tiffBits;
                saveFormat.tiffUncompressed = tiffUncompressed != 0;
                saveFormat.saveParams = saveParams != 0;
                entry->forceFormatOpts = forceFormatOpts != 0;
            } else {
                entry->forceFormatOpts = false;
            }

            fd.push_back (entry);
        }
    }

    redraw ();
    notifyListener (false);

    return !fd.empty ();
}

Glib::ustring BatchQueue::getTempFilenameForParams( const Glib::ustring filename )
{
    time_t rawtime;
    struct tm *timeinfo;
    char stringTimestamp [80];
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (stringTimestamp, sizeof(stringTimestamp), "_%Y%m%d%H%M%S_", timeinfo);
    Glib::ustring savedParamPath;
    savedParamPath = options.rtdir + "/batch/";
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

void BatchQueue::cancelItems (const std::vector<ThumbBrowserEntryBase*>& items)
{
    {
        MYWRITERLOCK(l, entryRW);

        for (const auto item : items) {

            const auto entry = static_cast<BatchQueueEntry*> (item);

            if (entry->processing)
                continue;

            const auto pos = std::find (fd.begin (), fd.end (), entry);

            if (pos == fd.end ())
                continue;

            fd.erase (pos);

            rtengine::ProcessingJob::destroy (entry->job);

            if (entry->thumbnail)
                entry->thumbnail->imageRemovedFromQueue ();

            g_idle_add (cancelItemUI, entry);
        }

        for (const auto entry : fd)
            entry->selected = false;

        lastClicked = nullptr;
        selected.clear ();
    }

    saveBatchQueue ();

    redraw ();
    notifyListener (false);
}

void BatchQueue::headItems (const std::vector<ThumbBrowserEntryBase*>& items)
{
    {
        MYWRITERLOCK(l, entryRW);

        for (auto item = items.rbegin(); item != items.rend(); ++item) {

            const auto entry = static_cast<BatchQueueEntry*> (*item);

            if (entry->processing)
                continue;

            const auto pos = std::find (fd.begin (), fd.end (), entry);

            if (pos == fd.end () || pos == fd.begin ())
                continue;

            fd.erase (pos);

            // find the first item that is not under processing
            const auto newPos = std::find_if (fd.begin (), fd.end (), [] (const ThumbBrowserEntryBase* fdEntry) { return !fdEntry->processing; });

            fd.insert (newPos, entry);
        }
    }

    saveBatchQueue ();

    redraw ();
}

void BatchQueue::tailItems (const std::vector<ThumbBrowserEntryBase*>& items)
{
    {
        MYWRITERLOCK(l, entryRW);

        for (const auto item : items) {

            const auto entry = static_cast<BatchQueueEntry*> (item);

            if (entry->processing)
                continue;

            const auto pos = std::find (fd.begin (), fd.end (), entry);

            if (pos == fd.end ())
                continue;

            fd.erase (pos);

            fd.push_back (entry);
        }
    }

    saveBatchQueue ();

    redraw ();
}

void BatchQueue::selectAll ()
{
    {
        MYWRITERLOCK(l, entryRW);

        lastClicked = NULL;
        selected.clear ();

        for (size_t i = 0; i < fd.size(); i++) {
            if (fd[i]->processing) {
                continue;
            }

            fd[i]->selected = true;
            selected.push_back (fd[i]);
        }
    }
    queue_draw ();
}

void BatchQueue::openLastSelectedItemInEditor()
{
    {
        MYREADERLOCK(l, entryRW);

        if (selected.size() > 0) {
            openItemInEditor(selected.back());
        }
    }
}

void BatchQueue::openItemInEditor(ThumbBrowserEntryBase* item)
{
    if (item) {
        std::vector< ::Thumbnail*> requestedItem;
        requestedItem.push_back(item->thumbnail);
        fileCatalog->openRequested(requestedItem);
    }
}


void BatchQueue::startProcessing ()
{

    if (!processing) {
        MYWRITERLOCK(l, entryRW);

        if (!fd.empty()) {
            BatchQueueEntry* next;

            next = static_cast<BatchQueueEntry*>(fd[0]);
            // tag it as processing and set sequence
            next->processing = true;
            next->sequence = sequence = 1;
            processing = next;

            // remove from selection
            if (processing->selected) {
                std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (selected.begin(), selected.end(), processing);

                if (pos != selected.end()) {
                    selected.erase (pos);
                }

                processing->selected = false;
            }

            MYWRITERLOCK_RELEASE(l);

            // remove button set
            next->removeButtonSet ();

            // start batch processing
            rtengine::startBatchProcessing (next->job, this, options.tunnelMetaData);
            queue_draw ();
        }
    }
}

rtengine::ProcessingJob* BatchQueue::imageReady (rtengine::IImage16* img)
{

    // save image img
    Glib::ustring fname;
    SaveFormat saveFormat;

    if (processing->outFileName == "") { // auto file name
        Glib::ustring s = calcAutoFileNameBase (processing->filename, processing->sequence);
        saveFormat = options.saveFormatBatch;
        fname = autoCompleteFileName (s, saveFormat.format);
    } else { // use the save-as filename with automatic completion for uniqueness
        if (processing->forceFormatOpts) {
            saveFormat = processing->saveFormat;
        } else {
            saveFormat = options.saveFormatBatch;
        }

        // The output filename's extension is forced to the current or selected output format,
        // despite what the user have set in the fielneame's field of the "Save as" dialgo box
        fname = autoCompleteFileName (removeExtension(processing->outFileName), saveFormat.format);
        //fname = autoCompleteFileName (removeExtension(processing->outFileName), getExtension(processing->outFileName));
    }

    //printf ("fname=%s, %s\n", fname.c_str(), removeExtension(fname).c_str());

    if (img && fname != "") {
        int err = 0;

        if (saveFormat.format == "tif") {
            err = img->saveAsTIFF (fname, saveFormat.tiffBits, saveFormat.tiffUncompressed);
        } else if (saveFormat.format == "png") {
            err = img->saveAsPNG (fname, saveFormat.pngCompression, saveFormat.pngBits);
        } else if (saveFormat.format == "jpg") {
            err = img->saveAsJPEG (fname, saveFormat.jpegQuality, saveFormat.jpegSubSamp);
        }

        img->free ();

        if (err) {
            throw Glib::FileError(Glib::FileError::FAILED, M("MAIN_MSG_CANNOTSAVE") + "\n" + fname);
        }

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
    bool queueEmptied = false;
    bool remove_button_set = false;

    {
        MYWRITERLOCK(l, entryRW);

        delete processing;
        processing = NULL;

        fd.erase (fd.begin());

        // return next job
        if (fd.empty()) {
            queueEmptied = true;
        } else if (listener && listener->canStartNext ()) {
            BatchQueueEntry* next = static_cast<BatchQueueEntry*>(fd[0]);
            // tag it as selected and set sequence
            next->processing = true;
            next->sequence = ++sequence;
            processing = next;

            // remove from selection
            if (processing->selected) {
                std::vector<ThumbBrowserEntryBase*>::iterator pos = std::find (selected.begin(), selected.end(), processing);

                if (pos != selected.end()) {
                    selected.erase (pos);
                }

                processing->selected = false;
            }

            // remove button set
            remove_button_set = true;
        }
    }

    if (remove_button_set) {
        // ButtonSet have Cairo::Surface which might be rendered while we're trying to delete them
        GThreadLock lock;
        processing->removeButtonSet ();
    }

    if (saveBatchQueue( )) {
        safe_g_remove( processedParams );
        // Delete all files in directory \batch when finished, just to be sure to remove zombies

        // Not sure that locking is necessary, but it should be safer
        MYREADERLOCK(l, entryRW);

        if( fd.empty() ) {
            MYREADERLOCK_RELEASE(l);

            std::vector<Glib::ustring> names;
            Glib::ustring batchdir = Glib::build_filename(options.rtdir, "batch");
            Glib::RefPtr<Gio::File> dir = Gio::File::create_for_path (batchdir);
            safe_build_file_list (dir, names, batchdir);

            for(std::vector<Glib::ustring>::iterator iter = names.begin(); iter != names.end(); iter++ ) {
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
Glib::ustring BatchQueue::calcAutoFileNameBase (const Glib::ustring& origFileName, int sequence)
{

    std::vector<Glib::ustring> pa;
    std::vector<Glib::ustring> da;

    for (size_t i = 0; i < origFileName.size(); i++) {
        while ((i < origFileName.size()) && (origFileName[i] == '\\' || origFileName[i] == '/')) {
            i++;
        }

        if (i >= origFileName.size()) {
            break;
        }

        Glib::ustring tok = "";

        while ((i < origFileName.size()) && !(origFileName[i] == '\\' || origFileName[i] == '/')) {
            tok = tok + origFileName[i++];
        }

        da.push_back (tok);
    }

    if (origFileName[0] == '/' || origFileName[0] == '\\') {
        pa.push_back ("/" + da[0]);
    } else {
        pa.push_back (da[0]);
    }

    for (size_t i = 1; i < da.size(); i++) {
        pa.push_back (pa[i - 1] + "/" + da[i]);
    }

//    for (int i=0; i<da.size(); i++)
//        printf ("da: %s\n", da[i].c_str());
//    for (int i=0; i<pa.size(); i++)
//        printf ("pa: %s\n", pa[i].c_str());

    // extracting filebase
    Glib::ustring filename;

    int extpos = origFileName.size() - 1;

    for (; extpos >= 0 && origFileName[extpos] != '.'; extpos--);

    for (int k = extpos - 1; k >= 0 && origFileName[k] != '/' && origFileName[k] != '\\'; k--) {
        filename = origFileName[k] + filename;
    }

//    printf ("%d, |%s|\n", extpos, filename.c_str());

    // constructing full output path
//    printf ("path=|%s|\n", options.savePath.c_str());

    Glib::ustring path = "";

    if (options.saveUsePathTemplate) {
        int ix = 0;

        while (options.savePathTemplate[ix] != 0) {
            if (options.savePathTemplate[ix] == '%') {
                ix++;

                if (options.savePathTemplate[ix] == 'p') {
                    ix++;
                    int i = options.savePathTemplate[ix] - '0';

                    if (i < pa.size()) {
                        path = path + pa[pa.size() - i - 1] + '/';
                    }

                    ix++;
                } else if (options.savePathTemplate[ix] == 'd') {
                    ix++;
                    int i = options.savePathTemplate[ix] - '0';

                    if (i < da.size()) {
                        path = path + da[da.size() - i - 1];
                    }
                } else if (options.savePathTemplate[ix] == 'f') {
                    path = path + filename;
                } else if (options.savePathTemplate[ix] == 'r') { // rank from pparams
                    char rank;
                    rtengine::procparams::ProcParams pparams;

                    if( pparams.load(origFileName + paramFileExtension) == 0 ) {
                        if (!pparams.inTrash) {
                            rank = pparams.rank + '0';
                        } else {
                            rank = 'x';
                        }
                    } else {
                        rank = '0';    // if param file not loaded (e.g. does not exist), default to rank=0
                    }

                    path += rank;
                } else if (options.savePathTemplate[ix] == 's') { // sequence
                    std::ostringstream seqstr;

                    int w = options.savePathTemplate[ix + 1] - '0';

                    if (w >= 1 && w <= 9) {
                        ix++;
                        seqstr << std::setw (w) << std::setfill ('0');
                    }

                    seqstr << sequence;
                    path += seqstr.str ();
                }
            }

            else {
                path = path + options.savePathTemplate[ix];
            }

            ix++;
        }
    } else {
        path = Glib::build_filename (options.savePathFolder, filename);
    }

    return path;
}

Glib::ustring BatchQueue::autoCompleteFileName (const Glib::ustring& fileName, const Glib::ustring& format)
{

    // separate filename and the path to the destination directory
    Glib::ustring dstdir = Glib::path_get_dirname (fileName);
    Glib::ustring dstfname = Glib::path_get_basename (fileName);
    Glib::ustring fname;

    // create directory, if does not exist
    if (safe_g_mkdir_with_parents (dstdir, 0755) ) {
        return "";
    }

    // In overwrite mode we TRY to delete the old file first.
    // if that's not possible (e.g. locked by viewer, R/O), we revert to the standard naming scheme
    bool inOverwriteMode = options.overwriteOutputFile;

    for (int tries = 0; tries < 100; tries++) {
        if (tries == 0) {
            fname = Glib::ustring::compose ("%1.%2", Glib::build_filename (dstdir,  dstfname), format);
        } else {
            fname = Glib::ustring::compose ("%1-%2.%3", Glib::build_filename (dstdir,  dstfname), tries, format);
        }

        int fileExists = safe_file_test (fname, Glib::FILE_TEST_EXISTS);

        if (inOverwriteMode && fileExists) {
            if (safe_g_remove(fname) == -1) {
                inOverwriteMode = false;    // failed to delete- revert to old naming scheme
            } else {
                fileExists = false;    // deleted now
            }
        }

        if (!fileExists) {
            return fname;
        }
    }

    return "";
}

int setProgressUI (void* p)
{
    (static_cast<BatchQueue*>(p))->redraw();
    return 0;
}

void BatchQueue::setProgress (double p)
{

    if (processing) {
        processing->progress = p;
    }

    // No need to acquire the GUI, setProgressUI will do it
    g_idle_add (setProgressUI, this);
}

void BatchQueue::buttonPressed (LWButton* button, int actionCode, void* actionData)
{

    std::vector<ThumbBrowserEntryBase*> bqe;
    bqe.push_back (static_cast<BatchQueueEntry*>(actionData));

    if (actionCode == 10) { // cancel
        cancelItems (bqe);
    } else if (actionCode == 8) { // to head
        headItems (bqe);
    } else if (actionCode == 9) { // to tail
        tailItems (bqe);
    }
}

struct NLParams {
    BatchQueueListener* listener;
    int qsize;
    bool queueEmptied;
    bool queueError;
    Glib::ustring queueErrorMessage;
};

int bqnotifylistenerUI (void* data)
{
    GThreadLock lock; // All GUI acces from idle_add callbacks or separate thread HAVE to be protected
    NLParams* params = static_cast<NLParams*>(data);
    params->listener->queueSizeChanged (params->qsize, params->queueEmptied, params->queueError, params->queueErrorMessage);
    delete params;
    return 0;
}

void BatchQueue::notifyListener (bool queueEmptied)
{

    if (listener) {
        NLParams* params = new NLParams;
        params->listener = listener;
        {
            MYREADERLOCK(l, entryRW);
            params->qsize = fd.size();
        }
        params->queueEmptied = queueEmptied;
        params->queueError = false;
        g_idle_add (bqnotifylistenerUI, params);
    }
}

void BatchQueue::redrawNeeded (LWButton* button)
{
    GThreadLock lock;
    queue_draw ();
}

void BatchQueue::error (Glib::ustring msg)
{

    if (processing && processing->processing) {
        // restore failed thumb
        BatchQueueButtonSet* bqbs = new BatchQueueButtonSet (processing);
        bqbs->setButtonListener (this);
        processing->addButtonSet (bqbs);
        processing->processing = false;
        processing->job = rtengine::ProcessingJob::create(processing->filename, processing->thumbnail->getType() == FT_Raw, processing->params);
        processing = NULL;
        redraw ();
    }

    if (listener) {
        NLParams* params = new NLParams;
        params->listener = listener;
        params->queueEmptied = false;
        params->queueError = true;
        params->queueErrorMessage = msg;
        g_idle_add (bqnotifylistenerUI, params);
    }
}
