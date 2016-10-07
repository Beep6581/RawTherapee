/*
 *  This file is part of RawTherapee.
 *
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
#ifndef _BATCHQUEUE_
#define _BATCHQUEUE_

#include <gtkmm.h>
#include "threadutils.h"
#include "batchqueueentry.h"
#include "../rtengine/rtengine.h"
#include "options.h"
#include "lwbuttonset.h"
#include "thumbbrowserbase.h"

class BatchQueueListener
{

public:
    virtual ~BatchQueueListener () {}
    virtual void queueSizeChanged     (int qsize, bool queueEmptied, bool queueError, Glib::ustring queueErrorMessage) = 0;
    virtual bool canStartNext         () = 0;
};

class FileCatalog;
class BatchQueue  : public ThumbBrowserBase,
    public rtengine::BatchProcessingListener,
    public LWButtonListener
{

protected:
    int getMaxThumbnailHeight() const;
    void saveThumbnailHeight (int height);
    int  getThumbnailHeight ();

    BatchQueueEntry* processing;  // holds the currently processed image
    FileCatalog* fileCatalog;
    int sequence; // holds the current sequence index

    Glib::ustring nameTemplate;

    Gtk::ImageMenuItem* cancel;
    Gtk::ImageMenuItem* head;
    Gtk::ImageMenuItem* tail;
    Gtk::MenuItem* selall;
    Gtk::MenuItem* open;
    Glib::RefPtr<Gtk::AccelGroup> pmaccelgroup;
    Gtk::Menu pmenu;

    BatchQueueListener* listener;

    Glib::ustring autoCompleteFileName (const Glib::ustring& fileName, const Glib::ustring& format);
    Glib::ustring getTempFilenameForParams( const Glib::ustring &filename );
    bool saveBatchQueue ();
    void notifyListener (bool queueEmptied);

public:
    BatchQueue (FileCatalog* aFileCatalog);
    ~BatchQueue ();

    void addEntries (const std::vector<BatchQueueEntry*>& entries, bool head = false, bool save = true);
    void cancelItems (const std::vector<ThumbBrowserEntryBase*>& items);
    void headItems (const std::vector<ThumbBrowserEntryBase *>& items);
    void tailItems (const std::vector<ThumbBrowserEntryBase *>& items);
    void selectAll ();
    void openItemInEditor(ThumbBrowserEntryBase* item);
    void openLastSelectedItemInEditor();

    void startProcessing ();

    bool hasJobs ()
    {
        MYREADERLOCK(l, entryRW);
        return (!fd.empty());
    }

    rtengine::ProcessingJob* imageReady (rtengine::IImage16* img);
    void error (Glib::ustring msg);
    void setProgress (double p);
    void rightClicked (ThumbBrowserEntryBase* entry);
    void doubleClicked (ThumbBrowserEntryBase* entry);
    bool keyPressed (GdkEventKey* event);
    void buttonPressed (LWButton* button, int actionCode, void* actionData);
    void redrawNeeded  (LWButton* button);

    void setBatchQueueListener (BatchQueueListener* l)
    {
        listener = l;
    }

    bool loadBatchQueue ();
    void resizeLoadedQueue();

    static Glib::ustring calcAutoFileNameBase (const Glib::ustring& origFileName, int sequence = 0);
    static int calcMaxThumbnailHeight();
};

#endif
