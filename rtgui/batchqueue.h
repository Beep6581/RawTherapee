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
#ifndef _BATCHQUEUE_
#define _BATCHQUEUE_

#include <gtkmm.h>
#include <batchqueueentry.h>
#include <rtengine.h>
#include <options.h>
#include <lwbuttonset.h>
#include <thumbbrowserbase.h>

class BatchQueueListener {

    public:
        virtual void queueSizeChanged     (int qsize) {}
        virtual void imageProcessingReady (Glib::ustring fname) {}
        virtual void queueEmpty           () {}
        virtual bool canStartNext         () {}
};

class FileCatalog;
class BatchQueue  : public ThumbBrowserBase, 
                    public rtengine::BatchProcessingListener, 
                    public LWButtonListener {  

  protected:

    BatchQueueEntry* processing;

    Glib::ustring nameTemplate;
    
    Gtk::MenuItem* cancel;
    Gtk::MenuItem* head;
    Gtk::MenuItem* tail;
    Gtk::MenuItem* selall;
    Gtk::Menu* pmenu;

    BatchQueueListener* listener;

    Glib::ustring obtainFileName (const Glib::ustring& origFileName);
    Glib::ustring autoCompleteFileName (const Glib::ustring& fileName, const Glib::ustring& format);

  public:
    BatchQueue ();

    void addEntry (BatchQueueEntry* entry, bool head=false);
    
    void cancelItems (std::vector<ThumbBrowserEntryBase*>* items);
    void headItems (std::vector<ThumbBrowserEntryBase*>* items);
    void tailItems (std::vector<ThumbBrowserEntryBase*>* items);
    void selectAll ();

    void startProcessing ();
    
    bool hasJobs () { return fd.size()>0; }

    rtengine::ProcessingJob* imageReady (rtengine::IImage16* img);
    void setProgress (double p);
    void rightClicked (ThumbBrowserEntryBase* entry);
    void buttonPressed (LWButton* button, int actionCode, void* actionData);
    void redrawNeeded  (LWButton* button);
    
    void setBatchQueueListener (BatchQueueListener* l) { listener = l; }
    void notifyListener ();
};

#endif
