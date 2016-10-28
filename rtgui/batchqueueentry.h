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
#ifndef _BATCHQUEUEENTRY_
#define _BATCHQUEUEENTRY_

#include <gtkmm.h>
#include "../rtengine/rtengine.h"
#include "thumbbrowserentrybase.h"
#include "thumbnail.h"
#include "bqentryupdater.h"

class BatchQueueEntry;
struct BatchQueueEntryIdleHelper {
    BatchQueueEntry* bqentry;
    bool destroyed;
    int pending;
};

class BatchQueueEntry : public ThumbBrowserEntryBase, public BQEntryUpdateListener
{

    guint8* opreview;
    int origpw, origph;
    BatchQueueEntryIdleHelper* bqih;
    bool opreviewDone;
    static bool iconsLoaded;

public:

    static Glib::RefPtr<Gdk::Pixbuf> savedAsIcon;

    rtengine::ProcessingJob* job;
    rtengine::procparams::ProcParams params;
    Glib::ustring savedParamsFile;
    double progress;
    Glib::ustring outFileName;
    int sequence;
    SaveFormat saveFormat;
    bool forceFormatOpts;

    BatchQueueEntry (rtengine::ProcessingJob* job, const rtengine::procparams::ProcParams& pparams, Glib::ustring fname, int prevw, int prevh, Thumbnail* thm = nullptr);
    ~BatchQueueEntry ();

    void refreshThumbnailImage ();
    void calcThumbnailSize ();

    void drawProgressBar (Glib::RefPtr<Gdk::Window> win, const Gdk::RGBA& foregr, const Gdk::RGBA& backgr, int x, int w, int y, int h);

    void removeButtonSet ();

    virtual std::vector<Glib::RefPtr<Gdk::Pixbuf> > getIconsOnImageArea ();
    virtual void getIconSize (int& w, int& h);
    virtual Glib::ustring getToolTip (int x, int y);

    // bqentryupdatelistener interface
    void updateImage (guint8* img, int w, int h, int origw, int origh, guint8* newOPreview);
    void _updateImage (guint8* img, int w, int h); // inside gtk thread
};



#endif
