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
#include <batchqueueentry.h>
#include <thumbbrowserbase.h>

BatchQueueEntry::BatchQueueEntry (rtengine::ProcessingJob* pjob, const rtengine::procparams::ProcParams& pparams, Glib::ustring fname, guint8* previmg, int prevw, int prevh, Thumbnail* thumbnail) 
    : job(pjob), ThumbBrowserEntryBase(fname), 
        opreview(previmg), origpw(prevw), origph(prevh), progress(0), thumbnail(thumbnail), 
        outFileName("") {

    params = pparams;

    bqih = new BatchQueueEntryIdleHelper;
    bqih->bqentry = this;
    bqih->destroyed = false;
    bqih->pending = 0;

    if (thumbnail)
        thumbnail->increaseRef ();       
}

BatchQueueEntry::~BatchQueueEntry () {

    batchQueueEntryUpdater.removeJobs (this);
    delete [] opreview;
    if (thumbnail)
        thumbnail->decreaseRef ();

    if (bqih->pending)
        bqih->destroyed = true;
    else
        delete bqih;
}

void BatchQueueEntry::refreshThumbnailImage () {

    if (!opreview)
        return;

    batchQueueEntryUpdater.add (opreview, origpw, origph, preh, this);
    batchQueueEntryUpdater.process ();
}

void BatchQueueEntry::calcThumbnailSize () {

    prew = preh * origpw / origph;
}


void BatchQueueEntry::drawProgressBar (Glib::RefPtr<Gdk::Window> win, Glib::RefPtr<Gdk::GC> gc, const Gdk::Color& foregr, const Gdk::Color& backgr, int x, int w, int y, int h) {

  if (processing) {
    Cairo::RefPtr<Cairo::Context> cr = win->create_cairo_context();
    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    double px = x + w/6.0;
    double pw = w*2.0/3.0;
    double py = y + h/4.0;
    double ph = h/2.0;
    cr->move_to (px, py);    
    cr->line_to (px+pw, py);    
    cr->set_line_width (ph);
    cr->set_line_cap (Cairo::LINE_CAP_ROUND);
    cr->set_source_rgb (foregr.get_red_p(), foregr.get_green_p(), foregr.get_blue_p());
    cr->stroke ();

    cr->move_to (px, py);    
    cr->line_to (px+pw, py);    
    cr->set_line_width (ph*3.0/4.0);
    cr->set_source_rgb (backgr.get_red_p(), backgr.get_green_p(), backgr.get_blue_p());
    cr->stroke ();

    cr->move_to (px, py);    
    cr->line_to (px+pw*progress, py);    
    cr->set_line_width (ph/2.0);
    cr->set_source_rgb (foregr.get_red_p(), foregr.get_green_p(), foregr.get_blue_p());
    cr->stroke ();
  }
}

void BatchQueueEntry::removeButtonSet () {

    delete buttonSet;
    buttonSet = NULL;
}
struct bqupdate {
    BatchQueueEntryIdleHelper* bqih;
    guint8* img;
    int w,h;
};

int bqeupdate (void* data) {

    gdk_threads_enter ();
    bqupdate* params = (bqupdate*)data;

    BatchQueueEntryIdleHelper* bqih = params->bqih;

    if (bqih->destroyed) {
        if (bqih->pending == 1)
            delete bqih;
        else    
            bqih->pending--;
        delete [] params->img;
        delete params;
        gdk_threads_leave ();
        return 0;
    }
    
    bqih->bqentry->_updateImage (params->img, params->w, params->h);
    bqih->pending--;
    
    gdk_threads_leave ();
    delete params;
    return 0;
}

void BatchQueueEntry::updateImage (guint8* img, int w, int h) {

    bqih->pending++;

    bqupdate* param = new bqupdate ();
    param->bqih = bqih;
    param->img = img;
    param->w = w;
    param->h = h;
    g_idle_add (bqeupdate, param);
}

void BatchQueueEntry::_updateImage (guint8* img, int w, int h) {

    if (preh == h) {
        prew = w;
        preview = new guint8 [prew*preh*3];
        memcpy (preview, img, prew*preh*3);
        if (parent)
            parent->redrawNeeded (this);
    }
    delete [] img;
}

