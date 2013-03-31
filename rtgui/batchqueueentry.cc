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
#include "batchqueueentry.h"
#include "thumbbrowserbase.h"

#include <cstring>
#include "guiutils.h"
#include "../rtengine/safegtk.h"
#include "multilangmgr.h"

bool BatchQueueEntry::iconsLoaded(false);
Glib::RefPtr<Gdk::Pixbuf> BatchQueueEntry::savedAsIcon;

BatchQueueEntry::BatchQueueEntry (rtengine::ProcessingJob* pjob, const rtengine::procparams::ProcParams& pparams, Glib::ustring fname, guint8* previmg, int prevw, int prevh, Thumbnail* thm) 
    : ThumbBrowserEntryBase(fname),
      opreview(previmg), origpw(prevw), origph(prevh),
      job(pjob), progress(0), outFileName(""), forceFormatOpts(false) {

    thumbnail=thm;
    params = pparams;

    #ifndef WIN32
	// The BatchQueueEntryIdleHelper tracks if an entry has been deleted while it was sitting wating for "idle"
    bqih = new BatchQueueEntryIdleHelper;
    bqih->bqentry = this;
    bqih->destroyed = false;
    bqih->pending = 0;
    #endif

    if (!iconsLoaded) {
        savedAsIcon = safe_create_from_file ("gtk-save.png");
        iconsLoaded = true;
    }

    if (thumbnail)
        thumbnail->increaseRef ();
}

BatchQueueEntry::~BatchQueueEntry () {

    batchQueueEntryUpdater.removeJobs (this);
    delete [] opreview; opreview=NULL;
    if (thumbnail)
        thumbnail->decreaseRef ();

    #ifndef WIN32
    if (bqih->pending)
        bqih->destroyed = true;
    else
        delete bqih;
    #endif
}

void BatchQueueEntry::refreshThumbnailImage () {

    if (!opreview)
        return;

    batchQueueEntryUpdater.process (opreview, origpw, origph, preh, this);  // this will asynchronously land at this.updateImage
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

std::vector<Glib::RefPtr<Gdk::Pixbuf> > BatchQueueEntry::getIconsOnImageArea () {

    std::vector<Glib::RefPtr<Gdk::Pixbuf> > ret;

    if (!outFileName.empty())
        ret.push_back (savedAsIcon);

   return ret;
}

void BatchQueueEntry::getIconSize (int& w, int& h) {

    w = savedAsIcon->get_width ();
    h = savedAsIcon->get_height ();
}


Glib::ustring BatchQueueEntry::getToolTip (int x, int y) {
    // get the parent class' tooltip first
    Glib::ustring tooltip = ThumbBrowserEntryBase::getToolTip(x, y);

    // add the saving param options
    if (!outFileName.empty()) {
        tooltip += Glib::ustring::compose("\n\n%1: %2", M("BATCHQUEUE_DESTFILENAME"), outFileName);
        if (forceFormatOpts) {
            tooltip += Glib::ustring::compose("\n\n%1: %2 (%3 bits)", M("SAVEDLG_FILEFORMAT"), saveFormat.format,
                       saveFormat.format == "png" ? saveFormat.pngBits :
                                                    saveFormat.format == "tif" ? saveFormat.tiffBits : 8);
            if (saveFormat.format == "jpg") {
                tooltip += Glib::ustring::compose("\n%1: %2\n%3: %4",
                           M("SAVEDLG_JPEGQUAL"), saveFormat.jpegQuality,
                           M("SAVEDLG_SUBSAMP"),
                           saveFormat.jpegSubSamp==1 ? M("SAVEDLG_SUBSAMP_1") :
                                                      saveFormat.jpegSubSamp==2 ? M("SAVEDLG_SUBSAMP_2") :
                                                                                  M("SAVEDLG_SUBSAMP_3"));
            }
            else if (saveFormat.format == "png")
                tooltip += Glib::ustring::compose("\n%1: %2", M("SAVEDLG_PNGCOMPR"), saveFormat.pngCompression);
            else if (saveFormat.format == "tif") {
                if (saveFormat.tiffUncompressed)
                tooltip += Glib::ustring::compose("\n%1", M("SAVEDLG_TIFFUNCOMPRESSED"));
            }
        }
    }

    return tooltip;

}

#ifndef WIN32

struct BQUpdateParam {
    BatchQueueEntryIdleHelper* bqih;
    guint8* img;
    int w,h;
};

int updateImageUIThread (void* data) {

    BQUpdateParam* params = static_cast<BQUpdateParam*>(data);

    BatchQueueEntryIdleHelper* bqih = params->bqih;

	// If the BQEntry was destroyed meanwhile, remove all the IdleHelper if all entries came through
    if (bqih->destroyed) {
        if (bqih->pending == 1)
            delete bqih;
        else    
            bqih->pending--;
        delete [] params->img;
        delete params;

        return 0;
    }
    
    bqih->bqentry->_updateImage (params->img, params->w, params->h);
    bqih->pending--;
    
    delete params;
    return 0;
}
#endif

// Starts a copy of img->preview via GTK thread
void BatchQueueEntry::updateImage (guint8* img, int w, int h) {
    // TODO: Check for Linux/Mac
#ifdef WIN32
    // since the update itself is already called in an async thread and there are problem with accessing opreview in thumbbrowserbase,
    // it's safer to do this synchrously
    {
        GThreadLock lock;

        _updateImage(img,w,h);
    }
#else
    bqih->pending++;

    BQUpdateParam* param = new BQUpdateParam ();
    param->bqih = bqih;
    param->img = img;
    param->w = w;
    param->h = h;
    g_idle_add (updateImageUIThread, param);
#endif
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

