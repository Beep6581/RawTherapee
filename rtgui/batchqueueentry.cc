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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "batchqueueentry.h"

#include <cstring>

#include "guiutils.h"
#include "threadutils.h"
#include "rtsurface.h"
#include "multilangmgr.h"
#include "thumbbrowserbase.h"
#include "thumbnail.h"
#include "rtsurface.h"

#include "rtengine/procparams.h"
#include "rtengine/rtengine.h"

bool BatchQueueEntry::iconsLoaded(false);
std::shared_ptr<RTSurface> BatchQueueEntry::savedAsIcon(std::shared_ptr<RTSurface>(nullptr));

BatchQueueEntry::BatchQueueEntry (rtengine::ProcessingJob* pjob, const rtengine::procparams::ProcParams& pparams, Glib::ustring fname, int prevw, int prevh, Thumbnail* thm, bool overwrite) :
    ThumbBrowserEntryBase(fname, thm),
    opreview(nullptr),
    origpw(prevw),
    origph(prevh),
    opreviewDone(false),
    job(pjob),
    params(new rtengine::procparams::ProcParams(pparams)),
    progress(0),
    sequence(0),
    forceFormatOpts(false),
    fast_pipeline(job->fastPipeline()),
    overwriteFile(overwrite)
{

    thumbnail = thm;

#if 1 //ndef _WIN32
    // The BatchQueueEntryIdleHelper tracks if an entry has been deleted while it was sitting waiting for "idle"
    bqih = new BatchQueueEntryIdleHelper;
    bqih->bqentry = this;
    bqih->destroyed = false;
    bqih->pending = 0;
#endif

    if (!iconsLoaded) {
        savedAsIcon = std::shared_ptr<RTSurface>(new RTSurface("save-small", Gtk::ICON_SIZE_SMALL_TOOLBAR));
        iconsLoaded = true;
    }

    if (thumbnail) {
        thumbnail->increaseRef ();
    }
}

BatchQueueEntry::~BatchQueueEntry ()
{

    batchQueueEntryUpdater.removeJobs (this);

    if (opreview) {
        delete [] opreview;
    }

    opreview = nullptr;

    if (thumbnail) {
        thumbnail->decreaseRef ();
    }

    if (bqih->pending) {
        bqih->destroyed = true;
    } else {
        delete bqih;
    }
}

void BatchQueueEntry::refreshThumbnailImage ()
{

    if (!opreviewDone) {
        // creating the image buffer first
        //if (!opreview) opreview = new guint8[(origpw+1) * origph * 3];
        // this will asynchronously compute the original preview and land at this.updateImage
        batchQueueEntryUpdater.process (nullptr, origpw, origph, preh, this, params.get(), thumbnail);
    } else {
        // this will asynchronously land at this.updateImage
        batchQueueEntryUpdater.process (opreview, origpw, origph, preh, this);
    }
}

void BatchQueueEntry::calcThumbnailSize ()
{
    prew = preh * origpw / origph;
    if (prew > options.maxThumbnailWidth) {
        const float s = static_cast<float>(options.maxThumbnailWidth) / prew;
        prew = options.maxThumbnailWidth;
        preh = std::max<int>(preh * s, 1);
    }
}


void BatchQueueEntry::drawProgressBar (Glib::RefPtr<Gdk::Window> win, const Gdk::RGBA& foregr, const Gdk::RGBA& backgr, int x, int w, int y, int h)
{

    if (processing) {
        Cairo::RefPtr<Cairo::Context> cr = win->create_cairo_context();
        cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        double px = x + w / 6.0;
        double pw = w * 2.0 / 3.0;
        double py = y + h / 4.0;
        double ph = h / 2.0;
        cr->move_to (px, py);
        cr->line_to (px + pw, py);
        cr->set_line_width (ph);
        cr->set_line_cap (Cairo::LINE_CAP_ROUND);
        cr->set_source_rgb (foregr.get_red(), foregr.get_green(), foregr.get_blue());
        cr->stroke ();

        cr->move_to (px, py);
        cr->line_to (px + pw, py);
        cr->set_line_width (ph * 3.0 / 4.0);
        cr->set_source_rgb (backgr.get_red(), backgr.get_green(), backgr.get_blue());
        cr->stroke ();

        cr->move_to (px, py);
        cr->line_to (px + pw * progress, py);
        cr->set_line_width (ph / 2.0);
        cr->set_source_rgb (foregr.get_red(), foregr.get_green(), foregr.get_blue());
        cr->stroke ();
    }
}

void BatchQueueEntry::removeButtonSet ()
{

    delete buttonSet;
    buttonSet = nullptr;
}

std::vector<std::shared_ptr<RTSurface>> BatchQueueEntry::getIconsOnImageArea ()
{

    std::vector<std::shared_ptr<RTSurface>> ret;

    if (!outFileName.empty()) {
        ret.push_back (savedAsIcon);
    }

    return ret;
}

void BatchQueueEntry::getIconSize (int& w, int& h) const
{

    w = savedAsIcon->getWidth ();
    h = savedAsIcon->getHeight ();
}


std::tuple<Glib::ustring, bool> BatchQueueEntry::getToolTip (int x, int y) const
{
    // get the parent class' tooltip first
    Glib::ustring tooltip;
    bool useMarkup;
    std::tie(tooltip, useMarkup) =  ThumbBrowserEntryBase::getToolTip(x, y);

    // add the saving param options
    if (!outFileName.empty()) {
        tooltip += Glib::ustring::compose("\n\n%1: %2", M("BATCHQUEUE_DESTFILENAME"), outFileName);

        if (forceFormatOpts) {
            tooltip += Glib::ustring::compose("\n\n%1: %2 (%3-bits%4)", M("SAVEDLG_FILEFORMAT"), saveFormat.format,
                                              saveFormat.format == "png" ? saveFormat.pngBits :
                                              saveFormat.format == "tif" ? saveFormat.tiffBits : 8,
                                              saveFormat.format == "tif" && saveFormat.tiffFloat ? M("SAVEDLG_FILEFORMAT_FLOAT") : "");

            if (saveFormat.format == "jpg") {
                tooltip += Glib::ustring::compose("\n%1: %2\n%3: %4",
                                                  M("SAVEDLG_JPEGQUAL"), saveFormat.jpegQuality,
                                                  M("SAVEDLG_SUBSAMP"),
                                                  saveFormat.jpegSubSamp == 1 ? M("SAVEDLG_SUBSAMP_1") :
                                                  saveFormat.jpegSubSamp == 2 ? M("SAVEDLG_SUBSAMP_2") :
                                                  M("SAVEDLG_SUBSAMP_3"));
            } else if (saveFormat.format == "tif") {
                if (saveFormat.tiffUncompressed) {
                    tooltip += Glib::ustring::compose("\n%1", M("SAVEDLG_TIFFUNCOMPRESSED"));
                }
                if (saveFormat.bigTiff) {
                    tooltip += Glib::ustring::compose("\n%1", M("SAVEDLG_BIGTIFF"));
                }
            }
        }
    }

    return std::make_tuple(std::move(tooltip), useMarkup);

}

struct BQUpdateParam {
    BatchQueueEntryIdleHelper* bqih;
    guint8* img;
    int w, h;
};

int updateImageUIThread (void* data)
{

    BQUpdateParam* params = static_cast<BQUpdateParam*>(data);

    BatchQueueEntryIdleHelper* bqih = params->bqih;

    GThreadLock tLock; // Acquire the GUI

    // If the BQEntry was destroyed meanwhile, remove all the IdleHelper if all entries came through
    if (bqih->destroyed) {
        if (bqih->pending == 1) {
            delete bqih;
        } else {
            bqih->pending--;
        }

        delete [] params->img;
        delete params;

        return 0;
    }

    bqih->bqentry->_updateImage (params->img, params->w, params->h);
    bqih->pending--;

    delete params;
    return 0;
}

// Starts a copy of img->preview via GTK thread
void BatchQueueEntry::updateImage (guint8* img, int w, int h, int origw, int origh, guint8* newOPreview)
{

    // since the update itself is already called in an async thread and there are problem with accessing opreview in thumbbrowserbase,
    // it's safer to do this synchronously
    {
        GThreadLock lock;

        _updateImage(img, w, h);
    }
}

void BatchQueueEntry::_updateImage (guint8* img, int w, int h)
{

    if (preh == h) {
        MYWRITERLOCK(l, lockRW);

        prew = w;
        preview.resize(prew * preh * 3);
        std::copy(img, img + preview.size(), preview.begin());

        if (parent) {
            parent->redrawNeeded (this);
        }
    }

    delete [] img;
}

