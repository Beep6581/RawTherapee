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
#include <thumbbrowserbase.h>
#include <glibmm.h>
#include <multilangmgr.h>
#include <options.h>
#include <mytime.h>

ThumbBrowserBase::ThumbBrowserBase () 
    : previewHeight(options.thumbSize), lastClicked(NULL) {

    inW = -1; inH = -1;

    Gtk::HBox* hb1 = new Gtk::HBox ();
    Gtk::HBox* hb2 = new Gtk::HBox ();
    Gtk::Frame* frame = new Gtk::Frame ();
    frame->add (internal);
    frame->set_shadow_type (Gtk::SHADOW_IN );
    hb1->pack_start (*frame);
    hb1->pack_end (vscroll, Gtk::PACK_SHRINK, 0);

    pack_start (*hb1);

    Gtk::HBox* tmp = new Gtk::HBox ();
    hb2->pack_start (hscroll);
    Gtk::Requisition vcr = vscroll.size_request ();
    tmp->set_size_request (vcr.width, vcr.width);
    hb2->pack_end (*tmp, Gtk::PACK_SHRINK, 0);

    pack_start (*hb2,Gtk::PACK_SHRINK, 0);

    internal.setParent (this);
    
    show_all ();
    
    hscroll.set_update_policy (Gtk::UPDATE_CONTINUOUS);
    vscroll.set_update_policy (Gtk::UPDATE_CONTINUOUS);

    vscroll.signal_value_changed().connect( sigc::mem_fun(*this, &ThumbBrowserBase::scrollChanged) );
    hscroll.signal_value_changed().connect( sigc::mem_fun(*this, &ThumbBrowserBase::scrollChanged) );

    internal.signal_size_allocate().connect( sigc::mem_fun(*this, &ThumbBrowserBase::internalAreaResized) );
    signal_style_changed().connect( sigc::mem_fun(*this, &ThumbBrowserBase::styleChanged) );
}

void ThumbBrowserBase::scrollChanged () {

    for (int i=0; i<fd.size(); i++)
        fd[i]->setOffset ((int)(hscroll.get_value()), (int)(vscroll.get_value()));

    internal.setPosition ((int)(hscroll.get_value()), (int)(vscroll.get_value()));

    if (!internal.isDirty()) {    
        internal.setDirty ();
        internal.queue_draw ();
//        gdk_window_process_updates (get_window()->gobj(), true);
    }
}

void ThumbBrowserBase::scroll (int direction) {

    if (arrangement==TB_Vertical)
        vscroll.set_value (vscroll.get_value() + (direction==GDK_SCROLL_DOWN ? +1 : -1) * vscroll.get_adjustment()->get_step_increment());
    else
        hscroll.set_value (hscroll.get_value() + (direction==GDK_SCROLL_DOWN ? +1 : -1) * hscroll.get_adjustment()->get_step_increment());
}

void ThumbBrowserBase::resizeThumbnailArea (int w, int h) {

    inW = w;
    inH = h;

    if (hscroll.get_value() + internal.get_width() > inW)
        hscroll.set_value (inW - internal.get_width());
    if (vscroll.get_value() + internal.get_height() > inH)
        vscroll.set_value (inH - internal.get_height());

    configScrollBars ();
}

void ThumbBrowserBase::internalAreaResized (Gtk::Allocation& req) {

  if (inW>0 && inH>0) {
    configScrollBars ();
    redraw ();
  }
}

void ThumbBrowserBase::configScrollBars () {

    if (inW>0 && inH>0) {
  
        int iw = internal.get_width ();
        int ih = internal.get_height ();

        hscroll.get_adjustment()->set_upper (inW);
        vscroll.get_adjustment()->set_upper (inH);
        hscroll.get_adjustment()->set_lower (0);
        vscroll.get_adjustment()->set_lower (0);
        hscroll.get_adjustment()->set_step_increment (32);
        vscroll.get_adjustment()->set_step_increment (32);
        hscroll.get_adjustment()->set_page_increment (iw);
        vscroll.get_adjustment()->set_page_increment (ih);
        hscroll.get_adjustment()->set_page_size (iw);
        vscroll.get_adjustment()->set_page_size (ih);
    }
}

void ThumbBrowserBase::arrangeFiles () {

    int N = fd.size ();
    // apply filter
    for (int i=0; i<N; i++) 
        fd[i]->filtered = !checkFilter (fd[i]);

    int rowHeight = 0;
    // compute size of the items
    for (int i=0; i<N; i++) 
        if (!fd[i]->filtered && fd[i]->getMinimalHeight() > rowHeight)
            rowHeight = fd[i]->getMinimalHeight ();
    
    if (arrangement==TB_Horizontal) {            

        int numOfRows = 1;
        if (rowHeight>0) {
            numOfRows = (internal.get_height()+rowHeight/2)/rowHeight;
            if (numOfRows<1)
                numOfRows = 1;
        }

        int ct = 0;
        int currx = 0; int curry = 0;
        while (ct<N) {
            // find widest item in the column
            int maxw = 0;
            for (int i=0; ct+i<N && i<numOfRows; i++)
                if (fd[ct+i]->getMinimalWidth() > maxw)
                    maxw = fd[ct+i]->getMinimalWidth ();

            // arrange items in the column
            curry = 0;
            for (int i=0; ct<N && i<numOfRows; i++, ct++) {
                while (ct<N && fd[ct]->filtered) 
                    fd[ct++]->drawable = false;
                if (ct<N) {
                    fd[ct]->setPosition (currx, curry, maxw, rowHeight);
                    fd[ct]->drawable = true;
                    curry += rowHeight;
                }
            }
            currx += maxw;
        }
        resizeThumbnailArea (currx, numOfRows*rowHeight);
    }
    else {
        int availWidth = internal.get_width();
        // initial number of columns
        int numOfCols = 0;
        int currColNum = 0;
        int colsWidth = 0;
        for (int i=0; i<N; i++)
            if (!fd[i]->filtered && colsWidth + fd[i]->getMinimalWidth() <= availWidth) {
                colsWidth += fd[numOfCols]->getMinimalWidth ();
                numOfCols++;
            }
        if (numOfCols<1)
            numOfCols = 1;
        std::vector<int> colWidths;
        for (; numOfCols>0; numOfCols--) {
            // compute column widths
            colWidths.resize (numOfCols);
            for (int i=0; i<numOfCols; i++)
                colWidths[i] = 0;
            for (int i=0, j=0; i<N; i++) {
                if (!fd[i]->filtered && fd[i]->getMinimalWidth() > colWidths[j%numOfCols])
                    colWidths[j%numOfCols] = fd[i]->getMinimalWidth ();
                if (!fd[i]->filtered)
                    j++;
            }
            // if not wider than the space available, arrange it and we are ready
            colsWidth = 0;
            for (int i=0; i<numOfCols; i++)
                colsWidth += colWidths[i];
            if (numOfCols==1 || colsWidth < availWidth)
                break;
        }
        // arrange files
        int ct = 0;
        int currx = 0; int curry = 0;
        while (ct<N) {
            // arrange items in the row
            currx = 0;
            for (int i=0; ct<N && i<numOfCols; i++, ct++) {
                while (ct<N && fd[ct]->filtered) 
                    fd[ct++]->drawable = false;
                if (ct<N) {
                    fd[ct]->setPosition (currx, curry, colWidths[i%numOfCols], rowHeight);
                    fd[ct]->drawable = true;
                    currx += colWidths[i%numOfCols];
                }
            }
            if (currx>0) // there were thumbnails placed in the row
                curry += rowHeight;
        }
        resizeThumbnailArea (colsWidth, curry);
    }
}


void ThumbBrowserBase::Internal::on_realize()
{
  Cairo::FontOptions cfo;
  cfo.set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
  get_pango_context()->set_cairo_font_options (cfo);

  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  set_flags (Gtk::CAN_FOCUS);
  add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::POINTER_MOTION_MASK | Gdk::SCROLL_MASK | Gdk::KEY_PRESS_MASK);
  gc_ = Gdk::GC::create(window);
  set_has_tooltip (true);
  signal_query_tooltip().connect( sigc::mem_fun(*this, &ThumbBrowserBase::Internal::on_query_tooltip) );
}

bool ThumbBrowserBase::Internal::on_query_tooltip (int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip) {

    Glib::ustring ttip = "";
    for (int i=0; i<parent->fd.size(); i++) 
        if (parent->fd[i]->drawable && parent->fd[i]->inside (x, y)) {
            ttip = parent->fd[i]->getToolTip (x, y);
            break;
        }
    if (ttip!="") {
        tooltip->set_text (ttip);
        return true;
    }
    else
        return false;
}

void ThumbBrowserBase::styleChanged (const Glib::RefPtr<Gtk::Style>& style) {

  refreshThumbImages ();
}

ThumbBrowserBase::Internal::Internal () : parent(NULL), ofsX(0), ofsY(0), dirty(true) {
}

void ThumbBrowserBase::Internal::setParent (ThumbBrowserBase* p) {

    parent = p;
}

void ThumbBrowserBase::Internal::setPosition (int x, int y) {

    ofsX = x;
    ofsY = y;
}

bool ThumbBrowserBase::Internal::on_key_press_event (GdkEventKey* event) {

    return parent->keyPressed (event);
}

bool ThumbBrowserBase::Internal::on_button_press_event (GdkEventButton* event) {

    grab_focus ();

    parent->eventTime = event->time;

    parent->buttonPressed ((int)event->x, (int)event->y, event->button, event->type, event->state, 0, 0, get_width(), get_height());
    Glib::RefPtr<Gdk::Window> window = get_window();

    GdkRectangle rect;
    rect.x = 0;
    rect.y = 0; 
    window->get_size (rect.width, rect.height);

    gdk_window_invalidate_rect (window->gobj(), &rect, true);
    gdk_window_process_updates (window->gobj(), true);
    return true;
}

void ThumbBrowserBase::buttonPressed (int x, int y, int button, GdkEventType type, int state, int clx, int cly, int clw, int clh) {

    ThumbBrowserEntryBase* fileDescr = NULL;
    bool handled = false;
    for (int i=0; i<fd.size(); i++) 
        if (fd[i]->drawable) {
            if (fd[i]->inside (x, y) && fd[i]->insideWindow (clx, cly, clw, clh))            
                fileDescr = fd[i];
            bool b = fd[i]->pressNotify (button, type, state, x, y);
            handled = handled || b;
        }   
        
    if (handled || (fileDescr && fileDescr->processing))
        return;
        
    if (selected.size()==1 && type==GDK_2BUTTON_PRESS && button==1) 
        doubleClicked (selected[0]);
    else if (button==1 && type==GDK_BUTTON_PRESS) {
        if (fileDescr && state & GDK_SHIFT_MASK) {
            if (selected.size()==0) {
                selected.push_back (fileDescr);
                fileDescr->selected = true;
                lastClicked = fileDescr;
                selectionChanged ();
            }
            else {
                // find the start and the end of the selection interval
                int startx = fd.size()-1;
                if (lastClicked) {
                    for (; startx>=0; startx--)
                        if (fd[startx]==lastClicked) 
                            break;
                }
                else {
                    for (; startx>=0; startx--)
                        if (fd[startx]==selected[0]) 
                            break;
                }
                int endx = 0;
                for (; endx<fd.size(); endx++)
                    if (fd[endx]==fileDescr) 
                        break;
                if (endx < startx) {
                    int tmp = endx;
                    endx = startx;
                    startx = tmp;
                }
                // clear current selection
                for (int i=0; i<selected.size(); i++)
                    selected[i]->selected = false;
                selected.clear ();
                // select thumbnails in the interval
                for (int i=startx; i<=endx; i++) {
                    if (!fd[i]->filtered) {
                        fd[i]->selected = true;
                        selected.push_back (fd[i]);
                    }
                }
                selectionChanged ();
            }
        }
        else if (fileDescr && state & GDK_CONTROL_MASK) {
            std::vector<ThumbBrowserEntryBase*>::iterator i = std::find (selected.begin(), selected.end(), fileDescr);
            if (i!=selected.end()) {
                (*i)->selected = false;
                selected.erase (i);
            }
            else {
                selected.push_back (fileDescr);
                fileDescr->selected = true;
            }
            lastClicked = fileDescr;
            selectionChanged ();
        }
        else {
            for (int i=0; i<selected.size(); i++)
                selected[i]->selected = false;
            selected.clear ();
            if (fileDescr) {
                selected.push_back (fileDescr);
                fileDescr->selected = true;
            }
            lastClicked = fileDescr;
            selectionChanged ();
        }
    }
    else if (fileDescr && button==3 && type==GDK_BUTTON_PRESS) {
        if (!fileDescr->selected) {
            for (int i=0; i<selected.size(); i++)
                selected[i]->selected = false;
            selected.clear ();
            fileDescr->selected = true;
            selected.push_back (fileDescr);
            lastClicked = fileDescr;
            selectionChanged ();
        }
        rightClicked (fileDescr);
    }
}

bool ThumbBrowserBase::Internal::on_expose_event(GdkEventExpose* event) {

    dirty = false;

    Glib::RefPtr<Gdk::Window> window = get_window();    

    int w = get_width();
    int h = get_height();

    window->clear();
    // draw thumbnails
    Glib::RefPtr<Pango::Context> context = get_pango_context ();
    context->set_font_description (get_style()->get_font());
    for (int i=0; i<parent->fd.size(); i++) {
        if (!parent->fd[i]->drawable || !parent->fd[i]->insideWindow (0, 0, w, h)) 
            parent->fd[i]->updatepriority = false;
        else {
            parent->fd[i]->updatepriority = true;
            parent->fd[i]->draw ();
        }
    }
    
    return true;
}

bool ThumbBrowserBase::Internal::on_button_release_event (GdkEventButton* event) {

    int w = get_width();
    int h = get_height();

    for (int i=0; i<parent->fd.size(); i++) 
        if (parent->fd[i]->drawable && parent->fd[i]->insideWindow (0, 0, w, h)) 
            parent->fd[i]->releaseNotify (event->button, event->type, event->state, (int)event->x, (int)event->y);
    return true;
}

bool ThumbBrowserBase::Internal::on_motion_notify_event (GdkEventMotion* event) {

    int w = get_width();
    int h = get_height();

    for (int i=0; i<parent->fd.size(); i++) 
        if (parent->fd[i]->drawable && parent->fd[i]->insideWindow (0, 0, w, h)) 
            parent->fd[i]->motionNotify ((int)event->x, (int)event->y);
    return true;
}

bool ThumbBrowserBase::Internal::on_scroll_event (GdkEventScroll* event) {

    parent->scroll (event->direction);
    return true;
}


void ThumbBrowserBase::redraw () {

    arrangeFiles ();
    queue_draw ();
}

void ThumbBrowserBase::zoomChanged (bool zoomIn) {

    int newHeight;
    int i=0;
    if (zoomIn)
        for (i=0; i<options.thumbnailZoomRatios.size(); i++) {
            newHeight = (int)(options.thumbnailZoomRatios[i] * options.maxThumbnailHeight);
            if (newHeight > options.thumbSize)
                break;
        }
    else
        for (i=options.thumbnailZoomRatios.size()-1; i>=0; i--) {
            newHeight = (int)(options.thumbnailZoomRatios[i] * options.maxThumbnailHeight);
            if (newHeight < options.thumbSize)
                break;
        }
    previewHeight = options.thumbSize = newHeight;
    for (int i=0; i<fd.size(); i++) 
        fd[i]->resize (previewHeight);
    redraw ();
#ifdef _WIN32
    gdk_window_process_updates (get_window()->gobj(), true);
#endif    
}
void ThumbBrowserBase::refreshThumbImages () {

    for (int i=0; i<fd.size(); i++) 
        fd[i]->refreshThumbnailImage ();

    redraw ();
}

void ThumbBrowserBase::refreshEditedState (const std::set<Glib::ustring>& efiles) {

    editedFiles = efiles;
    for (int i=0; i<fd.size(); i++)
        fd[i]->framed = editedFiles.find (fd[i]->filename)!=editedFiles.end();
    queue_draw ();
}
    
void ThumbBrowserBase::setArrangement (Arrangement a) {

    arrangement = a;
    redraw ();
}

void ThumbBrowserBase::initEntry (ThumbBrowserEntryBase* entry) {

        entry->setOffset ((int)(hscroll.get_value()), (int)(vscroll.get_value()));
}
void ThumbBrowserBase::getScrollPosition (double& h, double& v) {

    h = hscroll.get_value ();
    v = vscroll.get_value ();
}

void ThumbBrowserBase::setScrollPosition (double h, double v) {

    hscroll.set_value (h>hscroll.get_adjustment()->get_upper() ? hscroll.get_adjustment()->get_upper() : h);
    vscroll.set_value (v>vscroll.get_adjustment()->get_upper() ? vscroll.get_adjustment()->get_upper() : v);
}

/*void PreviewImgUpdater::processCustomOrder () {

    // find first filtered entry, if any
    std::list<ThumbBrowserEntryBase*>::iterator i;
    for (i=jqueue.begin (); i!=jqueue.end(); i++)
        if (!(*i)->filtered)
            break;
    if (i==jqueue.end())
        i = jqueue.begin();

    ThumbBrowserEntryBase* current = *i;
    jqueue.erase (i);

    current->updateImg ();
    if (parent) {
        gdk_threads_enter ();
        parent->queue_draw ();
        if (parent->get_window())
            gdk_window_process_updates (parent->get_window()->gobj(), true);
        gdk_threads_leave ();
    }
}
*/

