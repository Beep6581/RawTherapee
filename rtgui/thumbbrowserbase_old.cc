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

ThumbBrowserBase::ThumbBrowserBase () 
    : previewHeight(options.thumbSize), lastClicked(NULL) {

    signal_style_changed().connect( sigc::mem_fun(*this, &ThumbBrowserBase::styleChanged) );
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
        if (get_parent () && rowHeight>0) {
            Gtk::Allocation alloc = get_parent ()->get_allocation ();
            numOfRows = (alloc.get_height()+rowHeight/2)/rowHeight;
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
        set_size_request (currx, numOfRows*rowHeight);
    }
    else {
        int availWidth = 0;
        if (get_parent ()) {
            Gtk::Allocation alloc = get_parent ()->get_allocation ();
            availWidth = alloc.get_width();
        }

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
            curry += rowHeight;
        }
        set_size_request (colsWidth, curry);
    }
}


void ThumbBrowserBase::on_realize()
{
  Cairo::FontOptions cfo;
  cfo.set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
  get_pango_context()->set_cairo_font_options (cfo);

  add_events(Gdk::LEAVE_NOTIFY_MASK);
  Gtk::DrawingArea::on_realize();
  Glib::RefPtr<Gdk::Window> window = get_window();
  add_events(Gdk::EXPOSURE_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::POINTER_MOTION_MASK);
 
  gc_ = Gdk::GC::create(window);
}

void ThumbBrowserBase::styleChanged (const Glib::RefPtr<Gtk::Style>& style) {

  refreshAll ();
}

bool ThumbBrowserBase::on_expose_event(GdkEventExpose* event) {

    Gtk::Viewport* vp = (Gtk::Viewport*) get_parent ();
    Gtk::ScrolledWindow* sw = (Gtk::ScrolledWindow*) vp->get_parent ();
  
    int px = (int)(sw->get_hscrollbar()->get_value ());
    int py = (int)(sw->get_vscrollbar()->get_value ());
    int pw = vp->get_width ();
    int ph = vp->get_height ();
    Glib::RefPtr<Gdk::Window> window = get_window();    
    
    window->clear();
    // draw thumbnails
    Glib::RefPtr<Pango::Context> context = get_pango_context ();
    context->set_font_description (get_style()->get_font());
    for (int i=0; i<fd.size(); i++) {
        if (!fd[i]->drawable || !fd[i]->insideWindow (px, py, pw, ph))
            continue;
        fd[i]->draw (this);
    }
    return true;
}

bool ThumbBrowserBase::on_button_press_event (GdkEventButton* event) {

    Gtk::Viewport* vp = (Gtk::Viewport*) get_parent ();
    Gtk::ScrolledWindow* sw = (Gtk::ScrolledWindow*) vp->get_parent ();
  
    int px = (int)(sw->get_hscrollbar()->get_value ());
    int py = (int)(sw->get_vscrollbar()->get_value ());
    int pw = vp->get_width ();
    int ph = vp->get_height ();

    ThumbBrowserEntryBase* fileDescr = NULL;
    bool handled = false;
    for (int i=0; i<fd.size(); i++) 
        if (fd[i]->drawable && fd[i]->insideWindow (px, py, pw, ph)) {
            if (fd[i]->inside ((int)event->x, (int)event->y))            
                fileDescr = fd[i];
            bool b = fd[i]->pressNotify (this, (int)event->x, (int)event->y);
            handled = handled || b;
        }   
    if (handled || (fileDescr && fileDescr->processing))
        return true;
        
    if (selected.size()==1 && event->type==GDK_2BUTTON_PRESS && event->button==1) 
        doubleClicked (selected[0]);
    else if (event->button==1 && event->type==GDK_BUTTON_PRESS) {
        if (fileDescr && event->state & GDK_SHIFT_MASK) {
            if (selected.size()==0) {
                selected.push_back (fileDescr);
                fileDescr->selected = true;
                lastClicked = fileDescr;
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
                        fd[i]->selected = true;
                        selected.push_back (fd[i]);
                }
            }
        }
        else if (fileDescr && event->state & GDK_CONTROL_MASK) {
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
        }
    }
    else if (fileDescr && event->button==3 && event->type==GDK_BUTTON_PRESS) {
        if (!fileDescr->selected) {
            for (int i=0; i<selected.size(); i++)
                selected[i]->selected = false;
            selected.clear ();
            fileDescr->selected = true;
            selected.push_back (fileDescr);
            lastClicked = fileDescr;
        }
        rightClicked (fileDescr);
    }
    Glib::RefPtr<Gdk::Window> window = get_window();

    GdkRectangle rect;
    rect.x = 0;
    rect.y = 0; 
    window->get_size (rect.width, rect.height);

    gdk_window_invalidate_rect (window->gobj(), &rect, true);
    gdk_window_process_updates (window->gobj(), true);
  
    return true;
}

bool ThumbBrowserBase::on_button_release_event (GdkEventButton* event) {

    Gtk::Viewport* vp = (Gtk::Viewport*) get_parent ();
    Gtk::ScrolledWindow* sw = (Gtk::ScrolledWindow*) vp->get_parent ();
  
    int px = (int)(sw->get_hscrollbar()->get_value ());
    int py = (int)(sw->get_vscrollbar()->get_value ());
    int pw = vp->get_width ();
    int ph = vp->get_height ();

    for (int i=0; i<fd.size(); i++) 
        if (fd[i]->drawable && fd[i]->insideWindow (px, py, pw, ph)) 
            fd[i]->releaseNotify (this, (int)event->x, (int)event->y);
    return true;
}

bool ThumbBrowserBase::on_motion_notify_event (GdkEventMotion* event) {

    Gtk::Viewport* vp = (Gtk::Viewport*) get_parent ();
    Gtk::ScrolledWindow* sw = (Gtk::ScrolledWindow*) vp->get_parent ();
  
    int px = (int)(sw->get_hscrollbar()->get_value ());
    int py = (int)(sw->get_vscrollbar()->get_value ());
    int pw = vp->get_width ();
    int ph = vp->get_height ();

    for (int i=0; i<fd.size(); i++) 
        if (fd[i]->drawable && fd[i]->insideWindow (px, py, pw, ph)) 
            fd[i]->motionNotify (this, (int)event->x, (int)event->y);
    return true;
}

void ThumbBrowserBase::resized (Gtk::Allocation& req) {
   arrangeFiles ();
}

void ThumbBrowserBase::redraw () {

    arrangeFiles ();
    queue_draw ();
}

void ThumbBrowserBase::setPreviewHeight (int h) {

    previewHeight = h;
    for (int i=0; i<fd.size(); i++)
        fd[i]->initSizes (this, previewHeight);
    arrangeFiles ();
}

void ThumbBrowserBase::refreshAll () {

    for (int i=0; i<fd.size(); i++) {
        fd[i]->forceHeight (options.thumbSize);
        fd[i]->updateImg ();
        fd[i]->initSizes (this, options.thumbSize);
    }   
    redraw ();
}

void ThumbBrowserBase::setOpenedFileName (const Glib::ustring& fname) {

    fileInEditor = fname;
    for (int i=0; i<fd.size(); i++)
        fd[i]->framed = fd[i]->filename==fileInEditor;
}
    
void ThumbBrowserBase::setArrangement (Arrangement a) {

    arrangement = a;
    redraw ();
}
