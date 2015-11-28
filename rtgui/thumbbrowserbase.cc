/*
 *  This file is part of RawTherapee.
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
#include "../rtengine/rt_math.h"

#include "thumbbrowserbase.h"
#include "multilangmgr.h"
#include "options.h"
#include "../rtengine/mytime.h"

using namespace std;

ThumbBrowserBase::ThumbBrowserBase ()
    : lastClicked(NULL), previewHeight(options.thumbSize), numOfCols(1), inspector(NULL), isInspectorActive(false)
{
    location = THLOC_FILEBROWSER;
    inW = -1;
    inH = -1;

    Gtk::HBox* hb1 = Gtk::manage( new Gtk::HBox () );
    Gtk::HBox* hb2 = Gtk::manage( new Gtk::HBox () );
    Gtk::Frame* frame = Gtk::manage( new Gtk::Frame () );
    frame->add (internal);
    frame->set_shadow_type (Gtk::SHADOW_IN );
    hb1->pack_start (*frame);
    hb1->pack_end (vscroll, Gtk::PACK_SHRINK, 0);

    pack_start (*hb1);

    hb2->pack_start (hscroll);

    pack_start (*hb2, Gtk::PACK_SHRINK, 0);

    internal.setParent (this);

    show_all ();

    hscroll.set_update_policy (Gtk::UPDATE_CONTINUOUS);
    vscroll.set_update_policy (Gtk::UPDATE_CONTINUOUS);

    vscroll.signal_value_changed().connect( sigc::mem_fun(*this, &ThumbBrowserBase::scrollChanged) );
    hscroll.signal_value_changed().connect( sigc::mem_fun(*this, &ThumbBrowserBase::scrollChanged) );

    internal.signal_size_allocate().connect( sigc::mem_fun(*this, &ThumbBrowserBase::internalAreaResized) );
}

void ThumbBrowserBase::scrollChanged ()
{
    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        for (size_t i = 0; i < fd.size(); i++) {
            fd[i]->setOffset ((int)(hscroll.get_value()), (int)(vscroll.get_value()));
        }
    }

    internal.setPosition ((int)(hscroll.get_value()), (int)(vscroll.get_value()));

    if (!internal.isDirty()) {
        internal.setDirty ();
        internal.queue_draw ();
    }
}

void ThumbBrowserBase::scroll (int direction)
{
    // GUI already acquired when here
    if (arrangement == TB_Vertical) {
        vscroll.set_value (vscroll.get_value() + (direction == GDK_SCROLL_DOWN ? +1 : -1) * vscroll.get_adjustment()->get_step_increment());
    } else {
        hscroll.set_value (hscroll.get_value() + (direction == GDK_SCROLL_DOWN ? +1 : -1) * hscroll.get_adjustment()->get_step_increment());
    }
}

void ThumbBrowserBase::scrollPage (int direction)
{
    // GUI already acquired when here
    if (arrangement == TB_Vertical) {
        vscroll.set_value (vscroll.get_value() + (direction == GDK_SCROLL_DOWN ? +1 : -1) * vscroll.get_adjustment()->get_page_increment());
    } else {
        hscroll.set_value (hscroll.get_value() + (direction == GDK_SCROLL_DOWN ? +1 : -1) * hscroll.get_adjustment()->get_page_increment());
    }
}

void ThumbBrowserBase::selectSingle (ThumbBrowserEntryBase* fileDescr)
{
    for (size_t i = 0; i < selected.size(); i++) {
        selected[i]->selected = false;
    }

    selected.clear ();

    if (fileDescr) {
        selected.push_back (fileDescr);
        fileDescr->selected = true;
    }
}

void ThumbBrowserBase::selectRange (ThumbBrowserEntryBase* fileDescr)
{
    if (selected.empty()) {
        selected.push_back (fileDescr);
        fileDescr->selected = true;
    } else {
        // find the start and the end of the selection interval
        size_t startx = fd.size() - 1;

        if (lastClicked) {
            for (; startx > 0; startx--)
                if (fd[startx] == lastClicked) {
                    break;
                }
        } else {
            for (; startx > 0; startx--)
                if (fd[startx] == selected[0]) {
                    break;
                }
        }

        size_t endx = 0;

        for (; endx < fd.size(); endx++)
            if (fd[endx] == fileDescr) {
                break;
            }

        if (endx < startx) {
            int tmp = endx;
            endx = startx;
            startx = tmp;
        }

        // clear current selection
        for (size_t i = 0; i < selected.size(); i++) {
            selected[i]->selected = false;
        }

        selected.clear ();

        // select thumbnails in the interval
        for (size_t i = startx; i <= endx; i++) {
            if (!fd[i]->filtered) {
                fd[i]->selected = true;
                selected.push_back (fd[i]);
            }
        }
    }
}

void ThumbBrowserBase::selectSet (ThumbBrowserEntryBase* fileDescr)
{
    std::vector<ThumbBrowserEntryBase*>::iterator i = std::find (selected.begin(), selected.end(), fileDescr);

    if (i != selected.end()) {
        (*i)->selected = false;
        selected.erase (i);
    } else {
        selected.push_back (fileDescr);
        fileDescr->selected = true;
    }
}

static void scrollToEntry (double& h, double& v, int iw, int ih, ThumbBrowserEntryBase* entry)
{
    const int hmin = entry->getX ();
    const int hmax = hmin + entry->getEffectiveWidth () - iw;
    const int vmin = entry->getY ();
    const int vmax = vmin + entry->getEffectiveHeight () - ih;

    if (hmin < 0) {
        h += hmin;
    } else if (hmax > 0) {
        h += hmax;
    }

    if(vmin < 0) {
        v += vmin;
    } else if (vmax > 0) {
        v += vmax;
    }
}

void ThumbBrowserBase::selectPrev (int distance, bool enlarge)
{
    double h, v;
    getScrollPosition (h, v);

    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        if (!selected.empty ()) {
            std::vector<ThumbBrowserEntryBase*>::iterator front = std::find (fd.begin (), fd.end (), selected.front ());
            std::vector<ThumbBrowserEntryBase*>::iterator back = std::find (fd.begin (), fd.end (), selected.back ());
            std::vector<ThumbBrowserEntryBase*>::iterator last = std::find (fd.begin (), fd.end (), lastClicked);

            if (front > back) {
                std::swap(front, back);
            }

            std::vector<ThumbBrowserEntryBase*>::iterator& curr = last == front ? front : back;

            // find next thumbnail at filtered distance before current
            for (; curr >= fd.begin (); --curr) {
                if (!(*curr)->filtered) {
                    if (distance-- == 0) {
                        // clear current selection
                        for (size_t i = 0; i < selected.size (); ++i) {
                            selected[i]->selected = false;
                            redrawNeeded (selected[i]);
                        }

                        selected.clear ();

                        // make sure the newly selected thumbnail is visible and make it current
                        scrollToEntry (h, v, internal.get_width (), internal.get_height (), *curr);
                        lastClicked = *curr;

                        // either enlarge current selection or set new selection
                        if(enlarge) {
                            // reverse direction if distance is too large
                            if(front > back) {
                                std::swap(front, back);
                            }

                            for (; front <= back; ++front) {
                                if (!(*front)->filtered) {
                                    (*front)->selected = true;
                                    redrawNeeded (*front);
                                    selected.push_back (*front);
                                }
                            }
                        } else {
                            (*curr)->selected = true;
                            redrawNeeded (*curr);
                            selected.push_back (*curr);
                        }

                        break;
                    }
                }
            }
        }

#if PROTECT_VECTORS
        MYWRITERLOCK_RELEASE(l);
#endif
        selectionChanged ();
    }

    setScrollPosition (h, v);
}

void ThumbBrowserBase::selectNext (int distance, bool enlarge)
{
    double h, v;
    getScrollPosition (h, v);

    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        if (!selected.empty ()) {
            std::vector<ThumbBrowserEntryBase*>::iterator front = std::find (fd.begin (), fd.end (), selected.front ());
            std::vector<ThumbBrowserEntryBase*>::iterator back = std::find (fd.begin (), fd.end (), selected.back ());
            std::vector<ThumbBrowserEntryBase*>::iterator last = std::find (fd.begin (), fd.end (), lastClicked);

            if (front > back) {
                std::swap(front, back);
            }

            std::vector<ThumbBrowserEntryBase*>::iterator& curr = last == back ? back : front;

            // find next thumbnail at filtered distance after current
            for (; curr < fd.end (); ++curr) {
                if (!(*curr)->filtered) {
                    if (distance-- == 0) {
                        // clear current selection
                        for (size_t i = 0; i < selected.size (); ++i) {
                            selected[i]->selected = false;
                            redrawNeeded (selected[i]);
                        }

                        selected.clear ();

                        // make sure the newly selected thumbnail is visible and make it current
                        scrollToEntry (h, v, internal.get_width (), internal.get_height (), *curr);
                        lastClicked = *curr;

                        // either enlarge current selection or set new selection
                        if(enlarge) {
                            // reverse direction if distance is too large
                            if(front > back) {
                                std::swap(front, back);
                            }

                            for (; front <= back; ++front) {
                                if (!(*front)->filtered) {
                                    (*front)->selected = true;
                                    redrawNeeded (*front);
                                    selected.push_back (*front);
                                }
                            }
                        } else {
                            (*curr)->selected = true;
                            redrawNeeded (*curr);
                            selected.push_back (*curr);
                        }

                        break;
                    }
                }
            }
        }

#if PROTECT_VECTORS
        MYWRITERLOCK_RELEASE(l);
#endif
        selectionChanged ();
    }

    setScrollPosition (h, v);
}

void ThumbBrowserBase::selectFirst (bool enlarge)
{
    double h, v;
    getScrollPosition (h, v);

    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        if (!fd.empty ()) {
            // find first unfiltered entry
            std::vector<ThumbBrowserEntryBase*>::iterator first = fd.begin ();

            for (; first < fd.end (); ++first) {
                if (!(*first)->filtered) {
                    break;
                }
            }

            scrollToEntry (h, v, internal.get_width (), internal.get_height (), *first);

            ThumbBrowserEntryBase* lastEntry = lastClicked;
            lastClicked = *first;

            if(selected.empty ()) {
                (*first)->selected = true;
                redrawNeeded (*first);
                selected.push_back (*first);
            } else {
                std::vector<ThumbBrowserEntryBase*>::iterator back = std::find (fd.begin (), fd.end (), lastEntry ? lastEntry : selected.back ());

                if (first > back) {
                    std::swap(first, back);
                }

                // clear current selection
                for (size_t i = 0; i < selected.size (); ++i) {
                    selected[i]->selected = false;
                    redrawNeeded (selected[i]);
                }

                selected.clear ();

                // either enlarge current selection or set new selection
                for (; first <= back; ++first) {
                    if (!(*first)->filtered) {
                        (*first)->selected = true;
                        redrawNeeded (*first);
                        selected.push_back (*first);
                    }

                    if (!enlarge) {
                        break;
                    }
                }
            }
        }

#if PROTECT_VECTORS
        MYWRITERLOCK_RELEASE(l);
#endif
        selectionChanged ();
    }

    setScrollPosition (h, v);
}

void ThumbBrowserBase::selectLast (bool enlarge)
{
    double h, v;
    getScrollPosition (h, v);

    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        if (!fd.empty ()) {
            // find last unfiltered entry
            std::vector<ThumbBrowserEntryBase*>::iterator last = fd.end () - 1;

            for (; last >= fd.begin (); --last) {
                if (!(*last)->filtered) {
                    break;
                }
            }

            scrollToEntry (h, v, internal.get_width (), internal.get_height (), *last);

            ThumbBrowserEntryBase* lastEntry = lastClicked;
            lastClicked = *last;

            if(selected.empty()) {
                (*last)->selected = true;
                redrawNeeded (*last);
                selected.push_back (*last);
            } else {
                std::vector<ThumbBrowserEntryBase*>::iterator front = std::find (fd.begin (), fd.end (), lastEntry ? lastEntry : selected.front ());

                if (last < front) {
                    std::swap(last, front);
                }

                // clear current selection
                for (size_t i = 0; i < selected.size (); ++i) {
                    selected[i]->selected = false;
                    redrawNeeded (selected[i]);
                }

                selected.clear ();

                // either enlarge current selection or set new selection
                for (; front <= last; --last) {
                    if (!(*last)->filtered) {
                        (*last)->selected = true;
                        redrawNeeded (*last);
                        selected.push_back (*last);
                    }

                    if (!enlarge) {
                        break;
                    }
                }

                std::reverse(selected.begin (), selected.end ());
            }
        }

#if PROTECT_VECTORS
        MYWRITERLOCK_RELEASE(l);
#endif
        selectionChanged ();
    }

    setScrollPosition (h, v);
}

void ThumbBrowserBase::resizeThumbnailArea (int w, int h)
{

    inW = w;
    inH = h;

    if (hscroll.get_value() + internal.get_width() > inW) {
        hscroll.set_value (inW - internal.get_width());
    }

    if (vscroll.get_value() + internal.get_height() > inH) {
        vscroll.set_value (inH - internal.get_height());
    }

    configScrollBars ();
}

void ThumbBrowserBase::internalAreaResized (Gtk::Allocation& req)
{

    if (inW > 0 && inH > 0) {
        configScrollBars ();
        redraw ();
    }
}

void ThumbBrowserBase::configScrollBars ()
{

    // HOMBRE:DELETE ME?
    GThreadLock tLock; // Acquire the GUI

    if (inW > 0 && inH > 0) {

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

        if(iw >= inW) {
            hscroll.hide();
        } else {
            hscroll.show();
        }

        if(ih >= inH) {
            vscroll.hide();
        } else {
            vscroll.show();
        }
    }
}

void ThumbBrowserBase::arrangeFiles ()
{

#if PROTECT_VECTORS
    MYREADERLOCK(l, entryRW);
#endif

    // GUI already locked by ::redraw, the only caller of this method for now.
    // We could lock it one more time, there's no harm excepted (negligible) speed penalty
    //GThreadLock lock;

    int N = fd.size ();

    // apply filter
    for (int i = 0; i < N; i++) {
        fd[i]->filtered = !checkFilter (fd[i]);
    }

    int rowHeight = 0;

    // compute size of the items
    for (int i = 0; i < N; i++)
        if (!fd[i]->filtered && fd[i]->getMinimalHeight() > rowHeight) {
            rowHeight = fd[i]->getMinimalHeight ();
        }

    if (arrangement == TB_Horizontal) {
        numOfCols = 1;
        int numOfRows = 1;
//        if (rowHeight>0) {
//            numOfRows = (internal.get_height()+rowHeight/2)/rowHeight;
//            if (numOfRows<1)
//                numOfRows = 1;
//        }

        int ct = 0;
        int currx = 0;
        int curry = 0;

        while (ct < N) {
            // find widest item in the column
            int maxw = 0;

            for (int i = 0; ct + i < N && i < numOfRows; i++)
                if (fd[ct + i]->getMinimalWidth() > maxw) {
                    maxw = fd[ct + i]->getMinimalWidth ();
                }

            // arrange items in the column
            curry = 0;

            for (int i = 0; ct < N && i < numOfRows; i++, ct++) {
                while (ct < N && fd[ct]->filtered) {
                    fd[ct++]->drawable = false;
                }

                if (ct < N) {
                    fd[ct]->setPosition (currx, curry, maxw, rowHeight);
                    fd[ct]->drawable = true;
                    curry += rowHeight;
                }
            }

            currx += maxw;
        }

#if PROTECT_VECTORS
        MYREADERLOCK_RELEASE(l);
#endif
        // This will require a Writer access
        resizeThumbnailArea (currx, numOfRows * rowHeight);
    } else {
        int availWidth = internal.get_width();
        // initial number of columns
        numOfCols = 0;
        int colsWidth = 0;

        for (int i = 0; i < N; i++)
            if (!fd[i]->filtered && colsWidth + fd[i]->getMinimalWidth() <= availWidth) {
                colsWidth += fd[numOfCols]->getMinimalWidth ();
                numOfCols++;
            }

        if (numOfCols < 1) {
            numOfCols = 1;
        }

        std::vector<int> colWidths;

        for (; numOfCols > 0; numOfCols--) {
            // compute column widths
            colWidths.resize (numOfCols);

            for (int i = 0; i < numOfCols; i++) {
                colWidths[i] = 0;
            }

            for (int i = 0, j = 0; i < N; i++) {
                if (!fd[i]->filtered && fd[i]->getMinimalWidth() > colWidths[j % numOfCols]) {
                    colWidths[j % numOfCols] = fd[i]->getMinimalWidth ();
                }

                if (!fd[i]->filtered) {
                    j++;
                }
            }

            // if not wider than the space available, arrange it and we are ready
            colsWidth = 0;

            for (int i = 0; i < numOfCols; i++) {
                colsWidth += colWidths[i];
            }

            if (numOfCols == 1 || colsWidth < availWidth) {
                break;
            }
        }

        // arrange files
        int ct = 0;
        int currx = 0;
        int curry = 0;

        while (ct < N) {
            // arrange items in the row
            currx = 0;

            for (int i = 0; ct < N && i < numOfCols; i++, ct++) {
                while (ct < N && fd[ct]->filtered) {
                    fd[ct++]->drawable = false;
                }

                if (ct < N) {
                    fd[ct]->setPosition (currx, curry, colWidths[i % numOfCols], rowHeight);
                    fd[ct]->drawable = true;
                    currx += colWidths[i % numOfCols];
                }
            }

            if (currx > 0) { // there were thumbnails placed in the row
                curry += rowHeight;
            }
        }

#if PROTECT_VECTORS
        MYREADERLOCK_RELEASE(l);
#endif
        // This will require a Writer access
        resizeThumbnailArea (colsWidth, curry);
    }
}

void ThumbBrowserBase::disableInspector()
{
    if (inspector) {
        inspector->setActive(false);
    }
}

void ThumbBrowserBase::enableInspector()
{
    if (inspector) {
        inspector->setActive(true);
    }
}

void ThumbBrowserBase::Internal::on_realize()
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)
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

bool ThumbBrowserBase::Internal::on_query_tooltip (int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)
    Glib::ustring ttip = "";

    {
#if PROTECT_VECTORS
        MYREADERLOCK(l, parent->entryRW);
#endif

        for (size_t i = 0; i < parent->fd.size(); i++)
            if (parent->fd[i]->drawable && parent->fd[i]->inside (x, y)) {
                ttip = parent->fd[i]->getToolTip (x, y);
                break;
            }
    }

    if (ttip != "") {
        tooltip->set_markup (ttip);
        return true;
    } else {
        return false;
    }
}

void ThumbBrowserBase::on_style_changed (const Glib::RefPtr<Gtk::Style>& style)
{
    // GUI will be acquired by refreshThumbImages
    refreshThumbImages ();
}

ThumbBrowserBase::Internal::Internal () : ofsX(0), ofsY(0), parent(NULL), dirty(true)
{
}

void ThumbBrowserBase::Internal::setParent (ThumbBrowserBase* p)
{
    parent = p;
}

void ThumbBrowserBase::Internal::setPosition (int x, int y)
{
    ofsX = x;
    ofsY = y;
}

bool ThumbBrowserBase::Internal::on_key_press_event (GdkEventKey* event)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)
    return parent->keyPressed (event);
}

bool ThumbBrowserBase::Internal::on_button_press_event (GdkEventButton* event)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)
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

void ThumbBrowserBase::buttonPressed (int x, int y, int button, GdkEventType type, int state, int clx, int cly, int clw, int clh)
{
    // GUI already acquired

    ThumbBrowserEntryBase* fileDescr = NULL;
    bool handled = false;

    {
        IFPV_MYREADERLOCK(l, entryRW);

        for (size_t i = 0; i < fd.size(); i++)
            if (fd[i]->drawable) {
                if (fd[i]->inside (x, y) && fd[i]->insideWindow (clx, cly, clw, clh)) {
                    fileDescr = fd[i];
                }

                bool b = fd[i]->pressNotify (button, type, state, x, y);
                handled = handled || b;
            }
    }

    if (handled || (fileDescr && fileDescr->processing)) {
        return;
    }

    {
        IFPV_MYWRITERLOCK(l, entryRW);

        if (selected.size() == 1 && type == GDK_2BUTTON_PRESS && button == 1) {
            doubleClicked (selected[0]);
        } else if (button == 1 && type == GDK_BUTTON_PRESS) {
            if (fileDescr && (state & GDK_SHIFT_MASK))
                selectRange (fileDescr);
            else if (fileDescr && (state & GDK_CONTROL_MASK))
                selectSet (fileDescr);
            else
                selectSingle (fileDescr);

            lastClicked = fileDescr;
            IFPV_MYWRITERLOCK_RELEASE(l);
            selectionChanged ();
        } else if (fileDescr && button == 3 && type == GDK_BUTTON_PRESS) {
            if (!fileDescr->selected) {
                selectSingle (fileDescr);

                lastClicked = fileDescr;
                IFPV_MYWRITERLOCK_RELEASE(l);
                selectionChanged ();
            }

            IFPV_MYWRITERLOCK_RELEASE(l);
            rightClicked (fileDescr);
        }
    } // end of MYWRITERLOCK(l, entryRW);

}

bool ThumbBrowserBase::Internal::on_expose_event(GdkEventExpose* event)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)

    dirty = false;

    Glib::RefPtr<Gdk::Window> window = get_window();

    int w = get_width();
    int h = get_height();

    window->clear();
    // draw thumbnails
    Glib::RefPtr<Pango::Context> context = get_pango_context ();
    context->set_font_description (get_style()->get_font());

    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, parent->entryRW);
#endif

        for (size_t i = 0; i < parent->fd.size() && !dirty; i++) { // if dirty meanwhile, cancel and wait for next redraw
            if (!parent->fd[i]->drawable || !parent->fd[i]->insideWindow (0, 0, w, h)) {
                parent->fd[i]->updatepriority = false;
            } else {
                parent->fd[i]->updatepriority = true;
                parent->fd[i]->draw ();
            }
        }
    }

    return true;
}

bool ThumbBrowserBase::Internal::on_button_release_event (GdkEventButton* event)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)
    int w = get_width();
    int h = get_height();

#if PROTECT_VECTORS
    MYREADERLOCK(l, parent->entryRW);
#endif

    for (size_t i = 0; i < parent->fd.size(); i++)
        if (parent->fd[i]->drawable && parent->fd[i]->insideWindow (0, 0, w, h)) {
            ThumbBrowserEntryBase* tbe = parent->fd[i];
#if PROTECT_VECTORS
            MYREADERLOCK_RELEASE(l);
#endif
            // This will require a Writer access...
            tbe->releaseNotify (event->button, event->type, event->state, (int)event->x, (int)event->y);
#if PROTECT_VECTORS
            MYREADERLOCK_ACQUIRE(l);
#endif
        }

    return true;
}

bool ThumbBrowserBase::Internal::on_motion_notify_event (GdkEventMotion* event)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)
    int w = get_width();
    int h = get_height();

#if PROTECT_VECTORS
    MYREADERLOCK(l, parent->entryRW);
#endif

    for (size_t i = 0; i < parent->fd.size(); i++)
        if (parent->fd[i]->drawable && parent->fd[i]->insideWindow (0, 0, w, h)) {
            /*#if PROTECT_VECTORS
            MYREADERLOCK_RELEASE(l); // motionNotify calls the queue, which locks
            #endif*/
            parent->fd[i]->motionNotify ((int)event->x, (int)event->y);
        }

    return true;
}

bool ThumbBrowserBase::Internal::on_scroll_event (GdkEventScroll* event)
{
    // Gtk signals automatically acquire the GUI (i.e. this method is enclosed by gdk_thread_enter and gdk_thread_leave)

    parent->scroll (event->direction);
    return true;
}


void ThumbBrowserBase::redraw ()
{

    GThreadLock lock;
    arrangeFiles ();
    queue_draw ();
}

void ThumbBrowserBase::zoomChanged (bool zoomIn)
{

    int newHeight = 0;
    int optThumbSize = getThumbnailHeight();

    if (zoomIn)
        for (size_t i = 0; i < options.thumbnailZoomRatios.size(); i++) {
            newHeight = (int)(options.thumbnailZoomRatios[i] * getMaxThumbnailHeight());

            if (newHeight > optThumbSize) {
                break;
            }
        }
    else
        for (size_t i = options.thumbnailZoomRatios.size() - 1; i > 0; i--) {
            newHeight = (int)(options.thumbnailZoomRatios[i] * getMaxThumbnailHeight());

            if (newHeight < optThumbSize) {
                break;
            }
        }

    previewHeight = newHeight;

    saveThumbnailHeight(newHeight);

    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        for (size_t i = 0; i < fd.size(); i++) {
            fd[i]->resize (previewHeight);
        }
    }

    redraw ();
#ifdef WIN32
    gdk_window_process_updates (get_window()->gobj(), true);
#endif
}

void ThumbBrowserBase::refreshThumbImages ()
{

    int previewHeight = getThumbnailHeight();
    {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        for (size_t i = 0; i < fd.size(); i++) {
            fd[i]->resize (previewHeight);
        }
    }

    redraw ();
}

void ThumbBrowserBase::refreshQuickThumbImages ()
{
#if PROTECT_VECTORS
    MYWRITERLOCK(l, entryRW);
#endif

    for (size_t i = 0; i < fd.size(); ++i) {
        fd[i]->refreshQuickThumbnailImage ();
    }
}

void ThumbBrowserBase::refreshEditedState (const std::set<Glib::ustring>& efiles)
{

    editedFiles = efiles;
    {
#if PROTECT_VECTORS
        MYREADERLOCK(l, entryRW);
#endif

        for (size_t i = 0; i < fd.size(); i++) {
            fd[i]->framed = editedFiles.find (fd[i]->filename) != editedFiles.end();
        }
    }

    queue_draw ();
}

void ThumbBrowserBase::setArrangement (Arrangement a)
{

    arrangement = a;
    redraw ();
}

void ThumbBrowserBase::enableTabMode(bool enable)
{
    location = enable ? THLOC_EDITOR : THLOC_FILEBROWSER;
    arrangement = enable ? ThumbBrowserBase::TB_Horizontal : ThumbBrowserBase::TB_Vertical;

    if ((!options.sameThumbSize && (options.thumbSizeTab != options.thumbSize)) || (options.showFileNames || options.filmStripShowFileNames)) {
#if PROTECT_VECTORS
        MYWRITERLOCK(l, entryRW);
#endif

        for (size_t i = 0; i < fd.size(); i++) {
            fd[i]->resize (getThumbnailHeight());
        }
    }

    redraw ();

    // Scroll to selected position if going into ribbon mode or back
    // Tab mode is horizontal, file browser is vertical
    {
#if PROTECT_VECTORS
        MYREADERLOCK(l, entryRW);
#endif

        if (!selected.empty()) {
            if (enable) {
                double h = selected[0]->getStartX();
#if PROTECT_VECTORS
                MYREADERLOCK_RELEASE(l);
#endif
                hscroll.set_value (min(h, hscroll.get_adjustment()->get_upper()));
            } else {
                double v = selected[0]->getStartY();
#if PROTECT_VECTORS
                MYREADERLOCK_RELEASE(l);
#endif
                vscroll.set_value (min(v, vscroll.get_adjustment()->get_upper()));
            }
        }
    }
}

void ThumbBrowserBase::initEntry (ThumbBrowserEntryBase* entry)
{
    entry->setOffset ((int)(hscroll.get_value()), (int)(vscroll.get_value()));
}

void ThumbBrowserBase::getScrollPosition (double& h, double& v)
{
    h = hscroll.get_value ();
    v = vscroll.get_value ();
}

void ThumbBrowserBase::setScrollPosition (double h, double v)
{
    hscroll.set_value (h > hscroll.get_adjustment()->get_upper() ? hscroll.get_adjustment()->get_upper() : h);
    vscroll.set_value (v > vscroll.get_adjustment()->get_upper() ? vscroll.get_adjustment()->get_upper() : v);
}

// needed for auto-height in single tab
int ThumbBrowserBase::getEffectiveHeight()
{
    int h = hscroll.get_height() + 2; // have 2 pixels rounding error for scroll bars to appear

#if PROTECT_VECTORS
    MYREADERLOCK(l, entryRW);
#endif

    // Filtered items do not change in size, so take a non-filtered
    for (size_t i = 0; i < fd.size(); i++)
        if (!fd[i]->filtered) {
            h += fd[i]->getEffectiveHeight();
            break;
        }

    return h;
}

void ThumbBrowserBase::redrawNeeded (ThumbBrowserEntryBase* entry)
{

    // HOMBRE:DELETE ME?
    GThreadLock tLock; // Acquire the GUI

    if (entry->insideWindow (0, 0, internal.get_width(), internal.get_height())) {
        if (!internal.isDirty ()) {
            internal.setDirty ();
            internal.queue_draw ();
        }
    }
}


