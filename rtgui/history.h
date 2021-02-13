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
#pragma once

#include <gtkmm.h>

#include "paramsedited.h"
#include "pparamschangelistener.h"
#include "profilechangelistener.h"

#include "../rtengine/rtengine.h"

class HistoryBeforeLineListener
{
public:
    virtual ~HistoryBeforeLineListener() = default;
    virtual void historyBeforeLineChanged(const rtengine::procparams::ProcParams& params) = 0;
};

class History : public Gtk::Box, public PParamsChangeListener
{

public:

    class HistoryColumns : public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring>  text;
        Gtk::TreeModelColumn<Glib::ustring>  value;
        Gtk::TreeModelColumn<rtengine::procparams::ProcParams>     params;
        Gtk::TreeModelColumn<rtengine::ProcEvent>    chev;
        Gtk::TreeModelColumn<ParamsEdited>     paramsEdited;
        HistoryColumns()
        {
            add(text);
            add(value);
            add(chev);
            add(params);
            add(paramsEdited);
        }
    };
    HistoryColumns historyColumns;
    class BookmarkColumns : public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring>  text;
        Gtk::TreeModelColumn<rtengine::procparams::ProcParams>     params;
        Gtk::TreeModelColumn<ParamsEdited>     paramsEdited;
        BookmarkColumns()
        {
            add(text);
            add(params);
            add(paramsEdited);
        }
    };
    BookmarkColumns bookmarkColumns;

protected:
    Gtk::Paned*            historyVPaned;
    Gtk::TreeView*          hTreeView;
    Glib::RefPtr<Gtk::ListStore> historyModel;

    Gtk::ScrolledWindow*    bscrollw;
    Gtk::TreeView*          bTreeView;
    Glib::RefPtr<Gtk::ListStore> bookmarkModel;

    Gtk::Button*            addBookmark;
    Gtk::Button*            delBookmark;

    sigc::connection        selchangehist;
    sigc::connection        selchangebm;

    HistoryBeforeLineListener * blistener;
    ProfileChangeListener* tpc;
    ParamsEdited defParamsEdited;
    int bmnum;

    bool on_query_tooltip(int x, int y, bool keyboard_tooltip, const Glib::RefPtr<Gtk::Tooltip>& tooltip);

public:

    explicit History (bool bookmarkSupport = true);

    void setProfileChangeListener     (ProfileChangeListener* tpc_)
    {
        tpc = tpc_;
    }
    void setHistoryBeforeLineListener (HistoryBeforeLineListener* bll)
    {
        blistener = bll;
    }

    // pparamschangelistener interface
    void procParamsChanged(
        const rtengine::procparams::ProcParams* params,
        const rtengine::ProcEvent& ev,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr
    ) override;
    void clearParamChanges () override;

    void historySelectionChanged ();
    void bookmarkSelectionChanged ();
    void initHistory ();

    bool getBeforeLineParams (rtengine::procparams::ProcParams& params);

    void addBookmarkWithText (Glib::ustring text);
    void addBookmarkPressed ();
    void delBookmarkPressed ();

    //void resized (Gtk::Allocation& req);

    void undo ();
    void redo ();

    bool blistenerLock;
    void resetSnapShotNumber()
    {
        bmnum = 1;
    };
};
