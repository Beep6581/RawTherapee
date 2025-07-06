#pragma once

#include <memory>

#include <gtkmm.h>

#include <glibmm/ustring.h>

#include "adjuster.h"
#include "guiutils.h"
#include "toolpanel.h"

class ClutComboBox final :
    public MyComboBox
{
public:
    explicit ClutComboBox(const Glib::ustring &path);
    //int fillFromDir (const Glib::ustring& path);
    int foundClutsCount() const;
    Glib::ustring getSelectedClut();
    void setSelectedClut( Glib::ustring filename );
    void setBatchMode(bool yes);

    static void cleanup();

private:
    void updateUnchangedEntry(); // in batchMode we need to add an extra entry "(Unchanged)". We do this whenever the widget is mapped (connecting to signal_map()), unless options.multiDisplayMode (see the comment below about cm2 in this case)

    class ClutColumns : public Gtk::TreeModel::ColumnRecord
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring> label;
        Gtk::TreeModelColumn<Glib::ustring> clutFilename;
        ClutColumns();
    };

    class ClutModel {
    public:
        Glib::RefPtr<Gtk::TreeStore> m_model;
        ClutColumns m_columns;
        int count;
        explicit ClutModel(const Glib::ustring &path);
        int parseDir (const Glib::ustring& path);
    };

    Glib::RefPtr<Gtk::TreeStore> &m_model();
    ClutColumns &m_columns();

    Gtk::TreeIter findRowByClutFilename(  Gtk::TreeModel::Children childs, Glib::ustring filename );

    static std::unique_ptr<ClutModel> cm; // we use a shared TreeModel for all the combo boxes, to save time (no need to reparse the clut dir multiple times)...
    static std::unique_ptr<ClutModel> cm2; // ... except when options.multiDisplayMode (i.e. editors in their own window), where we need two. This is because we might have two combo boxes displayed at the same time in this case
    bool batchMode;
};

class FilmSimulation : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{
public:
    static const Glib::ustring TOOL_NAME;

    FilmSimulation();

    void adjusterChanged(Adjuster* a, double newval) override;
    void setBatchMode(bool batchMode) override;
    void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setAdjusterBehavior(bool strength);
    void trimValues(rtengine::procparams::ProcParams* pp) override;

private:
    void onClutSelected();
    void enabledChanged() override;

    void updateDisable( bool value );

    ClutComboBox *m_clutComboBox;
    sigc::connection m_clutComboBoxConn;
    Glib::ustring m_oldClutFilename;

    Adjuster *m_strength;
};
