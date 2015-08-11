#ifndef FILM_SIMULATION_INCLUDED
#define FILM_SIMULATION_INCLUDED

#include <gtkmm.h>
#include <glibmm.h>
#include <memory>
#include "toolpanel.h"
#include "guiutils.h"
#include "adjuster.h"

class ClutComboBox : public MyComboBox
{
public:
    int fillFromDir( Glib::ustring path );
    Glib::ustring getSelectedClut();
    void setSelectedClut( Glib::ustring filename );
    void addUnchangedEntry();

private:
    class ClutColumns : public Gtk::TreeModel::ColumnRecord 
    {
    public:
        Gtk::TreeModelColumn<Glib::ustring> label;
        Gtk::TreeModelColumn<Glib::ustring> clutFilename;
        ClutColumns();
    };

    int parseDir( Glib::ustring path, Gtk::TreeModel::Row *parentRow );
    Gtk::TreeIter findRowByClutFilename(  Gtk::TreeModel::Children childs, Glib::ustring filename );

    Glib::RefPtr<Gtk::TreeStore> m_model;
    ClutColumns m_columns;
};

class FilmSimulation : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{
public:
    FilmSimulation();

    void adjusterChanged( Adjuster* a, double newval );
    void setBatchMode( bool batchMode );
    void read( const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = NULL );
    void write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = NULL );
    void setAdjusterBehavior( bool strength );
    void trimValues( rtengine::procparams::ProcParams* pp );

private:
    void onClutSelected();
    void enabledChanged();

    void updateDisable( bool value );

    ClutComboBox *m_clutComboBox;
    sigc::connection m_clutComboBoxConn;
    Glib::ustring m_oldClutFilename;

    Adjuster *m_strength;
};

#endif
