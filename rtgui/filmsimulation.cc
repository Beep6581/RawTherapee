#include "filmsimulation.h"
#include "options.h"
#include "../rtengine/clutstore.h"
#include "../rtengine/safegtk.h"

using namespace rtengine;
using namespace rtengine::procparams;

typedef std::vector<Glib::ustring> Strings;

FilmSimulation::FilmSimulation()
    :   FoldableToolPanel( this, "filmsimulation", M("TP_FILMSIMULATION_LABEL"), false, true )
{
    m_clutComboBox = Gtk::manage( new ClutComboBox() );
    int foundClutsCount = m_clutComboBox->fillFromDir( options.clutsDir );

    if ( foundClutsCount == 0 ) {
        pack_start( *Gtk::manage( new Gtk::Label( M("TP_FILMSIMULATION_ZEROCLUTSFOUND") ) ) );
    }

    m_clutComboBoxConn = m_clutComboBox->signal_changed().connect( sigc::mem_fun( *this, &FilmSimulation::onClutSelected ) );
    pack_start( *m_clutComboBox );

    m_strength = Gtk::manage( new Adjuster( M("TP_FILMSIMULATION_STRENGTH"), 0., 100, 1., 100 ) );
    m_strength->setAdjusterListener( this );

    pack_start( *m_strength, Gtk::PACK_SHRINK, 0 );

}

void FilmSimulation::onClutSelected()
{
    Glib::ustring currentClutFilename = m_clutComboBox->getSelectedClut();

    if ( getEnabled() && !currentClutFilename.empty() && listener && currentClutFilename != m_oldClutFilename ) {
        Glib::ustring clutName, dummy;
        splitClutFilename( currentClutFilename, clutName, dummy, dummy );
        listener->panelChanged( EvFilmSimulationFilename, clutName );

        m_oldClutFilename = currentClutFilename;
    }
}

void FilmSimulation::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvFilmSimulationEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvFilmSimulationEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvFilmSimulationEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void FilmSimulation::adjusterChanged( Adjuster* a, double newval )
{
    if (listener && (multiImage || getEnabled()) ) {
        Glib::ustring value = a->getTextValue();
        listener->panelChanged ( EvFilmSimulationStrength, value );
    }
}

void FilmSimulation::setBatchMode( bool batchMode )
{
    ToolPanel::setBatchMode( batchMode );
    m_clutComboBox->addUnchangedEntry();
}

void FilmSimulation::read( const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited )
{
    //copypasted from lensprofile.cc & sharpening.cc
    disableListener();
    updateDisable( true );

    setEnabled (pp->filmSimulation.enabled);

    if ( !pp->filmSimulation.clutFilename.empty() ) {
        m_clutComboBox->setSelectedClut( pp->filmSimulation.clutFilename );
    }

    m_strength->setValue( pp->filmSimulation.strength );

    if (pedited) {
        set_inconsistent (multiImage && !pedited->filmSimulation.enabled);
        m_strength->setEditedState (pedited->filmSimulation.strength ? Edited : UnEdited);

        if ( !pedited->filmSimulation.clutFilename ) {
            m_clutComboBox->setSelectedClut("NULL");
        }
    }

    if ( !get_inconsistent() && !pp->filmSimulation.enabled ) {
        if (options.clutCacheSize == 1) {
            clutStore.clearCache();
        }
    }

    updateDisable( false );
    enableListener();
}

void FilmSimulation::updateDisable( bool value )
{
    m_clutComboBoxConn.block( value );
}

void FilmSimulation::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited )
{
    if ( pedited ) {
        pedited->filmSimulation.enabled  = !get_inconsistent();
        pedited->filmSimulation.strength = m_strength->getEditedState ();
        pedited->filmSimulation.clutFilename = m_clutComboBox->getSelectedClut() != "NULL";
    }

    pp->filmSimulation.enabled = getEnabled();
    Glib::ustring clutFName = m_clutComboBox->getSelectedClut();

    if ( clutFName != "NULL" ) { // We do not want to set "NULL" in clutFilename, even if "unedited"
        pp->filmSimulation.clutFilename = clutFName;
    }

    pp->filmSimulation.strength = m_strength->getValue();
}

void FilmSimulation::setAdjusterBehavior( bool strength )
{
    m_strength->setAddMode( strength );
}

void FilmSimulation::trimValues( rtengine::procparams::ProcParams* pp )
{
    pp->filmSimulation.strength = m_strength->trimValue( pp->filmSimulation.strength );
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ClutComboBox::ClutColumns::ClutColumns()
{
    add( label );
    add( clutFilename );
}

int ClutComboBox::fillFromDir( Glib::ustring path )
{
    int result = 0;

    if ( !path.empty() ) {
        m_model.clear();
        m_model = Gtk::TreeStore::create( m_columns );
        set_model( m_model );
        result = parseDir( path, 0 );
        pack_start( m_columns.label, false );
    }

    return result;
}

Gtk::TreeIter appendToModel( Glib::RefPtr<Gtk::TreeStore> model, Gtk::TreeModel::Row *parent )
{
    Gtk::TreeIter result;

    if ( parent ) {
        result = model->append( parent->children() );

    } else {
        result = model->append();
    }

    return result;
}

int ClutComboBox::parseDir( Glib::ustring path, Gtk::TreeModel::Row *parentRow )
{
    int result = 0;

    if ( path.empty() || !safe_file_test( path, Glib::FILE_TEST_EXISTS ) || !safe_file_test ( path, Glib::FILE_TEST_IS_DIR ) ) {
        return result;
    }

    Glib::Dir* dir = new Glib::Dir( path );

    Strings names;

    for( Glib::DirIterator it = dir->begin(); it != dir->end(); ++it ) {
        Glib::ustring current = *it;

        if ( current != "." && current != ".." ) {
            names.push_back( current );
        }
    }

    std::sort( names.begin(), names.end() );

    for ( Strings::iterator it = names.begin(); it != names.end(); ++it ) {
        Glib::ustring current = *it;
        Glib::ustring fullname = Glib::build_filename( path, current );

        if ( safe_file_test( fullname, Glib::FILE_TEST_IS_DIR ) ) {

            Gtk::TreeModel::Row newFolderMenu = *appendToModel( m_model, parentRow );
            newFolderMenu[ m_columns.label ] = current;
            result += parseDir( fullname, &newFolderMenu );
        } else {
            Glib::ustring name, extension, profileName;
            splitClutFilename( current, name, extension, profileName );

            if ( extension == "tif" ||
                    extension == "TIF" ||
                    extension == "png" ||
                    extension == "PNG" ) {
                Gtk::TreeModel::Row newClut = *appendToModel( m_model, parentRow );
                newClut[ m_columns.label ] = name;
                newClut[ m_columns.clutFilename ] = fullname;
                ++result;
            }
        }
    }

    return result;
}

Glib::ustring ClutComboBox::getSelectedClut()
{
    Glib::ustring result;
    Gtk::TreeModel::iterator current = get_active();
    Gtk::TreeModel::Row row = *current;

    if ( row ) {
        result = row[ m_columns.clutFilename ];
    }

    return result;
}

void ClutComboBox::setSelectedClut( Glib::ustring filename )
{
    if ( !filename.empty() ) {
        Gtk::TreeIter found = findRowByClutFilename( m_model->children(), filename );

        if ( found ) {
            set_active( found );
        }
    }
}

Gtk::TreeIter ClutComboBox::findRowByClutFilename( Gtk::TreeModel::Children childs, Glib::ustring filename )
{
    Gtk::TreeIter result = childs.end();

    for( Gtk::TreeModel::Children::iterator it = childs.begin(); !result && it != childs.end(); ++it ) {
        Gtk::TreeModel::Row row = *it;

        if ( row[ m_columns.clutFilename ] == filename ) {
            result = it;
        } else {
            result = findRowByClutFilename( it->children(), filename );
        }
    }

    return result;
}

void ClutComboBox::addUnchangedEntry()
{
    Gtk::TreeModel::Row row = *(m_model->append());
    row[m_columns.label] = M("GENERAL_UNCHANGED");
    row[m_columns.clutFilename] = "NULL";
}
