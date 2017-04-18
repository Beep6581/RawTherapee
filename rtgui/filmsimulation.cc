#include "filmsimulation.h"

#include <chrono>

#include "options.h"
#include "../rtengine/clutstore.h"

using namespace rtengine;
using namespace rtengine::procparams;

namespace
{

Glib::ustring stripPrefixDir(const Glib::ustring& filename, const Glib::ustring& dir)
{
    const Glib::ustring full_dir =
        !Glib::str_has_suffix(dir, G_DIR_SEPARATOR_S)
            ? dir + G_DIR_SEPARATOR_S
            : dir;
    return
        Glib::str_has_prefix(filename, full_dir)
            ? filename.substr(full_dir.size())
            : filename;
}

bool notifySlowParseDir (const std::chrono::system_clock::time_point& startedAt)
{
    enum Decision {
        UNDECIDED,
        CANCEL,
        CONTINUE
    };

    static Decision decision = UNDECIDED;

    if (decision == CANCEL) {
        return false;
    } else if (decision == CONTINUE) {
        return true;
    }

    const auto now = std::chrono::system_clock::now();
    if (now - startedAt < std::chrono::seconds(10)) {
        return true;
    }

    Gtk::MessageDialog dialog(M("TP_FILMSIMULATION_SLOWPARSEDIR"), false, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
    if (dialog.run() == Gtk::RESPONSE_YES) {
        decision = CANCEL;
        return false;
    } else {
        decision = CONTINUE;
        return true;
    }
}

}

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
        HaldCLUT::splitClutFilename( currentClutFilename, clutName, dummy, dummy );
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
    updateDisable(true);

    setEnabled(pp->filmSimulation.enabled);

    if (!pp->filmSimulation.clutFilename.empty()) {
        m_clutComboBox->setSelectedClut(
            !Glib::path_is_absolute(pp->filmSimulation.clutFilename)
                ? Glib::ustring(Glib::build_filename(options.clutsDir, pp->filmSimulation.clutFilename))
                : pp->filmSimulation.clutFilename
        );
        m_oldClutFilename = m_clutComboBox->getSelectedClut();
    }

    m_strength->setValue(pp->filmSimulation.strength);

    if (pedited) {
        set_inconsistent (multiImage && !pedited->filmSimulation.enabled);
        m_strength->setEditedState(
            pedited->filmSimulation.strength
                ? Edited
                : UnEdited
        );

        if (!pedited->filmSimulation.clutFilename) {
            m_clutComboBox->setSelectedClut("NULL");
        }
    }

    if (!get_inconsistent() && !pp->filmSimulation.enabled) {
        if (options.clutCacheSize == 1) {
            CLUTStore::getInstance().clearCache();
        }
    }

    updateDisable(false);
    enableListener();
}

void FilmSimulation::updateDisable( bool value )
{
    m_clutComboBoxConn.block( value );
}

void FilmSimulation::write( rtengine::procparams::ProcParams* pp, ParamsEdited* pedited )
{
    if (pedited) {
        pedited->filmSimulation.enabled = !get_inconsistent();
        pedited->filmSimulation.strength = m_strength->getEditedState();
        pedited->filmSimulation.clutFilename = m_clutComboBox->getSelectedClut() != "NULL";
    }

    pp->filmSimulation.enabled = getEnabled();
    const Glib::ustring clutFName = m_clutComboBox->getSelectedClut();

    if (clutFName != "NULL") { // We do not want to set "NULL" in clutFilename, even if "unedited"
        pp->filmSimulation.clutFilename = stripPrefixDir(clutFName, options.clutsDir);
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

int ClutComboBox::fillFromDir (const Glib::ustring& path)
{
    m_model = Gtk::TreeStore::create (m_columns);
    set_model (m_model);

    const auto result = parseDir (path);

    if (result > 0) {
        pack_start (m_columns.label, false);
    }

    return result;
}

int ClutComboBox::parseDir (const Glib::ustring& path)
{
    if (path.empty () || !Glib::file_test (path, Glib::FILE_TEST_IS_DIR)) {
        return 0;
    }

    const auto startedAt = std::chrono::system_clock::now ();

    // Build menu of limited directory structure using breadth-first search
    using Dirs = std::vector<std::pair<Glib::ustring, Gtk::TreeModel::Row>>;
    Dirs dirs;

    {
        Dirs currDirs;
        Dirs nextDirs;

        currDirs.emplace_back (path, Gtk::TreeModel::Row ());

        while (!currDirs.empty ()) {

            for (auto& dir : currDirs) {

                const auto& path = dir.first;
                const auto& row = dir.second;

                try {
                    for (const auto& entry : Glib::Dir (path)) {

                        const auto entryPath = Glib::build_filename (path, entry);

                        if (!Glib::file_test (entryPath, Glib::FILE_TEST_IS_DIR)) {
                            continue;
                        }

                        auto newRow = row ? *m_model->append (row.children ()) : *m_model->append ();
                        newRow[m_columns.label] = entry;

                        nextDirs.emplace_back (entryPath, newRow);
                    }
                } catch (Glib::Exception&) {}

                dirs.push_back (std::move (dir));

                if (!notifySlowParseDir (startedAt)) {
                    m_model->clear ();
                    return 0;
                }
            }

            currDirs.clear ();
            currDirs.swap (nextDirs);
        }
    }

    // Fill menu structure with CLUT files
    Strings entries;

    auto fileCount = 0;

    for (const auto& dir : dirs) {

        const auto& path = dir.first;
        const auto& row = dir.second;

        entries.clear ();

        try {
            for (const auto& entry : Glib::Dir (path)) {

                const auto entryPath = Glib::build_filename (path, entry);

                if (!Glib::file_test (entryPath, Glib::FILE_TEST_IS_REGULAR)) {
                    continue;
                }

                entries.push_back (entryPath);
            }
        } catch (Glib::Exception&) {}

        std::sort (entries.begin (), entries.end ());

        for (const auto& entry : entries) {

            Glib::ustring name, extension, profileName;
            HaldCLUT::splitClutFilename (entry, name, extension, profileName);

            extension = extension.casefold ();
            if (extension.compare ("tif") != 0 && extension.compare ("png") != 0) {
                continue;
            }

            auto newRow = row ? *m_model->append (row.children ()) : *m_model->append ();
            newRow[m_columns.label] = name;
            newRow[m_columns.clutFilename] = entry;

            ++fileCount;

            if (!notifySlowParseDir (startedAt)) {
                m_model->clear ();
                return 0;
            }
        }
    }

    return fileCount;
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
