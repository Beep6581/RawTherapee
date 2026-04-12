/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2025 Daniel Gao <daniel.gao.work@gmail.com>
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

#include "cropguide.h"

#include "eventmapper.h"
#include "guiutils.h"
#include "multilangmgr.h"
#include "paramsedited.h"
#include "rtimage.h"
#include "widgets/basic/colorpreview.h"

#include <algorithm>

namespace {

using CropParams = rtengine::procparams::CropParams;
using CropGuideParams = rtengine::procparams::CropGuideParams;

// clang-format off
constexpr std::array<const char*, 9> GUIDE_TYPE_OPTIONS = {
    "TP_CROP_GTRULETHIRDS",
    "TP_CROP_GTDIAGONALS",
    "TP_CROP_GTHARMMEANS",
    "TP_CROP_GTCROSSHAIR",
    "TP_CROP_GTGRID",
    "TP_CROP_GTTRIANGLE",
    "TP_CROP_GTGOLDENRATIO",
    "TP_CROP_GTEPASSPORT",
    "TP_CROP_GTCENTEREDSQUARE"
};
// clang-format on

// Bleed basis combo box data reusing text from resize tool
constexpr int EMPTY_COMBO_INDEX = -1;
constexpr int INDEX_BASIS_SCALE = 0;
constexpr int INDEX_BASIS_WIDTH = 1;
constexpr int INDEX_BASIS_HEIGHT = 2;
constexpr int INDEX_BASIS_LONG = 3;
constexpr int INDEX_BASIS_SHORT = 4;
constexpr int INDEX_BASIS_UNCHANGED = 5;
constexpr std::array<const char*, 5> BLEED_BASIS = {
    "TP_RESIZE_SCALE",
    "TP_RESIZE_WIDTH",
    "TP_RESIZE_HEIGHT",
    "TP_RESIZE_LONG",
    "TP_RESIZE_SHORT"
};

int mapBasis(CropGuideParams::Basis basis)
{
    using Basis = CropGuideParams::Basis;
    switch (basis) {
        case Basis::SCALE:
            return INDEX_BASIS_SCALE;
        case Basis::WIDTH:
            return INDEX_BASIS_WIDTH;
        case Basis::HEIGHT:
            return INDEX_BASIS_HEIGHT;
        case Basis::LONG:
            return INDEX_BASIS_LONG;
        case Basis::SHORT:
            return INDEX_BASIS_SHORT;
        default:
            return INDEX_BASIS_SCALE;
    }
}

CropGuideParams::Basis mapBasis(int comboIndex)
{
    using Basis = CropGuideParams::Basis;
    switch (comboIndex) {
        case INDEX_BASIS_SCALE:
            return Basis::SCALE;
        case INDEX_BASIS_WIDTH:
            return Basis::WIDTH;
        case INDEX_BASIS_HEIGHT:
            return Basis::HEIGHT;
        case INDEX_BASIS_LONG:
            return Basis::LONG;
        case INDEX_BASIS_SHORT:
            return Basis::SHORT;
        default:
            return Basis::SCALE;
    }
}

void updateImage(Gtk::ToggleButton* button, bool is_enabled)
{
    if (is_enabled) {
        button->set_image(*Gtk::manage(
            new RTImage("power-on-small", Gtk::ICON_SIZE_BUTTON)));
    } else {
        button->set_image(*Gtk::manage(
            new RTImage("power-off-small", Gtk::ICON_SIZE_BUTTON)));
    }
}

} // namespace

const Glib::ustring CropGuide::TOOL_NAME = "cropguide";

CropGuide::CropGuide()
    : FoldableToolPanel(this, TOOL_NAME, M("TP_CROP_GUIDE_LABEL"), false, true),
      m_aspect_ratio_store(Gio::ListStore<AspectRatioModel>::create()),
      m_presets{},
      m_mirror_golden_triangle(false),
      m_dirty_mirror_golden_triangle(false),
      m_rotate_golden_ratio(false),
      m_dirty_rotate_golden_ratio(false),
      m_mirror_golden_ratio(false),
      m_dirty_mirror_golden_ratio(false),
      m_dirty_aspect_ratios(false)
{
    setupEvents();

    setEnabledTooltipMarkup(M("TP_CROP_GUIDE_TOOLTIP"));

    setupPresets();
    pack_start(*Gtk::manage(new Gtk::Separator()));
    setupAspectRatioGuides();

    pack_start(*Gtk::manage(new Gtk::Separator()));
    m_bleed = Gtk::manage(new Adjuster(M("TP_CROP_GUIDE_BLEED"), 0, 10, 1, 0));
    m_bleed->setAdjusterListener(this);
    pack_start(*m_bleed);

    Gtk::Box* basis_box = Gtk::manage(new Gtk::Box());
    // Reuse text from resize tool
    auto basis_label = Gtk::manage(new Gtk::Label(M("TP_RESIZE_APPLIESTO")));
    basis_box->pack_start(*basis_label, false, false);
    m_basis = Gtk::manage(new MyComboBoxText());
    for (auto label : BLEED_BASIS) {
        m_basis->append(M(label));
    }
    m_basis->set_active(INDEX_BASIS_SCALE);
    basis_box->pack_start(*m_basis, true, true);
    pack_start(*basis_box);

    m_basis_changed = m_basis->signal_changed().connect(
        sigc::mem_fun(*this, &CropGuide::onBasisChanged));

    show_all();
}

void CropGuide::setupEvents()
{
    auto m = ProcEventMapper::getInstance();

    EvCropGuideEnabled = m->newEvent(MINUPDATE, "HISTORY_MSG_CROP_GUIDE_ENABLED");
    EvCropGuidePresetChanged =
        m->newEvent(MINUPDATE, "HISTORY_MSG_CROP_GUIDE_PRESET_CHANGED");
    EvCropGuideAspectRatioPresetChanged =
        m->newEvent(MINUPDATE, "HISTORY_MSG_CROP_GUIDE_ASPECT_RATIO_PRESET_CHANGED");
    EvCropGuideBleedChanged =
        m->newEvent(MINUPDATE, "HISTORY_MSG_CROP_GUIDE_BLEED_CHANGED");
    EvCropGuideBasisChanged =
        m->newEvent(MINUPDATE, "HISTORY_MSG_CROP_GUIDE_BLEED_BASIS_CHANGED");
}

void CropGuide::setupPresets()
{
    // N by 7 grid with columns:
    //
    // | Toggle | Label | [Rotate] | [Mirror] | [Undo] | Color Preview | Color Picker |
    //
    // where Label may span multiple columns to fill space
    auto grid = Gtk::manage(new Gtk::Grid());
    grid->set_row_homogeneous(false);
    grid->set_column_homogeneous(false);
    grid->set_hexpand(true);
    grid->set_halign(Gtk::ALIGN_FILL);
    pack_start(*grid, true, true);

    for (size_t i = 0; i < m_presets.size(); i++) {
        const auto type_index = static_cast<CropGuideParams::PresetIndex>(i);
        auto& preset = m_presets[i];

        auto visible_button = Gtk::manage(new Gtk::ToggleButton());
        preset.visibility_button = visible_button;
        updateImage(visible_button, false);
        visible_button->set_relief(Gtk::RELIEF_NONE);
        preset.visibility_conn = visible_button->signal_button_release_event().connect(
            sigc::bind(sigc::mem_fun(*this, &CropGuide::onPresetToggled), i));

        grid->attach(*visible_button, 0, i);

        auto label = Gtk::manage(new Gtk::Label(M(GUIDE_TYPE_OPTIONS.at(i))));
        label->set_hexpand(true);
        label->set_halign(Gtk::ALIGN_START);
        label->set_line_wrap(true);
        // Label is attached depending on preset type later

        auto add_button = [&](const char* icon, int column) {
            auto button = Gtk::manage(new Gtk::Button());
            button->set_image(*Gtk::manage(
                    new RTImage(icon, Gtk::ICON_SIZE_BUTTON)));
            button->set_relief(Gtk::RELIEF_NONE);
            button->set_hexpand(false);
            grid->attach(*button, column, i);

            return button;
        };

        auto color_event_box = Gtk::manage(new Gtk::EventBox());
        preset.preview = Gtk::manage(new ColorPreview());
        preset.preview->set_hexpand(false);
        preset.preview->set_halign(Gtk::ALIGN_CENTER);
        preset.preview->set_vexpand(false);
        preset.preview->set_valign(Gtk::ALIGN_CENTER);
        preset.preview->set_margin_start(4);
        preset.preview->set_margin_end(4);
        preset.preview->setRgb(preset.color.get_red(), preset.color.get_green(),
                               preset.color.get_blue());
        color_event_box->add(*(preset.preview));
        grid->attach(*color_event_box, 5, i);

        auto color_button = add_button("color-picker-small", 6);

        color_event_box->signal_button_release_event().connect(
            [this, i](GdkEventButton* ev) -> bool {
                this->onPresetPickColor(i);
                return true;
            });
        color_button->signal_clicked().connect(sigc::bind(
            sigc::mem_fun(this, &CropGuide::onPresetPickColor), i));

        if (type_index == CropGuideParams::PresetIndex::GOLDEN_TRIANGLE) {
            grid->attach(*label, 1, i, 2);

            auto mirror_button = add_button("flip-horizontal", 3);
            auto undo_button = add_button("undo-small", 4);

            mirror_button->signal_clicked().connect(
                sigc::mem_fun(this, &CropGuide::onGoldenTriangleMirror));
            undo_button->signal_clicked().connect(
                sigc::mem_fun(this, &CropGuide::onGoldenTriangleReset));
        } else if (type_index == CropGuideParams::PresetIndex::GOLDEN_RATIO) {
            grid->attach(*label, 1, i, 1);

            auto rotate_button = add_button("rotate-right-small", 2);
            auto mirror_button = add_button("flip-horizontal", 3);
            auto undo_button = add_button("undo-small", 4);

            rotate_button->signal_clicked().connect(
                sigc::mem_fun(this, &CropGuide::onGoldenRatioRotate));
            mirror_button->signal_clicked().connect(
                sigc::mem_fun(this, &CropGuide::onGoldenRatioMirror));
            undo_button->signal_clicked().connect(
                sigc::mem_fun(this, &CropGuide::onGoldenRatioReset));
        } else {
            grid->attach(*label, 1, i, 6);
        }
    }
}

void CropGuide::setupAspectRatioGuides()
{
    auto box = Gtk::manage(new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    auto label = Gtk::manage(new Gtk::Label(M("TP_CROP_GUIDE_ASPECT_RATIO")));
    label->set_halign(Gtk::ALIGN_START);
    box->pack_start(*label);
    m_available_aspect_ratios_combo = Gtk::manage(new MyComboBoxText());
    m_available_aspect_ratios_conn = m_available_aspect_ratios_combo->signal_changed()
        .connect(sigc::mem_fun(this, &CropGuide::onAspectRatioComboChanged));
    box->pack_start(*m_available_aspect_ratios_combo);
    pack_start(*box);

    {
        std::vector<AspectRatio> ratios = getAspectRatios();

        m_aspect_ratio_presets.reserve(ratios.size());
        size_t index = 0;
        for (auto& r : ratios) {
            m_available_aspect_ratios_combo->append(r.label);

            Glib::RefPtr<AspectRatioModel> preset = AspectRatioModel::create();
            preset->aspect_ratio = std::move(r);
            preset->index = index++;
            m_aspect_ratio_presets.push_back(preset);
        }
    }

    m_aspect_ratio_listbox = Gtk::manage(new Gtk::ListBox());
    m_aspect_ratio_listbox->set_selection_mode(Gtk::SELECTION_NONE);

    auto placeholder = Gtk::manage(new Gtk::Label(
        M("TP_CROP_GUIDE_NO_ASPECT_RATIO_GUIDES")));
    placeholder->show();
    m_aspect_ratio_listbox->set_placeholder(*placeholder);

    m_aspect_ratio_listbox->bind_list_store(
        m_aspect_ratio_store,
        sigc::mem_fun(this, &CropGuide::createAspectRatioModelControls));
    pack_start(*m_aspect_ratio_listbox);
}

Gtk::Widget* CropGuide::createAspectRatioModelControls(
    const Glib::RefPtr<AspectRatioModel>& item)
{
    // N by 6 grid with columns:
    //
    // | Toggle | Label | Rotate | Color | Color Select | Remove |
    //
    // where Label may span multiple columns to fill space
    auto grid = Gtk::manage(new Gtk::Grid());
    grid->set_row_homogeneous(false);
    grid->set_column_homogeneous(false);
    grid->set_hexpand(true);
    grid->set_halign(Gtk::ALIGN_FILL);

    auto visible_button = Gtk::manage(new Gtk::ToggleButton());
    visible_button->set_active(item->visible);
    updateImage(visible_button, item->visible);
    visible_button->set_relief(Gtk::RELIEF_NONE);
    grid->attach(*visible_button, 0, 0);

    visible_button->signal_clicked().connect(sigc::bind(
        sigc::mem_fun(this, &CropGuide::onAspectRatioPresetToggled),
        visible_button, item->index));

    auto label = Gtk::manage(new Gtk::Label(item->aspect_ratio.label));
    label->set_hexpand(true);
    label->set_halign(Gtk::ALIGN_START);
    label->set_line_wrap(true);
    grid->attach(*label, 1, 0);

    auto rotate_button = Gtk::manage(new Gtk::Button());
    rotate_button->set_image(*Gtk::manage(
        new RTImage("rotate-right-small", Gtk::ICON_SIZE_BUTTON)));
    rotate_button->set_relief(Gtk::RELIEF_NONE);
    grid->attach(*rotate_button, 2, 0);

    rotate_button->signal_clicked().connect(sigc::bind(
        sigc::mem_fun(this, &CropGuide::onAspectRatioPresetRotate),
        item->index));

    auto color_event_box = Gtk::manage(new Gtk::EventBox());
    auto color_preview = Gtk::manage(new ColorPreview());
    color_preview->set_hexpand(false);
    color_preview->set_halign(Gtk::ALIGN_CENTER);
    color_preview->set_vexpand(false);
    color_preview->set_valign(Gtk::ALIGN_CENTER);
    color_preview->set_margin_start(4);
    color_preview->set_margin_end(4);
    color_preview->setRgb(item->color.get_red(), item->color.get_green(),
                          item->color.get_blue());
    color_event_box->add(*color_preview);
    grid->attach(*color_event_box, 3, 0);

    auto color_button = Gtk::manage(new Gtk::Button());
    color_button->set_image(*Gtk::manage(
        new RTImage("color-picker-small", Gtk::ICON_SIZE_BUTTON)));
    color_button->set_relief(Gtk::RELIEF_NONE);
    grid->attach(*color_button, 4, 0);

    const size_t preset_index = item->index;
    color_event_box->signal_button_release_event().connect(
        [this, preset_index, color_preview](GdkEventButton* ev) -> bool {
            this->onAspectRatioPresetPickColor(preset_index, color_preview);
            return true;
        });
    color_button->signal_clicked().connect(sigc::bind(
        sigc::mem_fun(this, &CropGuide::onAspectRatioPresetPickColor),
        item->index,
        color_preview));

    auto remove_button = Gtk::manage(new Gtk::Button());
    remove_button->set_image(*Gtk::manage(
        new RTImage("cancel-small", Gtk::ICON_SIZE_BUTTON)));
    remove_button->set_relief(Gtk::RELIEF_NONE);
    grid->attach(*remove_button, 5, 0);

    remove_button->signal_clicked().connect(sigc::bind(
        sigc::mem_fun(this, &CropGuide::onAspectRatioPresetRemoved),
        item->index));

    grid->show_all();

    return grid;
}

void CropGuide::read(const rtengine::procparams::ProcParams* pp,
                     const ParamsEdited* pedited)
{
    DisableListener disable_listener(this);
    BlockAdjusterEvents block_bleed(m_bleed);
    ConnectionBlocker block_basis(m_basis_changed);

    setEnabled(pp->cropGuide.enabled);
    m_mirror_golden_triangle = pp->cropGuide.mirror_golden_triangle;
    m_rotate_golden_ratio = pp->cropGuide.rotate_golden_ratio;
    m_mirror_golden_ratio = pp->cropGuide.mirror_golden_ratio;
    m_bleed->setValue(pp->cropGuide.bleed);

    if (pedited) {
        set_inconsistent(multiImage && pedited->cropGuide.enabled);
        m_dirty_mirror_golden_triangle = pedited->cropGuide.mirror_golden_triangle;
        m_dirty_rotate_golden_ratio = pedited->cropGuide.rotate_golden_ratio;
        m_dirty_mirror_golden_ratio = pedited->cropGuide.mirror_golden_ratio;
        m_bleed->setEditedState(pedited->cropGuide.bleed ? Edited : UnEdited);
    } else {
        m_dirty_mirror_golden_triangle = false;
        m_dirty_rotate_golden_ratio = false;
        m_dirty_mirror_golden_ratio = false;
    }

    for (size_t i = 0; i < m_presets.size(); i++) {
        const auto& pp_preset = pp->cropGuide.presets[i];
        auto& preset = m_presets[i];
        bool is_enabled = pp_preset.enabled;

        if (preset.visibility_button->get_active() != is_enabled) {
            ConnectionBlocker block(preset.visibility_conn);
            preset.visibility_button->set_active(is_enabled);
            updateImage(preset.visibility_button, is_enabled);
        }

        Gdk::RGBA color;
        color.set_rgba(pp_preset.red, pp_preset.green, pp_preset.blue, pp_preset.alpha);
        preset.color = color;
        preset.preview->setRgb(color.get_red(), color.get_green(), color.get_blue());

        if (pedited) {
            preset.is_dirty = pedited->cropGuide.presets[i];
        } else {
            preset.is_dirty = false;
        }
    }

    m_aspect_ratio_store->remove_all();
    for (const auto& model : m_aspect_ratio_presets) {
        model->active = false;
        model->visible = true;

        Gdk::RGBA color;
        color.set_rgba(1.f, 1.f, 1.f);
        model->color = color;
    }

    for (const auto& param : pp->cropGuide.aspect_ratios) {
        auto& preset = m_aspect_ratio_presets.at(param.preset_index);
        preset->active = true;
        preset->is_portrait = param.is_portrait;
        preset->visible = param.enabled;

        Gdk::RGBA color;
        color.set_rgba(param.red, param.green, param.blue, param.alpha);
        preset->color = color;

        m_aspect_ratio_store->insert_sorted(
            preset, sigc::mem_fun(this, &CropGuide::compareAspectRatioModels));
    }

    refreshAvailableAspectRatios();

    if (pedited) {
        m_dirty_aspect_ratios = pedited->cropGuide.aspect_ratios;
    } else {
        m_dirty_aspect_ratios = false;
    }

    m_basis->set_active(mapBasis(pp->cropGuide.basis));
    if (pedited && !pedited->cropGuide.basis) {
        m_basis->set_active(EMPTY_COMBO_INDEX);
    }
}

void CropGuide::write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited)
{
    pp->cropGuide.enabled = getEnabled();
    pp->cropGuide.mirror_golden_triangle = m_mirror_golden_triangle;
    pp->cropGuide.rotate_golden_ratio = m_rotate_golden_ratio;
    pp->cropGuide.mirror_golden_ratio = m_mirror_golden_ratio;
    pp->cropGuide.bleed = m_bleed->getIntValue();

    if (pedited) {
        pedited->cropGuide.enabled = !get_inconsistent();
        pedited->cropGuide.mirror_golden_triangle = m_dirty_mirror_golden_triangle;
        pedited->cropGuide.rotate_golden_ratio = m_dirty_rotate_golden_ratio;
        pedited->cropGuide.mirror_golden_ratio = m_dirty_mirror_golden_ratio;
        pedited->cropGuide.bleed = m_bleed->getEditedState();
    }

    for (size_t i = 0; i < m_presets.size(); i++) {
        const auto& preset = m_presets[i];
        auto& pp_preset = pp->cropGuide.presets[i];

        pp_preset.enabled = preset.visibility_button->get_active();
        pp_preset.red = preset.color.get_red();
        pp_preset.green = preset.color.get_green();
        pp_preset.blue = preset.color.get_blue();
        pp_preset.alpha = preset.color.get_alpha();

        if (pedited) {
            pedited->cropGuide.presets[i] = preset.is_dirty;
        }
    }

    if (pedited) {
        pedited->cropGuide.aspect_ratios = m_dirty_aspect_ratios;
    }

    pp->cropGuide.aspect_ratios.clear();
    for (const auto& model : m_aspect_ratio_presets) {
        if (!model->active) continue;

        CropGuideParams::AspectRatioParams params(model->index);
        params.enabled = model->visible;
        params.is_portrait = model->is_portrait;
        params.red = model->color.get_red();
        params.green = model->color.get_green();
        params.blue = model->color.get_blue();
        params.alpha = model->color.get_alpha();

        pp->cropGuide.aspect_ratios.push_back(std::move(params));
    }

    const int active_basis = m_basis->get_active_row_number();
    const bool valid_basis = (active_basis > EMPTY_COMBO_INDEX)
                             && (active_basis < INDEX_BASIS_UNCHANGED);
    if (valid_basis) {
        pp->cropGuide.basis = mapBasis(active_basis);
    }
    if (pedited) {
        pedited->cropGuide.basis = valid_basis;
    }
}

void CropGuide::trimValues(rtengine::procparams::ProcParams* pp)
{
    m_bleed->trimValue(pp->cropGuide.bleed);
}

void CropGuide::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    m_basis->append(M("GENERAL_UNCHANGED"));
}

void CropGuide::setAdjusterBehavior(bool bleed)
{
    m_bleed->setAddMode(bleed);
}

void CropGuide::adjusterChanged(Adjuster* adj, double newVal)
{
    if (listener && (getEnabled() || batchMode)) {
        if (adj == m_bleed) {
            listener->panelChanged(EvCropGuideBleedChanged,
                                   Glib::ustring::format(adj->getIntValue()));
        }
    }
}

void CropGuide::enabledChanged()
{
    if (!listener) return;

    if (get_inconsistent()) {
        listener->panelChanged(EvCropGuideEnabled, M("GENERAL_UNCHANGED"));
    } else if (getEnabled()) {
        listener->panelChanged(EvCropGuideEnabled, M("GENERAL_ENABLED"));
    } else {
        listener->panelChanged(EvCropGuideEnabled, M("GENERAL_DISABLED"));
    }
}

bool CropGuide::onPresetToggled(GdkEventButton* event, size_t index)
{
    Preset& toggled_preset = m_presets.at(index);
    const bool is_enabled = toggled_preset.visibility_button->get_active();

    // Disable presets one at a time.
    // Enable presets one at a time if CTRL is held down.
    if (!is_enabled || (event->state & GDK_CONTROL_MASK)) {
        updateImage(toggled_preset.visibility_button, is_enabled);
        toggled_preset.is_dirty = true;
    } else {
        for (size_t i = 0; i < m_presets.size(); i++) {
            Preset& preset = m_presets[i];

            if (i == index) {
                bool is_visible = preset.visibility_button->get_active();
                updateImage(preset.visibility_button, is_visible);
                preset.is_dirty = true;
            } else {
                if (preset.visibility_button->get_active()) {
                    ConnectionBlocker block(preset.visibility_conn);
                    preset.visibility_button->set_active(false);
                    updateImage(preset.visibility_button, false);
                    preset.is_dirty = true;
                }
            }
        }
    }

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged, M(GUIDE_TYPE_OPTIONS.at(index)));
    }

    return true;
}

void CropGuide::onGoldenTriangleMirror()
{
    m_mirror_golden_triangle = !m_mirror_golden_triangle;
    m_dirty_mirror_golden_triangle = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged, M(GUIDE_TYPE_OPTIONS.at(
            CropGuideParams::PresetIndex::GOLDEN_TRIANGLE)));
    }
}

void CropGuide::onGoldenTriangleReset()
{
    m_mirror_golden_triangle = false;
    m_dirty_mirror_golden_triangle = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged, M(GUIDE_TYPE_OPTIONS.at(
            CropGuideParams::PresetIndex::GOLDEN_TRIANGLE)));
    }
}

void CropGuide::onGoldenRatioRotate()
{
    m_rotate_golden_ratio = !m_rotate_golden_ratio;
    m_dirty_rotate_golden_ratio = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged, M(GUIDE_TYPE_OPTIONS.at(
            CropGuideParams::PresetIndex::GOLDEN_RATIO)));
    }
}

void CropGuide::onGoldenRatioMirror()
{
    m_mirror_golden_ratio = !m_mirror_golden_ratio;
    m_dirty_mirror_golden_ratio = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged, M(GUIDE_TYPE_OPTIONS.at(
            CropGuideParams::PresetIndex::GOLDEN_RATIO)));
    }
}

void CropGuide::onGoldenRatioReset()
{
    m_mirror_golden_ratio = false;
    m_dirty_mirror_golden_ratio = true;
    m_rotate_golden_ratio = false;
    m_dirty_rotate_golden_ratio = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged, M(GUIDE_TYPE_OPTIONS.at(
            CropGuideParams::PresetIndex::GOLDEN_RATIO)));
    }
}

void CropGuide::onPresetPickColor(size_t index)
{
    auto& preset = m_presets.at(index);

    Gtk::ColorChooserDialog dialog;
    dialog.set_use_alpha(true);
    dialog.set_rgba(preset.color);

    int result = dialog.run();
    if (result != Gtk::RESPONSE_OK) return;

    Gdk::RGBA color = dialog.get_rgba();
    preset.color = color;
    preset.is_dirty = true;

    preset.preview->setRgb(color.get_red(), color.get_green(), color.get_blue());

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuidePresetChanged,
                               M(GUIDE_TYPE_OPTIONS.at(index)));
    }
}

void CropGuide::onAspectRatioComboChanged()
{
    ConnectionBlocker block(m_available_aspect_ratios_conn);

    int active_row = m_available_aspect_ratios_combo->get_active_row_number();
    Glib::ustring active_text = m_available_aspect_ratios_combo->get_active_text();
    m_available_aspect_ratios_combo->unset_active();
    m_available_aspect_ratios_combo->remove_text(active_row);

    auto is_matching_preset = [&](const Glib::RefPtr<AspectRatioModel>& preset) {
        return active_text == preset->aspect_ratio.label;
    };

    auto it = std::find_if(m_aspect_ratio_presets.begin(), m_aspect_ratio_presets.end(),
                           is_matching_preset);
    if (it == m_aspect_ratio_presets.end()) return;

    Glib::RefPtr<AspectRatioModel> model = *it;
    model->active = true;
    model->visible = true;
    m_dirty_aspect_ratios = true;

    m_aspect_ratio_store->insert_sorted(
        model, sigc::mem_fun(this, &CropGuide::compareAspectRatioModels));

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuideAspectRatioPresetChanged,
                               model->aspect_ratio.label);
    }
}

void CropGuide::onAspectRatioPresetToggled(Gtk::ToggleButton* button, size_t index)
{
    auto& preset = m_aspect_ratio_presets.at(index);

    bool is_visible = button->get_active();
    if (preset->visible == is_visible) return;

    updateImage(button, is_visible);
    preset->visible = is_visible;
    m_dirty_aspect_ratios = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuideAspectRatioPresetChanged,
                               preset->aspect_ratio.label);
    }
}

void CropGuide::onAspectRatioPresetPickColor(size_t index, ColorPreview* preview)
{
    auto& preset = m_aspect_ratio_presets.at(index);

    Gtk::ColorChooserDialog dialog;
    dialog.set_use_alpha(true);
    dialog.set_rgba(preset->color);

    int result = dialog.run();
    if (result != Gtk::RESPONSE_OK) return;

    Gdk::RGBA color = dialog.get_rgba();
    preset->color = color;
    m_dirty_aspect_ratios = true;

    preview->setRgb(color.get_red(), color.get_green(), color.get_blue());

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuideAspectRatioPresetChanged,
                               preset->aspect_ratio.label);
    }
}

void CropGuide::onAspectRatioPresetRotate(size_t index)
{
    auto& preset = m_aspect_ratio_presets.at(index);
    preset->is_portrait = !preset->is_portrait;
    m_dirty_aspect_ratios = true;

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuideAspectRatioPresetChanged,
                               preset->aspect_ratio.label);
    }
}

void CropGuide::onAspectRatioPresetRemoved(size_t index)
{
    auto& preset = m_aspect_ratio_presets.at(index);
    preset->active = false;
    preset->visible = false;
    m_dirty_aspect_ratios = true;

    auto size = m_aspect_ratio_store->get_n_items();
    for (guint i = 0; i < size; i++) {
        if (m_aspect_ratio_store->get_object(i) == preset) {
            m_aspect_ratio_store->remove(i);
            break;
        }
    }

    refreshAvailableAspectRatios();

    if (listener && getEnabled()) {
        listener->panelChanged(EvCropGuideAspectRatioPresetChanged,
                               preset->aspect_ratio.label);
    }
}

void CropGuide::onBasisChanged()
{
    if (listener && (getEnabled() || batchMode)) {
        listener->panelChanged(EvCropGuideBasisChanged, m_basis->get_active_text());
    }
}

void CropGuide::refreshAvailableAspectRatios()
{
    ConnectionBlocker block(m_available_aspect_ratios_conn);
    m_available_aspect_ratios_combo->remove_all();
    for (const auto& model : m_aspect_ratio_presets) {
        if (!model->active) {
            m_available_aspect_ratios_combo->append(model->aspect_ratio.label);
        }
    }
}

int CropGuide::compareAspectRatioModels(
    const Glib::RefPtr<const AspectRatioModel>& lhs,
    const Glib::RefPtr<const AspectRatioModel>& rhs) const
{
    if (lhs->index < rhs->index) {
        return -1;
    } else if (lhs->index > rhs->index) {
        return 1;
    } else {
        return 0;
    }
}
