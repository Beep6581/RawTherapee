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

#pragma once

#error "gettext-based l10n is not implemented yet!"

#include <libintl.h>

#include <utility>

// Note: The macros defined below are similar to those provided by
//       #include <glib/gi18n.h>.

// xgettext version must be greater than 0.15 and must be passed the following
// arguments to extract all translatable strings correctly:
//
// xgettext --keyword=_ --keyword=M --keyword=N_ --keyword=SP_:1,2
//   --keyword=C_:1c,2 --keyword=NC_:1c,2 --keyword=CSP_:1c,2,3
//   --keyword=EV_:3

#define gettext_noop(id) id
#define pgettext_noop(ctx, id) \
    std::pair<const char* const, const char* const>{ctx, id}

/**
 * Marks a string for translation.
 *
 * The string will be replaced by its translation at run time if the
 * translation exists. Otherwise, it will be passed as is.
 *
 * Example: _("Hello")
 */
#define _(id) gettext(id)
// For backwards compatibility with old multilangmgr implementation
// #define M(id) _(id)

/**
 * Marks a string for deferred translation.
 *
 * Unlike the macro `_(id)`, this macro does not replace the string with its
 * translation for special cases of translation strings like when the string
 * is inside a struct or array declaration.
 *
 * It is the programmer's responsibility to call gettext() for the deferred
 * translation string.
 */
#define N_(id) gettext_noop(id)

/**
 * Marks a singular-form and plural-form string for translation.
 *
 * Example: SP_("%d file removed", "%d files removed", num_files)
 */
#define SP_(id_single, id_plural, n) ngettext(id_single, id_plural, n)

/**
 * Marks a context-aware string for translation.
 *
 * Example: C_("Menu|File", "Open")
 */
#define C_(ctx, id) pgettext(ctx, id)

/**
 * Marks a context-aware string for deferred translation.
 *
 * It is the programmer's responsibility to call pgettext() for the deferred
 * translation string.
 */
#define NC_(ctx, id) pgettext_noop(ctx, id)

/**
 * Marks a context-aware singular-form and plural-form string for translation.
 *
 * Example: CSP_("editor", "%d file removed", "%d files removed", num_files)
 */
#define CSP_(ctx, id_single, id_plural, n) \
    cngettext(ctx, id_single, id_plural, n)

// Utility mapping for ProcEventMapper::newEvent() with translation strings
#define EV_(mapper, action, msg) \
    mapper->newEvent(action, msg)

namespace rt {
namespace i18n {

const char* const DOMAIN = "rawtherapee";

// TODO: Setup locales and gettext domain

} // i18n
} // rt
