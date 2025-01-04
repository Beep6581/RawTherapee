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

#include <iomanip>

#include "resize.h"

#include "eventmapper.h"
#include "guiutils.h"

#include "rtengine/procparams.h"

using namespace rtengine;
using namespace rtengine::procparams;

const Glib::ustring Resize::TOOL_NAME = "resize";

Resize::Resize () : FoldableToolPanel(this, TOOL_NAME, M("TP_RESIZE_LABEL"), false, true), maxw(100000), maxh(100000)
{
    auto m = ProcEventMapper::getInstance();
    EvResizeAllowUpscaling = m->newEvent(RESIZE, "HISTORY_MSG_RESIZE_ALLOWUPSCALING");
    EvResizeLongedge = m->newEvent (RESIZE, "HISTORY_MSG_RESIZE_LONGEDGE");
    EvResizeShortedge = m->newEvent (RESIZE, "HISTORY_MSG_RESIZE_SHORTEDGE");   

    cropw = 0;
    croph = 0;

    Gtk::Grid* combos = Gtk::manage (new Gtk::Grid());
    combos->set_row_spacing(4);

    appliesTo = Gtk::manage (new MyComboBoxText ());
    appliesTo->append (M("TP_RESIZE_CROPPEDAREA"));
    appliesTo->append (M("TP_RESIZE_FULLIMAGE"));
    appliesTo->set_active (0);
    appliesTo->set_hexpand();
    appliesTo->set_halign(Gtk::ALIGN_FILL);

    Gtk::Label *label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_APPLIESTO"), Gtk::ALIGN_START));
    
    combos->attach(*label, 0, 0, 1, 1);
    combos->attach(*appliesTo, 1, 0, 1, 1);

    // See Resize::methodChanged() when adding a new method.
    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_RESIZE_LANCZOS"));
    method->append (M("TP_RESIZE_NEAREST"));
    method->set_active (0);
    method->set_hexpand();
    method->set_halign(Gtk::ALIGN_FILL);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_METHOD"), Gtk::ALIGN_START));
    
    combos->attach(*label, 0, 1, 1, 1);
    combos->attach(*method, 1, 1, 1, 1);

    spec = Gtk::manage (new MyComboBoxText ());
    spec->append (M("TP_RESIZE_SCALE"));
    spec->append (M("TP_RESIZE_WIDTH"));
    spec->append (M("TP_RESIZE_HEIGHT"));
    spec->append (M("TP_RESIZE_FITBOX"));
    spec->append (M("TP_RESIZE_LONG"));
    spec->append (M("TP_RESIZE_SHORT"));
    spec->set_active (0);
    spec->set_hexpand();
    spec->set_halign(Gtk::ALIGN_FILL);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_SPECIFY"), Gtk::ALIGN_START));

    combos->attach(*label, 0, 2, 1, 1);
    combos->attach(*spec, 1, 2, 1, 1);

    pack_start (*combos, Gtk::PACK_SHRINK, 4);

    scale = new Adjuster (M("TP_RESIZE_SCALE"), 0.01, MAX_SCALE, 0.01, 1.);
    scale->setAdjusterListener (this);

    pack_start (*scale, Gtk::PACK_SHRINK, 4);

    sizeBox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));

    Gtk::Box* sbox = Gtk::manage (new Gtk::Box ());
    Gtk::Box* wbox = Gtk::manage (new Gtk::Box ());
    Gtk::Box* hbox = Gtk::manage (new Gtk::Box ());
    Gtk::Box* ebox = Gtk::manage (new Gtk::Box ());
    Gtk::Box* lebox = Gtk::manage (new Gtk::Box ());
    Gtk::Box* sebox = Gtk::manage (new Gtk::Box ());
    
    w = Gtk::manage (new MySpinButton ());
    w->set_width_chars(5);
    setExpandAlignProperties(w, false, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    h = Gtk::manage (new MySpinButton ());
    h->set_width_chars(5);
    setExpandAlignProperties(h, false, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    le = Gtk::manage (new MySpinButton ());
    le->set_width_chars(5);
    setExpandAlignProperties(le, false, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);
    se = Gtk::manage (new MySpinButton ());
    se->set_width_chars(5);
    setExpandAlignProperties(se, false, false, Gtk::ALIGN_END, Gtk::ALIGN_CENTER);

    wbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_W"))), Gtk::PACK_SHRINK, 0);
    wbox->pack_start (*w);
    hbox->set_spacing(3);
    hbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_H"))), Gtk::PACK_SHRINK, 0);
    hbox->pack_start (*h);
    lebox->set_spacing(3);
    lebox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_LE"))), Gtk::PACK_SHRINK, 0);
    lebox->pack_start (*le);
    sebox->set_spacing(3);
    sebox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_SE"))), Gtk::PACK_SHRINK, 0);
    sebox->pack_start (*se);

    sbox->set_spacing(4);
    sbox->pack_start (*wbox);
    sbox->pack_start (*hbox);
    sbox->set_homogeneous();
    ebox->set_spacing(4);
    ebox->pack_start (*lebox);
    ebox->pack_start (*sebox);
    ebox->set_homogeneous();
    
    sizeBox->pack_start (*sbox, Gtk::PACK_SHRINK, 0);
    sizeBox->pack_start (*ebox, Gtk::PACK_SHRINK, 0);
    sizeBox->show_all ();
    sizeBox->reference ();

    allowUpscaling = Gtk::manage(new Gtk::CheckButton(M("TP_RESIZE_ALLOW_UPSCALING")));
    pack_start(*allowUpscaling);
    allowUpscaling->signal_toggled().connect(sigc::mem_fun(*this, &Resize::allowUpscalingChanged));

    w->set_digits (0);
    w->set_increments (1, 100);
    w->set_range (MIN_SIZE, MAX_SCALE * maxw);
    w->set_value (800);           // Doesn't seem to have any effect (overwritten in Resize::read)

    h->set_digits (0);
    h->set_increments (1, 100);
    h->set_range (MIN_SIZE, MAX_SCALE * maxh);
    h->set_value (600);           // Doesn't seem to have any effect (overwritten in Resize::read)

    le->set_digits (0);
    le->set_increments (1, 100);
    le->set_range (MIN_SIZE, MAX_SCALE * maxw);
    le->set_value (900);

    se->set_digits (0);
    se->set_increments (1, 100);
    se->set_range (MIN_SIZE, MAX_SCALE * maxh);
    se->set_value (900);

    wconn = w->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryWChanged), true);
    hconn = h->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryHChanged), true);
    leconn = le->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryLEChanged), true);
    seconn = se->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entrySEChanged), true);
    aconn = appliesTo->signal_changed().connect ( sigc::mem_fun(*this, &Resize::appliesToChanged) );
    method->signal_changed().connect ( sigc::mem_fun(*this, &Resize::methodChanged) );
    sconn = spec->signal_changed().connect ( sigc::mem_fun(*this, &Resize::specChanged) );

    getSubToolsContainer()->hide();

    show_all();
}

Resize::~Resize ()
{
    idle_register.destroy();
    delete scale;
    delete sizeBox;
}

void Resize::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    disableListener ();
    aconn.block (true);
    wconn.block (true);
    hconn.block (true);
    leconn.block (true);
    seconn.block (true);
    sconn.block (true);
    scale->block(true);

    scale->setValue (pp->resize.scale);
    w->set_value (pp->resize.width);
    h->set_value (pp->resize.height);
    le->set_value (pp->resize.longedge);
    se->set_value (pp->resize.shortedge);
    setEnabled (pp->resize.enabled);
    spec->set_active (pp->resize.dataspec);
    allowUpscaling->set_active(pp->resize.allowUpscaling);
    setDimensions();    // Sets Width/Height in the GUI according to value of Specify after loading a .pp3 profile (same behavior as if changed manually)
    updateGUI();

    appliesTo->set_active (0);

    if (pp->resize.appliesTo == "Cropped area") {
        appliesTo->set_active (0);
    } else if (pp->resize.appliesTo == "Full image") {
        appliesTo->set_active (1);
    }

    if (pp->resize.method == "Lanczos") {
        method->set_active (0);
    } else if (pp->resize.method == "Nearest") {
        method->set_active (1);
    } else {
        method->set_active (0);
    }

    wDirty = false;
    hDirty = false;
    leDirty = false;
    seDirty = false;

    if (pedited) {
        wDirty = pedited->resize.width;
        hDirty = pedited->resize.height;
        leDirty = pedited->resize.longedge;
        seDirty = pedited->resize.shortedge;
        scale->setEditedState (pedited->resize.scale ? Edited : UnEdited);

        if (!pedited->resize.appliesTo) {
            appliesTo->set_active (2);
        }

        if (!pedited->resize.method) {
            method->set_active (3);
        }

        if (!pedited->resize.dataspec) {
            spec->set_active (6);
        }

        allowUpscaling->set_inconsistent(!pedited->resize.allowUpscaling);
        set_inconsistent (multiImage && !pedited->resize.enabled);
    }

    setDimensions(); // fixes the issue that values in GUI are not recomputed when loading profile
    
    scale->block(false);
    sconn.block (false);
    wconn.block (false);
    hconn.block (false);
    leconn.block (false);
    seconn.block (false);
    aconn.block (false);
    enableListener ();
}

void Resize::write (ProcParams* pp, ParamsEdited* pedited)
{
    int dataSpec = spec->get_active_row_number();
    
    pp->resize.scale  = scale->getValue();

    pp->resize.appliesTo = "Cropped area";

    if (appliesTo->get_active_row_number() == 0) {
        pp->resize.appliesTo = "Cropped area";
    } else if (appliesTo->get_active_row_number() == 1) {
        pp->resize.appliesTo = "Full image";
    }

    pp->resize.method = "Lanczos";

    if (method->get_active_row_number() == 0) {
        pp->resize.method = "Lanczos";
    } else if (method->get_active_row_number() == 1) {
        pp->resize.method = "Nearest";
    }

    pp->resize.dataspec = dataSpec;
    pp->resize.width = w->get_value_as_int ();
    pp->resize.height = h->get_value_as_int ();
    pp->resize.longedge = le->get_value_as_int ();
    pp->resize.shortedge = se->get_value_as_int ();
    pp->resize.enabled = getEnabled ();
    //printf("  L:%d   H:%d\n", pp->resize.width, pp->resize.height);

    pp->resize.allowUpscaling = allowUpscaling->get_active();

    if (pedited) {
        pedited->resize.enabled   = !get_inconsistent();
        pedited->resize.dataspec  = dataSpec != 6;
        pedited->resize.appliesTo = appliesTo->get_active_row_number() != 2;
        pedited->resize.method    = method->get_active_row_number() != 3;

        if (pedited->resize.dataspec) {
            pedited->resize.scale     = scale->getEditedState ();
            pedited->resize.width     = wDirty;
            pedited->resize.height    = hDirty;
            pedited->resize.longedge  = leDirty;
            pedited->resize.shortedge = seDirty;
        } else {
            pedited->resize.scale     = false;
            pedited->resize.width     = false;
            pedited->resize.height    = false;
            pedited->resize.longedge  = false;
            pedited->resize.shortedge = false;
        }
        pedited->resize.allowUpscaling = !allowUpscaling->get_inconsistent();
    }
}

void Resize::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    scale->setDefault (defParams->resize.scale);

    if (pedited) {
        scale->setDefaultEditedState (pedited->resize.scale ? Edited : UnEdited);
    } else {
        scale->setDefaultEditedState (Irrelevant);
    }
}

void Resize::adjusterChanged(Adjuster* a, double newval)
{
    if (!batchMode) {
        wconn.block (true);
        hconn.block (true);
        h->set_value ((croph && appliesTo->get_active_row_number() == 0 ? croph : maxh) * a->getValue ());
        w->set_value ((cropw && appliesTo->get_active_row_number() == 0 ? cropw : maxw) * a->getValue ());
        wconn.block (false);
        hconn.block (false);
    }

    if (listener && (getEnabled () || batchMode)) {
        listener->panelChanged (EvResizeScale, Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(2), scale->getValue()));
    }
}

int Resize::getComputedWidth(double height)
{

    if (cropw && appliesTo->get_active_row_number() == 0)
        // we use the crop dimensions
    {
        return (int)((double)(cropw) * (height / (double)(croph)) + 0.5);
    } else
        // we use the image dimensions
    {
        return (int)((double)(maxw) * (height / (double)(maxh)) + 0.5);
    }
}

int Resize::getComputedHeight(double width)
{

    if (croph && appliesTo->get_active_row_number() == 0)
        // we use the crop dimensions
    {
        return (int)((double)(croph) * (width / (double)(cropw)) + 0.5);
    } else
        // we use the image dimensions
    {
        return (int)((double)(maxh) * (width / (double)(maxw)) + 0.5);
    }
}

void Resize::appliesToChanged ()
{

    //printf("\nPASSAGE EN MODE \"%s\"\n\n", appliesTo->get_active_text().c_str());
    setDimensions();

    if (listener && (getEnabled () || batchMode)) {
        //printf("Appel du listener\n");
        listener->panelChanged (EvResizeAppliesTo, appliesTo->get_active_text());
    }
}

void Resize::methodChanged ()
{

    if (listener && (getEnabled () || batchMode)) {
        listener->panelChanged (EvResizeMethod, method->get_active_text());
    }

    // Post-resize Sharpening assumes the image is in Lab space, and currently Lanczos is the only method which uses that space, and Lanczos is on row 0.
    if (method->get_active_row_number() == 0) {
        getSubToolsContainer()->set_sensitive(true);
    } else {
        getSubToolsContainer()->set_sensitive(false);
    }
}

void Resize::update (bool isCropped, int cw, int ch, int ow, int oh)
{

    // updating crop values now
    if (isCropped) {
        cropw = cw;
        croph = ch;
    } else {
        cropw = 0;
        croph = 0;
    }

    // updating the full image dimensions
    if (ow && oh) {
        maxw = ow;
        maxh = oh;
    }

    // updating the GUI synchronously
    setDimensions();
}

void Resize::sizeChanged(int mw, int mh, int ow, int oh)
{
    // updating max values now
    maxw = ow;
    maxh = oh;

    // updating the GUI synchronously
    setDimensions();
}

void Resize::setDimensions ()
{
    idle_register.add(
        [this]() -> bool
        {
            wconn.block(true);
            hconn.block(true);
            leconn.block(true);
            seconn.block(true);
            scale->block(true);
            
            int refw, refh;

            if (appliesTo->get_active_row_number() == 0 && cropw) {
                // Applies to Cropped area
                refw = cropw;
                refh = croph;
            } else {
                // Applies to Full image or crop is disabled
                refw = maxw;
                refh = maxh;
            }

            w->set_range(32, MAX_SCALE * refw);
            h->set_range(32, MAX_SCALE * refh);

            switch (spec->get_active_row_number()) {
                case 0: {
                    // Scale mode
                    w->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * scale->getValue() + 0.5)));
                    h->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * scale->getValue() + 0.5)));
                    break;
                }

                case 1: {
                    // Width mode
                    const double tmp_scale = w->get_value() / static_cast<double>(refw);
                    scale->setValue(tmp_scale);
                    h->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * tmp_scale + 0.5)));
                    break;
                }

                case 2: {
                    // Height mode
                    const double tmp_scale = h->get_value() / static_cast<double>(refh);
                    scale->setValue(tmp_scale);
                    w->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * tmp_scale + 0.5)));
                    break;
                }

                case 3: {
                    // Bounding box mode
                    const double tmp_scale =
                        w->get_value() / h->get_value() < static_cast<double>(refw) / static_cast<double>(refh)
                            ? w->get_value() / static_cast<double>(refw)
                            : h->get_value() / static_cast<double>(refh);

                    scale->setValue(tmp_scale);
                    break;
                }

                case 4: {
                    // Long edge mode
                    if (refw > refh) {
                        const double tmp_scale = le->get_value() / static_cast<double>(refw);
                        scale->setValue(tmp_scale);
                        se->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * tmp_scale + 0.5)));
                    } else {
                        const double tmp_scale = le->get_value() / static_cast<double>(refh);
                        scale->setValue(tmp_scale);
                        se->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * tmp_scale + 0.5)));
                    }
                    break;
                }

                case 5: {
                    // Short edge mode
                    if (refw > refh) {
                        const double tmp_scale = se->get_value() / static_cast<double>(refh);
                        scale->setValue(tmp_scale);
                        le->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * tmp_scale + 0.5)));
                    } else {
                        const double tmp_scale = se->get_value() / static_cast<double>(refw);
                        scale->setValue(tmp_scale);
                        le->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * tmp_scale + 0.5)));
                    }
                    break;
                }

                default: {
                    break;
                }
            }

            scale->block(false);
            wconn.block(false);
            hconn.block(false);
            leconn.block(false);
            seconn.block(false);

            return false;
        }
    );
}

void Resize::fitBoxScale()
{
    double tmpScale;
    double neww = w->get_value ();
    double newh = h->get_value ();

    if (cropw && appliesTo->get_active_row_number() == 0) {
        // we use the crop dimensions
        if (((double)(cropw) / (double)(croph)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(cropw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(croph);
        }
    } else {
        // we use the image dimensions
        if (((double)(maxw) / (double)(maxh)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(maxw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(maxh);
        }
    }

    scale->setValue (tmpScale);
}

void Resize::entryWChanged ()
{

    wDirty = true;

    // updating width
    if (!batchMode) {
        if (spec->get_active_row_number() == 3) {
            // Fit box mode
            fitBoxScale();
        } else {
            // Other modes
            hconn.block (true);
            scale->block (true);

            h->set_value ((double)(getComputedHeight(w->get_value())));
            scale->setValue (w->get_value () / (cropw && appliesTo->get_active_row_number() == 0 ? (double)cropw : (double)maxw));

            scale->block (false);
            hconn.block (false);
        }
    }

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            notifyBBox();
        } else {
            if (getEnabled () || batchMode) {
                listener->panelChanged (EvResizeWidth, Glib::ustring::format (w->get_value_as_int()));
            }
        }
    }
}

void Resize::entryHChanged ()
{

    hDirty = true;

    if (!batchMode && listener) {
        if (spec->get_active_row_number() == 3) {
            // Fit box mode
            fitBoxScale();
        } else {
            // Other modes
            wconn.block (true);
            scale->block (true);

            w->set_value ((double)(getComputedWidth(h->get_value())));
            scale->setValue (h->get_value () / (croph && appliesTo->get_active_row_number() == 0 ? (double)croph : (double)maxh));

            scale->block (false);
            wconn.block (false);
        }
    }

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            notifyBBox();
        } else {
            if (getEnabled () || batchMode) {
                listener->panelChanged (EvResizeHeight, Glib::ustring::format (h->get_value_as_int()));
            }
        }
    }
}

void Resize::entryLEChanged ()
{
    
    leDirty = true;

    // updating long edge
    if (!batchMode && listener) {
        int refw, refh;
        
        seconn.block (true);
        scale->block (true);
        
        if (cropw && appliesTo->get_active_row_number() == 0) {
            // we use the crop dimensions
            refw = cropw;
            refh = croph;
        } else {
            // we use the image dimensions
            refw = maxw;
            refh = maxh;
        } 

        if (refw > refh) {
            se->set_value ((double) (getComputedHeight(le->get_value())));
            scale->setValue (le->get_value () / (cropw && appliesTo->get_active_row_number() == 0 ? (double)cropw : (double)maxw));
        } else {
            se->set_value ((double)(getComputedWidth(le->get_value())));
            scale->setValue (le->get_value () / (croph && appliesTo->get_active_row_number() == 0 ? (double)croph : (double)maxh));
        }
         
        scale->block (false);
        seconn.block (false);
    }

    if (listener) {
        if (getEnabled () || batchMode) {
           listener->panelChanged (EvResizeLongedge, Glib::ustring::format (le->get_value_as_int()));
       }
    }
}

void Resize::entrySEChanged ()
{

    seDirty = true;

    // updating short edge
    if (!batchMode && listener) {
        int refw, refh;
        
        leconn.block (true);
        scale->block (true);
        
        if (cropw && appliesTo->get_active_row_number() == 0) {
            // we use the crop dimensions
            refw = cropw;
            refh = croph;
        } else {
            // we use the image dimensions
            refw = maxw;
            refh = maxh;
        } 

        if (refw > refh) {
            le->set_value ((double)(getComputedWidth(se->get_value())));
            scale->setValue (se->get_value () / (croph && appliesTo->get_active_row_number() == 0 ? (double)croph : (double)maxh));
        } else {
            le->set_value ((double)(getComputedHeight(se->get_value())));
            scale->setValue (se->get_value () / (cropw && appliesTo->get_active_row_number() == 0 ? (double)cropw : (double)maxw));
        }

        scale->block (false);
        leconn.block (false);
    }

    if (listener) {
        if (getEnabled () || batchMode) {
           listener->panelChanged (EvResizeShortedge, Glib::ustring::format (se->get_value_as_int()));
       }
    }
}

void Resize::specChanged ()
{

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        scale->sliderChanged();
        break;

    case (1):
        // Width mode
        w->set_value((double)(getComputedWidth(h->get_value())));
        entryWChanged();
        break;

    case (2):
        // Height mode
        h->set_value((double)(getComputedHeight(w->get_value())));
        entryHChanged();
        break;

    case (3):
        // Bounding box mode
        notifyBBox();
        break;

    case (4):
        // Long edge mode
        entryLEChanged();
        break;

    case (5):
        // Short edge mode
        entrySEChanged();
        break;

    default:
        break;
    }

    updateGUI();
}

void Resize::updateGUI ()
{

    removeIfThere (this, scale, false);
    removeIfThere (this, sizeBox, false);

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        pack_start (*scale, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        break;

    case (1):
        // Width mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        w->set_sensitive (true);
        h->set_sensitive (false);
        w->get_parent()->get_parent()->show();
        le->get_parent()->get_parent()->hide();
        break;

    case (2):
        // Height mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        w->set_sensitive (false);
        h->set_sensitive (true);
        w->get_parent()->get_parent()->show();
        le->get_parent()->get_parent()->hide();
        break;

    case (3):
        // Bounding box mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        w->set_sensitive (true);
        h->set_sensitive (true);
        w->get_parent()->get_parent()->show();
        le->get_parent()->get_parent()->hide();
        break;

    case (4):
        // Long edge mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        le->set_sensitive (true);
        se->set_sensitive (false);
        w->get_parent()->get_parent()->hide();
        le->get_parent()->get_parent()->show();
        break;

    case (5):
        // Short edge mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        le->set_sensitive (false);
        se->set_sensitive (true);
        w->get_parent()->get_parent()->hide();
        le->get_parent()->get_parent()->show();
        break;

    default:
        break;
    }
}

void Resize::notifyBBox()
{
    if (listener && (getEnabled () || batchMode)) {
        listener->panelChanged (EvResizeBoundingBox, Glib::ustring::compose("(%1x%2)", (int)w->get_value(), (int)h->get_value() ));
    }
}

void Resize::setBatchMode (bool batchMode)
{

    method->append (M("GENERAL_UNCHANGED"));
    spec->append (M("GENERAL_UNCHANGED"));
    ToolPanel::setBatchMode (batchMode);
    scale->showEditedCB ();
}

void Resize::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Resize::allowUpscalingChanged()
{

    if (listener) {
        if (allowUpscaling->get_inconsistent()) {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_UNCHANGED"));
        } else if (allowUpscaling->get_active()) {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_DISABLED"));
        }
    }
}


void Resize::setAdjusterBehavior (bool scaleadd)
{

    scale->setAddMode(scaleadd);
}

void Resize::trimValues (rtengine::procparams::ProcParams* pp)
{

    scale->trimValue(pp->resize.scale);
}
