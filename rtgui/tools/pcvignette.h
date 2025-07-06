/*
 *  This file is part of RawTherapee.
 */
#pragma once

#include <gtkmm.h>

#include "adjuster.h"
#include "toolpanel.h"

class PCVignette final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel
{

protected:
    Adjuster* strength;
    Adjuster* feather;
    Adjuster* roundness;

public:
    static const Glib::ustring TOOL_NAME;

    PCVignette ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) override;
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) override;
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr) override;
    void setBatchMode   (bool batchMode) override;

    void adjusterChanged (Adjuster* a, double newval) override;
    void enabledChanged  () override;
    void setAdjusterBehavior (bool strengthadd, bool featheradd, bool roundnessadd);
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
};
