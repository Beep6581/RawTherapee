# Image Processing Tools

This file describes the basics of adding/modifying an image processing tool
to RawTherapee. Examples of image processing tools are the cropping tool,
and vignetting tool.

This file does not give fully functional example code. Instead, it provides
an overview of the steps and important functions. You'll need to grep and use
example code from existing tools.

## Processing Parameters

All the parameters that control how an image is processed are stored in a
massive `ProcParams` struct. These parameters are what get written/read to
RawTherapee's .pp3 files. The struct is found in `rtengine/procparams.h`.

To add a new tool, you should create a struct holding values for your tool.
Afterwards, update `ProcParams` to support your new struct.

```cpp
// rtengine/procparams.h

struct FakeToolParams {
    double value;

    // Assign default values in constructor
    FakeToolParams();

    bool operator==(const FakeToolParams& other) const;
    bool operator!=(const FakeToolParams& other) const { return !(*this == other); }
};

struct ProcParams {
    // ...
    FakeToolParams fake_tool;
    // ...

    // Update setDefaults() and operator==() to include your new tool.
    // Update load() and save() to serialize/deserialize to the pp3 file.
};
```

### Edited Parameters

To support batch processing and partial profiles, we need to track which
parameters are "active" and should have edits applied. This is accomplished
using the `ParamsEdited` struct. It simply holds booleans for each parameter
based on if it needs to be updated or not.

The struct is found in `rtgui/paramsedited.h`. As with `ProcParams`, create a
struct for your tool's parameters and update `ParamsEdited`'s member variables
and functions to support your tool.

```cpp
struct FakeToolParamsEdited {
    bool slider;
};

struct ParamsEdited {
    // ...
    FakeToolParamsEdited fake_tool;
    // ...
};
```

## GUI

RawTherapee uses the GTK 3 UI framework through gtkmm (C++ wrapper library).

The GUI code for each processing tool is available under `rtgui/tools/`.
Existing code should be used as a reference for creating new tool UIs. A
basic boilerplate header example is given below which should help guide you in
identifying important functions to implement.

If you need to have slider controls (a.k.a. adjusters), you'll need to inherit
from `AdjusterListener` and add the `setAdjusterBehavior()` function.

```cpp
#include "adjuster.h"
#include "guiutils.h"
#include "toolpanel.h"

#include "rtengine/procevents.h"

class FakeTool final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel {
public:
    static const Glib::ustring TOOL_NAME;

    FakeTool();

    // FoldableToolPanel
    void read(const rtengine::procparams::ProcParams* pp,
              const ParamsEdited* pedited = nullptr) override;
    void write(rtengine::procparams::ProcParams* pp,
               ParamsEdited* pedited = nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams* defaults,
                     const ParamsEdited* pedited = nullptr) override;
    void trimValues(rtengine::procparams::ProcParams* pp) override;
    void setBatchMode(bool batch_mode) override;
    void enabledChanged() override;

    void setAdjusterBehavior(bool bool_for_each_adjuster);
    // AdjusterListener
    void adjusterChanged(Adjuster* adj, double new_val) override;

    // ...
};
```

Once the boilerplate for your tool is set up, you can add it to the editor
by updating the code for `ToolPanelCoordinator`. You need to update:

- Update `ToolPanelCoordinator`
  - Add pointer member variable for tool
  - Add enum entry to `ToolPanelCoordinator::Tool`
  - Add case to `getToolName()`
  - Add case to `getFoldableToolPanel()`
  - Place your tool in the desired tool panel group (`std::vector<ToolTree>`)
    - e.g. `const std::vector<ToolTree> COLOR_PANEL_TOOLS = { /* ... */ };`
- Add case to `toollocationpref.cc:getToolTitleKey()`

### Localization

Currently, an in-house localization method is used. Text to be localized are
stored in key-text pairs. These can be found in `rtdata/languages/default`.
Translators will take the values in `default` and provide translations for them
in their respective language files.

In the code, the localized text can be found using the `M()` function. For
example, to obtain the text corresponding to key `TP_OK`, call `M("TP_OK")`.
When creating your GUI, make sure to add all strings to the `default`
localization file and use `M()` to load them.

### Trigering Preview/History Updates

In order for your tool to trigger image preview updates and be displayed in
the edit history, you need to connect changes to the event system. Example
code to do so is below.

Events are categorized such that they trigger only part of the image processing
pipeline when emitted. See `rtengine/refreshmap.h` for more details.

```cpp
// rtgui/tools/faketool.h

class FakeTool /* ... */ {
    // ...

private:
    void setupEvents();

    rtengine::ProcEvent EvFakeToolEnabled;
};
```

```cpp
// rtgui/tools/faketool.cc

// This function should be called in the constructor.
void FakeTool::setupEvents()
{
    auto m = ProcEventMapper::getInstance();
    // ALL is one possible stage when regenerating image preview.
    // See refreshmap.h for more options.
    //
    // HISTORY_MSG_FAKE_TOOL_ENABLED is a localization string key.
    EvFakeToolEnabled = m->newEvent(ALL, "HISTORY_MSG_FAKE_TOOL_ENABLED");
}
```

When registering a new event, the localization key is for the general action
that is displayed in the history. Later, when emitting the event, you can
include additional details in the `listener->panelChanged()` call. The string
should be added to `rtdata/languages/default`.

## Maintaining Feature Compatibility

### Batch Edit Mode Support

Your `read()` and `write()` functions in the tool class should be implemented
while accounting for the `ParamsEdited` values.

#### Adjusters

If you have adjusters in your UI, you'll need to allow choosing between setting
vs adding to existing value as the adjuster behavior. To do so, add an entry
in `addsetids.h` for each new adjuster. Then, call your tool's
`setAdjusterBehavior()` from `BatchToolPanelCoordinator::initSession()`.

You'll also need to zero out the parameter later in `initSession()` based on
adjuster behavior preferences using code like the following:

```cpp
if (options.baBehav[ADDSET_FAKE_TOOL_PARAM_1]) { pparams.fake_tool.param = 0; }
```

In `ParamsEdited::combine()`, make sure to consider the add/set option  when
combining parameters (grep for `ADDSET` in `paramsedited.cc` for examples).

The user should be able to change the setting from preferences. To do so, edit
`Preferences::getBatchProcPanel()` to call `appendBehavList()` for your new
tool's adjusters.

### Partial Profile Copy/Paste

Like for batch mode support, you'll need to implement `read()` and `write()`
with proper `ParamsEdited` handling. Afterwards, you need to update
`rtgui/partialpastedlg.*` to have an entry for your tool's parameters.
