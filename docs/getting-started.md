# Getting started

[Guide contents](README.md) · [Editing and masking →](editing.md)

- [Workspace tabs](#workspace-tabs)
- [Command ribbon](#command-ribbon)
- [EDI import and nearby periods](#edi-import-and-nearby-periods)
- [Computed responses](#computed-responses)
- [Transfer-function view](#transfer-function-view)
- [Plot hotkeys](#plot-hotkeys)

## Workspace tabs

**TF view** is the default, pinned first tab. Open a profile, period map,
statistics, neighbouring-station comparison, Rho+ analysis, survey average,
polarity review, period layout or detailed auto-mask review from **Views**
or **Auto mask**. The analysis opens beside TF view in the same project window.
Its plots and tables occupy the workspace; its controls appear in the ribbon.

A view command returns to an existing tab. **Views → New view** opens another
instance; **Duplicate view** copies the current analysis and its settings.
Double-click a tab or use **Rename tab** to change its name; **F2** also works
when the tab bar has keyboard focus. Drag tabs to reorder them. TF view stays
first and cannot be closed. Close an analysis with its cross or **Ctrl+W**;
**Ctrl+Shift+T** reopens the most recently closed tab, including its settings.

Each analysis has its own data selection, display settings and view extent.
Survey data, masks and Undo/Redo are shared by the project. A profile's corridor
selection and an average's checked stations do not change the survey's enabled
stations. Selecting a station in a profile updates the current TF station while
keeping the profile tab active. **Undo vertex** edits the line; the main **Undo**
command restores survey edits. Inactive analyses refresh when selected.

Saving a project stores the open tabs, their names, order, settings and active
tab. Profile vertices are stored geographically, so changing the UTM origin
preserves the line. Version 10 projects include this workspace state; older
projects open in TF view. Calculated previews are recomputed from the loaded
data. Closed-tab recovery lasts for the current session.

**TF tools → Station navigator** collapses the station/map sidebar. Analyses
use their own full-width workspace. **Open as survey** and **Create survey**
produce an independent project window with its own data and edit history.

## Command ribbon

**Project**, **Edit**, **Views** and **Export** stay available. Additional pages
follow the active workspace: **TF tools / TF appearance / Auto mask / Mask
settings**, **Profile tools / Profile appearance / Strike tools**, **Map tools /
Map appearance**, and the corresponding analysis tools. Open project, Save,
Undo and Redo remain in the upper left. Component visibility controls remain
beside their plots.

Controls appear when their settings apply. For example, rotation reveals strike
options; Manual reveals the angle; automatic strike or PT filtering reveals
quality thresholds. Arrow settings appear when an arrow layer is enabled, fixed
limits appear when Auto is off, and resampling shows the fields for its chosen
grid. Hidden settings keep their values. Commands that apply but lack data or a
valid preview are disabled. The ribbon height stays fixed; wide groups scroll
horizontally on smaller screens.

Use the arrow or **Ctrl+F1** to collapse the commands; click a ribbon tab to
expand them. **Interface → Classic menus** restores the menu bar. Existing
keyboard shortcuts work in both layouts. **Plot → Configure TF view** opens
**TF appearance**. Opening an analysis exposes its contextual ribbon controls.
Hover over a command for its description.

**Export → Active view PNG** captures the active view. **Active view PDF / CSV**
uses that view's export where available. **Survey PDF**, **Native MT** and
**GoFEM** operate on the survey and have their own selection controls.

Each station-list entry shows its name, number of periods and coordinates.
Coordinates follow the map's **Lat/Lon** or **UTM** selection, including
**Center UTM**. Hover over a station for its full name, coordinates, UTM zone
and origin. Period counts include masked periods and refresh after decimation
or undo/redo. Station checkboxes still enable or disable the station.

## EDI import and nearby periods

Both impedance (`MTSECT`) and seven-channel remote-reference spectra
(`SPECTRASECT`) are processed natively in C++. No Fortran compiler or runtime is
required. Spectral matrices use the channel order Hx, Hy, Hz, Ex, Ey, Rx, Ry.
The inline `// 49` count may be omitted; the reader validates the matrix size
against `NCHAN` and rejects missing or extra values.
The estimator calculates impedance, tipper and their standard errors from the
packed spectral powers and `AVGT`; derived quantities use the same routines as
other input formats. Invalid reference powers, nonpositive averaging counts,
singular systems and truncated matrices are rejected. Missing output channels
remain missing without discarding valid outputs from other channels.

After selecting EDI files, the import preview shows station/sample counts, the
period range, the distinct-period count before and after merging, and a table
of all periods involving new stations, including exact matches already shared
across stations. The preview states when no merging is needed; click **Import**
to keep those periods unchanged. The period selector offers **Keep original**
(default), **Merge close periods**, and **Interpolate to existing**.
For merging, the initial **Tolerance (%)** is **0.1**. Changing it updates the preview.

A group must satisfy `(maximum period - minimum period) / minimum period <
tolerance / 100`; adjacent matches cannot form a chain spanning a wider range.
New stations receive the group's lower median distinct period, or an existing
survey period when there is exactly one. Existing stations remain unchanged.
Groups containing multiple samples from one station or multiple distinct
existing periods are left unchanged and identified in the table.

Merging assigns common frequencies without averaging or interpolation. It
preserves impedance, tipper, phase tensor, their errors and masks. Apparent
resistivity and its error are adjusted for the assigned frequency. Original EDI
files remain unchanged. Duplicate station names are skipped and listed in the
preview help. Canceling the preview or a file-read failure leaves the current
survey and project association unchanged.

When adding stations to an existing project, **Interpolate to existing** puts
only the new stations on the union of valid periods already stored at existing
sites, including disabled sites and masked periods. Existing station values,
errors, masks and period layouts remain unchanged. This option is unavailable
without an existing layout or when no usable target values can be produced.

Interpolation uses the same engine as **Edit → Layout / resample**: real and
imaginary transfer functions are interpolated against log(period), with no
extrapolation and a default **Max gap ratio** of **2**. **Respect masks** is on
by default. Errors use weighted endpoint standard errors; derived quantities
are recalculated, and exact samples retain their values, errors and masks.
The preview lists retained, interpolated, derived and unavailable values per
target period; hover over unavailable counts for the reasons. Counts represent
four complex Z components, two complex tippers and four real PT entries per
new station, rather than independent scalar observations.

Incoming samples outside the chosen layout are omitted from the imported
curves. The project retains the complete pre-interpolation survey and an audit
of the new stations. **Open source** and **Saved audit…** in the Period Layout
window expose these after saving/reopening the project. The existing survey's
earlier resampling history remains in that source snapshot. Canceling or
switching back to **Keep original** discards the interpolation preview;
**Undo** can reverse the accepted import.

## Computed responses

Load the observed EDI data or an existing project, then use **Survey → Responses**. Select any mix of GoFEM and native MT inversion response files;
the format is detected from their contents. Select a file in **Calculated
responses** to overlay its curves on the selected station.
The last selected filename in sorted order is displayed immediately.

Both formats support the same plots, visibility controls, fit statistics, RMS
maps, period maps and project storage.
Reloading a file replaces its data and marks open analysis tabs for refresh.

| Input format | Six columns |
| --- | --- |
| Native MT inversion | `frequency_hz receiver observable component value error` |
| GoFEM | `type frequency_hz source receiver value error` |

All 24 supported impedance, tipper, resistivity, phase and phase-tensor scalar
types are accepted. Rows may be unordered and sparse. Missing values leave
gaps, valid zeros are preserved, and absent stations show no computed curves.
Impedances produce resistivity, phase and phase-tensor curves where the necessary
components exist; explicitly supplied derived values take precedence. Station
names must match the observed survey for overlays and maps. In the main TF view,
tipper responses appear as dashed arrows and phase tensors as unfilled ellipse
outlines. Missing vector components or incomplete tensors are omitted; projects
with saved tipper projections retain those projections.

The original scalar values and separate real/imaginary errors are retained for
fit statistics and project saving. GoFEM permits zero response errors for plotting;
those rows can be compared using **Errors: Observed**, but are excluded
from **Errors: Response** normalization. Blank lines and `#` comments
are accepted; GoFEM also accepts trailing `#` comments. Malformed or duplicate
scalar rows are reported with a line number.

Reload GoFEM response files imported by older versions once to recover their
original missing-value and error information, then save the project.

## Transfer-function view

The main plots always use the transfer-function (TF) layout: XY/YX on the left,
XX/YY on the right, then full-width induction arrows and phase-tensor ellipses.
This also applies when opening older projects.
The current station name appears above the plots and in exported plot images.

Configure the display directly in **TF appearance**:

- **Axis ranges**: choose the shared period axis or a resistivity/phase column.
  Enter minimum and maximum, then press **Apply** or Enter to fix the range.
  Scientific notation is accepted; period and resistivity limits must be positive.
  **Auto** fits that axis when changing stations; **Fit all TF axes** restores
  automatic ranges throughout the view. The two impedance columns can have
  different Y ranges.
- **Colors**: choose a component, then select a swatch or enter a color name or
  `#RRGGBB` value. Points, errors, responses, legends and rho+ guides use the
  corresponding colors. Masked Z/PT and tipper colors are configurable too.
- **Sizes (pixels)**: adjust marker diameter, line width, tipper arrow scale
  (pixels per unit magnitude) and phase-tensor ellipse diameter.
- **Phase tensor · φmin**: choose a colormap, reverse it, and set its phase limits.
- **Display**: toggle error bars and phase wrapping, fit the axes, or reset
  colors, sizes and colormap with **Reset TF appearance**.

View settings are saved with the project. Appearance changes preserve the
current zoom, observations and masks. Invalid limits or colors are marked in
the ribbon and leave the applied setting unchanged.

Component checkboxes control observations, errors, calculated responses and
rho+ guides. The period axes pan and zoom together. Diagonal phases remain signed
(−180° to 180° with automatic ranges); **Wrap phase** applies to YX.
Masking or unmasking XY also edits XX at the same periods; YX similarly edits YY, even if
the paired diagonal is hidden. Selecting a diagonal leaves other impedance
components unchanged. **Link Z / PT** independently controls the corresponding PT masks.

Drag a rectangle around observations in either impedance column, then use the
plot's context menu to **Mask**, **Unmask**, or apply an **Inverted Mask**. Masked
points remain selectable. Hidden observations and calculated responses cannot
be selected. For inverted masks, a mask takes precedence when selections give
conflicting instructions for a paired diagonal or linked tensor entry.
**Undo / Redo** restores the whole edit, including its paired and linked masks.

Ellipses have equal major diameters, with geographic north up and east right.
Their color shows the minimum principal phase, with a default 0–90° scale (values
outside the chosen interval use the end colors). Masked ellipses use the
**Masked Z / PT** color (grey by default); unfilled
outlines show the selected calculated response. Missing or undefined tensors
are omitted. Select ellipse centers to mask/unmask all four PT entries at those
periods. **Link Z / PT** controls the corresponding impedance masks; tipper masks
remain independent. Use **Edit → Group masking** to inspect or mask individual
PT entries. **Export → Active view PNG** exports the displayed layout; the survey PDF
retains its own plot settings.

## Plot hotkeys

The following hotkeys are also useful when you click on a transfer function plot:

- **Scroll**: zoom in/out both axes.
- **X + Scroll**: zoom in/out only the *Y* axis.
- **D + Click**: pan the plot.

---

[Guide contents](README.md) · [Editing and masking →](editing.md)
