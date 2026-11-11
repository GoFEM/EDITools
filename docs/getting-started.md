# Getting started

[Guide contents](README.md) · [Editing and masking →](editing.md)

- [Command ribbon](#command-ribbon)
- [EDI import and nearby periods](#edi-import-and-nearby-periods)
- [Computed responses](#computed-responses)
- [Plot hotkeys](#plot-hotkeys)

## Command ribbon

The compact ribbon groups commands into **Survey**, **Edit**, **Analysis**,
**Maps** and **Report**. Open project, Save, Undo and Redo stay available in the
upper left. Component visibility controls remain beside their plots.

Use the arrow or **Ctrl+F1** to collapse the commands; click a tab to expand them.
**Interface → Classic menus** restores the menu bar. The layout, collapse state
and selected tab are remembered between sessions. Existing keyboard shortcuts
work in both layouts. Hover over a command for its description; tools that need
survey data are disabled until a survey is loaded.

The instructions below use ribbon labels. The classic menus expose the same
commands and share their checked and enabled states.

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
Reloading a file replaces its data and updates every open analysis window.

| Input format | Six columns |
| --- | --- |
| Native MT inversion | `frequency_hz receiver observable component value error` |
| GoFEM | `type frequency_hz source receiver value error` |

All 24 supported impedance, tipper, resistivity, phase and phase-tensor scalar
types are accepted. Rows may be unordered and sparse. Missing values leave
gaps, valid zeros are preserved, and absent stations show no computed curves.
Impedances produce resistivity, phase and phase-tensor curves where the necessary
components exist; explicitly supplied derived values take precedence. Station
names must match the observed survey for overlays and maps. Tippers follow the
selected scalar/arrow display mode.

The original scalar values and separate real/imaginary errors are retained for
fit statistics and project saving. GoFEM permits zero response errors for plotting;
those rows can be compared using **Errors: Observed**, but are excluded
from **Errors: Response** normalization. Blank lines and `#` comments
are accepted; GoFEM also accepts trailing `#` comments. Malformed or duplicate
scalar rows are reported with a line number.

Reload GoFEM response files imported by older versions once to recover their
original missing-value and error information, then save the project.

## Plot hotkeys

The following hotkeys are also useful when you click on a transfer function plot:

- **Scroll**: zoom in/out both axes.
- **X + Scroll**: zoom in/out only the *Y* axis.
- **D + Click**: pan the plot.

---

[Guide contents](README.md) · [Editing and masking →](editing.md)
