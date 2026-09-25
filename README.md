# EDI Tools

A simple GUI tool to work with large collections of [Magnetitelluric](https://en.wikipedia.org/wiki/Magnetotellurics) data in the Electrical Data Interchange (EDI) format, which is [SEG standard for MT and EMAP data](https://pubs.usgs.gov/of/2003/of03-056/Data/Edistd.doc).

![](EDITools_demo.gif)

## Features

- Visualize all transfer functions 
- Assign error floors (with a proper error propagation from impedances to derived quantities such as apparent resitivity, phase or phase tensor)
- Decimate data
- Save a PDF graphic of all plots
- Mask any data
- Save/Load a project with all masking and settings preserved
- Export a subset of data to a [GoFEM](https://github.com/GoFEM/pyGoFEM) format for the subsequent inversion.
- Visualize responses calculated by GoFEM or MT inversion (Tellurion3D) as lines over the observed data.

## Dependencies

- Qt >= 5.0
- Eigen >= 3
- CMake >= 3.10
- Boost >= 1.84
- QCustomPlot >= 2.0
- PROJ (UTM coordinate conversion)

## Building

The code was tested for linux using GNU Toolchain. Provided dependencies above are installed (for instance using conda) and CMake knows where to find them, the configuring and building is rather simple

> mkdir build; cd build

> cmake -DQCUSTOMPLOT_PATH:PATH=/path/to/qcustomplot -DEIGEN_PATH:PATH=/path/to/eigen ../

> make

## Usage

### Computed responses

Load the observed EDI data or an existing project, then use **File → Load
responses…**. Select any mix of GoFEM and native MT inversion response files;
the format is detected from their contents. For Tellurion3D, select one or more
`inversion_predicted_iterNNNN.txt` files in the run's `output` folder. Select a
file in **Calculated responses** to overlay its curves on the selected station.
The last selected filename in sorted order is displayed immediately.

Only the input readers differ. Both formats use the same response model, curves,
visibility controls, fit statistics, RMS maps, period maps and project storage.
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

### Data fit analysis

Open **Data → Data fit statistics…** for a separate comparison window. Check or
uncheck any loaded response to include it in the overall nRMS, period, component,
station and normalized-residual plots. Each response keeps the same color across
the comparison plots; the list shows its nRMS and matched/available scalar count.
All residuals, including outliers, contribute to the histograms and statistics.

**Maps and heatmaps** compares two checked responses side by side, with shared
axes and color scales. Stations are ordered north to south, and grey cells mark
missing comparisons. Select either response from the dropdowns. Pan and zoom
the plots, or use **Save PDF…** to export the overview, station plot and selected
spatial comparison.

Use **Plot ranges…** to set minimum and maximum values for each statistics
plot's X and Y axes. Leave **Auto** checked to fit the enabled responses, or
choose **Autoscale all** to reset the ranges. **Colors** changes the
comparison palette. In **Maps and heatmaps**, choose the **Colormap** and
optionally **Reverse** it; uncheck **Auto** beside either nRMS range to edit
its limits. Map and heatmap ranges are independent and shared by their two
panels. Values beyond the limits use the endpoint colors. These display
settings survive response selection and refresh, apply to the PDF export,
and leave the underlying statistics unchanged. The **Station names** checkbox
shows or hides labels on both RMS maps.

The default normalization uses the errors stored with the inversion responses,
as in an inversion data report. Choose **Errors: Observed** to use
the survey's current error floors instead. Masked data and disabled stations
are excluded in both modes. Statistics update after imports, masks, error-floor
changes and decimation; **Refresh** recomputes them on demand. Imported response
values and individual scalar errors are preserved in saved projects. Different
matched counts mean the responses cover different subsets of the observations.

### Component visibility

Each plot's legend has one checkbox per component. Uncheck **XX** and **YY**
to keep only off-diagonal impedances, or leave just one component checked.
Tipper legends offer **Re / Im Tzx / Tzy** in scalar mode and **Real / Imaginary**
in arrow mode. Hidden components remain in the legend so they can be restored.
Click the checkbox or label, or use Tab and Space with the keyboard.

Each checkbox controls that plot's observed points, error bars and computed
curves. Plots have independent visibility settings, saved with the project;
data masks and fit statistics are unchanged. Partial tipper arrows are labeled
with their selected direction, and toggling them preserves that projection.

### Station, induction-vector and phase-tensor maps

Use **Station names** below the small station map to show or hide its labels.
The same option is available under **Plot → Show station names** and is saved
with the project. New sessions start with names hidden.

Click **Maps…** beside that checkbox, or **Data → Induction / phase tensor
maps…**, to open a separate map window. Choose **Observed data** or a loaded
computed response, then select a **Period (s)**; **Previous** and **Next** step
through available periods. Layer checkboxes above the plot independently control
phase tensors, real induction arrows, imaginary induction arrows and station
names. Hover near a station for its actual period and numerical values.

Each station uses its closest available period within **Tolerance (%)**
(default 0.1%). Values are not interpolated. Disabled stations and unavailable,
incomplete or masked components are omitted; valid zero vectors have no arrow.
The summary reports coverage. Response maps use the observed station locations
and enabled stations, with the selected response's available components.

To mask observations from this map, enable **Select stations**, then click a
station, arrow or ellipse, or drag a box around stations. **Ctrl-click** toggles
individual stations; **Ctrl-drag** adds stations. Blue rings mark the selection.
Choose **Tippers**, **Phase tensor**, or **Impedance + tensor**, then click
**Mask** or **Unmask**. This changes all components in that group at the closest
observed period within the tolerance; other periods are unaffected. Tipper masks
apply to both real and imaginary parts. Impedance masks also apply to derived
resistivity and phase. Masked stations retain their center dots for selection
and unmasking. Turn off **Select stations** to pan again.

Masking updates the data plots and fit statistics, and is saved with the project.
It always edits observations, including when viewing a computed response;
computed response values and curves remain unchanged. Disabled stations and
stations without an observed period within the tolerance are skipped.

Induction arrows use **Tzx = north**, **Tzy = east**. **Parkinson (−T)** reverses
the stored components; **Wiese (+T)** keeps their signs. The chosen sign applies
to both real and imaginary arrows. Set their length in kilometres for **|T| = 1**;
the plot includes a **|T| = 0.5** reference arrow.

Phase tensors use **Φ = Re(Z)⁻¹ Im(Z)** (or directly loaded tensor components).
Ellipses have equal major diameters, minor/major ratio **|Φmin / Φmax|**, and
azimuth **α − β**, clockwise from geographic north. Color represents maximum
phase **atan(Φmax)**, minimum phase **atan(Φmin)**, or skew **β**, in degrees.
Choose the diameter, colormap, reversal and color limits above the plot.
Undefined tensors and singular real impedance matrices are omitted.

Maps use WGS84 / UTM with equal horizontal and vertical scales, incorporating
the survey's chosen zone and origin. Glyph directions account for the difference
between true north and grid north. Pan and zoom, use **Fit view** to frame the
visible layers, and **Save PDF…** to export. Map controls survive data refresh
while the window is open.

The tensor geometry follows
[Caldwell et al. (2004)](https://doi.org/10.1111/j.1365-246X.2004.02281.X).

### Survey PDF report

Use **File → Export survey report…** to create a landscape PDF with a survey
overview followed by one page per station. The overview includes a station map,
period coverage, full/partial/masked/missing period counts and, when selected,
response fit statistics. Each station page places resistivity, phase, tipper and
phase-tensor plots beside a map highlighting that station.

Each stored period counts once for each data type (impedance, tipper, phase
tensor). **Full** means all components are finite and enabled: four complex
impedance components, two complex tipper components, or four tensor entries.
**Part.** means at least one scalar is usable, but the set is incomplete or
partly masked. **Mask** means finite values exist but none are enabled (including
disabled stations). **Miss.** means no finite values exist for that type.
These categories do not check uncertainty quality; zero-valued data are valid.
Derived phase tensors count as stored data and follow their own masks.

Every station-table row sums to the station's stored period count. The overview
sums these **station-periods**, rather than multiplying periods by components;
periods absent from a station's layout are not added to its missing count.
The coverage plot shows, separately for each data type, stations with full or
partial data at each exact period. Periods with zero coverage remain visible.
Counts use all components regardless of the plot visibility option.

Choose all stations or enabled stations, an optional loaded response overlay,
error bars, and whether to include only the currently visible components.
GoFEM and native inversion responses use the same report rendering. nRMS uses
observed errors and all matched, enabled scalar data, independently of plot
visibility. Fixed Y ranges from **Plot → Axis ranges** are retained on every
station page, including response overlays and empty panels. Y axes set to
autoscale and period axes scale to each station. Phase wrapping follows the
main window. Maps include all valid survey locations. Missing coordinates are
noted on the station page. The **?** button explains the options.

Export leaves the survey and current plots unchanged. Canceling or a failed
export preserves an existing destination PDF.

### Period layout and resampling

The period-layout, period-map and fit-statistics windows use compact labels.
Hover over a **?** button, click it, or focus it and press **Space** for the
full explanation. Table headers also have tooltips. In the period-layout window,
**Coverage ?** explains that **Usable on grid** includes estimates, while
**Original on grid** counts retained source samples at the exact target periods;
increased grid coverage does not mean additional measurements.

Open **Data → Period layout / Resampling…** to inspect station periods on a
logarithmic axis. Choose an impedance, tipper or phase-tensor component to view
its source samples, masks and missing values. The **Coverage** tab counts enabled
stations at each period. **Station summary** lists coverage and gaps larger than
the chosen bounding-period ratio. Component filters affect the preview only;
creating a survey resamples all components.

Choose a **Logarithmic grid** (limits in seconds and points per decade), a
**Reference station**, or **Custom periods** separated by spaces, commas or
semicolons. Scientific notation is accepted. Custom periods are sorted and
deduplicated. A logarithmic grid includes both requested limits; its final
interval may be shorter to include the upper limit. Click **Preview** to calculate
coverage before creating anything. Changes to data or grid settings invalidate
the preview. **Target preview** distinguishes existing, interpolated, derived,
masked and unavailable values; hover for bounding source periods and skip reasons.

Interpolation rules:

- Interpolate complete complex impedance and tipper components linearly in
  **log(period)**, handling real and imaginary values separately with the same
  weights. Valid zeros are retained.
- Preserve matching source samples, errors and masks (including isolated and
  partially populated samples). Only floating-point roundoff counts as a match.
- Require valid, unmasked endpoints for each complex component. **Never
  extrapolate**, including where a component has a narrower range than its station.
- **Max gap ratio** defaults to **2**. This tests the entire
  bounding interval, not just proximity to its nearest endpoint. For example,
  1–2 s is allowed, while 1–10 s is blocked even for a target very near 1 s.
- **Respect masks** blocks interpolation across masked samples by default. Turning this off allows
  interpolation across masks subject to the same gap limit; exact masked samples
  remain masked. Missing endpoints are skipped, and the resulting larger bracket
  must still pass the gap limit.
- Recompute apparent resistivity, phase and phase tensors from interpolated
  impedance, respecting independent tensor masks and rejecting singular tensors.
  Stations with no impedance use direct component-wise tensor interpolation.
  Resistivity/phase-only inputs are preserved at matching periods; their wrapped
  phase values are not directly interpolated.
- Interpolate the currently effective endpoint standard errors using the same
  positive weights. This avoids an artificial reduction under an independence
  assumption. Invalid endpoint errors block interpolation. Derived errors use
  the existing propagation formulas. Cross-period covariance is not modeled:
  a denser grid does **not** add independent information for inversion.

**Create survey** opens an independent window, with the current loaded
responses available for comparison. Save it to a new `.mtd` project. Every station
gets the same target grid; unsupported values remain missing, and disabled
stations remain disabled. The original stays open. The new project also embeds
its source survey and an audit of source periods and outcomes for each component.
Reopen the tool to **Open source** or **Saved audit…**; **Preview CSV…**
exports the current preview before applying. The saved audit
describes the resampling operation, before subsequent edits or masking.

#### MTPy-v2 review

The implementation was reviewed against MTPy-v2 revision
[`c57f20e`](https://github.com/MTgeophysics/mtpy-v2/tree/c57f20ee8d77fe5c6e979a747eb5f52b50f2bdcc).
Its [transfer-function accessor](https://github.com/MTgeophysics/mtpy-v2/blob/c57f20ee8d77fe5c6e979a747eb5f52b50f2bdcc/mtpy/core/transfer_function/accessor.py)
interpolates real and imaginary parts separately, handles each channel pair's
finite data range and disables extrapolation by default. It removes non-finite
samples before interpolation; its interpolation routine has no bounding-gap
limit or mask-barrier control. The default method is `slinear` against period
itself. Its [survey interpolation](https://github.com/MTgeophysics/mtpy-v2/blob/c57f20ee8d77fe5c6e979a747eb5f52b50f2bdcc/mtpy/core/mt_data.py)
supports a shared target grid and clips it to station ranges when requested.
EDITools adopts the component-level range checks and separate real/imaginary
treatment, with log-period weighting, explicit gap limits, mask barriers,
retention of exact isolated samples, and an auditable source copy.

### Coordinates and export

Use **Coordinates → Convert latitude/longitude to UTM…** to project the loaded
station locations using WGS84. The tool suggests one zone and hemisphere for
the survey; both can be changed. Check **Center coordinates on the survey**
to subtract the midpoint of the stations' projected bounds. The origin easting
and northing are shown in selectable fields, with a **Copy origin** button.
They also remain visible below the map and are recorded in receiver-file comments.
Both receiver exporters always write UTM coordinates, choosing an uncentered
UTM projection automatically if needed. Each receiver file includes a single
`# UTM origin (easting northing, m): ...` line and the zone/hemisphere.
When centering is enabled, receiver coordinates are offsets from that origin;
otherwise the origin is zero and the coordinates are absolute UTM values.

The coordinates apply to the map and data exports, and the zone/origin are
saved with the project. All stations, including disabled ones, define the
center; changing masks does not move the origin. Original latitude/longitude
are retained. Use the **Map coordinates: Lat/Lon / UTM** switch below the map
to change the view without changing export coordinates or the saved UTM origin.
Station selection and map labels follow the selected view. **Coordinates →
Show latitude/longitude on map** is a shortcut to the geographic view.
Elevations and response tensors are unchanged by this
coordinate operation. Projection uses [PROJ](https://proj.org/en/stable/development/quickstart.html).

**File → Export to MT Inversion…** writes the format described in
[MT_DATA_EXPORT_SPEC.md](MT_DATA_EXPORT_SPEC.md). Select the observations and
periods, then export; no coordinate file or convention confirmations are
required. It uses the survey's UTM coordinates, automatically choosing an
uncentered WGS84/UTM projection if none has been applied. The native receiver
file uses `x = northing`, `y = easting`, `z = -elevation` (metres).
Optional sign changes and a distortion-compatibility check are under **Advanced**.

All 24 scalar types are supported, including apparent resistivity and phase.
Phases retain their impedance quadrant, independent of plot wrapping. Missing
or masked data are omitted, valid zeros are preserved, and invalid errors stop
export with a diagnostic. Resistivity/phase data require `Distortion estimate = false`.

**Every data export writes a plain frequency file**: sorted distinct frequencies
in Hz, one number per line, with no header. It contains only frequencies for
which observations were actually exported. Native output `survey.data` has
companions `survey.recvs` and `survey.freqs`; GoFEM output retains its appended
companion naming, `<output>.recvs` and `<output>.freqs`.

Run the format, coordinate and GUI workflow tests after building with
`ctest --test-dir build/Release --output-on-failure`.

The animation above illustrating main features. 

The following hotkeys are also useful when you click on a transfer function plot:

- **Scroll**: zoom in/out both axes.
- **X + Scroll**: zoom in/out only the *Y* axis.
- **D + Click**: pan the plot.

## Contributing

Feel free to create an issue or pull request in case you want to report a problem or contribute/fix the code.
