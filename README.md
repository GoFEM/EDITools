# EDI Tools

A simple GUI tool to work with large collections of [Magnetitelluric](https://en.wikipedia.org/wiki/Magnetotellurics) data in the Electrical Data Interchange (EDI) format, which is [SEG standard for MT and EMAP data](https://pubs.usgs.gov/of/2003/of03-056/Data/Edistd.doc).


## Features

- Visualize all transfer functions
- Assign error floors (with a proper error propagation from impedances to derived quantities such as apparent resistivity, phase or phase tensor)
- Decimate data
- Export a PDF survey report or a PNG image of the current plots
- Mask any data
- Save/Load a project with all masking and settings preserved
- Export selected data in supported inversion formats.
- Overlay computed inversion responses on observed data.
- Load, display and export topography, including optional DEM downloads.

## Dependencies

- A C++14 compiler (C++17 when using GDAL >= 3.9)
- Qt >= 5.12 (Widgets, PrintSupport and Network)
- Eigen >= 3
- CMake >= 3.10
- Boost >= 1.84
- QCustomPlot >= 2.0
- PROJ (UTM coordinate conversion)
- GDAL >= 3.7 (GeoTIFF input, raster reprojection and export)

## Building

The app has been tested on Linux with the GNU toolchain. Install the dependencies
above and configure their locations before building:

```sh
mkdir -p build/Release
cd build/Release
cmake -DQCUSTOMPLOT_PATH:PATH=/path/to/qcustomplot \
  -DEIGEN_PATH:PATH=/path/to/eigen ../..
cmake --build .
cd ../..
```

The build shares the data and plotting libraries between the app and tests.
Run checks with `ctest --test-dir build/Release --output-on-failure`.

## Usage

### Command ribbon

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

### EDI import and nearby periods

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

### Computed responses

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

### Data fit analysis

Open **Analysis → Statistics / RMS** for a separate comparison window. Check or
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

Component colors use two families: **XX / XY** are light / dark blue, and
**YX / YY** are dark / light orange. Tippers retain their blue/light-blue **Tzx**
and orange/yellow **Tzy** palette for real/imaginary parts. The same colors identify observed
points, error bars, response lines, legends and survey PDF plots.

Each plot's legend has one checkbox per component. Uncheck **XX** and **YY**
to keep only off-diagonal impedances, or leave just one component checked.
Tipper legends offer **Re / Im Tzx / Tzy** in scalar mode and **Real / Imaginary**
in arrow mode. Hidden components remain in the legend so they can be restored.
Click the checkbox or label, or use Tab and Space with the keyboard.

Each checkbox controls that plot's observed points, error bars and computed
curves. Plots have independent visibility settings, saved with the project;
data masks and fit statistics are unchanged. Partial tipper arrows are labeled
with their selected direction, and toggling them preserves that projection.

### Linked impedance and phase-tensor masks

**Edit → Link Z / PT** is on by default and saved with the project.
It links mask/unmask operations from plots (including inverted masking),
station component controls, and period maps. Turning it off restores independent
editing; toggling the option does not change existing masks. Tippers are independent.

The link follows the general phase tensor definition `Phi = Re(Z)^-1 Im(Z)`.
Each PT entry depends on all four real impedance entries, so there is no general
one-to-one mapping between Z and PT components. For example, only in the special
2D case do `Zxy` and `Zyx` correspond to `PTyy` and `PTxx`, respectively.
See [Caldwell et al. (2004), equations 13–15](https://doi.org/10.1111/j.1365-246X.2004.02281.X).

Because each complex impedance component shares a real/imaginary mask, masking
one Z component (or its apparent resistivity/phase) also masks **all four PT
entries at that period**. Masking one PT entry also masks **all four Z components
at that period**. This applies only to the opposite tensor: it does not recursively
change the other entries in the starting tensor. Unmasking applies the same links.
For inverted selections, a mask takes precedence if multiple selected components
give conflicting instructions for a shared dependent entry.

### Edit history

**Edit → Undo / Redo** shares one history across station plots, map masking and
rho+ masking. It also covers error floors, station enable/disable, rename/delete,
EDI additions, decimation, UTM settings and topography/elevation edits. A batch
operation is one step; edits that leave the data unchanged add no step.

**Edit → History** lists the latest 20 operations. Select a step to move back
or forward; a new edit after undo discards the redo branch. **Open original
survey** opens an independent copy of the survey as it was loaded in this
session, for comparison. History is kept in memory and starts afresh when a
project is loaded. Project files save the current survey, not the undo history.
Resampling opens a separate survey with its own history and retains its source.

### Neighbouring-station comparison

**Analysis → Neighbours** overlays the stations selected in the
main station list. With one station selected, it chooses that reference and its
nearest enabled neighbours. Change **Reference** and click **Nearest** to choose
a different group; the count specifies neighbours in addition to the reference.
Distances use WGS84 geodesics and handle the date line.

Each station has one visibility checkbox; clicking its marker on the comparison
map toggles the same selection. Choose a component of resistivity/phase, tipper,
or phase tensor. Colours identify stations, **Errors** shows their effective
uncertainties, and phases retain their stored convention. Masked or missing
samples remain gaps. Disabled stations are excluded; stations without locations
can still be checked manually. This comparison does not change data or masks.

### Impedance validity and rho+

Open **Analysis → Rho+ analysis** to compare a station's measured
apparent resistivity and phase with fitted rho+ responses. Select the station;
**XY / YX** shows both off-diagonal components together, with independent fits
drawn as thin dashed lines in their component colors. The **XY** and **YX**
checkboxes control the plots, table rows and masking scope; **Rho+** toggles the
fitted lines. Individual diagonal components are not offered. Invariant responses
remain available in the response selector. The table and CSV identify each
component explicitly. **YX**
is shown as `−Zyx` to use the same sign convention as XY. Phase is not folded
into the first quadrant, so anomalous phases remain visible.

Use **Analysis → Rho+ guides** to overlay thin, dashed XY/YX fits on the main
resistivity and phase plots in their component colors. They follow the existing
component visibility and phase-wrapping controls and refresh after masking or
unmasking, preserving the current axes. The guides cannot be selected for masking
and can be shown alongside loaded responses. They use the main data's propagated
errors. Grid density comes from the validity window; this setting and guide visibility are saved in
the project. Unconverged or grid-sensitive fits are omitted from the main plots.

The native C++ implementation uses the bounded spectral representation of
log admittance from [Parker and Booker (1996)](https://doi.org/10.1016/S0031-9201(96)03191-3).
It minimizes weighted residuals in log resistivity and phase jointly, using
exact spectral-interval integrals, a free real offset, and both spectral tails.
It is a **finite-grid approximation**, with a second solve at twice the selected
**Grid/decade**. The window reports solver convergence and the change in χ²;
a grid-sensitive result should be rechecked at higher density. Curves are drawn
only within the usable period range.

Rho+ uses the **main data's propagated resistivity and phase errors**, including
the error floor assigned in the main interface. It has no separate error-floor
setting and never changes these errors. Fit weights use `σρ / ρ` for log
resistivity and `σφ` for phase. The same errors are used for the displayed error
bars, residuals, guides and masking. Invariant errors are propagated from the
main impedance errors, assuming independent tensor entries. Missing, zero or
invalid errors are counted and excluded from fitting and masking; no fallback
floor is supplied. At least three usable complex-impedance samples with positive
errors spanning distinct periods are required. Masked, missing and disabled
data are omitted. Resistivity/phase-only inputs without complex impedance are
not used. Obsolete rho+ floor settings in older projects are ignored.

The table gives observed and fitted values and normalized residuals. Residuals
beyond ±2 are highlighted. `nRMS = sqrt(χ² / (2N))` uses both measurements at
N periods; it is not reduced χ² or a formal significance test. A poor fit may
reflect multidimensional structure, a convention mismatch, or underestimated
errors. In particular, diagonal components and tensor averages are not generally
1-D responses. This is a consistency diagnostic, not an automatic bad-data test;
fits alone do not change observations or masks. **CSV…** exports the displayed values.

For the selected station's checked **XY / YX** components, the **Preview masks**
checkbox highlights
periods where either normalized log-resistivity or phase residual exceeds the
**Residual threshold (σ)**. The default is 3σ; use 2σ for closer review. This uses
the main data's propagated measurement errors, not a confidence band around
the fitted curve. Only converged, grid-stable fits can provide candidates.
Crosses in component colors and red table rows identify candidates; the preview
reports each component's count and how many paired diagonal and linked PT masks would change,
counting shared PT changes once. Components without a stable fit are skipped
and identified in the preview. **Apply masks** performs one pass,
masking both resistivity and phase for those component/periods. A **Zxy** candidate
also masks **Zxx**, and a **Zyx** candidate also masks **Zyy**, at the same period.
This polarization pairing applies even when **Edit → Link Z / PT** is off;
that option independently controls PT masking. Tippers remain unchanged.
Refitting does not trigger further rejection. **Undo** uses the shared survey
history to restore the latest edit, including paired diagonals and linked PT masks. Preview is enabled by default and automatically refreshes
when navigating stations or changing the threshold, components or fit settings.
Turning it off hides candidates and disables single-station Apply.

**Apply all stations** fits the checked components at every enabled station and
applies one masking pass with the current threshold, main-data errors, grid density
and Z/PT linking setting. Disabled stations and unavailable or unstable fits are
skipped and counted in the summary. All fits finish before any masks change;
**Cancel** leaves the survey unchanged. Each completed batch is one step in
the shared history. Navigating stations preserves undo; other edits add their own steps. Automatic preview never applies additional masking passes.

### Survey impedance averages

Open **Analysis → Survey averages** to display station responses and
their **Regional geometric mean**. Check stations in the list to include them;
the two checkboxes above the plots control curve visibility. Disabled survey
stations and masked or incomplete required components cannot contribute.

The response selector offers these explicitly defined quantities:

| Response | Impedance |
| --- | --- |
| SSQ | `sqrt((Zxx² + Zxy² + Zyx² + Zyy²) / 2)` |
| Berdichevsky | `(Zxy − Zyx) / 2` |
| Determinant | `sqrt(Zxx Zyy − Zxy Zyx)` |

SSQ uses **complex squares**, not squared magnitudes. The square-root quantities
use the branch with nonnegative imaginary part, giving phases from 0° to 180°.
Missing diagonal entries are not assumed to be zero. The determinant-based
regional geometric mean is also called the Berdichevsky regional average;
see [Rung-Arunwan et al. (2017)](https://doi.org/10.1186/s40623-017-0665-8)
for its distinction from regional SSQ averaging.

Each checked station has equal weight. Mean resistivity is
`exp(mean(ln(rho)))`; phases are aligned to a common branch before averaging.
If phase directions cancel, mean phase is left undefined. At an exact common
period these are the apparent resistivity and phase of the complex geometric
mean impedance. **SD lnρ** is the between-station sample spread, not an error
bar on the mean; it is undefined for a single station.

**Tolerance (%)** defaults to zero (exact periods). An optional tolerance groups
nearby periods only when the entire group fits the relative span; it cannot
chain across a larger gap. A station with multiple samples in a group is omitted
from that group. Resistivities retain their measured periods during averaging;
the group is plotted at the lower median contributing period. There is no
interpolation or extrapolation. **Min. stations** suppresses underpopulated means.
The table and CSV include the contributor count at each period, since changing
station coverage can itself change a regional curve. The original survey and
its station enable flags remain unchanged.

### Polarity review

Open **Analysis → Polarity review** (or **Data → Polarity review…** in classic
menus) to flag likely station polarity reversals relative to nearby stations.
The table shows reversed / matched period counts for XY and YX; hover for the
usable-period count. **No flag** means the test did not identify a sustained
reversal, not that the station is necessarily correct.

Comparison uses the nearest five enabled stations within 100 km by default.
Only unmasked samples with impedance SNR ≥ 3, or equivalent phase uncertainty,
are used. Periods must agree within 0.1%; ambiguous duplicates are excluded,
with no interpolation or extrapolation. Each period needs at least three
neighbours and 75% agreement within the angular tolerance (default 30°).
A candidate requires at least three matched periods, coverage of at least 60%
of its usable periods and a near-180° reversal at 80% of matched periods.
Distance, neighbour count and angular tolerance are adjustable.

The plots show raw XY/YX phases without 180° folding, the neighbour consensus
and a correction preview. **Suggest Ex / Ey** selects an electric-row hypothesis
for flagged components; it does not modify data. Select the actual channels and
click **Apply to station** after reviewing the preview and acquisition records.

| Reversed channel | Quantities that change sign |
| --- | --- |
| Ex | Zxx, Zxy |
| Ey | Zyx, Zyy |
| Hx | Zxx, Zyx, Tzx, PTxy, PTyx |
| Hy | Zxy, Zyy, Tzy, PTxy, PTyx |
| Hz | Tzx, Tzy |

Channel choices combine algebraically in the stored x/y basis:
`Z' = De Z Dh`, `T' = sz T Dh`, and `PT' = Dh PT Dh`, with diagonal sign
matrices `De` and `Dh`. Electric reversals therefore leave tippers and the
phase tensor unchanged. Sign-reversed impedance phases shift by 180°;
apparent resistivities, standard errors, error floors and masks are preserved.
The correction applies at every stored period, including masked samples, so
unmasking cannot bring back an uncorrected value. **Undo** restores the edit.
Save the project to retain corrections; reports and data exports use the
corrected observations. Original input files and loaded calculated responses
are unchanged.

This is a relative consistency test, not a diagnosis of wiring faults. It
assumes comparable station orientations and conventions. It cannot identify
a survey-wide common reversal, distinguish electric from magnetic channel
reversal using off-diagonal phases alone, or detect Hz reversals automatically.
Stations with insufficient overlap or no geographic positions remain available
for manual correction. Real structure can also cause phase reversals; see
[the dispersion-relation study](https://doi.org/10.1186/s40623-020-01273-4).
The phase-tensor transformation follows its definition
[PT = Re(Z)⁻¹ Im(Z)](https://doi.org/10.1111/j.1365-246X.2004.02281.X).

### Station, induction-vector and phase-tensor maps

Use **Station names** below the small station map to show or hide its labels.
The same option is available under **Maps → Station names** and is saved
with the project. New sessions start with names hidden.

Click **Maps…** beside that checkbox, or **Maps → PT / induction**, to open a separate map window. Choose **Observed data** or a loaded
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

### Geographic map layers

**Map layers…** adds optional coastlines, country boundaries, state/province
boundaries, and rivers/lakes to station maps, RMS maps and period maps. Each layer
has one checkbox. Report export has the same menu and applies its choices to all
maps in the PDF. Its initial selection follows the main station map.

Natural Earth data are **not bundled in this repository or embedded in the app**.
They are optional external files. Without them, geographic layer checkboxes are
disabled; station plots, maps and report export continue to work normally. No map
data are downloaded automatically.

#### Optional data setup

Download these Natural Earth **1:10-million** ZIP files into a directory outside
the repository, retaining their filenames:

| Layer | Download |
| --- | --- |
| Coastlines | [ne_10m_coastline.zip](https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_coastline.zip) |
| Country boundaries | [ne_10m_admin_0_boundary_lines_land.zip](https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_0_boundary_lines_land.zip) |
| State/province boundaries | [ne_10m_admin_1_states_provinces_lines.zip](https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces_lines.zip) |
| Rivers | [ne_10m_rivers_lake_centerlines.zip](https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_rivers_lake_centerlines.zip) |
| Lakes | [ne_10m_lakes.zip](https://naturalearth.s3.amazonaws.com/10m_physical/ne_10m_lakes.zip) |

Convert them with the provided helper (Python standard library only). The app
automatically looks in the user's data directory:

- Linux: `~/.local/share/EDITools/natural-earth` (or under `XDG_DATA_HOME`).
- Windows: `%LOCALAPPDATA%\EDITools\natural-earth`.
- macOS: `~/Library/Application Support/EDITools/natural-earth`.

For another location, set `EDITOOLS_MAP_DATA_DIR` before starting the app:

```sh
python3 tools/pack_natural_earth.py /path/to/downloads /path/to/map-data
export EDITOOLS_MAP_DATA_DIR=/path/to/map-data
./build/Release/main
```

On Windows, set `$env:EDITOOLS_MAP_DATA_DIR = 'D:\Maps\natural-earth'` in
PowerShell before launching the executable. Restart the app after installing or
changing data. The helper produces five compressed `.bin` files (about 8 MB) plus
`manifest.json` recording source URLs, versions and checksums. Keep these files
outside the repository. It retains line parts and lake shoreline rings without
additional simplification, converts longitude/latitude to float32, and uses
Qt-compatible compression.

Once installed, layers work offline and export as vector lines. They use the
same latitude/longitude or WGS84/UTM projection, origin and units as the survey
map. Layers do not alter station selection, masking or automatic axis extents.
Zoom out when the survey is far from a regional boundary or coastline. UTM
context is limited to ±30° of the zone's central meridian and latitudes −80° to 84°.

Natural Earth is **public domain**, permitting use, modification and redistribution,
including commercial use; attribution is optional. See the
[Natural Earth terms](https://www.naturalearthdata.com/about/terms-of-use/).
These are generalized regional reference features; small rivers and lakes may
be absent, and political boundaries follow Natural Earth's de facto representation.

#### OpenStreetMap online basemap

Enable **Map layers… → OpenStreetMap (online)** on a station, RMS or period map.
It starts off. **OSM opacity** adjusts the background in all maps and exports.
Tiles are aligned with the current geographic or UTM coordinates, including the
survey origin. Station symbols, arrows and ellipses remain on top.

Only the visible map in the active window requests tiles, after panning or zooming
settles. Hidden windows, background windows and exports do not download tiles.
The app shares its cache between maps, identifies itself to the server, respects
HTTP cache lifetimes, and revalidates expired tiles with server validators. If
expiry information is absent, tiles are cached for seven days. At most two
requests run together; rate limits and access-denied responses stop further
requests. There is no area-download, prefetch or offline-download feature.

**PDF and image exports use available tiles only**, including in survey reports.
View the desired area online before exporting. Higher-resolution export does not
fetch extra tiles; incomplete exports carry a short “Partial basemap” note.
On screen, a compact **© OpenStreetMap contributors** link opens the licence;
hover over it for tile availability and connection details. Exports include the
printed licence URL. **OSM attribution / licence…** also opens the licence information.

The default service is `https://tile.openstreetmap.org/{z}/{x}/{y}.png`, used under
the [tile policy](https://operations.osmfoundation.org/policies/tiles/).
Map credits follow the [attribution guidelines](https://osmfoundation.org/wiki/Licence/Attribution_Guidelines).
Tile requests disclose the viewed map area to the service; station names and
measurements are not transmitted. Availability depends on the service and the
network. Natural Earth layers remain independent of the online basemap.

The cache is stored outside the repository, under the user's cache directory
(`~/.cache/EDITools/osm` on Linux, respecting `XDG_CACHE_HOME`). It targets 512 MB
by removing old, expired entries while preserving fresh tiles and at least seven
days of storage. To use another compatible OSM-derived service, set
`EDITOOLS_OSM_TILE_URL` to its HTTPS `{z}/{x}/{y}` template before launching.
Additional provider credits can be supplied in `EDITOOLS_OSM_ATTRIBUTION`;
use the provider's required attribution and access terms. Restart after changes.

### Survey PDF report

Use **Report → Survey PDF** to create a landscape PDF with a survey
overview, station pages, phase-tensor ellipse maps, and/or induction-vector maps.
Each section can be selected independently. The overview includes a station map
with station names, period coverage, full/partial/masked/missing period counts and, when selected,
response fit statistics. Each station page places resistivity, phase, tipper and
phase-tensor plots beside a map highlighting that station. Station tippers can
be exported as lines or arrows; masked samples interrupt observation lines.

Local project paths are omitted from the report. Response labels use only the
filename, without directories or machine names.

Map pages use a comma-separated list of periods in seconds (**All** selects all
survey periods), with one page per period and dataset. Selected phase-tensor
ellipses and induction vectors appear together on the same map. Choose
observations, the selected response, or both; choose real and/or imaginary
induction vectors. The dialog previews the total page count. Nearest-period
matching, tolerance, vector convention, sizes, colors and station labels follow
the period map window when open; maps never interpolate or extrapolate. Only
valid, enabled data appear. Map extents fit each exported page. The UTM zone and
**origin easting/northing in metres** accompany the overview, station and period
maps so plotted offsets can be converted back to absolute coordinates.

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
visibility. Fixed Y ranges from **Report → Axis ranges** are retained on every
station page, including response overlays and empty panels. Y axes set to
autoscale and period axes scale to each station. Phase wrapping follows the
main window. Station locator maps include all valid survey locations. Missing coordinates are
noted on the station page. Hover over an option to read its tooltip.

Export leaves the survey and current plots unchanged. Canceling or a failed
export preserves an existing destination PDF.

### Period layout and resampling

The period-layout, period-map and fit-statistics windows use compact labels.
Hover over controls, their labels, or table headers for explanations. In the
period-layout window, the **Coverage** tab and plot explain that **Usable on
grid** includes estimates, while
**Original on grid** counts retained source samples at the exact target periods;
increased grid coverage does not mean additional measurements.

Open **Edit → Layout / resample** to inspect station periods on a
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

### Coordinates and export

#### Topography

Open **Maps → Topography**, or **Topography…** beside the station map.
Load one or more single-band **GeoTIFFs**, or choose **Copernicus GLO-30** or
**USGS 3DEP** and click **Download**. Downloads occur only on request. The default
area covers the survey plus a 5 km margin; change the margin and click
**Survey extent**, or edit the absolute UTM bounds directly. **Cell size** sets
the output spacing, not the source resolution. Grids are limited to 16 million
cells. Large areas need coarser spacing.

- [Copernicus GLO-30 Public](https://registry.opendata.aws/copernicus-dem/)
  provides approximately 30 m surface elevations, including vegetation and
  buildings. The public cloud collection uses the 2021 release; some tiles are
  unavailable. The app downloads at most 16 one-degree tiles per request and
  asks you to accept the linked dataset licence. Heights use EGM2008. Required
  attribution and licence notices accompany numeric exports and PDF maps.
- [USGS 3DEP](https://www.usgs.gov/3d-elevation-program/about-3dep-products-services)
  provides US bare-earth elevations from sources of varying resolution. The
  app requests numeric elevation GeoTIFFs from the National Map service, with
  a default output spacing of 10 m. Service requests are limited to 8000 rows
  or columns. Products are free without use restrictions. The service combines
  source datasets: verify the applicable vertical datum before changing station
  heights or using elevations in a model.
- Local files must have a coordinate reference system and geotransform. Confirm
  **Input units** and **Vertical datum**; file units can also be read automatically.
  Unknown units require an explicit choice. The app converts heights to metres
  and reprojects horizontally into the survey's WGS84/UTM zone. It performs no
  vertical datum conversion. Missing coverage remains missing. Overlapping files
  use the first valid raster value in the selected input order.

**Show**, **Colours**, **Hillshade** and **Contours** control the background on the
station, period and fit-statistics maps. Choose a palette, elevation range,
opacity and contour interval. **Auto** restores the loaded grid's elevation
range. Hover over controls for details. Survey PDF export has a separate
**Topography** checkbox and uses these display settings on its maps, including
combined phase-tensor/induction maps. Rendering and PDF export never download
terrain data.

**GeoTIFF…** exports a clipped, projected elevation raster with absolute UTM
coordinates, survey-origin metadata, positive-up heights in metres, NoData, source attribution and
vertical-datum metadata. **XYZ…** exports valid cell centres as
`northing_offset easting_offset height`, subtracting the survey's UTM origin.
Choose positive-up elevation or positive-down depth. Its header records the
column order, units, UTM EPSG code, origin, vertical datum and source. Export
uses the current bounds and cell size; it preserves missing coverage. Failed
or canceled exports preserve existing destination files.

Loading a DEM does not change station elevations. **Station heights…** previews
the current elevation, sampled DEM height and difference for every station.
Uncheck any stations to preserve, confirm the vertical datum, then **Apply**.
Stations outside usable coverage stay unchanged. Response tensors and masks
are unaffected.

The clipped grid and its display settings are embedded in `.mtd` projects, so
reopening and exporting work offline without the original files. Downloads are
cached in the application's standard local cache directory under `dem`; no DEM
tiles are bundled with the repository. Automatic exports contain dataset
metadata, not local source/cache paths. The editable **Credit** field is public
text included in exports.

Use **Maps → UTM / origin** to project the loaded
station locations using WGS84. The default zone and hemisphere come from the
midpoint of the geographic survey extent, using the shortest longitude arc
across the date line. Existing choices are retained when reopening the dialog.
The preview overlays nominal zone boundaries, including Norway and Svalbard,
and identifies stations outside the selected zone or hemisphere. Choose a zone
with the controls, a comparison-table row, or a double-click on the map.

The comparison table shows maximum absolute local scale error at stations and
across a sampled 33 × 33 grid covering the survey extent. The extent result is
an estimate, not a guaranteed bound. The adjustable **Reference** highlights
larger errors without preventing a valid choice; 0.1% equals 1 metre per
kilometre. The station table also shows signed scale error and grid convergence.
These describe projection effects on the reference ellipsoid; terrain-height
effects are not included. A common projection keeps coordinates continuous
across zone boundaries, but wide surveys can still have substantial distortion.

Check **Center UTM** in the dialog or beside the main map coordinate selector
to subtract the midpoint of the stations' projected bounds. This changes the
shared survey origin, without changing the zone or its distortion. The origin easting
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
Station selection and map labels follow the selected view. **Maps → Geographic** is a shortcut to the geographic view.
Elevations and response tensors are unchanged by this
coordinate operation. Projection uses [PROJ](https://proj.org/en/stable/development/quickstart.html).

**Survey → Native MT** exports selected observations and periods
to native data, receiver and frequency files. It uses the survey's UTM coordinates,
automatically choosing an uncentered WGS84/UTM projection if none has been applied. The native receiver
file uses `x = northing`, `y = easting`, `z = -elevation` (metres).
Optional sign changes and a distortion-compatibility check are under **Advanced**.

All 24 scalar types are supported, including apparent resistivity and phase.
Phases retain their impedance quadrant, independent of plot wrapping. Missing
or masked data are omitted, valid zeros are preserved, and invalid errors stop
export with a diagnostic. Resistivity/phase data require `Distortion estimate = false`.

Both export dialogs show a live validation preview before writing. Tabs list
exactly the stations, periods and scalar components that will be exported,
including effective uncertainty ranges. The summary counts omitted data and
shows units, sign choices, projection and origin. Invalid selections, coordinates
or errors prevent export; hover for the full coordinate description. Previewing
and exporting do not change the survey's saved coordinates.

Station names may start with digits. Names containing spaces or unsupported
characters can be exported using **Use compatible names** in either export
dialog. This optional conversion prefixes names such as `001` with `S_`, replaces
unsupported characters with underscores, and adds suffixes to prevent collisions.
The Stations tab shows original and exported names. Project names stay unchanged;
a `.names.csv` companion records changed names. Use this mapping to restore
original receiver names in calculated responses before loading them into the
project for comparison.

Both formats use the same data preparation and validation. GoFEM keeps complex
real/imaginary pairs together and uses positive-up receiver elevation; native
output uses positive-down elevation. Valid zero tippers are retained, and tipper
errors retain measured uncertainty when it exceeds the configured floor.
Period selection is exact, close-frequency collisions are diagnosed, and receiver
companions contain only stations with exported observations.

**Every data export writes a plain frequency file**: sorted distinct frequencies
in Hz, one number per line, with no header. It contains only frequencies for
which observations were actually exported. Native output `survey.data` has
companions `survey.recvs` and `survey.freqs`; GoFEM output retains its appended
companion naming, `<output>.recvs` and `<output>.freqs`.

Run the format, coordinate and GUI workflow tests after building with
`ctest --test-dir build/Release --output-on-failure`.

The following hotkeys are also useful when you click on a transfer function plot:

- **Scroll**: zoom in/out both axes.
- **X + Scroll**: zoom in/out only the *Y* axis.
- **D + Click**: pan the plot.

## Contributing

Feel free to create an issue or pull request in case you want to report a problem or contribute/fix the code.
