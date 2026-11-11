# Maps and geographic layers

[Guide contents](README.md) · [← Response analysis](analysis.md) · [Profile analysis →](profiles.md)

- [Station, induction-vector and phase-tensor maps](#station-induction-vector-and-phase-tensor-maps)
- [Geographic map layers](#geographic-map-layers)
  - [Data setup for source builds](#data-setup-for-source-builds)
  - [OpenStreetMap online basemap](#openstreetmap-online-basemap)

## Station, induction-vector and phase-tensor maps

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
**Stroke** adjusts line and arrowhead thickness; the default is a thin 0.7-pixel
line for dense surveys.

Phase tensors use **Φ = Re(Z)⁻¹ Im(Z)** (or directly loaded tensor components).
Ellipses have equal major diameters, minor/major ratio **|Φmin / Φmax|**, and
azimuth **α − β**, clockwise from geographic north. Color represents maximum
phase **atan(Φmax)**, minimum phase **atan(Φmin)**, or skew **β**, in degrees.
Choose the diameter, colormap, reversal and color limits above the plot.
Undefined tensors and singular real impedance matrices are omitted.

**Color by** also offers apparent resistivity and phase for **XY**, **YX**,
**Berdichevsky** and **SSQ**. These appear as colored station markers; phase
tensors can remain visible as outlines and induction arrows can be overlaid.
Resistivity uses a logarithmic color scale in Ω m; phase uses a linear scale in
degrees. Limits and colormaps are adjustable and apply to PDF exports too.
YX uses −Zyx to align its phase convention with XY. The local Berdichevsky
impedance is `(Zxy − Zyx)/2`; SSQ is
`sqrt((Zxx² + Zxy² + Zyx² + Zyy²)/2)`, taking the complex square root with
nonnegative imaginary part. These use the same calculations as the impedance
analysis window. Required components must be finite and unmasked: SSQ requires
all four components, while the off-diagonal average requires only XY and YX.
Hover over a station to read its value. No spatial averaging or interpolation
is performed.

Maps use WGS84 / UTM with equal horizontal and vertical scales, incorporating
the survey's chosen zone and origin. Glyph directions account for the difference
between true north and grid north. Pan and zoom, use **Fit view** to frame the
station locations, and **Save PDF…** to export. Map controls survive data refresh
while the window is open.

The tensor geometry follows
[Caldwell et al. (2004)](https://doi.org/10.1111/j.1365-246X.2004.02281.X).

## Geographic map layers

**Map layers…** adds optional coastlines, country boundaries, state/province
boundaries, and rivers/lakes to station maps, RMS maps and period maps. Each layer
has one checkbox. Report export has the same menu and applies its choices to all
maps in the PDF. Its initial selection follows the main station map.
Period maps inherit the main map's layer and station-name choices when first
opened; their **Map layers…** menu then controls backgrounds independently.

Portable packages include Natural Earth coastlines, country and state/province
boundaries, rivers and lakes in `share/natural-earth` (about 8 MB). They work
offline without setup. Their public-domain notice is included in `licenses.zip`.
The geometry is kept outside the executable and source repository. No map data
are downloaded automatically. Without local data, the corresponding checkboxes
are disabled; station plots, maps and report export continue to work normally.

### Data setup for source builds

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
first honors `EDITOOLS_MAP_DATA_DIR`, then checks `share/natural-earth` beside
the executable (or one directory above `bin`), then the user's data directory:

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
context covers the forward hemisphere of the selected UTM projection and
latitudes −80° to 84°. This includes surveys spanning several UTM zones;
the UTM dialog reports projection distortion before a zone is chosen.

Natural Earth is **public domain**, permitting use, modification and redistribution,
including commercial use; attribution is optional. See the
[Natural Earth terms](https://www.naturalearthdata.com/about/terms-of-use/).
These are generalized regional reference features; small rivers and lakes may
be absent, and political boundaries follow Natural Earth's de facto representation.

### OpenStreetMap online basemap

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

---

[Guide contents](README.md) · [← Response analysis](analysis.md) · [Profile analysis →](profiles.md)
