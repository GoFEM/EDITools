# Coordinates, topography and data export

[Guide contents](README.md) · [← Period layout and resampling](periods.md)

- [Topography](#topography)
- [UTM and survey origin](#utm-and-survey-origin)
- [Data export](#data-export)

## Topography

Open **TF tools → Topography**, or **Topography…** beside the station map.
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

## UTM and survey origin

Use **TF tools → UTM / origin** to project the loaded
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
Station selection and map labels follow the selected view. **TF tools → Geographic** is a shortcut to the geographic view.
Elevations and response tensors are unchanged by this
coordinate operation. Projection uses [PROJ](https://proj.org/en/stable/development/quickstart.html).

## Data export

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

---

[Guide contents](README.md) · [← Period layout and resampling](periods.md)
