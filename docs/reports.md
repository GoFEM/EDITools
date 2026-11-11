# Survey reports

[Guide contents](README.md) · [← Profile analysis](profiles.md) · [Period layout and resampling →](periods.md)

- [Current plot image](#current-plot-image)
- [Survey PDF report](#survey-pdf-report)

## Current plot image

Use **Report → Plot image** to save a PNG of the selected station's TF view:
the station heading, two impedance columns, induction arrows and phase-tensor ellipses. The
image includes the current zoom, ribbon colors and sizes, component visibility, masks, error bars and
calculated response overlays. This exports the displayed layout without
changing the survey.

## Survey PDF report

Use **Report → Survey PDF** to create a landscape PDF with a survey
overview, station pages, phase-tensor ellipse maps, and/or induction-vector maps.
Each section can be selected independently. The overview includes a station map
with station names, period coverage, full/partial/masked/missing period counts and, when selected,
response fit statistics. Each station page places resistivity, phase, tipper and
phase-tensor plots beside a map highlighting that station. Station tippers can
be exported as lines or arrows; masked samples interrupt observation lines.
PDF station pages retain their own four-plot layout and scalar PT components;
use **Plot image** when you want the main TF layout.

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
visibility. Fixed XY/YX resistivity and phase ranges from **View → Axis ranges**
are retained on every station page, including response overlays and empty panels.
The PDF uses one axis per quantity, so separate XX/YY ranges and the TF view's
colors, glyph sizes and ellipse colormap apply only to the on-screen view and
its PNG export. Y axes set to
autoscale and period axes scale to each station. Phase wrapping follows the
main window. Station locator maps include all valid survey locations. Missing coordinates are
noted on the station page. Hover over an option to read its tooltip.

Export leaves the survey and current plots unchanged. Canceling or a failed
export preserves an existing destination PDF.

---

[Guide contents](README.md) · [← Profile analysis](profiles.md) · [Period layout and resampling →](periods.md)
