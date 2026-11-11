# Profile analysis

[Guide contents](README.md) · [← Maps and geographic layers](maps.md) · [Survey reports →](reports.md)

- [Draw and select](#draw-and-select)
- [Resistivity and phase sections](#resistivity-and-phase-sections)
- [Rotate to strike](#rotate-to-strike)
- [Phase-tensor section](#phase-tensor-section)
- [Strike roses](#strike-roses)
- [Extract and export](#extract-and-export)

Open **Views → 2D profile**, or click **2D profile…** beside the station map,
to draw a corridor and inspect the stations along it.
The source survey, coordinates, errors and masks are unchanged.

## Draw and select

1. Leave **Draw line** checked and click successive vertices on the map.
2. Right-click to finish. Turn off **Draw line** to pan; the wheel zooms.
3. Set **Distance ± (km)** to the corridor distance on each side of the line.
   **Auto** restores the median nonzero nearest-neighbour spacing of enabled
   survey stations (1 km when spacing cannot be determined).
4. Use **Undo vertex**, **Reverse** or **Clear** to edit the line. **Fit** restores the
   survey map extent. Reverse changes the distance origin and station order.

Orange points and the table show the selected stations. Distance is measured
along the polyline to the nearest segment; offset is positive to the left of
its drawing direction. Endpoints have round corridor caps: beyond an endpoint,
the offset magnitude is distance to that endpoint. A station occurs only once
at a bend or crossing; equal-distance ties use the first segment.

Distances use one common WGS84/UTM projection. The tool uses the survey's saved
UTM zone and origin, or suggests a zone from its enabled stations. The map shows
the zone and origin. For wide surveys, review projection distortion under
[UTM / origin](coordinates-export.md#utm-and-survey-origin) before choosing a corridor.
Station names, offline map layers, online basemap and loaded topography use the
same map controls as the other analysis tabs. Line and corridor controls are in
**Profile tools**, section colors are in **Profile appearance**, and rotation
and rose settings are in **Strike tools**.

## Resistivity and phase sections

Select observed data or a loaded calculated response under **Data**. Calculated
stations are matched by name to the corridor's observed stations. Choose
**XY / YX**, **Berdichevsky / SSQ**, or **Determinant**. These use the same complex
impedance definitions as [survey averages](analysis.md#survey-impedance-averages).
YX is sign-reversed. Components use their original axes unless **Rotate to strike** is checked.

Horizontal position is cumulative profile distance. **Vertical** switches between
period and frequency; longer periods always lie lower. The **Periods (s)** bounds
filter all sections and roses. Section axes pan and zoom together, including
the tensor section. Small vertical station labels above each section follow the distance axis; zoom
in to separate crowded labels.
Click a column to highlight its station on the map and select and scroll to its
row in the table. Selecting a table row also updates the map, column highlight
and current TF station without switching workspace tabs.
Dragging pans without changing the selected station. Hover over a sample to
identify its station, period and value.

**Cells** paints bounded cells around measured samples; uncheck it for markers.
There is no interpolation. Masked and missing samples stay blank. **Gap ratio**
limits a cell's width in logarithmic period, and station cells are capped by the
survey's median nearest-neighbour spacing. Large gaps therefore remain blank.
Stations with identical along-profile distances overlap. Their label names the
frontmost station and counts the additional sites; select each site separately
in the table, or narrow the corridor. Long labels are shortened to fit; hover over
a sample or use the table for the full station name.

Choose the palette and separate resistivity/phase color limits. Resistivity
colors are logarithmic, phase colors linear; values beyond the limits use the
endpoint colors. Changing color limits never changes masks.

## Rotate to strike

Check **Rotate to strike** to rotate the displayed transfer-function axes while
keeping the drawn line, station positions, offsets and distances fixed. Angles
are clockwise from north, assuming the stored x/y channels already use common
north/east axes. This is a display transformation; the original survey, masks,
error floors and calculated responses remain unchanged.

**Auto (observed)** is the default angle source. It uses usable observed phase
tensors from the corridor, within the selected period range and rose band. The
skew and ellipticity thresholds always apply to this estimate, even if **Filter
PT** is unchecked for the rose. Nearly circular, singular, masked and missing
tensors are excluded. **Equal stations** gives each station total weight one;
otherwise each eligible station-period sample has equal weight.

The estimate averages principal-axis directions modulo 90° using the fourth
circular moment. Its station count and concentration **R₄** indicate agreement,
not statistical confidence; values below 0.25 are flagged as weak agreement.
The two possible strike axes retain a 90° ambiguity: **+90°** chooses the other
axis. Use **Manual** to enter a preferred angle. If no automatic estimate is
available, rotated plots stay blank and PDF export is disabled; choose Manual
or uncheck rotation to restore the original-frame view.

With rotation enabled, primed XY/YX quantities use the selected strike axes:
`Z′ = R Z Rᵀ`, with the new x axis along the chosen strike. Rotated components
require every contributing input component to be usable; a quarter-turn only
needs the corresponding original component. Current effective impedance errors
propagate through the rotation without resetting error floors. Berdichevsky,
SSQ and determinant averages are rotation-invariant and retain their original
values and uncertainties.

Tensor orientations and tipper-derived directions rotate into the same frame.
Ellipse up and rose 0° mean the chosen x/strike axis; right is the y axis at
+90°. Rose cardinal labels change to x/y. Tipper validity is independent of
impedance masks. The selected frame is recorded in plot/PDF headings and the
station CSV header. **Open as survey** still extracts the original data.

## Phase-tensor section

The **Phase tensors / strike** tab shows the tensor section and strike roses together.
Drag the divider to resize them. Its **Colormap** selector changes the tensor
colors and the rose colors; the resistivity/phase tab uses the same palette.
The tensor section places an ellipse at each usable station-period sample.
Choose maximum phase, minimum phase or skew for its color; adjust the color limits
and ellipse size. Ellipse geometry remains independent of the plot aspect ratio.
Within each glyph, up is the data x axis (north unless rotated) and right is the y axis. The vertical axis
represents period or frequency, not depth. Masked, missing and singular tensors
are omitted using the same validity checks as the period maps.

## Strike roses

Enable **PT azimuth**, **Real tipper normal** and/or **Imaginary tipper normal**.
North is 0° before rotation; after rotation the chosen strike/x axis is 0°. Clockwise is positive. Each rose bins axial directions over 0–180°
and mirrors them into the opposite half of the circle. The outer-ring scale is
shared across visible roses; percentages refer to the original axial bins.

PT azimuth has an additional 90° strike ambiguity. Tipper normals are perpendicular
to induction vectors and are indicators, not an independently determined geological
strike. Exactly circular tensors and zero tipper vectors are excluded. Optional
**Filter PT** restricts the rose by absolute skew and ellipticity
`(1 − axis ratio) / (1 + axis ratio)`; it does not alter sections or masks.

Choose all selected periods or a period decade, and a bin width. **Equal stations**
gives each contributing station total weight one, so dense period sampling does
not dominate. Uncheck it to weight every usable station-period sample equally.
The sample/station counts, axial circular mean and concentration `R` appear below
each rose. `R` near zero indicates no preferred axial direction; an undefined mean
is not reported as a strike.

## Extract and export

**Open as survey** opens an independent copy containing the corridor stations,
with their original data, masks and errors. It retains their original components
and coordinates, rather than rotating or projecting the transfer functions.

**Station CSV** exports station names, along-profile distances, signed offsets
and original latitude/longitude. Its header records the polyline vertices,
corridor distance, UTM zone and origin. **Save PDF** exports the corridor map,
resistivity/phase panels, tensor section and enabled roses using the current
selection, color limits and plot ranges. Map attributions are retained. Neither
export includes local file paths.

The project saves the line, corridor, rotation and display settings for each named
profile tab. The line is stored in geographic coordinates and follows changes
to the survey's UTM zone or origin. Survey mask edits and Undo/Redo refresh the
profile; masked or missing samples remain blank.

---

[Guide contents](README.md) · [← Maps and geographic layers](maps.md) · [Survey reports →](reports.md)
