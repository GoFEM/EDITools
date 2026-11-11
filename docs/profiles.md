# Profile analysis

[Guide contents](README.md) · [← Maps and geographic layers](maps.md) · [Survey reports →](reports.md)

- [Draw and select](#draw-and-select)
- [Resistivity and phase sections](#resistivity-and-phase-sections)
- [Phase-tensor section](#phase-tensor-section)
- [Strike roses](#strike-roses)
- [Extract and export](#extract-and-export)

Open **Maps → Profile** to draw a corridor and inspect the stations along it.
The source survey, coordinates, errors and masks are unchanged.

## Draw and select

1. Leave **Draw line** checked and click successive vertices on the map.
2. Right-click to finish. Turn off **Draw line** to pan; the wheel zooms.
3. Set **Distance ± (km)** to the corridor distance on each side of the line.
   **Auto** restores the median nonzero nearest-neighbour spacing of enabled
   survey stations (1 km when spacing cannot be determined).
4. Use **Undo**, **Reverse** or **Clear** to edit the line. **Fit** restores the
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
same map controls as the other analysis windows.

## Resistivity and phase sections

Select observed data or a loaded calculated response under **Data**. Calculated
stations are matched by name to the corridor's observed stations. Choose
**XY / YX**, **Berdichevsky / SSQ**, or **Determinant**. These use the same complex
impedance definitions as [survey averages](analysis.md#survey-impedance-averages).
YX is sign-reversed; components are not rotated into the profile direction.

Horizontal position is cumulative profile distance. **Vertical** switches between
period and frequency; longer periods always lie lower. The **Periods (s)** bounds
filter all sections and roses. Section axes pan and zoom together, including
the tensor section. Hover over a sample to identify its station, period and value.

**Cells** paints bounded cells around measured samples; uncheck it for markers.
There is no interpolation. Masked and missing samples stay blank. **Gap ratio**
limits a cell's width in logarithmic period, and station cells are capped by the
survey's median nearest-neighbour spacing. Large gaps therefore remain blank.
Stations with identical along-profile distances overlap; narrow the corridor
or use the table to identify them.

Choose the palette and separate resistivity/phase color limits. Resistivity
colors are logarithmic, phase colors linear; values beyond the limits use the
endpoint colors. Changing color limits never changes masks.

## Phase-tensor section

The **Phase tensors** tab places an ellipse at each usable station-period sample.
Choose maximum phase, minimum phase or skew for its color; adjust the color limits
and ellipse size. Ellipse geometry remains independent of the plot aspect ratio.
Within each glyph, up is geographic north and right is east. The vertical axis
represents period or frequency, not depth. Masked, missing and singular tensors
are omitted using the same validity checks as the period maps.

## Strike roses

Enable **PT azimuth**, **Real tipper normal** and/or **Imaginary tipper normal**.
North is 0°, clockwise is positive. Each rose bins axial directions over 0–180°
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

The line remains available while the window is open and updates after data/mask
edits. Loading another survey or changing its projection clears the line. Profile
geometry is not stored in project files; keep the CSV to record the selection.

---

[Guide contents](README.md) · [← Maps and geographic layers](maps.md) · [Survey reports →](reports.md)
