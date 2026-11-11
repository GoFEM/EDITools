# Response analysis

[Guide contents](README.md) · [← Editing and masking](editing.md) · [Maps and geographic layers →](maps.md)

- [Data fit analysis](#data-fit-analysis)
- [Neighbouring-station comparison](#neighbouring-station-comparison)
- [Impedance validity and rho+](#impedance-validity-and-rho)
- [Survey impedance averages](#survey-impedance-averages)
- [Polarity review](#polarity-review)

## Data fit analysis

Open **Views → Statistics / RMS** to open a comparison tab. Check or
uncheck any loaded response to include it in the overall nRMS, period, component,
station and normalized-residual plots. Each response keeps the same color across
the comparison plots; the list shows its nRMS and matched/available scalar count.
All residuals, including outliers, contribute to the histograms and statistics.

**Maps and heatmaps** compares two checked responses. The station map appears
in the shared area at the left; use **Fit appearance → Shared map** to choose
response A or B. Only one geographic map is visible. Heatmaps remain side by
side in the analysis tab, with shared axes and color scales. Stations are ordered north to south, and grey cells mark
missing comparisons. Select either response from the dropdowns. Pan and zoom
the plots, or use **Save PDF…** to export the overview, station plot and selected
spatial comparison.

Use **Fit appearance → Axis ranges** to choose a statistics plot and its X or Y
axis. Leave **Auto range** checked to fit the enabled responses, or turn it off
to set the minimum and maximum. **Fit all axes** resets the ranges. **Colors** changes the
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

## Neighbouring-station comparison

**Views → Neighbours** uses the current station in the main station table as its
reference and overlays its nearest enabled neighbours. Select another row in
the main table to change the reference and choose its nearest neighbours.
The ribbon shows the reference name; **Nearest** reapplies the neighbour count,
which specifies neighbours in addition to the reference.
Distances use WGS84 geodesics and handle the date line.

Each station has one visibility checkbox; clicking its marker on the shared
map at the left toggles the same selection. Select a row in the comparison list
to highlight its visible curves with thicker lines and larger filled markers,
without changing the reference or plot zoom. Coloured neighbour markers appear
only while that comparison tab is active; the tab has no duplicate map. Choose a component of resistivity/phase, tipper,
or phase tensor. Colours identify stations, **Errors** shows their effective
uncertainties, and phases retain their stored convention. Masked or missing
samples remain gaps. Disabled stations are excluded; stations without locations
can still be checked manually. This comparison does not change data or masks.

## Impedance validity and rho+

Open **Views → Rho+ analysis** to compare a station's measured
apparent resistivity and phase with fitted rho+ responses. Select the station;
**XY / YX** shows both off-diagonal components together, with independent fits
drawn as thin dashed lines in their component colors. The **XY** and **YX**
checkboxes control the plots, table rows and masking scope; **Rho+** toggles the
fitted lines. Individual diagonal components are not offered. Invariant responses
remain available in the response selector. The table and CSV identify each
component explicitly. **YX**
is shown as `−Zyx` to use the same sign convention as XY. Phase is not folded
into the first quadrant, so anomalous phases remain visible.

Use **TF tools → Rho+ guides** to overlay thin, dashed XY/YX fits on the main
resistivity and phase plots in their component colors. They follow the existing
component visibility and phase-wrapping controls and refresh after masking or
unmasking, preserving the current axes. The guides cannot be selected for masking
and can be shown alongside loaded responses. They use the main data's propagated
errors. Grid density comes from the validity tab; this setting and guide visibility are saved in
the project. Every fit uses only the remaining unmasked, enabled samples with
valid propagated errors. Finite provisional fits remain visible as dotted lines
when convergence or grid refinement checks fail; hover over **Rho+ guides** for
each component's status and fitted sample count. Such fits cannot drive automatic
masking. Fewer than three usable samples, invalid errors, or no finite response
leave no guide; the tooltip explains why. No previous fit is reused after masks change.

The native C++ implementation uses the bounded spectral representation of
log admittance from [Parker and Booker (1996)](https://doi.org/10.1016/S0031-9201(96)03191-3).
It minimizes weighted residuals in log resistivity and phase jointly, using
exact spectral-interval integrals, a free real offset, and both spectral tails.
It is a **finite-grid approximation**, with a second solve at twice the selected
**Grid/decade**. The view reports solver convergence and the change in χ²;
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

## Survey impedance averages

Open **Views → Survey averages** to display station responses and
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

## Polarity review

Open **Views → Polarity review** (or **Data → Polarity review…** in classic
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

---

[Guide contents](README.md) · [← Editing and masking](editing.md) · [Maps and geographic layers →](maps.md)
