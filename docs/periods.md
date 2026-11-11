# Period layout and resampling

[Guide contents](README.md) · [← Survey reports](reports.md) · [Coordinates, topography and data export →](coordinates-export.md)

- [Period layout and resampling](#period-layout-and-resampling)

## Period layout and resampling

Open the view from **Edit → Layout / resample**. Grid, interpolation, preview
and output controls are in **Period tools**; source, target, coverage and
station summaries remain in the workspace. Grid fields appear for the selected
grid type and retain their values when hidden.

The period-layout, period-map and fit-statistics tabs use compact labels.
Hover over controls, their labels, or table headers for explanations. In the
period-layout tab, the **Coverage** tab and plot explain that **Usable on
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

---

[Guide contents](README.md) · [← Survey reports](reports.md) · [Coordinates, topography and data export →](coordinates-export.md)
