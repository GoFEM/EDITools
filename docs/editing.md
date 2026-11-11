# Editing and masking

[Guide contents](README.md) · [← Getting started](getting-started.md) · [Response analysis →](analysis.md)

- [Component visibility](#component-visibility)
- [Linked impedance and phase-tensor masks](#linked-impedance-and-phase-tensor-masks)
- [Group masking](#group-masking)
- [Auto masking](#auto-masking)
- [Edit history](#edit-history)

## Component visibility

Component colors use two families: **XX / XY** are light / dark blue, and
**YX / YY** are dark / light orange. The same colors identify observed impedance
points, error bars, response lines, legends and survey PDF plots. Induction arrows
use black for real and grey for imaginary values, with pale grey for masked data.
Scalar tipper plots in analysis windows and reports retain the blue/orange palette.

Each impedance plot's legend has one checkbox per component. Uncheck **XX** and **YY**
to keep only off-diagonal impedances, or leave just one component checked.
Tipper legends offer **Real / Imaginary** arrows, and the phase-tensor strip has
one **PT** checkbox. Hidden components remain in the legend so they can be restored.
Click the checkbox or label, or use Tab and Space with the keyboard.

Each checkbox controls that plot's observed points, error bars and computed
curves. Plots have independent visibility settings, saved with the project;
data masks and fit statistics are unchanged. Partial tipper arrows are labeled
with their selected direction, and toggling them preserves that projection.

## Linked impedance and phase-tensor masks

Masking or unmasking **XY** also edits **XX**, and **YX** also edits **YY**, at
the same periods. This applies to resistivity and phase plots, inverted masks,
station component controls and group masking, including hidden diagonals.
Selecting a diagonal does not edit its off-diagonal partner. Rho+ uses the same
pairs; auto masking pairs them by default through **Pair diagonals**.

**Edit → Link Z / PT** is on by default and saved with the project.
It links mask/unmask operations from plots (including inverted masking),
station component controls, and period maps. Turning it off separates Z from PT;
diagonal pairing still applies. Toggling the option does not change existing masks.
Tippers are independent.

The link follows the general phase tensor definition `Phi = Re(Z)^-1 Im(Z)`.
Each PT entry depends on all four real impedance entries, so there is no general
one-to-one mapping between Z and PT components. For example, only in the special
2D case do `Zxy` and `Zyx` correspond to `PTyy` and `PTxx`, respectively.
See [Caldwell et al. (2004), equations 13–15](https://doi.org/10.1111/j.1365-246X.2004.02281.X).

Because each complex impedance component shares a real/imaginary mask, masking
one Z component (or its apparent resistivity/phase) also masks **all four PT
entries at that period**. Masking one PT entry also masks **all four Z components
at that period**. These Z/PT links do not propagate recursively; diagonal pairing
is applied separately. Unmasking applies the same links.
For inverted selections, a mask takes precedence if multiple selected components
give conflicting instructions for a shared dependent entry.

## Group masking

Select stations with Ctrl-click or Shift-click, then open **Edit → Group masking**
or **Group masking…** in the station context menu. Choose **Mask** or **Unmask**,
components, and all periods, an inclusive range, periods outside that range, or
individual stored periods. Ctrl+A selects all entries in the period list.
Individual periods match exactly; nearby periods are not merged or interpolated.

The preview shows matching periods and changes to Z, tipper and PT masks per
station, including linked components. **Link Z / PT** starts with the main setting
and can be overridden for this operation. Zxy also edits Zxx, and Zyx also edits
Zyy, regardless of that setting; tippers stay independent. Component choices are
independent of main-plot visibility. Selected disabled stations are included without enabling them.
Counts describe mask flags, including flags for missing values; unmasking does not
create data. **Apply** makes one undoable edit; **Cancel** leaves masks unchanged.

Four plots overlay all selected stations: apparent resistivity, stored phase
(without wrapping), real/imaginary tippers and phase tensor. Component colours
match the main plots. Grey crosses mark existing masks; large red circled crosses
show points to mask, and green circles show points to unmask. These markers use
the same edit plan as the table, including linked components. Missing values
remain gaps and have no marker; nonpositive resistivity cannot appear on a log axis.
Hover over points to identify stations. Each plot has component visibility
checkboxes, independent of the masking choices. Drag or scroll to inspect the
curves; period axes stay synchronized and **Fit view** restores the overview.

## Auto masking

Use the **Auto mask** ribbon tab (or the **Auto mask** classic menu). **Preview**
marks candidates with red circles directly on the main station plots, including
tipper lines or arrow anchors. **Previous / Next** switches stations without
changing the selected batch. **Apply** masks the current preview; **Apply + next**
also advances. **Apply to selected** evaluates the selected enabled sites with
the same settings and applies them in one Undo step. Cancelling that calculation
applies nothing. Unselected sites can provide neighbour references but are not edited.

The ribbon has independent check and component switches for **XY**, **YX**, **Tzx**
and **Tzy**. **Settings** opens the shared threshold editor; **Detailed review**
provides per-point checkboxes, scores, plots and CSV export for selected or all
enabled sites. Checks are **Rho+ / dispersion**, **Period outliers**, **Neighbours**,
and **Uncertainty**. Disabled checks are not evaluated.
**Minimum flags** defaults to one: any enabled check can nominate points for review.
Increase it to require agreement; this can exclude noisy points that only fail the
uncertainty check. Main-plot previews refresh after station, mask or setting changes.
Detailed review clears old suggestions until **Analyze** runs again.

- Rho+ uses the existing joint resistivity/phase fit and propagated main-data
  errors. Only converged, grid-stable fits contribute. Dispersion consistency is
  included in this check, rather than counted again as an independent test.
  Rho+ does not apply to tippers and never supplies a tipper vote.
- Period outliers use a robust trend in log period, excluding the tested point,
  with two reference samples on each side. Neither endpoints nor masked/missing
  reference samples or excessive period gaps are bridged. Resistivity and circular
  phase residuals share one vote, with measured errors and reference scatter.
  Tippers use real/imaginary residuals without logarithms or phase wrapping.
- Neighbours compares uniquely matched periods among nearby enabled sites,
  including sites outside the target selection. It requires at least three
  neighbours and 75% phase agreement. Distance, period tolerance, angular departure
  and agreement limits are adjustable. Impedance amplitude differences do not
  contribute. Tippers compare complex values with an adjustable minimum offset
  (default 0.05) and agreement distance (0.1).
- Uncertainty flags poorly resolved data, including endpoints, when the main-data
  one-sigma resistivity error exceeds 50% of resistivity or the phase error exceeds
  15°. Both limits are adjustable in **Settings**. Its score is the larger error /
  limit ratio, and its threshold defaults to 1×. This check catches noisy data whose
  large errors suppress normalized residuals; it does not change those errors.
  Tippers use an absolute error limit (default 0.05), so small or zero amplitudes
  with small errors are not rejected merely for their relative error.

The three consistency checks default to 3 sigma. These are standardized residuals for screening,
not calibrated probabilities: genuine spatial or multidimensional structure can
also disagree with these models. In particular, amplitude–phase dispersion
relations require additional assumptions beyond causality; see the
[dispersion-relation study](https://link.springer.com/article/10.1186/s40623-020-01273-4).
This tool does not implement a general 3-D Kramers–Kronig validity test.

Review scores and hover explanations in the table; check or uncheck individual
suggestions while viewing the plotted mask plan. **Flagged only** includes every
point flagged by any check, even if it falls below **Minimum flags**. **CSV…** saves
settings, scores, reasons and choices without source paths. **Pair diagonals** adds
Zxx for checked Zxy points and Zyy for Zyx; **Link Z / PT** adds PT masks. Both can be
disabled. Impedance/PT and tipper edits remain independent. A tipper candidate
masks both real and imaginary parts of that component only. **Apply checked** changes only the reviewed masks
in one Undo step; preview and analysis never alter data, errors or error floors.

## Edit history

**Edit → Undo / Redo** shares one history across station plots, map masking and
rho+ masking. It also covers error floors, station enable/disable, rename/delete,
EDI additions, decimation, UTM settings and topography/elevation edits. A batch
operation is one step; edits that leave the data unchanged add no step.
Select several stations and choose **Delete** from their context menu to remove
them together. One **Undo** restores the whole selection.

**Edit → History** lists the latest 20 operations. Select a step to move back
or forward; a new edit after undo discards the redo branch. **Open original
survey** opens an independent copy of the survey as it was loaded in this
session, for comparison. History is kept in memory and starts afresh when a
project is loaded. Project files save the current survey, not the undo history.
Resampling opens a separate survey with its own history and retains its source.

---

[Guide contents](README.md) · [← Getting started](getting-started.md) · [Response analysis →](analysis.md)
