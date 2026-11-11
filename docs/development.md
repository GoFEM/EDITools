# Workspace implementation

The application uses one `MainWindow` per survey. The station table and shared
map host sit outside `WorkspaceTabs`; the pinned TF view and analysis views sit
inside it. `MainWindowWorkspace.cpp` creates and connects views, switches the
active ribbon context, and saves the tab arrangement.

## Adding or changing an analysis

Derive from `WorkspaceView` and add the view to `createWorkspaceView` and the
workspace kind list. Keep existing kind identifiers stable: project files use
them to restore tabs. Use local pointers during construction; each tab owns its
analysis instance, so opening multiple views must not rely on a singleton window.

Override `refreshData` to receive the shared survey and calculated responses.
Survey edits mark views dirty; the active view refreshes immediately and other
views refresh on activation. `WorkspaceView::refreshSurvey` retains plot ranges
when the coordinate frame is unchanged. It discards pending ranges when the
frame changes. Views that compute plots later override `hasDeferredPlots` and
call `restorePendingPlotRanges` when the refreshed plots are ready.

Override `selectStation` when the view follows the main station table. Emit
`stationActivated` for navigation originating inside the view; avoid emitting it
from `selectStation`, which would feed the same selection back to the main window.

## Shared map and ribbon ownership

Use `setCentralMaps` for geographic plots that replace the survey map in the
shared map area. The view retains ownership even when the main window reparents
those widgets. `centralMapCoordinates` supplies the station-table coordinate
frame. Multiple registered maps use the ribbon selector, with only one displayed.

Profile and neighbour views instead draw overlays on the survey map. Their
`setSharedMapActive` implementations enable only the active tab's overlay.
`SharedMapInteraction` captures and restores the map's cursor and interaction
flags. Deactivate the previous view before activating the next one, so the next
view captures the original settings. Overlay callbacks must check whether their
map interaction is active. The helper also tolerates the map being destroyed
first during shutdown.

`addToolGroup` moves controls into view-owned ribbon pages. The ribbon borrows
these pages while the tab is active; it does not become their lifetime owner.
Use `control<T>` to locate editors or maps that have been reparented, and
`setToolVisible` to show or hide an editor together with its label.

## Saved state and checks

`src/WorkspaceState.cc` implements the common state format for editors, plot
ranges, splitters, internal tabs, and list/tree checkboxes. Keep widget object
names and station/response `Qt::UserRole` identities stable across refactors.
The plot traversal includes maps displayed outside the tab. Restoration blocks
editor signals, then runs the view's state-refresh callback. Domain-specific
state, such as profile vertices, remains in the derived view.

Run the full CTest suite after changing workspace lifecycle or ownership.
`command_ribbon_test` covers dirty refresh, coordinate changes, deferred zoom,
selection restoration, and shared-map lifetime. `response_plot_test` exercises
the integrated tab, station, map, masking, and project-save workflows. Profile,
survey-workflow, and auto-mask tests cover the individual analyses. Build and
portable-runtime checks should pass on both Linux and Windows before updating
the release packages.
