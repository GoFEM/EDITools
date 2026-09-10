# EDI Tools

A simple GUI tool to work with large collections of [Magnetitelluric](https://en.wikipedia.org/wiki/Magnetotellurics) data in the Electrical Data Interchange (EDI) format, which is [SEG standard for MT and EMAP data](https://pubs.usgs.gov/of/2003/of03-056/Data/Edistd.doc).

![](EDITools_demo.gif)

## Features

- Visualize all transfer functions 
- Assign error floors (with a proper error propagation from impedances to derived quantities such as apparent resitivity, phase or phase tensor)
- Decimate data
- Save a PDF graphic of all plots
- Mask any data
- Save/Load a project with all masking and settings preserved
- Export a subset of data to a [GoFEM](https://github.com/GoFEM/pyGoFEM) format for the subsequent inversion.
- Visualize reponses calculated by the [GoFEM](https://github.com/GoFEM/pyGoFEM)

## Dependencies

- Qt >= 5.0
- Eigen >= 3
- CMake >= 3.10
- Boost >= 1.84
- QCustomPlot >= 2.0
- PROJ (UTM coordinate conversion)

## Building

The code was tested for linux using GNU Toolchain. Provided dependencies above are installed (for instance using conda) and CMake knows where to find them, the configuring and building is rather simple

> mkdir build; cd build

> cmake -DQCUSTOMPLOT_PATH:PATH=/path/to/qcustomplot -DEIGEN_PATH:PATH=/path/to/eigen ../

> make

## Usage

### Coordinates and export

Use **Coordinates → Convert latitude/longitude to UTM…** to project the loaded
station locations using WGS84. The tool suggests one zone and hemisphere for
the survey; both can be changed. Check **Center coordinates on the survey**
to subtract the midpoint of the stations' projected bounds. The origin easting
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
Station selection and map labels follow the selected view. **Coordinates →
Show latitude/longitude on map** is a shortcut to the geographic view.
Elevations and response tensors are unchanged by this
coordinate operation. Projection uses [PROJ](https://proj.org/en/stable/development/quickstart.html).

**File → Export to MT Inversion…** writes the format described in
[MT_DATA_EXPORT_SPEC.md](MT_DATA_EXPORT_SPEC.md). Select the observations and
periods, then export; no coordinate file or convention confirmations are
required. It uses the survey's UTM coordinates, automatically choosing an
uncentered WGS84/UTM projection if none has been applied. The native receiver
file uses `x = northing`, `y = easting`, `z = -elevation` (metres).
Optional sign changes and a distortion-compatibility check are under **Advanced**.

All 24 scalar types are supported, including apparent resistivity and phase.
Phases retain their impedance quadrant, independent of plot wrapping. Missing
or masked data are omitted, valid zeros are preserved, and invalid errors stop
export with a diagnostic. Resistivity/phase data require `Distortion estimate = false`.

**Every data export writes a plain frequency file**: sorted distinct frequencies
in Hz, one number per line, with no header. It contains only frequencies for
which observations were actually exported. Native output `survey.data` has
companions `survey.recvs` and `survey.freqs`; GoFEM output retains its appended
companion naming, `<output>.recvs` and `<output>.freqs`.

Run the format, coordinate and GUI workflow tests after building with
`ctest --test-dir build/Release --output-on-failure`.

The animation above illustrating main features. 

The following hotkeys are also useful when you click on a transfer function plot:

- **Scroll**: zoom in/out both axes.
- **X + Scroll**: zoom in/out only the *Y* axis.
- **D + Click**: pan the plot.

## Contributing

Feel free to create an issue or pull request in case you want to report a problem or contribute/fix the code.
