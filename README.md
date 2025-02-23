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

## Building

The code was tested for linux using GNU Toolchain. Provided dependencies above are installed (for instance using conda) and CMake knows where to find them, the configuring and building is rather simple

> mkdir build; cd build

> cmake -DQCUSTOMPLOT_PATH:PATH=/path/to/qcustomplot -DEIGEN_PATH:PATH=/path/to/eigen ../

> make

## Usage

The animation above illustrating main features. 

The following hotkeys are also useful when you click on a transfer function plot:

- **Scroll**: zoom in/out both axes.
- **X + Scroll**: zoom in/out only the *Y* axis.
- **D + Click**: pan the plot.

## Contributing

Feel free to create an issue or pull request in case you want to report a problem or contribute/fix the code.

