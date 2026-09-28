#include "include/Topography.h"
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <ogr_spatialref.h>
#include <cpl_string.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <locale>
#include <stdexcept>
#include <sstream>
#include <cctype>

namespace Topography {
namespace {
constexpr size_t maxCells = 16000000;
using Dataset = std::unique_ptr<GDALDataset, decltype(&GDALClose)>;
void require(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
void initialize() {
  static const bool done = [] {
    GDALAllRegister();
    return true;
  }();
  (void)done;
}
std::string number(double value) {
  std::ostringstream s;
  s.imbue(std::locale::classic());
  s << std::setprecision(17) << value;
  return s.str();
}
OGRSpatialReference reference(int epsg) {
  OGRSpatialReference s;
  require(s.importFromEPSG(epsg) == OGRERR_NONE, "Unknown raster coordinate system.");
  s.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
  return s;
}
void checkUtm(int epsg) {
  require((epsg >= 32601 && epsg <= 32660) || (epsg >= 32701 && epsg <= 32760), "Choose a WGS84 UTM grid.");
}
Dataset memory(const Grid &grid) {
  initialize();
  grid.validate();
  Dataset ds(GetGDALDriverManager()->GetDriverByName("MEM")->Create("", grid.width, grid.height, 1,
                                                                    GDT_Float32, nullptr),
             GDALClose);
  require(bool(ds), "Cannot allocate terrain raster.");
  double gt[]{grid.west, grid.spacing, 0., grid.north, 0., -grid.spacing};
  auto crs = reference(grid.epsg);
  ds->SetSpatialRef(&crs);
  ds->SetGeoTransform(gt);
  auto *band = ds->GetRasterBand(1);
  band->SetNoDataValue(NAN);
  band->SetUnitType("m");
  require(band->RasterIO(GF_Write, 0, 0, grid.width, grid.height, const_cast<float *>(grid.values.data()),
                         grid.width, grid.height, GDT_Float32, 0, 0) == CE_None,
          "Cannot write terrain pixels.");
  ds->SetMetadataItem("DEM_SOURCE", grid.source.c_str());
  ds->SetMetadataItem("VERTICAL_DATUM", grid.datum.c_str());
  ds->SetMetadataItem("DEM_LICENCE_NOTICE", grid.notice.c_str());
  ds->SetMetadataItem("HEIGHT_CONVENTION", "Elevation in metres, positive up");
  return ds;
}
int progressCallback(double value, const char *, void *context) {
  const auto *progress = static_cast<const Progress *>(context);
  return !*progress || (*progress)(value);
}
std::shared_ptr<Grid> warp(GDALDataset *input, const Request &r, const Progress &progress) {
  r.validate();
  require(input->GetRasterCount() == 1, "A DEM must have one elevation band.");
  double gt[6];
  require(input->GetGeoTransform(gt) == CE_None && input->GetSpatialRef(),
          "DEM needs a geotransform and coordinate reference system.");
  const int width = int(std::ceil((r.east - r.west) / r.spacing)),
            height = int(std::ceil((r.north - r.south) / r.spacing));
  // Use an explicit output grid. GDAL performs horizontal reprojection only;
  // vertical datum changes are deliberately not inferred.
  char **args = nullptr;
  for (const auto &s : std::vector<std::string>{
           "-of", "MEM", "-t_srs", "EPSG:" + std::to_string(r.epsg), "-te", number(r.west),
           number(r.north - height * r.spacing), number(r.west + width * r.spacing), number(r.north), "-ts",
           std::to_string(width), std::to_string(height), "-r", "bilinear", "-ot", "Float32", "-dstnodata",
           "nan", "-novshift"})
    args = CSLAddString(args, s.c_str());
  auto *options = GDALWarpAppOptionsNew(args, nullptr);
  CSLDestroy(args);
  require(options != nullptr, "Cannot configure DEM reprojection.");
  GDALWarpAppOptionsSetProgress(options, progressCallback, const_cast<Progress *>(&progress));
  GDALDatasetH source = input;
  int usageError = 0;
  Dataset output(static_cast<GDALDataset *>(GDALWarp("", nullptr, 1, &source, options, &usageError)),
                 GDALClose);
  GDALWarpAppOptionsFree(options);
  require(bool(output) && !usageError, "DEM reprojection failed or was canceled.");
  auto result = std::make_shared<Grid>();
  result->width = width;
  result->height = height;
  result->epsg = r.epsg;
  result->west = r.west;
  result->north = r.north;
  result->spacing = r.spacing;
  result->source = r.source;
  result->datum = r.datum;
  result->notice = r.notice;
  if (r.source == "User-provided DEM")
    if (const auto *value = input->GetMetadataItem("DEM_SOURCE"))
      result->source = value;
  if (r.datum == "Unknown")
    if (const auto *value = input->GetMetadataItem("VERTICAL_DATUM"))
      result->datum = value;
  if (r.notice.empty())
    if (const auto *value = input->GetMetadataItem("DEM_LICENCE_NOTICE"))
      result->notice = value;
  result->values.resize(size_t(width) * height);
  require(output->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, width, height, result->values.data(), width,
                                             height, GDT_Float32, 0, 0) == CE_None,
          "Cannot read terrain pixels.");
  auto *band = input->GetRasterBand(1);
  double factor = r.unitScale;
  if (factor == 0.) {
    std::string units = band->GetUnitType();
    std::transform(units.begin(), units.end(), units.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (units == "m" || units == "metre" || units == "meter" || units == "metres" || units == "meters")
      factor = 1.;
    else if (units == "ft" || units == "foot" || units == "feet")
      factor = .3048;
    else if (units == "us_survey_foot" || units == "us survey foot")
      factor = 1200. / 3937.;
    else
      throw std::runtime_error("DEM elevation units are unknown. Choose metres or feet explicitly.");
  }
  int present = 0;
  const double scale = band->GetScale(&present);
  const double multiplier = present ? scale : 1.;
  const double offset = band->GetOffset(&present);
  const double shift = present ? offset : 0.;
  for (auto &v : result->values)
    if (std::isfinite(v))
      v = float((v * multiplier + shift) * factor);
  result->validate();
  return result;
}
} // namespace
int utmEpsg(const SurveyCoordinates &c) {
  return (c.north ? 32600 : 32700) + c.zone;
}
void Request::validate() const {
  checkUtm(epsg);
  require(std::isfinite(west) && std::isfinite(east) && std::isfinite(south) && std::isfinite(north) &&
              east > west && north > south && std::isfinite(spacing) && spacing > 0.,
          "Invalid DEM bounds or cell size.");
  const double w = std::ceil((east - west) / spacing), h = std::ceil((north - south) / spacing);
  require(
      w >= 2 && h >= 2 && w * h <= maxCells,
      "Choose a smaller area or larger cell size (at most 16 million cells, at least two rows and columns).");
  require(std::isfinite(unitScale) && unitScale >= 0., "Invalid elevation unit conversion.");
}
void Grid::validate() const {
  checkUtm(epsg);
  require(width >= 2 && height >= 2 && size_t(width) * height <= maxCells &&
              values.size() == size_t(width) * height && std::isfinite(west) && std::isfinite(north) &&
              std::isfinite(spacing) && spacing > 0.,
          "Invalid stored DEM grid.");
}
std::array<double, 2> Grid::range() const {
  double lo = INFINITY, hi = -INFINITY;
  for (float v : values)
    if (std::isfinite(v)) {
      lo = std::min(lo, double(v));
      hi = std::max(hi, double(v));
    }
  return {{lo, hi}};
}
double Grid::sample(double x, double y) const {
  const double col = (x - west) / spacing - .5, row = (north - y) / spacing - .5;
  if (!std::isfinite(col) || !std::isfinite(row) || col < 0 || row < 0 || col > width - 1 || row > height - 1)
    return NAN;
  const int c = std::min(width - 2, int(col)), r = std::min(height - 2, int(row));
  const double dx = col - c, dy = row - r;
  double result = 0.;
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 2; ++i) {
      const double weight = (i ? dx : 1 - dx) * (j ? dy : 1 - dy);
      if (weight <= 1e-12)
        continue;
      const auto value = values[size_t(r + j) * width + c + i];
      if (!std::isfinite(value))
        return NAN;
      result += weight * value;
    }
  return result;
}
void transform(int from, int to, std::vector<double> &x, std::vector<double> &y) {
  require(x.size() == y.size(), "Mismatched coordinates.");
  if (from == to || x.empty())
    return;
  auto a = reference(from), b = reference(to);
  std::unique_ptr<OGRCoordinateTransformation, decltype(&OCTDestroyCoordinateTransformation)> t(
      OGRCreateCoordinateTransformation(&a, &b), OCTDestroyCoordinateTransformation);
  require(bool(t), "Cannot transform DEM coordinates.");
  std::vector<int> ok(x.size());
  t->Transform(int(x.size()), x.data(), y.data(), nullptr, ok.data());
  for (size_t i = 0; i < x.size(); ++i)
    if (!ok[i])
      x[i] = y[i] = NAN;
}
Request surveyArea(const std::vector<std::array<double, 3>> &locations, const SurveyCoordinates &coordinates,
                   double margin, double spacing) {
  require(std::isfinite(margin) && margin >= 0., "Invalid survey margin.");
  auto c = coordinates.utm ? coordinates : SurveyCoordinates::suggested(locations);
  c.origin_easting = c.origin_northing = 0.;
  const auto points = c.transform(locations);
  require(!points.empty(), "Load stations with valid coordinates first.");
  Request r;
  r.epsg = utmEpsg(c);
  r.spacing = spacing;
  r.west = r.south = INFINITY;
  r.east = r.north = -INFINITY;
  for (const auto &p : points) {
    r.west = std::min(r.west, p[1]);
    r.east = std::max(r.east, p[1]);
    r.south = std::min(r.south, p[0]);
    r.north = std::max(r.north, p[0]);
  }
  r.west -= std::max(margin, spacing);
  r.east += std::max(margin, spacing);
  r.south -= std::max(margin, spacing);
  r.north += std::max(margin, spacing);
  r.validate();
  return r;
}
std::array<double, 4> geographicBounds(const Request &r) {
  r.validate();
  std::vector<double> x, y;
  for (int i = 0; i <= 32; ++i) {
    const double t = i / 32.;
    for (double e : {r.west, r.east}) {
      x.push_back(e);
      y.push_back(r.south + t * (r.north - r.south));
    }
    for (double n : {r.south, r.north}) {
      x.push_back(r.west + t * (r.east - r.west));
      y.push_back(n);
    }
  }
  transform(r.epsg, 4326, x, y);
  std::array<double, 4> b{{INFINITY, INFINITY, -INFINITY, -INFINITY}};
  for (size_t i = 0; i < x.size(); ++i) {
    require(std::isfinite(x[i]) && std::isfinite(y[i]), "DEM bounds cannot be projected.");
    b[0] = std::min(b[0], x[i]);
    b[1] = std::min(b[1], y[i]);
    b[2] = std::max(b[2], x[i]);
    b[3] = std::max(b[3], y[i]);
  }
  require(b[2] - b[0] < 180., "Split DEM downloads at the date line.");
  return b;
}
std::shared_ptr<Grid> load(const std::vector<std::string> &files, const Request &r,
                           const Progress &progress) {
  initialize();
  r.validate();
  require(!files.empty(), "Choose a DEM file.");
  std::shared_ptr<Grid> result;
  for (size_t i = 0; i < files.size(); ++i) {
    const char *drivers[]{"GTiff", nullptr};
    Dataset ds(static_cast<GDALDataset *>(GDALOpenEx(files[i].c_str(), GDAL_OF_RASTER | GDAL_OF_READONLY,
                                                     drivers, nullptr, nullptr)),
               GDALClose);
    require(bool(ds), "Cannot open DEM GeoTIFF.");
    auto tile = warp(ds.get(), r, [&](double f) { return !progress || progress((i + f) / files.size()); });
    if (!result)
      result = tile;
    else {
      require(result->datum == tile->datum,
              "DEM tiles have different vertical datums. Load consistent inputs.");
      if (result->source != tile->source && result->source.find(tile->source) == std::string::npos)
        result->source += "; " + tile->source;
      if (!tile->notice.empty() && result->notice.find(tile->notice) == std::string::npos)
        result->notice += " " + tile->notice;
      for (size_t p = 0; p < result->values.size(); ++p)
        if (!std::isfinite(result->values[p]))
          result->values[p] = tile->values[p];
    }
  }
  require(std::isfinite(result->range()[0]), "No valid DEM coverage within the requested area.");
  return result;
}
std::shared_ptr<Grid> reproject(const Grid &g, const Request &r, const Progress &p) {
  auto ds = memory(g);
  auto request = r;
  request.unitScale = 1.;
  request.source = g.source;
  request.datum = g.datum;
  request.notice = g.notice;
  return warp(ds.get(), request, p);
}
void writeGeoTiff(const std::string &path, const Grid &grid, const SurveyCoordinates &coordinates) {
  auto ds = memory(grid);
  require(!coordinates.utm || utmEpsg(coordinates) == grid.epsg,
          "GeoTIFF grid must match the survey UTM zone.");
  ds->SetMetadataItem("SURVEY_ORIGIN_EASTING_M", number(coordinates.origin_easting).c_str());
  ds->SetMetadataItem("SURVEY_ORIGIN_NORTHING_M", number(coordinates.origin_northing).c_str());
  char **options = nullptr;
  options = CSLSetNameValue(options, "COMPRESS", "DEFLATE");
  options = CSLSetNameValue(options, "PREDICTOR", "3");
  options = CSLSetNameValue(options, "TILED", "YES");
  Dataset out(GetGDALDriverManager()->GetDriverByName("GTiff")->CreateCopy(path.c_str(), ds.get(), false,
                                                                           options, nullptr, nullptr),
              GDALClose);
  CSLDestroy(options);
  require(bool(out), "Could not export DEM GeoTIFF.");
  require(out->FlushCache() == CE_None, "Could not finish DEM GeoTIFF.");
}
void writeXYZ(const std::string &path, const Grid &g, const SurveyCoordinates &c, bool down) {
  g.validate();
  require(c.utm && utmEpsg(c) == g.epsg, "XYZ grid must match the survey UTM zone.");
  std::ofstream f(path);
  f.imbue(std::locale::classic());
  f << std::setprecision(12);
  require(bool(f), "Cannot create DEM XYZ file.");
  auto line = [](std::string s) {
    for (auto &c : s)
      if (c == '\n' || c == '\r')
        c = ' ';
    return s;
  };
  f << "# DEM source: " << line(g.source) << "\n# Vertical datum: " << line(g.datum)
    << "\n# Horizontal CRS: EPSG:" << g.epsg << "\n# UTM origin easting northing (m): " << c.origin_easting
    << ' ' << c.origin_northing << "\n# Columns: northing_offset easting_offset "
    << (down ? "depth_positive_down" : "elevation_positive_up")
    << " (m)\n# Missing cells omitted; cell size (m): " << g.spacing << '\n';
  if (!g.notice.empty())
    f << "# " << line(g.notice) << '\n';
  for (int row = 0; row < g.height; ++row)
    for (int col = 0; col < g.width; ++col) {
      const float z = g.values[size_t(row) * g.width + col];
      if (std::isfinite(z))
        f << g.north - (row + .5) * g.spacing - c.origin_northing << ' '
          << g.west + (col + .5) * g.spacing - c.origin_easting << ' ' << (down ? -z : z) << '\n';
    }
  f.close();
  require(bool(f), "Could not finish DEM XYZ export.");
}
} // namespace Topography
