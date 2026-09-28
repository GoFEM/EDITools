#include "include/SurveyExport.h"
#include "include/MTSurveyData.h"
#include <cmath>
#include <iomanip>
#include <sstream>

NativeMT::Export SurveyExport::prepare(const MTSurveyData &survey, NativeMT::Options options, Format format) {
  std::vector<std::array<double, 3>> valid;
  for (const auto &p : survey.geographic_locations())
    if (std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]) && p[0] >= -80. && p[0] <= 84. &&
        std::abs(p[1]) <= 180.)
      valid.push_back(p);
  const auto coordinates =
      survey.coordinates().utm ? survey.coordinates() : SurveyCoordinates::suggested(valid);
  options.horizontal_axes_confirmed = options.conventions_confirmed = true;
  options.coordinate_description =
      coordinates.description() + (format == Format::Native
                                       ? "\nVertical coordinate = -elevation (m, positive down); origin "
                                         "elevation = 0; source datum retained."
                                       : "\nVertical coordinate = elevation (m, positive up); origin "
                                         "elevation = 0; source datum retained.");
  std::vector<const MTStationData *> stations;
  for (const auto &name : survey.get_stations_names()) {
    const auto &station = survey.get_station_data(name);
    stations.push_back(&station);
    if (!station.active())
      continue;
    std::array<double, 3> p{{NAN, NAN, NAN}};
    try {
      p = coordinates.transform({station.position()}).front();
    } catch (const std::exception &) {
    }
    if (format == Format::Native)
      p[2] = -p[2];
    options.receivers[name] = p;
  }
  if (format == Format::GoFEM) {
    std::set<RealDataType> types(options.types.begin(), options.types.end());
    for (auto type : options.types) {
      const auto other = Datum::get_complex_complement(type);
      if (other != InvalidType)
        types.insert(other);
    }
    options.types.assign(types.begin(), types.end());
  }
  auto result = NativeMT::prepare(stations, options);
  if (format == Format::GoFEM) {
    std::ostringstream data;
    data.imbue(std::locale::classic());
    data << std::setprecision(17);
    data << "# DataType Frequency SourceName ReceiverName Value Error\n"
         << "# SI units; x=North y=East; receiver z=elevation (positive up); phase in degrees\n";
    for (const auto &row : result.observations)
      data << Datum::convert_type_to_string(row.type) << ' ' << row.frequency << " Plane_wave "
           << row.receiver << ' ' << row.value << ' ' << row.error << '\n';
    result.data_text = data.str();
  }
  return result;
}
