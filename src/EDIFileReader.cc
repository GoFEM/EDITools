/*
 * The EDI Tools application.
 *
 * Copyright (C) 2024 Alexander Grayver <agrayver.geophysics@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "include/EDIFileReader.h"

#include <fstream>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

namespace {
using Matrix2cd = Eigen::Matrix<std::complex<double>, 2, 2, Eigen::RowMajor>;
using Matrix2d = Eigen::Matrix<double, 2, 2, Eigen::RowMajor>;

// Some writers emit KEY= VALUE. Normalize once for both measurement and
// spectra headers, then extract their inline options.
StringMap inlineOptions(std::string line)
{
    for(size_t p = 0; (p = line.find('=', p)) != std::string::npos; ++p) {
        size_t q = p + 1;
        while(q < line.size() && std::isspace(static_cast<unsigned char>(line[q]))) ++q;
        line.erase(p + 1, q - p - 1);
    }
    std::istringstream tokens(line);
    StringMap options;
    for(std::string token; tokens >> token;) {
        const auto pos = token.find('=');
        if(pos != std::string::npos) options.emplace(token.substr(0, pos), token.substr(pos + 1));
    }
    return options;
}
}

EDIFileReader::EDIFileReader(const std::string &edi_file_name):
    m_edi_file_name(edi_file_name), empty_value(1.0E+32)
{
    read_edi();
    if(is_data_spectra)
        calculate_data_from_spectra();
    else
        calculate_data_from_mtsect();

    fill_stations_location();
}

void EDIFileReader::read_edi()
{
    std::ifstream ifs(m_edi_file_name);

    if(!ifs.is_open())
        throw std::ios_base::failure("Cannot open file " + m_edi_file_name);

    station_data.file_name = m_edi_file_name;

    std::string clean_content;

    std::string line;
    while (std::getline(ifs, line))
    {
        // Trim string
        line = trim(line, " \t\r\n");

        // Skip empty lines and comments
        if (line.empty() || line.compare(0, 2, ">!") == 0)
            continue;

        clean_content += line + '\n';
    }

    std::istringstream ss(clean_content);

    while(ss.good())
    {
        std::getline (ss, line);

        if( ss.eof() ) break;

        const BlockID block_id = get_block_id(line);

        if(block_id == EDI_END)
            break;

        switch (block_id) {
        case EDI_HEAD:
            read_head_block(ss);
            break;
        case EDI_INFO:
            read_info_block(ss);
            break;
        case EDI_DEFINEMEAS:
            read_definemeas_section(ss);
            break;
        case EDI_MTSECT:
            read_mtsect_section(ss);
            break;
        case EDI_SPECTRASECT:
            read_spectrasect_section(ss);
            break;
        default:
            break;
        }
    }
}

BlockID EDIFileReader::get_block_id(const std::string &line) const
{
    if(line.empty() || line[0] != '>')
        return EDI_INVALID;

    for(auto &p: block_names)
        if(line.find(p.first) != std::string::npos)
            return p.second;

    return EDI_INVALID;
}

void EDIFileReader::read_head_block(std::istringstream &ss)
{
    std::string line;

    while(ss.good())
    {
        auto pos = ss.tellg();
        std::getline (ss, line);

        if(!line.empty() && line[0] == '>')
        {
            ss.seekg(pos);
            break;
        }

        auto pos_token = line.find_first_of('=');

        if(pos_token == std::string::npos)
            continue;

        head_options.insert({line.substr(0, pos_token),
                             line.substr(pos_token + 1, line.length())});
    }

    if(head_options.count("EMPTY"))
        empty_value = get_option_value<double>(head_options, "EMPTY");
    
    boost::filesystem::path p(m_edi_file_name);
    station_data.station_name = p.stem().string();
}

void EDIFileReader::read_info_block(std::istringstream &ss)
{
    std::string line;

    while(ss.good())
    {
        auto pos = ss.tellg();
        std::getline (ss, line);

        if(!line.empty() && line[0] == '>')
        {
            ss.seekg(pos);
            break;
        }

        edi_info += line + "\n";
    }
}

void EDIFileReader::read_definemeas_section(std::istringstream &ss)
{
    std::string line;

    while(ss.good())
    {
        auto pos = ss.tellg();
        std::getline (ss, line);

        if(!line.empty() && line[0] == '>')
        {
            if(line.find("HMEAS") != std::string::npos ||
                line.find("EMEAS") != std::string::npos)
            {
                definemeas.MEAS.push_back(inlineOptions(line));

                continue;
            }
            else
            {
                ss.seekg(pos);
                break;
            }
        }

        auto pos_token = line.find_first_of('=');

        if(pos_token == std::string::npos)
            continue;

        definemeas.options.insert({line.substr(0, pos_token),
                                   line.substr(pos_token + 1, line.length())});
    }
}

void EDIFileReader::read_spectrasect_section(std::istringstream &ss)
{
    is_data_spectra = true;
    unsigned channels = 0;
    for(std::string line; ss.good();) {
        const auto position = ss.tellg();
        if(!std::getline(ss, line)) break;
        if(line.empty()) continue;
        if(line[0] == '>') {
            if(line.find("SPECTRA") == std::string::npos) { ss.seekg(position); break; }
            if(channels != 7) throw std::runtime_error("Spectra require seven channels: Hx, Hy, Hz, Ex, Ey, Rx, Ry.");
            const auto separator = line.find("//");
            std::istringstream count(separator == std::string::npos ? "" : line.substr(separator + 2));
            unsigned values = 0;
            if(!(count >> values) || values != 49) throw std::runtime_error("Spectral matrix must contain 49 values.");
            SPECTRA_DATA spectra; spectra.options = inlineOptions(line);
            for(auto &row: spectra.data) for(auto &value: row) {
                std::string token;
                if(!(ss >> token)) throw std::runtime_error("Truncated spectral matrix.");
                std::replace(token.begin(), token.end(), 'D', 'E');
                std::replace(token.begin(), token.end(), 'd', 'e');
                size_t consumed = 0;
                try { value = std::stod(token, &consumed); }
                catch(const std::exception &) { throw std::runtime_error("Invalid value in spectral matrix."); }
                if(consumed != token.size()) throw std::runtime_error("Invalid value in spectral matrix.");
                if(value == empty_value) value = std::numeric_limits<double>::quiet_NaN();
            }
            spectrasect.spectra_data.push_back(std::move(spectra));
        } else if(line.find('=') != std::string::npos) {
            const auto separator = line.find('=');
            const auto key = trim(line.substr(0, separator), " \t");
            spectrasect.options[key] = trim(line.substr(separator + 1), " \t");
            if(key == "NCHAN") channels = get_option_value<unsigned>(spectrasect.options, key);
        } else if(line.find("//") != std::string::npos) {
            std::istringstream count(line.substr(line.find("//") + 2));
            if(!(count >> channels) || channels != 7) throw std::runtime_error("Spectra require seven channel identifiers.");
            for(unsigned i = 0; i < channels; ++i) {
                std::string id;
                if(!(ss >> id) || id[0] == '>') throw std::runtime_error("Truncated spectral channel list.");
                spectrasect.channel_ids.push_back(id);
            }
        }
    }
    if(spectrasect.spectra_data.empty()) throw std::runtime_error("No spectral matrices found.");
}

void EDIFileReader::read_mtsect_section(std::istringstream &ss)
{
    is_data_spectra = false;

    mtsect.Zdata.resize(8);
    mtsect.Zvar.resize(4);
    mtsect.Tdata.resize(4);
    mtsect.Tvar.resize(2);

    const std::map<std::string, unsigned> Ztype2column = {{"ZXXR", 0}, {"ZXXI", 1},
                                                          {"ZXYR", 2}, {"ZXYI", 3},
                                                          {"ZYXR", 4}, {"ZYXI", 5},
                                                          {"ZYYR", 6}, {"ZYYI", 7}};

    const std::map<std::string, unsigned> ZEtype2column = {{"ZXX.VAR", 0}, {"ZXY.VAR", 1},
                                                           {"ZYX.VAR", 2}, {"ZYY.VAR", 3}};

    const std::map<std::string, unsigned> Ttype2column = {{"TXR.EXP", 0}, {"TXI.EXP", 1},
                                                          {"TYR.EXP", 2}, {"TYI.EXP", 3}};

    const std::map<std::string, unsigned> TEtype2column = {{"TXVAR.EXP", 0}, {"TYVAR.EXP", 1}};

    std::string line;
    unsigned n_freq = 0;

    while(ss.good())
    {
        std::getline (ss, line);

        // std::cout << line << std::endl;
        // if(line.find(">END") != std::string::npos)
        //     break;

        if(!line.empty() && line[0] == '>')
        {
            if(line.find("FREQ") != std::string::npos &&
                line.find("NFREQ") == std::string::npos)
            //      if(line.find("FREQ") != std::string::npos)
            {
                auto pos_token = line.find_first_of("//");
                if(pos_token != std::string::npos)
                {
                    const std::string str = trim(line.substr(pos_token + 2, line.length()), " ");
                    n_freq = boost::lexical_cast<unsigned>(str);
                    for(unsigned i = 0; i < n_freq; ++i)
                    {
                        double val;
                        ss >> val;
                        mtsect.frequencies.push_back(val);
                    }
                }

                continue;
            }

            bool found = false;
            for(auto &p: Ztype2column)
            {
                if(line.find(p.first) != std::string::npos)
                {
                    for(unsigned i = 0; i < n_freq; ++i)
                    {
                        double val;
                        ss >> val;

                        val = (val == empty_value) ? std::numeric_limits<double>::quiet_NaN() : val;

                        mtsect.Zdata[p.second].push_back(val);
                    }

                    found = true;
                    break;
                }
            }

            if(found)
                continue;

            for(auto &p: ZEtype2column)
            {
                if(line.find(p.first) != std::string::npos)
                {
                    for(unsigned i = 0; i < n_freq; ++i)
                    {
                        double val;
                        ss >> val;

                        val = (val == empty_value) ? std::numeric_limits<double>::quiet_NaN() : val;

                        mtsect.Zvar[p.second].push_back(val);
                    }

                    found = true;
                    break;
                }
            }

            if(found)
                continue;

            for(auto &p: Ttype2column)
            {
                if(line.find(p.first) != std::string::npos)
                {
                    for(unsigned i = 0; i < n_freq; ++i)
                    {
                        double val;
                        ss >> val;

                        val = (val == empty_value) ? std::numeric_limits<double>::quiet_NaN() : val;

                        mtsect.Tdata[p.second].push_back(val);
                    }

                    found = true;
                    break;
                }
            }

            if(found)
                continue;

            for(auto &p: TEtype2column)
            {
                if(line.find(p.first) != std::string::npos)
                {
                    for(unsigned i = 0; i < n_freq; ++i)
                    {
                        double val;
                        ss >> val;

                        val = (val == empty_value) ? std::numeric_limits<double>::quiet_NaN() : val;

                        mtsect.Tvar[p.second].push_back(val);
                    }

                    found = true;
                    break;
                }
            }
        }
    }
}

std::string EDIFileReader::trim(const std::string& str, const std::string& whitespace) const
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    if(strBegin+strRange > str.length())
        throw std::runtime_error("trim: out of range");

    return str.substr(strBegin, strRange);
}

void EDIFileReader::calculate_data_from_spectra()
{
    // Convert EDI field impedance (mV/km)/nT to SI E/H; tippers are dimensionless.
    const double factor = 4. * std::acos(-1.) * 1.e-4;
    station_data.set_size(spectrasect.spectra_data.size(), true);
    for(unsigned f = 0; f < spectrasect.spectra_data.size(); ++f) {
        const auto &spectra = spectrasect.spectra_data[f];
        const double frequency = get_option_value<double>(spectra.options, "FREQ");
        if(!std::isfinite(frequency) || frequency <= 0.) throw std::runtime_error("Spectral FREQ must be finite and positive.");
        const auto transfer = MTSpectra::estimate(spectra.data, get_option_value<double>(spectra.options, "AVGT"));
        station_data.freqs[f] = frequency;
        for(unsigned c = 0; c < 4; ++c) {
            station_data.Z[c][f] = factor * transfer.impedance[c];
            station_data.Z_err[c][f] = factor * transfer.impedanceError[c];
        }
        for(unsigned c = 0; c < 2; ++c) {
            station_data.T[c][f] = transfer.tipper[c];
            station_data.T_err[c][f] = transfer.tipperError[c];
        }
    }
    station_data.Z_err_floor = station_data.Z_err;
    station_data.T_err_floor = station_data.T_err;
    station_data.calculate_apparent_resistivity();
    station_data.calculate_phase();
    station_data.calculate_phase_tensor();
    station_data.propagate_rho_phase_error();
    station_data.propagate_phase_tensor_error();
}

void EDIFileReader::calculate_data_from_mtsect()
{
    // factor to convert impedances from field units to S.I.
    // This is mu0 (to convert from B in nT to H in gamma)
    // divided by 10^-3 (to convert from mV/km to V/m)
    const double factor = 4.*M_PI*1.e-4;
    const double mu0 = 4.*M_PI*1e-7;

    Matrix2cd Z;
    Matrix2d PT;

    station_data.set_size(mtsect.frequencies.size());
    station_data.file_name = m_edi_file_name;

    station_data.freqs = mtsect.frequencies;

    for(unsigned fidx = 0; fidx < mtsect.frequencies.size(); ++fidx)
    {
        if(mtsect.Zdata[0].size() == mtsect.frequencies.size())
        {
            Z(0, 0) = std::complex<double>(mtsect.Zdata[0][fidx], mtsect.Zdata[1][fidx])*factor;
            Z(0, 1) = std::complex<double>(mtsect.Zdata[2][fidx], mtsect.Zdata[3][fidx])*factor;
            Z(1, 0) = std::complex<double>(mtsect.Zdata[4][fidx], mtsect.Zdata[5][fidx])*factor;
            Z(1, 1) = std::complex<double>(mtsect.Zdata[6][fidx], mtsect.Zdata[7][fidx])*factor;

            Matrix2d X, Y;

            X(0,0) = Z(0,0).real(); X(0,1) = Z(0,1).real();
            X(1,0) = Z(1,0).real(); X(1,1) = Z(1,1).real();

            Y(0,0) = Z(0,0).imag(); Y(0,1) = Z(0,1).imag();
            Y(1,0) = Z(1,0).imag(); Y(1,1) = Z(1,1).imag();

            PT = X.inverse() * Y;

            for(unsigned i = 0; i < 2; ++i)
            {
                for(unsigned j = 0; j < 2; ++j)
                {
                    station_data.Z[i*2+j][fidx] = Z(i, j);
                    station_data.PT[i*2+j][fidx] = PT(i, j);
                    station_data.Rho[i*2+j][fidx] = 1. / (mu0 * 2. * M_PI * station_data.freqs[fidx]) * pow(std::abs(Z(i, j)), 2.);
                    station_data.Phs[i*2+j][fidx] = 180. / M_PI * atan2(Z(i, j).imag(), Z(i, j).real());

                    station_data.Z_err[i*2+j][fidx] = sqrt(mtsect.Zvar[i*2+j][fidx])*factor;
                    station_data.PT_err[i*2+j][fidx] = 0.;
                    station_data.Rho_err[i*2+j][fidx] = 0;
                    station_data.Phs_err[i*2+j][fidx] = 0;
                }
            }
        }

        if(mtsect.Tdata[0].size() == mtsect.frequencies.size())
        {
            for(unsigned i = 0; i < 2; ++i)
            {
                station_data.T[i][fidx] = std::complex<double>(mtsect.Tdata[i*2][fidx], mtsect.Tdata[i*2+1][fidx]);
                station_data.T_err[i][fidx] = sqrt(mtsect.Tvar[i][fidx]);
            }
        }
    }
}

void EDIFileReader::fill_stations_location()
{
    auto lat_str = get_option_value<std::string>(head_options, "LAT");
    if(lat_str.empty())
        lat_str = get_option_value<std::string>(definemeas.options, "REFLAT");

    auto long_str = get_option_value<std::string>(head_options, "LONG");
    if(long_str.empty())
        long_str = get_option_value<std::string>(head_options, "LON");
    if(long_str.empty())
        long_str = get_option_value<std::string>(definemeas.options, "REFLONG");

    std::vector<double> lat_dms, long_dms;

    boost::char_separator<char> sep(":");
    {
        boost::tokenizer<boost::char_separator<char>> tokens(lat_str, sep);
        for(auto &token: tokens)
            lat_dms.push_back(boost::lexical_cast<double>(trim(token, " \t\r\n")));

        if(lat_dms.size() != 3)
        {
            station_data.location[0] = boost::lexical_cast<double>(trim(lat_str, " \t\r\n"));
        }
        else
        {
            const double multiplier = (lat_dms[0] < 0 ? -1 : 1);
            station_data.location[0] = multiplier * (fabs(lat_dms[0]) + (lat_dms[1] / 60) + (lat_dms[2] / 3600));
        }
    }

    {
        boost::tokenizer<boost::char_separator<char>> tokens(long_str, sep);
        for(auto &token: tokens)
            long_dms.push_back(boost::lexical_cast<double>(trim(token, " \t\r\n")));

        //if(long_dms.size() != 3)
        //    throw std::ios_base::failure("Wrong coordinate format.");

        if(long_dms.size() != 3)
        {
            station_data.location[1] = boost::lexical_cast<double>(trim(long_str, " \t\r\n"));
        }
        else
        {
            const double multiplier = (long_dms[0] < 0 ? -1 : 1);
            station_data.location[1] = multiplier * (fabs(long_dms[0]) + (long_dms[1] / 60) + (long_dms[2] / 3600));
        }

    }

    station_data.location[2] = get_option_value<double>(head_options, "ELEV");
    if(station_data.location[2] == 0.)
        station_data.location[2] = get_option_value<double>(definemeas.options, "REFELEV");
}

template<typename T>
T EDIFileReader::get_option_value(const StringMap &options,
                                  const std::string &option_name) const
{
    const auto it = options.find(option_name);
    if(it == options.end())
        return T();

    try
    {
        return boost::lexical_cast<T>(trim(it->second, " \t\r\n\""));
    }
    catch(const std::exception &)
    {
        return T();
    }
}

const MTStationData &EDIFileReader::get_mt_data() const
{
    return station_data;
}
