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

#include "EDIFileReader.h"

#include <fstream>
#include <sstream>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>

// Declaration of procedures from mtcomp_interface.f90
extern "C"
{
void call_mtcomp(const double *spectra, const double frequency,
                 const double avgt,
                 std::complex<double> *Z, std::complex<double> *T,
                 double *rho, double *phase,
                 double *Ze, double *Te,
                 double *rho_err, double *phase_err);
}

EDIFileReader::EDIFileReader(const std::string edi_file_name):
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

    std::string clean_content;

    std::string line;
    while (!ifs.eof())
    {
        std::getline (ifs, line);
        // Trim string
        line = trim(line, " \t\r\n");

        // Skip empty lines and comments
        if (line.length() < 1 || (line[0] == '>' && line[1] == '!'))
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
    if(line[0] != '>')
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

        if(line[0] == '>')
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

    empty_value = get_option_value<double>(head_options, "EMPTY");
}

void EDIFileReader::read_info_block(std::istringstream &ss)
{
    std::string line;

    while(ss.good())
    {
        auto pos = ss.tellg();
        std::getline (ss, line);

        if(line[0] == '>')
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

        if(line[0] == '>')
        {
            if(line.find("HMEAS") != std::string::npos ||
                line.find("EMEAS") != std::string::npos)
            {
                StringMap meas_options;

                boost::char_separator<char> sep(" ");
                boost::tokenizer<boost::char_separator<char>> tokens(line, sep);

                for(auto &token: tokens)
                {
                    auto pos_token = token.find_first_of('=');

                    if(pos_token == std::string::npos)
                        continue;

                    meas_options.insert({token.substr(0, pos_token),
                                         token.substr(pos_token + 1, token.length())});
                }

                definemeas.MEAS.push_back(meas_options);

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

    std::string line;
    unsigned n_chanells;

    while(ss.good())
    {
        auto pos = ss.tellg();
        std::getline (ss, line);

        if(line[0] == '>')
        {
            if(line.find("SPECTRA") != std::string::npos)
            {
                SPECTRA_DATA spectra;

                boost::char_separator<char> sep(" ");
                boost::tokenizer<boost::char_separator<char>> tokens(line, sep);

                std::vector<std::string> str_tokens;
                for(auto &token: tokens)
                {
                    str_tokens.push_back(token);

                    auto pos_token = token.find_first_of('=');

                    if(pos_token == std::string::npos)
                        continue;

                    spectra.options.insert({token.substr(0, pos_token),
                                            token.substr(pos_token + 1, token.length())});

                }

                unsigned n_data = boost::lexical_cast<unsigned>(str_tokens.back());
                if(n_data != n_chanells*n_chanells)
                    throw std::ios_base::failure("Error reading spectral data. Number of channels does not match number of data");

                spectra.data.resize(n_chanells, n_chanells);

                for(unsigned i = 0; i < n_chanells; ++i)
                {
                    for(unsigned j = 0; j < n_chanells; ++j)
                    {
                        double val;
                        ss >> val;
                        spectra.data(i, j) = val;
                    }
                }

                spectrasect.spectra_data.push_back(spectra);

                continue;
            }
            else
            {
                ss.seekg(pos);
                break;
            }
        }

        auto pos_token = line.find_first_of('=');
        if(pos_token != std::string::npos)
        {
            spectrasect.options.insert({line.substr(0, pos_token),
                                        line.substr(pos_token + 1, line.length())});

            continue;
        }

        pos_token = line.find_first_of("//");
        if(pos_token != std::string::npos)
        {
            const std::string str = trim(line.substr(pos_token + 2, line.length()), " ");
            n_chanells = boost::lexical_cast<unsigned>(str);
            for(unsigned i = 0; i < n_chanells; ++i)
            {
                std::string ch_id;
                ss >> ch_id;
                //        std::getline (ss, line);
                spectrasect.channel_ids.push_back(ch_id);
            }

            continue;
        }
    }
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

        if(line[0] == '>')
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

std::string EDIFileReader::trim(const std::string& str, const std::string& whitespace)
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
    // factor to convert impedances from field units to S.I.
    // This is mu0 (to convert from B in nT to H in gamma)
    // divided by 10^-3 (to convert from mV/km to V/m)
    const double factor = 4.*M_PI*1.e-4;

    Matrix2cd Z;
    Matrix2d Ze, Rho, Phs, Rho_e, Phs_e, PT;
    Eigen::Vector2cd T;
    Eigen::Vector2d Te;

    station_data.set_size(spectrasect.spectra_data.size());

    //  station_data.station_name = get_option_value<std::string>(head_options, "DATAID");
    boost::filesystem::path p(m_edi_file_name);
    station_data.station_name = p.stem().string();

    station_data.file_name = m_edi_file_name;

    for(unsigned fidx = 0; fidx < spectrasect.spectra_data.size(); ++fidx)
    {
        auto &spectra = spectrasect.spectra_data[fidx];

        const double frequency = get_option_value<double>(spectra.options, "FREQ");
        const double avgt = get_option_value<double>(spectra.options, "AVGT");

        call_mtcomp(spectra.data.data(), frequency, avgt,
                    Z.data(), T.data(),
                    Rho.data(), Phs.data(),
                    Ze.data(), Te.data(),
                    Rho_e.data(), Phs_e.data());

        Z = factor*Z;
        Ze = factor*Ze;

        Matrix2d X, Y;

        X(0,0) = Z(0,0).real(); X(0,1) = Z(0,1).real();
        X(1,0) = Z(1,0).real(); X(1,1) = Z(1,1).real();

        Y(0,0) = Z(0,0).imag(); Y(0,1) = Z(0,1).imag();
        Y(1,0) = Z(1,0).imag(); Y(1,1) = Z(1,1).imag();

        PT = X.inverse() * Y;

        station_data.freqs[fidx] = frequency;
        for(unsigned i = 0; i < 2; ++i)
        {
            station_data.T[i][fidx] = T(i);
            station_data.T_err[i][fidx] = Te(i);

            for(unsigned j = 0; j < 2; ++j)
            {
                station_data.Z[i*2+j][fidx] = Z(i, j);
                station_data.PT[i*2+j][fidx] = PT(i, j);
                station_data.Rho[i*2+j][fidx] = Rho(i, j);
                station_data.Phs[i*2+j][fidx] = Phs(i, j);

                station_data.Z_err[i*2+j][fidx] = Ze(i, j);
                station_data.PT_err[i*2+j][fidx] = 0.;
                station_data.Rho_err[i*2+j][fidx] = Rho_e(i, j);
                station_data.Phs_err[i*2+j][fidx] = Phs_e(i, j);
            }
        }
    }
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

    //  station_data.station_name = get_option_value<std::string>(head_options, "DATAID");

    boost::filesystem::path p(m_edi_file_name);
    station_data.station_name = p.stem().string();

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
    auto long_str = get_option_value<std::string>(head_options, "LONG");
    if(long_str.empty())
        long_str = get_option_value<std::string>(head_options, "LON");

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
        return boost::lexical_cast<T>(it->second);
    }
    catch(const boost::bad_lexical_cast &)
    {
        return T();
    }
}

const MTStationData &EDIFileReader::get_mt_data() const
{
    return station_data;
}
