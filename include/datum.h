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

#ifndef _DATUM_H_
#define _DATUM_H_

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <set>
#include <ostream>

// Define real data types
enum RealDataType {// CSEM Data
    RealEx = 1, ImagEx = 2,
    RealEy = 3, ImagEy = 4,
    RealEz = 5, ImagEz = 6,
    RealHx = 11, ImagHx = 12,
    RealHy = 13, ImagHy = 14,
    RealHz = 15, ImagHz = 16,
    AmpEx = 21, PhsEx = 22,
    AmpEy = 23, PhsEy = 24,
    AmpEz = 25, PhsEz = 26,
    log10AmpEx = 27,
    log10AmpEy = 28,
    log10AmpEz = 29,
    AmpHx = 31, PhsHx = 32,
    AmpHy = 33, PhsHy = 34,
    AmpHz = 35, PhsHz = 36,
    log10AmpHx = 37,
    log10AmpHy = 38,
    log10AmpHz = 39,
    // MT Data
    RhoZxx = 101, PhsZxx = 102,
    RhoZxy = 103, PhsZxy = 104,
    RhoZyx = 105, PhsZyx = 106,
    RhoZyy = 107, PhsZyy = 108,
    RealZxx = 111, ImagZxx = 112,
    RealZxy = 113, ImagZxy = 114,
    RealZyx = 115, ImagZyx = 116,
    RealZyy = 117, ImagZyy = 118,
    log10RhoZxx = 122, log10RhoZxy = 123,
    log10RhoZyx = 124, log10RhoZyy = 125,
    // Phase tensor
    PTxx = 126, PTxy = 127,
    PTyx = 128, PTyy = 129,
    // Vertical magnetic TFs
    RealTzx = 133, ImagTzx = 134,
    RealTzy = 135, ImagTzy = 136,
    // Magnetic fields for two polarizations (3D MT inversion specific)
    RealHx1 = 41, ImagHx1 = 42,
    RealHx2 = 43, ImagHx2 = 44,
    RealHy1 = 45, ImagHy1 = 46,
    RealHy2 = 47, ImagHy2 = 48,
    // Global Data
    RealCResponse = 201, ImagCResponse = 202,
    RealQResponse = 203, ImagQResponse = 204,
    // Geoelectric data
    dU = 301, RhoApp = 302,
    // Receiver Functions
    RFValue = 501, InvalidType = 999
};

// Define complex data types
enum ComplexDataType {
    Ex = 501, Ey = 502, Ez = 503,
    Hx = 504, Hy = 505, Hz = 506,
    Zxx = 507, Zxy = 508,
    Zyx = 509, Zyy = 510,
    Tzx = 511, Tzy = 512,
    Hx1 = 513, Hx2 = 514,
    Hy1 = 515, Hy2 = 516,
    CResponse = 601,
    InvalidComplexType = 998
};

static std::map<ComplexDataType, std::pair<RealDataType, RealDataType>> cmplx_to_real_table =
    {
        {Ex, {RealEx, ImagEx}},
        {Ey, {RealEy, ImagEy}},
        {Ez, {RealEz, ImagEz}},
        {Hx, {RealHx, ImagHx}},
        {Hy, {RealHy, ImagHy}},
        {Hz, {RealHz, ImagHz}},
        {Zxx, {RealZxx, ImagZxx}},
        {Zxy, {RealZxy, ImagZxy}},
        {Zyx, {RealZyx, ImagZyx}},
        {Zyy, {RealZyy, ImagZyy}},
        {Tzx, {RealTzx, ImagTzx}},
        {Tzy, {RealTzy, ImagTzy}},
        {Hx1, {RealHx1, ImagHx1}},
        {Hx2, {RealHx2, ImagHx2}},
        {Hy1, {RealHy1, ImagHy1}},
        {Hy2, {RealHy2, ImagHy2}},
        {CResponse, {RealCResponse, ImagCResponse}}
};

static std::map<RealDataType, RealDataType> cmplx_counterparts_table =
    {
        {RealEx, ImagEx},
        {RealEy, ImagEy},
        {RealEz, ImagEz},
        {RealHx, ImagHx},
        {RealHy, ImagHy},
        {RealHz, ImagHz},
        {RealZxx, ImagZxx},
        {RealZxy, ImagZxy},
        {RealZyx, ImagZyx},
        {RealZyy, ImagZyy},
        {RealTzx, ImagTzx},
        {RealTzy, ImagTzy},
        {RealCResponse, ImagCResponse},
        {RealQResponse, ImagQResponse},
        {ImagEx, RealEx},
        {ImagEy, RealEy},
        {ImagEz, RealEz},
        {ImagHx, RealHx},
        {ImagHy, RealHy},
        {ImagHz, RealHz},
        {ImagZxx, RealZxx},
        {ImagZxy, RealZxy},
        {ImagZyx, RealZyx},
        {ImagZyy, RealZyy},
        {ImagTzx, RealTzx},
        {ImagTzy, RealTzy},
        {ImagCResponse, RealCResponse},
        {ImagQResponse, RealQResponse}
};

static std::vector<std::pair<std::string, RealDataType>> data_type_conversion =
    {
        {"RealEx", RealEx}, {"ImagEx", ImagEx},
        {"RealEy", RealEy}, {"ImagEy", ImagEy},
        {"RealEz", RealEz}, {"ImagEz", ImagEz},
        {"RealHx", RealHx}, {"ImagHx", ImagHx},
        {"RealHy", RealHy}, {"ImagHy", ImagHy},
        {"RealHz", RealHz}, {"ImagHz", ImagHz},
        {"AmpEx", AmpEx}, {"PhsEx", PhsEx},
        {"AmpEy", AmpEy}, {"PhsEy", PhsEy},
        {"AmpEz", AmpEz}, {"PhsEz", PhsEz},
        {"AmpHx", AmpHx}, {"PhsHx", PhsHx},
        {"AmpHy", AmpHy}, {"PhsHy", PhsHy},
        {"AmpHz", AmpHz}, {"PhsHz", PhsHz},
        {"log10AmpEx", log10AmpEx}, {"log10AmpEy", log10AmpEy},
        {"log10AmpEz", log10AmpEz}, {"log10AmpHx", log10AmpHx},
        {"log10AmpHy", log10AmpHy}, {"log10AmpHz", log10AmpHz},
        {"RhoZxx", RhoZxx}, {"PhsZxx", PhsZxx},
        {"RhoZxy", RhoZxy}, {"PhsZxy", PhsZxy},
        {"RhoZyx", RhoZyx}, {"PhsZyx", PhsZyx},
        {"RhoZyy", RhoZyy}, {"PhsZyy", PhsZyy},
        {"RealZxx", RealZxx}, {"ImagZxx", ImagZxx},
        {"RealZxy", RealZxy}, {"ImagZxy", ImagZxy},
        {"RealZyx", RealZyx}, {"ImagZyx", ImagZyx},
        {"RealZyy", RealZyy}, {"ImagZyy", ImagZyy},
        {"RealTzy", RealTzy}, {"ImagTzy", ImagTzy},
        {"RealTzx", RealTzx}, {"ImagTzx", ImagTzx},
        {"log10RhoZxx", log10RhoZxx}, {"log10RhoZxy", log10RhoZxy},
        {"log10RhoZyx", log10RhoZyx}, {"log10RhoZyy", log10RhoZyy},
        {"PTxx", PTxx}, {"PTxy", PTxy},
        {"PTyx", PTyx}, {"PTyy", PTyy},
        {"dU", dU}, {"RhoApp", RhoApp},
        {"RealCResponse", RealCResponse}, {"ImagCResponse", ImagCResponse},
        {"RealQResponse", RealQResponse}, {"ImagQResponse", ImagQResponse},
        {"RFValue", RFValue},
        {"InvalidType", InvalidType}
};

static std::vector<std::pair<std::string, ComplexDataType>> complex_data_type_conversion =
    {
        {"Ex", Ex}, {"Ey", Ey},
        {"Ez", Ez}, {"Hx", Hx},
        {"Hy", Hy}, {"Hz", Hz},
        {"Zxx", Zxx}, {"Zxy", Zxy},
        {"Zyx", Zyx}, {"Zyy", Zyy},
        {"Tzx", Tzx}, {"Tzy", Tzy},
        {"CResponse", CResponse},
        {"InvalidComplexType", InvalidComplexType}
};

// Data specific. Note that survey method is different from the
// forward calculator type in a sense that several codes can be
// used to model the same method. For example, depending on the
// model dimension, MT can be calculator with 1D analytical code
// or numerically using GoFEM, therefore both method and forward
// operator type need to be specified
enum SurveyMethod {Geoelectric, CSEM, MT, TEM, GlobalEM, MTSphere, ReceiverFunction};
enum ForwardType {MT3D_GOFEM, MT2D_GOFEM, MTLayered, C1D, CSEMLayered, X3DG, GEM1D, RF_AMMON, RF_REFLECTIVITY};
enum SourceType {PlaneWave = 0, Dipole,
                  SphericalHarmonic, DCDipole,
                  RadialSheet, CurrentShell,
                  ExternalSourceFile, PointEarthquake};

static std::map<SurveyMethod, std::set<RealDataType>> method_data_types_table =
    {
        {MT, {RhoZxx, PhsZxx, RhoZxy, PhsZxy, RhoZyx, PhsZyx, RhoZyy, PhsZyy,
              RealZxx, ImagZxx, RealZxy, ImagZxy, RealZyx, ImagZyx, RealZyy,
              ImagZyy, log10RhoZxx, log10RhoZxy, log10RhoZyx, log10RhoZyy,
              RealTzy, ImagTzy, RealTzx, ImagTzx, PTxx, PTxy, PTyx, PTyy}},
        {CSEM, {RealEx, ImagEx, RealEy, ImagEy, RealEz, ImagEz, RealHx, ImagHx, RealHy,
                ImagHy, RealHz, ImagHz, AmpEx, PhsEx, AmpEy, PhsEy, AmpEz, PhsEz,
                log10AmpEx, log10AmpEy, log10AmpEz, AmpHx, PhsHx, AmpHy, PhsHy,
                AmpHz, PhsHz, log10AmpHx, log10AmpHy, log10AmpHz}},
        {Geoelectric, {dU, RhoApp}},
        {GlobalEM, {RealEx, ImagEx, RealEy, ImagEy, RealEz, ImagEz, RealHx, ImagHx, RealHy,
                    ImagHy, RealHz, ImagHz, AmpEx, PhsEx, AmpEy, PhsEy, AmpEz, PhsEz,
                    log10AmpEx, log10AmpEy, log10AmpEz, AmpHx, PhsHx, AmpHy, PhsHy,
                    AmpHz, PhsHz, log10AmpHx, log10AmpHy, log10AmpHz,
                    RealCResponse, ImagCResponse, RealQResponse, ImagQResponse}},
        {TEM, {RealEx, RealEy, RealEz, RealHx, RealHy, RealHz}},
        {ReceiverFunction, {RFValue}}
};

/*
 * This error floor suggests that absolute data
 * values should not be smaller than that.
 */
static const double minimum_impedance = 1e-12;
static const double minimum_tipper = 1e-12;
static const double minimum_potential = 1e-7;
static const double minimum_rho = 1e-6;
static const double minimum_electric_field = 1e-15;
static const double minimum_magnetic_field = 1e-16;

enum DatumStatus { DatumObserved, DatumPredicted, DatumEmpty };

/*
 * This structure describes a single measurement, namely
 * it store its type, value and standard error
 */
struct Datum
{
    SurveyMethod method;
    ForwardType calculator_type;
    RealDataType type;
    double value, stderr;
    DatumStatus status;
    unsigned linear_index;

    // Convert data type string code to the internal type
    static RealDataType convert_string_to_type(const std::string &string_code)
    {
        try
        {
            unsigned data_code = std::stoul(string_code);
            return static_cast<RealDataType>(data_code);
        }
        catch(std::invalid_argument &e)
        {
            for(const auto &d: data_type_conversion)
                if(Datum::istrcompare(string_code, d.first))
                    return d.second;

            throw std::runtime_error("Unknown data type " + string_code);
        }
    }

    static std::string convert_type_to_string(RealDataType data_type)
    {
        std::string string_code;

        for(const auto &d: data_type_conversion)
            if(data_type == d.second)
                return d.first;

        return string_code;
    }

    static std::string convert_type_to_string(ComplexDataType data_type)
    {
        std::string string_code;

        for(const auto &d: complex_data_type_conversion)
            if(data_type == d.second)
                return d.first;

        return string_code;
    }

    /*
   * Returns true if data type is naturally real or
   * represents real part of a complex valued measurement.
   */
    static bool is_real(RealDataType data_type)
    {
        for(auto t: cmplx_to_real_table)
        {
            if(data_type == t.second.second)
                return false;
        }

        return true;
    }

    /*
   * If given type is part of a complex valued measurement,
   * return correponding complex data type, otherwise return
   * InvalidComplexType
   */
    static ComplexDataType get_parent_type(RealDataType type)
    {
        for(auto t: cmplx_to_real_table)
        {
            if(type == t.second.second || type == t.second.first)
                return t.first;
        }

        return InvalidComplexType;
    }

    /*
   * Given complex valued measurement, return correponding
   * real pair ofl data types, otherwise return InvalidComplexType
   */
    static std::pair<RealDataType, RealDataType> get_real_pair(ComplexDataType type)
    {
        auto it = cmplx_to_real_table.find(type);
        if(it != cmplx_to_real_table.end())
            return it->second;

        return std::make_pair(InvalidType, InvalidType);
    }

    static RealDataType get_complex_complement(RealDataType type)
    {
        auto it = cmplx_counterparts_table.find(type);
        if(it == cmplx_counterparts_table.end())
            return InvalidType;
        else
            return it->second;
    }

    static bool icompare_pred(unsigned char a, unsigned char b)
    {
        return std::tolower(a) == std::tolower(b);
    }

    static bool istrcompare(const std::string& a, const std::string& b)
    {
        if(a.length() != b.length())
            return false;
        else
            return std::equal(a.begin(), a.end(), b.begin(), icompare_pred);
    }
};

/*
 * Triple index to locate data easily using
 * transmitter#, receiver#, component
 */
template<typename DataType>
struct DatumKey
{
    DatumKey(const std::string &source, const std::string &receiver, DataType component):
        _source(source), _receiver(receiver), _component(component)
    {}

    bool operator==(const DatumKey &other) const
    {
        return (_source == other._source &&
                _receiver == other._receiver &&
                _component == other._component);
    }

    std::ostream& operator<< (std::ostream& stream)
    {
        stream << _source << " " << _receiver << " " << Datum::convert_type_to_string(_component);
        return stream;
    }

    std::string _source, _receiver;
    DataType _component;
};

/*
 * To be able to use index above in unordered map,
 * we need to define a hash function for it
 */
template<typename DataType>
struct DatumHash
{
    std::size_t operator()(const DatumKey<DataType>& k) const
    {
        return ((std::hash<std::string>()(k._source)
                 ^(std::hash<std::string>()(k._receiver) << 1)) >> 1)
               ^(std::hash<unsigned>()(k._component) << 1);
    }
};

#endif
