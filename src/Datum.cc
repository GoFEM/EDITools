#include "include/datum.h"

const std::map<ComplexDataType, std::pair<RealDataType, RealDataType>> cmplx_to_real_table =
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

const std::map<RealDataType, RealDataType> cmplx_counterparts_table =
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

const std::vector<std::pair<std::string, RealDataType>> data_type_conversion =
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

const std::vector<std::pair<std::string, ComplexDataType>> complex_data_type_conversion =
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

const std::map<SurveyMethod, std::set<RealDataType>> method_data_types_table =
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
