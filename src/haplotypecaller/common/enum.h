#ifndef ROVACA_HC_ENUM_H_
#define ROVACA_HC_ENUM_H_
#include <cstdint>
#include <string>
#include <unordered_map>

enum class ReferenceConfidenceMode : uint8_t { NONE = 0, BP_RESOLUTION, GVCF };
enum class PcrIndelModel : uint8_t { NONE = 0, HOSTILE, AGGRESSIVE, CONSERVATIVE };

const std::unordered_map<std::string, PcrIndelModel> s_str2pcr{{"NONE", PcrIndelModel::NONE},
                                                               {"HOSTILE", PcrIndelModel::HOSTILE},
                                                               {"AGGRESSIVE", PcrIndelModel::AGGRESSIVE},
                                                               {"CONSERVATIVE", PcrIndelModel::CONSERVATIVE}};
const std::unordered_map<std::string, ReferenceConfidenceMode> s_str2erc{{"NONE", ReferenceConfidenceMode::NONE},
                                                                         {"BP_RESOLUTION", ReferenceConfidenceMode::BP_RESOLUTION},
                                                                         {"GVCF", ReferenceConfidenceMode::GVCF}};

static inline std::string pcr2str(PcrIndelModel pim)
{
    switch (pim) {
        case PcrIndelModel::NONE: return "NONE";
        case PcrIndelModel::HOSTILE: return "HOSTILE";
        case PcrIndelModel::AGGRESSIVE: return "AGGRESSIVE";
        case PcrIndelModel::CONSERVATIVE: return "CONSERVATIVE";
    }
    return "NULL";
}

static inline std::string erc2str(ReferenceConfidenceMode erc)
{
    switch (erc) {
        case ReferenceConfidenceMode::NONE: return "NONE";
        case ReferenceConfidenceMode::BP_RESOLUTION: return "BP_RESOLUTION";
        case ReferenceConfidenceMode::GVCF: return "GVCF";
    }
    return "NULL";
}

#endif  // ROVACA_HC_ENUM_H_
