#ifndef CONST_H_
#define CONST_H_
#include <cstdint>

constexpr uint32_t s_ambig_char = 4;
constexpr uint32_t s_bits_per_byte = 8;
constexpr float s_min_accepted = 1e-28f;
constexpr uint32_t s_num_distinct_chars = 5;

constexpr uint32_t s_mask_all_ones_s = 0xFFFFFFFF;
constexpr uint64_t s_mask_all_ones_d = 0xFFFFFFFFFFFFFFFF;

// read base index
constexpr uint32_t s_rbi_0 = 0x0;
constexpr uint32_t s_rbi_1 = 0x1;
constexpr uint32_t s_rbi_2 = 0x2;
constexpr uint32_t s_rbi_3 = 0x3;
constexpr uint32_t s_rbi_4 = 0x4;
constexpr uint32_t s_rbi_5 = 0x5;
constexpr uint32_t s_rbi_6 = 0x6;
constexpr uint32_t s_rbi_7 = 0x7;
constexpr uint32_t s_rbi_8 = 0x8;
constexpr uint32_t s_rbi_9 = 0x9;
constexpr uint32_t s_rbi_A = 0xA;
constexpr uint32_t s_rbi_B = 0xB;
constexpr uint32_t s_rbi_C = 0xC;
constexpr uint32_t s_rbi_D = 0xD;
constexpr uint32_t s_rbi_E = 0xE;
constexpr uint32_t s_rbi_F = 0xF;

// move count
constexpr int32_t s_mc_0 = 0x0;
constexpr int32_t s_mc_1 = 0x1;
constexpr int32_t s_mc_2 = 0x2;
constexpr int32_t s_mc_3 = 0x3;
constexpr int32_t s_mc_4 = 0x4;
constexpr int32_t s_mc_5 = 0x5;
constexpr int32_t s_mc_6 = 0x6;
constexpr int32_t s_mc_7 = 0x7;
constexpr int32_t s_mc_8 = 0x8;
constexpr int32_t s_mc_9 = 0x9;
constexpr int32_t s_mc_A = 0xA;
constexpr int32_t s_mc_B = 0xB;
constexpr int32_t s_mc_C = 0xC;
constexpr int32_t s_mc_D = 0xD;
constexpr int32_t s_mc_E = 0xE;
constexpr int32_t s_mc_F = 0xF;

// shift count
constexpr int32_t s_sc_40 = 0x40;
constexpr int32_t s_sc_3F = 0x3F;
constexpr int32_t s_sc_3E = 0x3E;
constexpr int32_t s_sc_3D = 0x3D;
constexpr int32_t s_sc_3C = 0x3C;
constexpr int32_t s_sc_3B = 0x3B;
constexpr int32_t s_sc_3A = 0x3A;
constexpr int32_t s_sc_39 = 0x39;
constexpr int32_t s_sc_20 = 0x20;
constexpr int32_t s_sc_1F = 0x1F;
constexpr int32_t s_sc_1E = 0x1E;
constexpr int32_t s_sc_1D = 0x1D;
constexpr int32_t s_sc_1C = 0x1C;
constexpr int32_t s_sc_1B = 0x1B;
constexpr int32_t s_sc_1A = 0x1A;
constexpr int32_t s_sc_19 = 0x19;
constexpr int32_t s_sc_18 = 0x18;
constexpr int32_t s_sc_17 = 0x17;
constexpr int32_t s_sc_16 = 0x16;
constexpr int32_t s_sc_15 = 0x15;
constexpr int32_t s_sc_14 = 0x14;
constexpr int32_t s_sc_13 = 0x13;
constexpr int32_t s_sc_12 = 0x12;
constexpr int32_t s_sc_11 = 0x11;

// reserved mask, 后缀代表保留位数
constexpr int32_t s_rm_0 = 0b0;
constexpr int32_t s_rm_1 = 0b1;
constexpr int32_t s_rm_2 = 0b11;
constexpr int32_t s_rm_3 = 0b111;
constexpr int32_t s_rm_4 = 0b1111;
constexpr int32_t s_rm_5 = 0b11111;
constexpr int32_t s_rm_6 = 0b111111;
constexpr int32_t s_rm_7 = 0b1111111;
constexpr int32_t s_rm_8 = 0b11111111;
constexpr int32_t s_rm_9 = 0b111111111;
constexpr int32_t s_rm_10 = 0b1111111111;
constexpr int32_t s_rm_11 = 0b11111111111;
constexpr int32_t s_rm_12 = 0b111111111111;
constexpr int32_t s_rm_13 = 0b1111111111111;
constexpr int32_t s_rm_14 = 0b11111111111111;
constexpr int32_t s_rm_15 = 0b111111111111111;

#endif  // CONST_H_
