#pragma once
// Rhombus Face Module - Quantization constants
// Decision: Global symmetric int8 quantization for both W and q.
// Confirmed: 2026-06-15

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>

namespace rhombus {
namespace face {

// Quantization bit width
constexpr uint32_t QUANT_BITS = 8;

// Quantization range: [-QUANT_MAX, +QUANT_MAX]
constexpr int8_t QUANT_MAX = 127;

// HyDia engine parameters derived from quantization
constexpr uint32_t MAT_BITS = 8;  // W matrix element bit width
constexpr uint32_t VEC_BITS = 8;  // q vector element bit width

// Default HE parameters for face pipeline
constexpr uint32_t DEFAULT_POLY_DEGREE = 8192;
constexpr uint32_t DEFAULT_SLOTS = DEFAULT_POLY_DEGREE / 2;  // N/2 = 4096

// Default threshold for face matching (cosine similarity)
constexpr double DEFAULT_THRESHOLD = 0.8;

// Quantize a float vector (assumed L2-normalized, values in [-1, 1])
// to symmetric int8 using global scale = QUANT_MAX.
inline std::vector<int8_t> quantize_int8(const std::vector<double>& v) {
    std::vector<int8_t> out(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        double clamped = std::max(-1.0, std::min(1.0, v[i]));
        int val = static_cast<int>(std::round(clamped * QUANT_MAX));
        // Clamp to int8 range
        val = std::max(-128, std::min(127, val));
        out[i] = static_cast<int8_t>(val);
    }
    return out;
}

// Dequantize int8 back to float (for reference comparison)
inline std::vector<double> dequantize_int8(const std::vector<int8_t>& v) {
    std::vector<double> out(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        out[i] = static_cast<double>(v[i]) / QUANT_MAX;
    }
    return out;
}

}  // namespace face
}  // namespace rhombus
