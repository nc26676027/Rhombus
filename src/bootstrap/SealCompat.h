#pragma once

/**
 * SealCompat.h
 *
 * Compatibility layer: maps RBOOT's modified-SEAL (3.6.6-modified) private API calls
 * to standard SEAL 4.1.1 public API.
 *
 * Strategy: We extend seal::Evaluator with the missing methods via a wrapper
 * class CompatEvaluator, which holds references to the necessary SEAL objects.
 *
 * This file must be included AFTER "seal/seal.h".
 */

#include "seal/seal.h"
#include <NTL/RR.h>
#include <vector>
#include <complex>
#include <algorithm>

using namespace seal;

// ─── iter() macro ───
// In seal-modified, iter() returns an iterator adapter.
// For coeff_modulus: iter(parms.coeff_modulus()) → parms.coeff_modulus()
// For RNS polynomial: iter(ct)[poly_idx] → ct[poly_idx] (data_ptr access)
#define iter(x) (x)

/**
 * CompatEvaluator
 *
 * Wraps seal::Evaluator and provides the missing "reduced error" methods.
 * In the modified SEAL, "reduced error" variants skip some noise growth
 * optimizations. In standard SEAL, we simply use the normal operations.
 * The noise impact is slightly larger but functionally correct.
 */
class CompatEvaluator {
public:
    Evaluator &evaluator;
    CKKSEncoder &encoder;
    const SEALContext &context;
    double default_scale;

    CompatEvaluator(Evaluator &ev, CKKSEncoder &enc, const SEALContext &ctx, double scale)
        : evaluator(ev), encoder(enc), context(ctx), default_scale(scale) {}

    // ─── safe_scale_before_op ───
    // Compute a safe scale to temporarily set on a ciphertext before calling
    // a SEAL multiply/square/rescale operation, so that the internal
    // "is_scale_within_bounds" check passes.
    //
    // For multiply_plain(ct, plain): SEAL checks log2(ct.scale * plain.scale) < total_bits.
    //   We need: temp_scale such that log2(temp_scale * plain.scale) < total_bits.
    //   temp_scale = 2^(total_bits - log2(plain.scale) - 1)
    //
    // For multiply(ct1, ct2) / square(ct): SEAL checks log2(ct.scale * ct.scale) < total_bits.
    //   We need: temp_scale such that log2(temp_scale * temp_scale) < total_bits.
    //   temp_scale = 2^(total_bits/2 - 1)
    //
    // For rescale_to_next(ct): SEAL checks log2(ct.scale) < total_bits at next level.
    //   We need: temp_scale such that log2(temp_scale) < next_total_bits.
    //   temp_scale = 2^(next_total_bits - 1)
    //
    // The caller must restore the original scale after the SEAL operation.
    double safe_scale_for_multiply_plain(const Ciphertext &ct, double plain_scale) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return 1.0;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        int plain_bits = (plain_scale > 0) ? (int)ceil(log2(plain_scale)) : 1;
        int safe_bits = total_bits - plain_bits - 2; // -2 for safety margin
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    double safe_scale_for_square(const Ciphertext &ct) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return 1.0;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        int safe_bits = total_bits / 2 - 2; // need log2(s^2) < total_bits → log2(s) < total_bits/2
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    double safe_scale_for_rescale(const Ciphertext &ct) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return 1.0;
        auto next_data = ctx_data->next_context_data();
        if (!next_data) return ct.scale(); // at last level, can't rescale anyway
        int next_bits = next_data->total_coeff_modulus_bit_count();
        int safe_bits = next_bits - 2;
        if (safe_bits < 10) safe_bits = 10;
        return pow(2.0, safe_bits);
    }

    // ─── multiply_vector_reduced_error ───
    // ct × vector<complex<double>>
    // Strategy: use safe_encode_scale for the plaintext, and temporarily set
    // ct.scale to a safe value before multiply_plain to avoid "scale out of bounds".
    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const std::vector<std::complex<double>> &coeff_vec,
        Ciphertext &destination)
    {
        if (coeff_vec.empty()) { destination = encrypted; return; }
        double ct_scale = encrypted.scale();
        double enc_scale = safe_encode_scale(encrypted);
        Plaintext plain;
        encoder.encode(coeff_vec, encrypted.parms_id(), enc_scale, plain);
        // Temporarily lower ct.scale so SEAL's internal scale check passes
        Ciphertext ct_tmp = encrypted;
        ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain.scale());
        evaluator.multiply_plain(ct_tmp, plain, destination);
        destination.scale() = ct_scale; // emulate seal-modified: don't update scale
    }

    // ct × single Plaintext
    // Temporarily lower ct.scale before multiply_plain, restore after.
    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const Plaintext &plain,
        Ciphertext &destination)
    {
        double ct_scale = encrypted.scale();
        Ciphertext ct_tmp = encrypted;
        ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain.scale());
        evaluator.multiply_plain(ct_tmp, plain, destination);
        destination.scale() = ct_scale; // emulate seal-modified: don't update scale
    }

    // ct × vector<Plaintext>
    // Temporarily lower ct.scale before each multiply_plain, restore after.
    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const std::vector<Plaintext> &plain_vec,
        Ciphertext &destination)
    {
        if (plain_vec.empty()) return;
        double ct_scale = encrypted.scale();
        double safe_ct_scale = safe_scale_for_multiply_plain(encrypted, plain_vec[0].scale());
        Ciphertext ct_tmp = encrypted;
        ct_tmp.scale() = safe_ct_scale;
        Ciphertext tmp;
        evaluator.multiply_plain(ct_tmp, plain_vec[0], destination);
        destination.scale() = ct_scale; // emulate seal-modified: don't update scale
        for (size_t i = 1; i < plain_vec.size(); ++i) {
            ct_tmp.scale() = safe_scale_for_multiply_plain(encrypted, plain_vec[i].scale());
            evaluator.multiply_plain(ct_tmp, plain_vec[i], tmp);
            tmp.scale() = ct_scale; // emulate seal-modified: don't update scale
            evaluator.add_inplace(destination, tmp);
        }
    }

    // ─── add_inplace_reduced_error ───
    // In seal-modified, multiply doesn't change scale or level, so add always
    // operates on same-level ciphertexts. Standard SEAL's multiply (ct×ct/square)
    // doesn't change level either, but rescale does. After our scale-override
    // trick, different computation paths may end up at different levels.
    // Fix: mod-switch the higher-level ciphertext down before add/sub.
    void align_levels(Ciphertext &ct1, Ciphertext &ct2) {
        if (ct1.parms_id() == ct2.parms_id()) return;
        auto lv1 = context.get_context_data(ct1.parms_id())->chain_index();
        auto lv2 = context.get_context_data(ct2.parms_id())->chain_index();
        if (lv1 > lv2) {
            evaluator.mod_switch_to_inplace(ct1, ct2.parms_id());
        } else {
            evaluator.mod_switch_to_inplace(ct2, ct1.parms_id());
        }
    }

    void add_inplace_reduced_error(Ciphertext &ct1, const Ciphertext &ct2) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct2_copy = ct2;
            align_levels(ct1, ct2_copy);
        }
        evaluator.add_inplace(ct1, ct2);
    }
    void add_inplace_reduced_error(Ciphertext &ct1, Ciphertext &ct2) {
        align_levels(ct1, ct2);
        evaluator.add_inplace(ct1, ct2);
    }

    // ─── add_reduced_error (3-arg: ct1 + ct2 → result) ───
    void add_reduced_error(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &result) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct1_copy = ct1, ct2_copy = ct2;
            align_levels(ct1_copy, ct2_copy);
            evaluator.add(ct1_copy, ct2_copy, result);
        } else {
            evaluator.add(ct1, ct2, result);
        }
    }
    void add_reduced_error(const Ciphertext &ct1, Ciphertext &ct2, Ciphertext &result) {
        if (ct1.parms_id() != ct2.parms_id()) {
            Ciphertext ct1_copy = ct1, ct2_copy = ct2;
            align_levels(ct1_copy, ct2_copy);
            evaluator.add(ct1_copy, ct2_copy, result);
        } else {
            evaluator.add(ct1, ct2, result);
        }
    }

    // ─── multiply_const: multiply ciphertext by RR constant ───
    // Strategy: encode at safe_encode_scale, temporarily lower ct.scale before
    // multiply_plain, then override result.scale back to ct.scale.
    void multiply_const(Ciphertext &ct, const NTL::RR &constant, Ciphertext &dest) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), enc_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain(ct_tmp, plain_const, dest);
        dest.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void multiply_const(Ciphertext &ct, double constant, Ciphertext &dest) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), enc_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain(ct_tmp, plain_const, dest);
        dest.scale() = ct_scale; // emulate seal-modified: don't update scale
    }

    // ─── multiply_reduced_error ───
    // ct1 × ct2 with relinearization. Temporarily lower scales before multiply.
    void multiply_reduced_error(
        const Ciphertext &ct1, const Ciphertext &ct2,
        const RelinKeys &rk, Ciphertext &dest)
    {
        double ct1_scale = ct1.scale();
        double safe_scale = safe_scale_for_square(ct1);
        Ciphertext ct1_tmp = ct1, ct2_tmp = ct2;
        ct1_tmp.scale() = safe_scale;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply(ct1_tmp, ct2_tmp, dest);
        dest.scale() = ct1_scale; // emulate seal-modified: don't update scale
        evaluator.relinearize_inplace(dest, rk);
    }
    void multiply_reduced_error(
        Ciphertext &ct1, Ciphertext &ct2,
        const RelinKeys &rk, Ciphertext &dest)
    {
        double ct1_scale = ct1.scale();
        double safe_scale = safe_scale_for_square(ct1);
        Ciphertext ct1_tmp = ct1, ct2_tmp = ct2;
        ct1_tmp.scale() = safe_scale;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply(ct1_tmp, ct2_tmp, dest);
        dest.scale() = ct1_scale; // emulate seal-modified: don't update scale
        evaluator.relinearize_inplace(dest, rk);
    }

    // ─── sub_inplace_reduced_error ───
    void sub_inplace_reduced_error(Ciphertext &ct1, const Ciphertext &ct2) {
        evaluator.sub_inplace(ct1, ct2);
    }

    // ─── Forward common Evaluator methods ───
    void rotate_vector(const Ciphertext &ct, int steps, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_vector(ct, steps, gk, dest);
    }
    void complex_conjugate(const Ciphertext &ct, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.complex_conjugate(ct, gk, dest);
    }
    void multiply_plain(const Ciphertext &ct, const Plaintext &pt, Ciphertext &dest) {
        double ct_scale = ct.scale();
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, pt.scale());
        evaluator.multiply_plain(ct_tmp, pt, dest);
        dest.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void add_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        if (ct1.parms_id() != ct2.parms_id()) {
            std::cerr << "[DIAG] fwd_add_inplace MISMATCH: lv1=" << ct1.coeff_modulus_size()
                 << " s1=" << log2(ct1.scale()) << " lv2=" << ct2.coeff_modulus_size()
                 << " s2=" << log2(ct2.scale()) << std::endl;
        }
        evaluator.add_inplace(ct1, ct2);
    }
    void add(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &dest) {
        if (ct1.parms_id() != ct2.parms_id()) {
            std::cerr << "[DIAG] fwd_add MISMATCH: lv1=" << ct1.coeff_modulus_size()
                 << " s1=" << log2(ct1.scale()) << " lv2=" << ct2.coeff_modulus_size()
                 << " s2=" << log2(ct2.scale()) << std::endl;
        }
        evaluator.add(ct1, ct2, dest);
    }
    void sub_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        if (ct1.parms_id() != ct2.parms_id()) {
            std::cerr << "[DIAG] fwd_sub_inplace MISMATCH: lv1=" << ct1.coeff_modulus_size()
                 << " s1=" << log2(ct1.scale()) << " lv2=" << ct2.coeff_modulus_size()
                 << " s2=" << log2(ct2.scale()) << std::endl;
        }
        evaluator.sub_inplace(ct1, ct2);
    }
    void rescale_to_next_inplace(Ciphertext &ct) {
        // In seal-modified lazy scaling, rescale does divide_and_round_q_last
        // (same as standard SEAL rescale) but does NOT update the scale field.
        // Standard SEAL rescale checks: log2(ct.scale) < next_level.total_bits
        // and sets: ct.scale /= q_last.
        // We temporarily lower ct.scale to pass the bounds check, call standard
        // rescale (which does the mathematically correct divide_and_round), then
        // restore the original scale (lazy scaling: don't update scale field).
        double ct_scale = ct.scale();
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return;
        auto next_data = ctx_data->next_context_data();
        if (!next_data) return; // at last level, can't rescale
        int next_bits = next_data->total_coeff_modulus_bit_count();
        // Set a safe scale: log2(safe_scale) < next_bits
        int safe_bits = next_bits - 2;
        if (safe_bits < 10) safe_bits = 10;
        ct.scale() = pow(2.0, safe_bits);
        evaluator.rescale_to_next_inplace(ct);
        ct.scale() = ct_scale; // lazy scaling: don't update scale field
    }
    void relinearize_inplace(Ciphertext &ct, const RelinKeys &rk) {
        evaluator.relinearize_inplace(ct, rk);
    }
    void mod_switch_to_next_inplace(Ciphertext &ct) {
        evaluator.mod_switch_to_next_inplace(ct);
    }
    void mod_switch_to_inplace(Ciphertext &ct, parms_id_type pid) {
        evaluator.mod_switch_to_inplace(ct, pid);
    }

    // Plaintext variant: standard SEAL 4.1.1 supports mod_switch_to_inplace on Plaintext
    void mod_switch_to_inplace(Plaintext &pt, parms_id_type pid) {
        evaluator.mod_switch_to_inplace(pt, pid);
    }
    void rotate_rows(Ciphertext &ct, int steps, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_rows(ct, steps, gk, dest);
    }
    void rotate_columns(Ciphertext &ct, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_columns(ct, gk, dest);
    }
    // square: ct² — temporarily lower scale before square to avoid bounds check.
    void square(const Ciphertext &ct, Ciphertext &dest) {
        double ct_scale = ct.scale();
        double safe_scale = safe_scale_for_square(ct);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale;
        evaluator.square(ct_tmp, dest);
        dest.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void square_inplace(Ciphertext &ct) {
        double ct_scale = ct.scale();
        double safe_scale = safe_scale_for_square(ct);
        ct.scale() = safe_scale;
        evaluator.square_inplace(ct);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    // multiply_plain_inplace: ct×plain — temporarily lower ct.scale before multiply.
    void multiply_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        double ct_scale = ct.scale();
        ct.scale() = safe_scale_for_multiply_plain(ct, pt.scale());
        evaluator.multiply_plain_inplace(ct, pt);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void add_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.add_plain_inplace(ct, pt);
    }
    void sub_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.sub_plain_inplace(ct, pt);
    }
    // multiply ct×ct: temporarily lower scales before multiply to avoid bounds check.
    void multiply(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &dest) {
        double ct1_scale = ct1.scale();
        double safe_scale = safe_scale_for_square(ct1);
        Ciphertext ct1_tmp = ct1, ct2_tmp = ct2;
        ct1_tmp.scale() = safe_scale;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply(ct1_tmp, ct2_tmp, dest);
        dest.scale() = ct1_scale; // emulate seal-modified: don't update scale
    }
    void multiply_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        double ct1_scale = ct1.scale();
        double safe_scale = safe_scale_for_square(ct1);
        ct1.scale() = safe_scale;
        Ciphertext ct2_tmp = ct2;
        ct2_tmp.scale() = safe_scale;
        evaluator.multiply_inplace(ct1, ct2_tmp);
        ct1.scale() = ct1_scale; // emulate seal-modified: don't update scale
    }
    void transform_to_ntt_inplace(Ciphertext &ct) {
        evaluator.transform_to_ntt_inplace(ct);
    }
    void transform_from_ntt_inplace(Ciphertext &ct) {
        evaluator.transform_from_ntt_inplace(ct);
    }
    void negate_inplace(Ciphertext &ct) {
        evaluator.negate_inplace(ct);
    }

    // ─── Methods from seal-modified that need standard SEAL equivalents ───

    // double_inplace: multiply ciphertext by 2 — temporarily lower scale before multiply.
    void double_inplace(Ciphertext &ct) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_two;
        encoder.encode(2.0, ct.parms_id(), enc_scale, plain_two);
        ct.scale() = safe_scale_for_multiply_plain(ct, plain_two.scale());
        evaluator.multiply_plain_inplace(ct, plain_two);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }

    // add_const: add a constant to ciphertext (double and NTL::RR overloads)
    void add_const(Ciphertext &ct, double constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain(ct, plain_const, dest);
    }
    void add_const(Ciphertext &ct, const NTL::RR &constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain(ct, plain_const, dest);
    }
    void add_const(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain_inplace(ct, plain_const);
    }

    // add_const_inplace: same as add_const above
    void add_const_inplace(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain_inplace(ct, plain_const);
    }
    void add_const_inplace(Ciphertext &ct, const NTL::RR &constant) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), ct.scale(), plain_const);
        evaluator.add_plain_inplace(ct, plain_const);
    }

    // multiply_const_inplace: multiply ciphertext by a constant (double and NTL::RR overloads)
    // In seal-modified, this uses a fixed internal scale and scale field doesn't change.
    // Strategy: encode at ct.scale(), multiply, then override result.scale back to ct.scale.
    // Using ct.scale() avoids "scale out of bounds" when ct.scale is large (e.g. 2^60):
    //   multiply_plain checks ct.scale * plain.scale = ct.scale * ct.scale must fit
    //   in the current level's coeff_modulus, which is more controlled than
    //   ct.scale * default_scale (=2^10) that could exceed the modulus at lower levels.
    // Helper: find the largest valid encoding scale for a given parms_id.
    // SEAL requires: log2(encode_scale) < total_coeff_modulus_bit_count at that level.
    // Also, after multiply_plain, ct.scale * encode_scale must be representable,
    // so we cap encode_scale to ensure ct.scale * encode_scale doesn't overflow.
    double safe_encode_scale(const Ciphertext &ct) const {
        auto ctx_data = context.get_context_data(ct.parms_id());
        if (!ctx_data) return default_scale;
        int total_bits = ctx_data->total_coeff_modulus_bit_count();
        double ct_scale = ct.scale();
        int ct_scale_bits = (ct_scale > 0) ? (int)ceil(log2(ct_scale)) : 0;
        // We need: log2(encode_scale) < total_bits AND ct_scale_bits + log2(encode_scale) < total_bits
        int max_encode_bits = total_bits - ct_scale_bits - 1; // -1 for safety margin
        if (max_encode_bits < 10) max_encode_bits = 10; // at least 2^10
        double max_scale = pow(2.0, max_encode_bits);
        return std::min(max_scale, ct_scale);
    }

    void multiply_const_inplace(Ciphertext &ct, double constant) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), enc_scale, plain_const);
        ct.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct, plain_const);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void multiply_const_inplace(Ciphertext &ct, const NTL::RR &constant) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), enc_scale, plain_const);
        ct.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct, plain_const);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void multiply_const_inplace(Ciphertext &ct, double constant, Ciphertext &dest) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(constant, ct.parms_id(), enc_scale, plain_const);
        Ciphertext ct_tmp = ct;
        ct_tmp.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain(ct_tmp, plain_const, dest);
        dest.scale() = ct_scale; // emulate seal-modified: don't update scale
    }

    // multiply_const_inplace_fit_to_scale: multiply by RR constant, then adjust scale
    // Strategy: use safe_encode_scale and safe_scale_for_multiply_plain.
    void multiply_const_inplace_fit_to_scale(Ciphertext &ct, const NTL::RR &constant, double /*scale_to_align*/) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), enc_scale, plain_const);
        ct.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct, plain_const);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }
    void multiply_const_inplace_fit_to_scale(Ciphertext &ct, const NTL::RR &constant, const NTL::RR &/*scale_to_align*/) {
        double ct_scale = ct.scale();
        double enc_scale = safe_encode_scale(ct);
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.parms_id(), enc_scale, plain_const);
        ct.scale() = safe_scale_for_multiply_plain(ct, plain_const.scale());
        evaluator.multiply_plain_inplace(ct, plain_const);
        ct.scale() = ct_scale; // emulate seal-modified: don't update scale
    }

    // ─── multiply_inplace_reduced_error ───
    // ct1 *= ct2 with relinearization. Temporarily lower scales before multiply.
    void multiply_inplace_reduced_error(Ciphertext &ct1, Ciphertext &ct2, const RelinKeys &rk) {
        double ct1_scale = ct1.scale();
        double safe_scale = safe_scale_for_square(ct1);
        ct1.scale() = safe_scale;
        ct2.scale() = safe_scale;
        evaluator.multiply_inplace(ct1, ct2);
        ct1.scale() = ct1_scale; // emulate seal-modified: don't update scale
        evaluator.relinearize_inplace(ct1, rk);
    }

    // ─── sub_reduced_error (3-arg: ct1 - ct2 → result) ───
    void sub_reduced_error(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &result) {
        evaluator.sub(ct1, ct2, result);
    }
    void sub_reduced_error(Ciphertext &ct1, Ciphertext &ct2, Ciphertext &result) {
        evaluator.sub(ct1, ct2, result);
    }
};
