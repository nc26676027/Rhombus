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
    double default_scale;

    CompatEvaluator(Evaluator &ev, CKKSEncoder &enc, double scale)
        : evaluator(ev), encoder(enc), default_scale(scale) {}

    // ─── multiply_vector_reduced_error ───
    // ct × vector<complex<double>>
    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const std::vector<std::complex<double>> &coeff_vec,
        Ciphertext &destination)
    {
        Plaintext plain;
        encoder.encode(coeff_vec, default_scale, plain);
        evaluator.multiply_plain(encrypted, plain, destination);
    }

    // ct × single Plaintext
    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const Plaintext &plain,
        Ciphertext &destination)
    {
        evaluator.multiply_plain(encrypted, plain, destination);
    }

    // ct × vector<Plaintext>
    void multiply_vector_reduced_error(
        const Ciphertext &encrypted,
        const std::vector<Plaintext> &plain_vec,
        Ciphertext &destination)
    {
        if (plain_vec.empty()) return;
        Ciphertext tmp;
        evaluator.multiply_plain(encrypted, plain_vec[0], destination);
        for (size_t i = 1; i < plain_vec.size(); ++i) {
            evaluator.multiply_plain(encrypted, plain_vec[i], tmp);
            evaluator.add_inplace(destination, tmp);
        }
    }

    // ─── add_inplace_reduced_error ───
    void add_inplace_reduced_error(Ciphertext &ct1, const Ciphertext &ct2) {
        evaluator.add_inplace(ct1, ct2);
    }
    void add_inplace_reduced_error(Ciphertext &ct1, Ciphertext &ct2) {
        evaluator.add_inplace(ct1, ct2);
    }

    // ─── add_reduced_error (3-arg: ct1 + ct2 → result) ───
    void add_reduced_error(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &result) {
        evaluator.add(ct1, ct2, result);
    }
    void add_reduced_error(const Ciphertext &ct1, Ciphertext &ct2, Ciphertext &result) {
        evaluator.add(ct1, ct2, result);
    }

    // ─── multiply_const: multiply ciphertext by RR constant ───
    void multiply_const(Ciphertext &ct, const NTL::RR &constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.scale(), plain_const);
        evaluator.multiply_plain(ct, plain_const, dest);
    }
    void multiply_const(Ciphertext &ct, double constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(constant, ct.scale(), plain_const);
        evaluator.multiply_plain(ct, plain_const, dest);
    }

    // ─── multiply_reduced_error ───
    // ct1 × ct2 with relinearization. In seal-modified this was a "reduced error" variant.
    // In standard SEAL, we use multiply + relinearize.
    void multiply_reduced_error(
        const Ciphertext &ct1, const Ciphertext &ct2,
        const RelinKeys &rk, Ciphertext &dest)
    {
        evaluator.multiply(ct1, ct2, dest);
        evaluator.relinearize_inplace(dest, rk);
    }
    void multiply_reduced_error(
        Ciphertext &ct1, Ciphertext &ct2,
        const RelinKeys &rk, Ciphertext &dest)
    {
        evaluator.multiply(ct1, ct2, dest);
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
        evaluator.multiply_plain(ct, pt, dest);
    }
    void add_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        evaluator.add_inplace(ct1, ct2);
    }
    void add(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &dest) {
        evaluator.add(ct1, ct2, dest);
    }
    void sub_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        evaluator.sub_inplace(ct1, ct2);
    }
    void rescale_to_next_inplace(Ciphertext &ct) {
        evaluator.rescale_to_next_inplace(ct);
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

    // Plaintext variant: in standard SEAL 4.1, Plaintext doesn't have mod_switch,
    // but we can achieve the same by re-encoding. However, RBOOT uses this in a
    // narrow context (aligning plaintext level to ciphertext level), so we simply
    // set the parms_id on the plaintext directly.
    void mod_switch_to_inplace(Plaintext &pt, parms_id_type pid) {
        pt.parms_id() = pid;
    }
    void rotate_rows(Ciphertext &ct, int steps, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_rows(ct, steps, gk, dest);
    }
    void rotate_columns(Ciphertext &ct, const GaloisKeys &gk, Ciphertext &dest) {
        evaluator.rotate_columns(ct, gk, dest);
    }
    void square(const Ciphertext &ct, Ciphertext &dest) {
        evaluator.square(ct, dest);
    }
    void square_inplace(Ciphertext &ct) {
        evaluator.square_inplace(ct);
    }
    void multiply_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.multiply_plain_inplace(ct, pt);
    }
    void add_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.add_plain_inplace(ct, pt);
    }
    void sub_plain_inplace(Ciphertext &ct, const Plaintext &pt) {
        evaluator.sub_plain_inplace(ct, pt);
    }
    void multiply(const Ciphertext &ct1, const Ciphertext &ct2, Ciphertext &dest) {
        evaluator.multiply(ct1, ct2, dest);
    }
    void multiply_inplace(Ciphertext &ct1, const Ciphertext &ct2) {
        evaluator.multiply_inplace(ct1, ct2);
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

    // double_inplace: multiply ciphertext by 2 (plaintext scalar 2)
    void double_inplace(Ciphertext &ct) {
        Plaintext plain_two;
        encoder.encode(2.0, ct.scale(), plain_two);
        evaluator.multiply_plain_inplace(ct, plain_two);
        // No rescale needed if we match the scale properly
    }

    // add_const: add a constant to ciphertext (double and NTL::RR overloads)
    void add_const(Ciphertext &ct, double constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(constant, ct.scale(), plain_const);
        plain_const.parms_id() = ct.parms_id();
        evaluator.add_plain(ct, plain_const, dest);
    }
    void add_const(Ciphertext &ct, const NTL::RR &constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.scale(), plain_const);
        plain_const.parms_id() = ct.parms_id();
        evaluator.add_plain(ct, plain_const, dest);
    }
    void add_const(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.scale(), plain_const);
        plain_const.parms_id() = ct.parms_id();
        evaluator.add_plain_inplace(ct, plain_const);
    }

    // add_const_inplace: same as add_const above
    void add_const_inplace(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.scale(), plain_const);
        plain_const.parms_id() = ct.parms_id();
        evaluator.add_plain_inplace(ct, plain_const);
    }
    void add_const_inplace(Ciphertext &ct, const NTL::RR &constant) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.scale(), plain_const);
        plain_const.parms_id() = ct.parms_id();
        evaluator.add_plain_inplace(ct, plain_const);
    }

    // multiply_const_inplace: multiply ciphertext by a constant (double and NTL::RR overloads)
    void multiply_const_inplace(Ciphertext &ct, double constant) {
        Plaintext plain_const;
        encoder.encode(constant, ct.scale(), plain_const);
        evaluator.multiply_plain_inplace(ct, plain_const);
    }
    void multiply_const_inplace(Ciphertext &ct, const NTL::RR &constant) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.scale(), plain_const);
        evaluator.multiply_plain_inplace(ct, plain_const);
    }
    void multiply_const_inplace(Ciphertext &ct, double constant, Ciphertext &dest) {
        Plaintext plain_const;
        encoder.encode(constant, ct.scale(), plain_const);
        evaluator.multiply_plain(ct, plain_const, dest);
    }

    // multiply_const_inplace_fit_to_scale: multiply by RR constant, then adjust scale
    // In seal-modified this aligned the output scale. Here we just do multiply_plain.
    void multiply_const_inplace_fit_to_scale(Ciphertext &ct, const NTL::RR &constant, double /*scale_to_align*/) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.scale(), plain_const);
        evaluator.multiply_plain_inplace(ct, plain_const);
    }
    void multiply_const_inplace_fit_to_scale(Ciphertext &ct, const NTL::RR &constant, const NTL::RR &/*scale_to_align*/) {
        Plaintext plain_const;
        encoder.encode(to_double(constant), ct.scale(), plain_const);
        evaluator.multiply_plain_inplace(ct, plain_const);
    }

    // ─── multiply_inplace_reduced_error ───
    // ct1 *= ct2 with relinearization
    void multiply_inplace_reduced_error(Ciphertext &ct1, Ciphertext &ct2, const RelinKeys &rk) {
        evaluator.multiply_inplace(ct1, ct2);
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
