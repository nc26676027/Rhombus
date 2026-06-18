#pragma once

/**
 * BootstrapHelper.h
 *
 * A thin wrapper around the RBOOT Bootstrapper, providing a clean interface
 * for use within the Rhombus project.
 *
 * Typical usage:
 *
 *   // 1. Build SEAL context + keys as usual (CKKS scheme)
 *   RhombusBootstrapParams params;
 *   params.loge      = 10;      // log2(scaling factor exponent)
 *   params.logn      = 13;      // log2(number of slots used in bootstrapping)
 *   params.logNh     = 13;      // log2(N/2), where N = poly_modulus_degree
 *   params.L         = 14;      // total number of moduli (levels)
 *   params.final_scale       = pow(2.0, 40);
 *   params.boundary_K        = 25;
 *   params.sin_cos_deg       = 59;
 *   params.double_angle_iter = 3;
 *   params.inverse_deg       = 7;
 *
 *   // 2. Collect the rotation steps required and add them to your GaloisKeys
 *   vector<int> boot_rot_steps;
 *   RhombusBootstrapHelper::collectRotationSteps(params, boot_rot_steps);
 *   // ... keygen.create_galois_keys(boot_rot_steps, gal_keys) ...
 *
 *   // 3. Construct the helper (expensive: precomputes FFT coefficients)
 *   RhombusBootstrapHelper helper(params,
 *       context, keygen, encoder, encryptor, decryptor, evaluator,
 *       relin_keys, gal_keys);
 *
 *   // 4. Bootstrap a ciphertext
 *   Ciphertext refreshed;
 *   helper.bootstrap(ct, refreshed);
 *   // or in-place:
 *   helper.bootstrap_inplace(ct);
 */

#include <vector>
#include <memory>
#include "seal/seal.h"
#include "Bootstrapper.h"

namespace rhombus {

/// Parameters that control the bootstrapping circuit.
struct RhombusBootstrapParams {
    long loge               = 10;    ///< log2(q/Delta): exponent of the encoding scale
    long logn               = 13;    ///< log2(slot count used inside bootstrapping)
    long logNh              = 13;    ///< log2(N/2)
    long L                  = 14;    ///< total number of coefficient moduli in the chain
    double final_scale      = 1LL << 40; ///< target scale after bootstrapping
    long boundary_K         = 25;    ///< range boundary K for modular reduction polynomial
    long sin_cos_deg        = 59;    ///< degree of the sin/cos polynomial approximation
    long double_angle_iter  = 3;     ///< number of double-angle iterations
    long inverse_deg        = 7;     ///< degree of the arcsine inverse polynomial
};

/**
 * RhombusBootstrapHelper
 *
 * Wraps the RBOOT Bootstrapper class, managing lifetime and providing
 * a minimal API surface.  The constructor is expensive (FFT coefficient
 * precomputation); construct once and reuse for many bootstrapping calls.
 */
class RhombusBootstrapHelper {
public:
    /**
     * Collect the rotation steps that must be present in the GaloisKeys
     * before constructing RhombusBootstrapHelper.
     *
     * Call this before key generation, merge the returned steps with any
     * other rotation steps your application needs, then create the keys.
     *
     * @param params        Bootstrapping parameters (same values used later).
     * @param steps_out     Output: rotation step integers to include in GaloisKeys.
     */
    static void collectRotationSteps(const RhombusBootstrapParams &params,
                                     std::vector<int> &steps_out,
                                     seal::SEALContext &context,
                                     seal::KeyGenerator &keygen,
                                     seal::CKKSEncoder &encoder,
                                     seal::Encryptor &encryptor,
                                     seal::Decryptor &decryptor,
                                     seal::Evaluator &evaluator,
                                     seal::RelinKeys &relin_keys,
                                     seal::GaloisKeys &gal_keys)
    {
        // Temporarily build a Bootstrapper to query which keys it needs.
        Bootstrapper tmp(
            params.loge, params.logn, params.logNh, params.L, params.final_scale,
            params.boundary_K, params.sin_cos_deg, params.double_angle_iter, params.inverse_deg,
            context, keygen, encoder, encryptor, decryptor, evaluator, relin_keys, gal_keys);
        tmp.addBootKeys_3_other_keys(gal_keys, steps_out);
    }

    /**
     * Construct the helper and precompute all FFT/polynomial coefficients.
     * This call is expensive; do it once during setup, not per-inference.
     */
    RhombusBootstrapHelper(const RhombusBootstrapParams &params,
                           seal::SEALContext &context,
                           seal::KeyGenerator &keygen,
                           seal::CKKSEncoder &encoder,
                           seal::Encryptor &encryptor,
                           seal::Decryptor &decryptor,
                           seal::Evaluator &evaluator,
                           seal::RelinKeys &relin_keys,
                           seal::GaloisKeys &gal_keys)
        : bootstrapper_(std::make_unique<Bootstrapper>(
              params.loge, params.logn, params.logNh, params.L, params.final_scale,
              params.boundary_K, params.sin_cos_deg, params.double_angle_iter, params.inverse_deg,
              context, keygen, encoder, encryptor, decryptor, evaluator, relin_keys, gal_keys))
    {
        bootstrapper_->prepare_mod_polynomial();
        bootstrapper_->generate_LT_coefficient_3();
    }

    /**
     * Bootstrap a ciphertext: refresh its modulus level.
     *
     * @param input     Ciphertext whose level has been consumed.
     * @param output    Freshly bootstrapped ciphertext at high level.
     */
    void bootstrap(const seal::Ciphertext &input, seal::Ciphertext &output) const {
        bootstrapper_->bootstrap_3(output, const_cast<seal::Ciphertext &>(input));
    }

    /**
     * Bootstrap in-place.
     */
    void bootstrap_inplace(seal::Ciphertext &ct) const {
        bootstrapper_->bootstrap_inplace_3(ct);
    }

    /**
     * "Slim" bootstrap variant (lower noise, fewer levels consumed).
     */
    void slim_bootstrap(const seal::Ciphertext &input, seal::Ciphertext &output) const {
        bootstrapper_->slim_bootstrap(output, const_cast<seal::Ciphertext &>(input));
    }

    void slim_bootstrap_inplace(seal::Ciphertext &ct) const {
        bootstrapper_->slim_bootstrap_inplace(ct);
    }

    /// Direct access to the underlying Bootstrapper for advanced use.
    Bootstrapper &raw() { return *bootstrapper_; }
    const Bootstrapper &raw() const { return *bootstrapper_; }

private:
    std::unique_ptr<Bootstrapper> bootstrapper_;
};

} // namespace rhombus
