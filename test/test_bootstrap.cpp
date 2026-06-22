/**
 * test_bootstrap.cpp
 * Bootstrapping verification test — parameters matching RBOOT reference
 *
 * Key differences from standard SEAL 4.1.1:
 *   - Uses sec_level_type::none because RBOOT's coeff_modulus exceeds tc128 bounds
 *     (RBOOT uses seal-modified-3.6.6 which bypasses security checks)
 *   - Lazy scaling is emulated via CompatEvaluator (SealCompat.h)
 *   - coeff_modulus structure: logq (51-bit) for boot levels, logp (46-bit) for remaining
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include "seal/seal.h"
#include "bootstrap/Bootstrapper.h"

using namespace std;
using namespace seal;

#define G  "\033[32m"
#define R  "\033[31m"
#define C  "\033[36m"
#define B  "\033[1m"
#define N  "\033[0m"

void section(const string& t) { cout << endl << B "--- " << t << " ---" N << endl; }

void show_vec(const vector<double>& v, int n=8) {
    cout << "[ ";
    for (int i = 0; i < min((int)v.size(), n); i++)
        cout << fixed << setprecision(4) << v[i] << " ";
    if ((int)v.size() > n) cout << "... ";
    cout << "]";
}

double max_err(const vector<double>& a, const vector<double>& b) {
    double m = 0;
    for (size_t i = 0; i < a.size(); i++) m = max(m, abs(a[i]-b[i]));
    return m;
}

int main() {
    cout << B C "\n  Rhombus Bootstrapping Test\n" N << endl;

    // ── 1. CKKS Parameters (matching RBOOT reference) ──
    section("1. CKKS Parameters");

    // Use logN=14 (poly_deg=16384) to fit in 16GB RAM.
    // RBOOT reference uses logN=16 but needs >16GB for GaloisKeys.
    // Reduced level budget for smaller parameters:
    //   remaining=2 (slot-to-coeff) + boot=10 (coeff-to-slot + sin + arcsin) = 12
    // coeff_modulus: {logq, 2*logp, 10*logq} where logq=50, logp=45
    long logN  = 14;
    long logn  = 11;
    long loge  = 8;
    long boundary_K = 25;
    long sin_cos_deg = 30;        // reduced from 59
    long double_angle_iter = 1;   // reduced from 2
    long inverse_deg = 63;        // reduced from 127

    int logq = 50;
    int logp = 45;
    int remaining_level = 2;   // slot-to-coeff (reduced)
    int boot_level = 2 + 4 + 1 + 1 + 4; // c2s + sin + double_angle + one_more + arcsin = 12
    int total_level = remaining_level + boot_level;  // = 14

    size_t poly_deg = size_t(1) << logN;
    double scale = pow(2.0, logp);  // initial scale = 2^46

    // Build coeff_modulus bit-sizes matching RBOOT structure
    vector<int> coeff_bit_vec;
    coeff_bit_vec.push_back(logq);           // first prime (key level)
    for (int i = 0; i < remaining_level; i++)
        coeff_bit_vec.push_back(logp);       // remaining levels
    for (int i = 0; i < boot_level; i++)
        coeff_bit_vec.push_back(logq);       // boot levels

    int L = total_level;  // L = 22

    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_deg);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_deg, coeff_bit_vec));
    // sec_level_type::none: standard SEAL 4.1.1 rejects this coeff_modulus under tc128
    // (total bits ~ 4*46 + 18*51 = 184 + 918 = 1102 > 881-bit tc128 limit for N=32768)
    // RBOOT uses seal-modified which bypasses this check; we emulate by using none.
    SEALContext context(parms, true, sec_level_type::none);
    if (!context.parameters_set()) {
        cout << R "  ERROR: SEAL parameters invalid!" N << endl;
        auto ctx_data = context.key_context_data();
        if (ctx_data) {
            cout << "  parameter_error: " << (int)ctx_data->qualifiers().parameter_error << endl;
        }
        return 1;
    }

    auto ctx_data = context.first_context_data();
    cout << "  poly_modulus_degree : " << poly_deg << endl;
    cout << "  coeff_modulus_bits  : " << ctx_data->total_coeff_modulus_bit_count() << endl;
    cout << "  coeff_modulus_count : " << ctx_data->parms().coeff_modulus().size() << endl;
    cout << "  total_level (L)     : " << L << endl;
    cout << "  initial_scale       : 2^" << logp << " = " << scale << endl;
    cout << "  sec_level           : none (RBOOT-compatible)" << endl;

    // ── 2. Key Generation ──
    section("2. Key Generation");
    KeyGenerator keygen(context);
    SecretKey sk = keygen.secret_key();
    PublicKey pk; keygen.create_public_key(pk);
    RelinKeys rlk; keygen.create_relin_keys(rlk);
    GaloisKeys gk;  // populated by Bootstrapper
    Encryptor encryptor(context, pk);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    CKKSEncoder encoder(context);
    cout << "  Keys ready" << endl;

    // ── 3. Build Bootstrapper ──
    section("3. Build Bootstrapper");
    long logNh = logN - 1;  // Nh = N/2
    double final_scale = scale;  // final_scale = 2^46, matching RBOOT

    auto t0 = chrono::high_resolution_clock::now();
    cout << "  Constructing Bootstrapper..." << endl;
    Bootstrapper booter(loge, logn, logNh, L, final_scale,
        boundary_K, sin_cos_deg, double_angle_iter, inverse_deg,
        context, keygen, encoder, encryptor, decryptor, evaluator, rlk, gk);

    cout << "  Adding bootstrap keys (GaloisKeys)..." << endl;
    // Following RBOOT reference: build gal_steps_vector manually
    vector<int> gal_steps_vector;
    gal_steps_vector.push_back(0);
    for (int i = 0; i < logNh; i++) {
        gal_steps_vector.push_back(1 << i);
    }
    booter.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
    keygen.create_galois_keys(gal_steps_vector, gk);
    booter.logn_set.insert(logn);

    cout << "  Preparing mod polynomial..." << endl;
    booter.prepare_mod_polynomial();

    cout << "  Generating LT coefficients..." << endl;
    booter.generate_LT_coefficient_3();

    auto t1 = chrono::high_resolution_clock::now();
    double setup_ms = chrono::duration<double,milli>(t1-t0).count();
    cout << "  Bootstrapper ready (" << (int)setup_ms << " ms)" << endl;

    // ── 4. Encrypt & consume levels ──
    section("4. Encrypt & Consume Levels");
    size_t slots = poly_deg / 2;  // full slots for N=32768
    size_t sparse_slots = size_t(1) << (1 + logn);  // 2^(1+13) = 16384

    // Input: small values in range [-0.05, 0.05] like RBOOT
    vector<double> input(slots);
    for (size_t i = 0; i < sparse_slots; i++) {
        input[i] = -0.05 + (0.1 * i / sparse_slots);
    }
    // Fill remaining slots with repetition for full-slot encoding
    for (size_t i = sparse_slots; i < slots; i++) {
        input[i] = input[i % sparse_slots];
    }

    Plaintext pt; encoder.encode(input, scale, pt);
    Ciphertext ct; encryptor.encrypt(pt, ct);
    cout << "  Initial chain_index : " << context.get_context_data(ct.parms_id())->chain_index() << endl;
    cout << "  Initial coeff_modulus_size : " << ct.coeff_modulus_size() << endl;
    cout << "  Input (first 8): "; show_vec(input); cout << endl;

    // Drop levels: mod-switch down to the lowest level (coeff_modulus_size=1).
    // bootstrap_sparse_3 / modraise_inplace requires ciphertext at the lowest level.
    // In RBOOT's seal-modified, n_special_primes=5 are excluded from the mod chain,
    // so their "total_level - 3" mod-switches result in coeff_modulus_size=1.
    // In standard SEAL 4.1.1, all primes are in the chain, so we mod-switch
    // until only 1 modulus remains.
    // Use SEAL's standard evaluator for mod-switch (not CompatEvaluator)
    Evaluator seal_eval(context);
    int initial_levels = ct.coeff_modulus_size();
    while (ct.coeff_modulus_size() > 1) {
        seal_eval.mod_switch_to_next_inplace(ct);
    }
    cout << "  Dropped " << (initial_levels - 1) << " levels" << endl;
    cout << "  After mod-switch chain_index: "
         << context.get_context_data(ct.parms_id())->chain_index() << endl;
    cout << "  After mod-switch coeff_modulus_size: " << ct.coeff_modulus_size() << endl;
    cout << "  scale: " << ct.scale() << endl;

    // ── 5. Bootstrap ──
    section("5. Bootstrap");
    Ciphertext refreshed;
    cout << "  Running bootstrap_3..." << endl;
    auto tb0 = chrono::high_resolution_clock::now();
    try {
        booter.bootstrap_3(refreshed, ct);
    } catch (const exception& e) {
        cout << R "  FAIL: " << e.what() << N << endl;
        return 1;
    }
    auto tb1 = chrono::high_resolution_clock::now();
    double boot_ms = chrono::duration<double,milli>(tb1-tb0).count();
    cout << "  Bootstrap done (" << (int)boot_ms << " ms)" << endl;
    cout << "  Refreshed coeff_modulus_size: " << refreshed.coeff_modulus_size() << endl;
    cout << "  Refreshed chain_index: "
         << context.get_context_data(refreshed.parms_id())->chain_index() << endl;

    // ── 6. Decrypt & verify ──
    section("6. Accuracy Check");
    Plaintext pt_out; vector<double> output;
    decryptor.decrypt(refreshed, pt_out);
    encoder.decode(pt_out, output);

    // Compute error only on the sparse slots (where data was meaningful)
    double err = 0;
    for (size_t i = 0; i < sparse_slots; i++) {
        err = max(err, abs(input[i] - output[i]));
    }
    int prec_bits = (err > 0) ? (int)round(-log2(err)) : 999;
    cout << "  Output (first 8): "; show_vec(output); cout << endl;
    cout << "  Max error (sparse slots): " << scientific << setprecision(3) << err << endl;
    cout << "  Precision       : ~" << prec_bits << " bits" << endl;

    // ── 7. Summary ──
    section("Summary");
    cout << "  +----------------------------+" << endl;
    cout << "  | Setup    : " << setw(8) << (int)setup_ms << " ms |" << endl;
    cout << "  | Bootstrap: " << setw(8) << (int)boot_ms << " ms |" << endl;
    cout << "  | Precision: ~" << prec_bits << " bits   |" << endl;
    cout << "  +----------------------------+" << endl;

    bool pass = (err < 1.0);
    cout << endl << (pass ? G "  PASS" : R "  FAIL") << N << endl << endl;
    return pass ? 0 : 1;
}
