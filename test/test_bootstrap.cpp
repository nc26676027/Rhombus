/**
 * test_bootstrap.cpp
 * Bootstrapping verification test - concise & visual
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

    // 1. CKKS params
    section("1. CKKS Parameters");
    size_t poly_deg = 32768;
    vector<int> bit_sizes = {60, 45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45};
    int L = (int)bit_sizes.size();
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_deg);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_deg, bit_sizes));
    SEALContext context(parms, true, sec_level_type::tc128);
    if (!context.parameters_set()) {
        cout << R "  ERROR: SEAL parameters invalid!" N << endl;
        return 1;
    }
    auto ctx_data = context.first_context_data();
    cout << "  poly_modulus_degree : " << poly_deg << endl;
    cout << "  coeff_modulus_bits  : " << ctx_data->total_coeff_modulus_bit_count() << endl;
    cout << "  coeff_modulus_count : " << ctx_data->parms().coeff_modulus().size() << endl;

    // 2. Keys
    section("2. Key Generation");
    KeyGenerator keygen(context);
    SecretKey sk = keygen.secret_key();
    PublicKey pk; keygen.create_public_key(pk);
    RelinKeys rlk; keygen.create_relin_keys(rlk);
    GaloisKeys gk; // will be populated by Bootstrapper
    Encryptor encryptor(context, pk);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    CKKSEncoder encoder(context);
    cout << "  Keys ready" << endl;

    // 3. Build Bootstrapper (this generates its own GaloisKeys internally)
    section("3. Build Bootstrapper");
    long loge = 10, logn = 14, logNh = 14;
    double final_scale = pow(2.0, 40);
    long boundary_K = 25, sin_cos_deg = 59, double_angle_iter = 3, inverse_deg = 7;

    auto t0 = chrono::high_resolution_clock::now();
    cout << "  Constructing Bootstrapper..." << endl;
    Bootstrapper booter(loge, logn, logNh, L, final_scale,
        boundary_K, sin_cos_deg, double_angle_iter, inverse_deg,
        context, keygen, encoder, encryptor, decryptor, evaluator, rlk, gk);

    cout << "  Adding bootstrap keys (GaloisKeys)..." << endl;
    vector<int> extra_keys; // no extra keys needed
    booter.addBootKeys_3_other_keys(gk, extra_keys);

    cout << "  Preparing mod polynomial..." << endl;
    booter.prepare_mod_polynomial();

    cout << "  Generating LT coefficients..." << endl;
    booter.generate_LT_coefficient_3();

    auto t1 = chrono::high_resolution_clock::now();
    double setup_ms = chrono::duration<double,milli>(t1-t0).count();
    cout << "  Bootstrapper ready (" << (int)setup_ms << " ms)" << endl;

    // 4. Encrypt & consume levels
    section("4. Encrypt & Consume Levels");
    size_t slots = poly_deg / 2;
    vector<double> input(slots);
    for (size_t i = 0; i < slots; i++) input[i] = 0.01 * (i % 100);

    double scale = pow(2.0, 40);
    Plaintext pt; encoder.encode(input, scale, pt);
    Ciphertext ct; encryptor.encrypt(pt, ct);
    cout << "  Initial level  : " << ct.coeff_modulus_size() << endl;
    cout << "  Input (first 8): "; show_vec(input); cout << endl;

    // Use mod-switch to drop levels (bootstrap requires lowest level)
    while (ct.coeff_modulus_size() > 1) {
        evaluator.mod_switch_to_next_inplace(ct);
    }
    cout << "  After mod-switch: " << ct.coeff_modulus_size() << " level(s)" << endl;

    // 5. Bootstrap
    section("5. Bootstrap");
    Ciphertext refreshed;
    cout << "  Running bootstrap..." << endl;
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
    cout << "  Refreshed level: " << refreshed.coeff_modulus_size() << endl;

    // 6. Decrypt & verify
    section("6. Accuracy Check");
    Plaintext pt_out; vector<double> output;
    decryptor.decrypt(refreshed, pt_out);
    encoder.decode(pt_out, output);

    double err = max_err(input, output);
    int prec_bits = (err > 0) ? (int)round(-log2(err)) : 999;
    cout << "  Output (first 8): "; show_vec(output); cout << endl;
    cout << "  Max error       : " << scientific << setprecision(3) << err << endl;
    cout << "  Precision       : ~" << prec_bits << " bits" << endl;

    // 7. Summary
    section("Summary");
    cout << "  +----------------------------+" << endl;
    cout << "  | Setup    : " << setw(8) << (int)setup_ms << " ms |" << endl;
    cout << "  | Bootstrap: " << setw(8) << (int)boot_ms << " ms |" << endl;
    cout << "  | Gained   : +" << (refreshed.coeff_modulus_size() - ct.coeff_modulus_size()) << " levels  |" << endl;
    cout << "  | Precision: ~" << prec_bits << " bits   |" << endl;
    cout << "  +----------------------------+" << endl;

    bool pass = (err < 1.0) && (refreshed.coeff_modulus_size() > ct.coeff_modulus_size());
    cout << endl << (pass ? G "  PASS" : R "  FAIL") << N << endl << endl;
    return pass ? 0 : 1;
}
