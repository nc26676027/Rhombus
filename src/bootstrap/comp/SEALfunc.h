#pragma	once

#include "seal/seal.h"
#include "../SealCompat.h"
#include "../func.h"
#include "MinicompFunc.h"
#include "PolyUpdate.h"
#include <memory>

using namespace seal;

namespace seal{
    void MultipleAdd_SEAL(CompatEvaluator &evaluator, Ciphertext &cipher, Ciphertext& result, long long n);
	void test_evaluation(CompatEvaluator &evaluator, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, const Ciphertext& cipher_in, Ciphertext& cipher_out);
    void geneT0T1(Encryptor &encryptor, CompatEvaluator &evaluator, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& T0, Ciphertext& T1, Ciphertext& cipher);
    void evalT(CompatEvaluator &evaluator, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& Tmplusn, const Ciphertext& Tm, const Ciphertext& Tn, const Ciphertext& Tmminusn);
	void eval_polynomial_integrate(Encryptor &encryptor, CompatEvaluator &evaluator, Decryptor &decryptor, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& res, Ciphertext& cipher, long deg, const vector<RR> &decomp_coeff, Tree &tree);
    long coeff_number(long deg, Tree& tree);
    void coeff_change(long comp_no, long deg[], double* coeff[], long type[], vector<Tree> &tree);
	long ShowFailure_ReLU(Decryptor &decryptor, CKKSEncoder &encoder, Ciphertext& cipher, vector<double>& x, long precision, long n);
}   