#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include "RemezCos.h"
#include "RemezArcsin.h"
#include "Polynomial.h"
#include "SealCompat.h"
// #include "ScaleInvEvaluator.h"
#include "PolyUpdate.h"

using namespace std;

class ModularReducer {
public:
	long boundary_K;
	double log_width;
	long deg;
	long num_double_formula;

	double inverse_log_width;
	long inverse_deg;

    double scale_inverse_coeff;

    SEALContext &context;
	CKKSEncoder &encoder;
    Encryptor &encryptor;
    CompatEvaluator &evaluator;
    RelinKeys &relin_keys;
    Decryptor &decryptor; 
	
	RemezParam rmparm;

	RemezSin* sin_generator;
	RemezSin* cos_generator;
	RemezArcsin* inverse_sin_generator;

	boot::Polynomial sin_polynomial;
	boot::Polynomial cos_polynomial;
	boot::Polynomial inverse_sin_polynomial_v1;
	boot::Polynomial inverse_sin_polynomial_v1_original;

	ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg,
    SEALContext &_context,
	CKKSEncoder &_encoder,
    Encryptor &_encryptor,
    CompatEvaluator &_evaluator,
    RelinKeys &_relin_keys,
    Decryptor &_decryptor
    );

	void double_angle_formula(Ciphertext &cipher);
	void double_angle_formula_scaled(Ciphertext &cipher, double scale_coeff);
	void generate_sin_cos_polynomial();
	void generate_inverse_sine_polynomial();
	void scaling_for_turn_back_q();
    void write_polynomials();
	void modular_reduction(Ciphertext &rtn, Ciphertext &cipher); 

	vector<RR> arcsin_decomp_coeff;
	vector<RR> arcsin_decomp_coeff_original;
	minicomp::Tree arcsin_tree;
	double abs_lift = 0.125;
};

static double scale_for_eval = 4.0 * 0.9; // tuning
extern double scale_for_boost_relu_range;
extern double dyn_scale_for_turn_back_q;
extern bool after_cnn_first_boot;
