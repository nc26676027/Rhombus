#pragma once

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<algorithm>
#include<NTL/RR.h>
#include"func.h"
#include"SealCompat.h"
// #include"ScaleInvEvaluator.h"

using namespace std;
using namespace NTL;
using namespace seal;

namespace boot
{
	class Polynomial {
	public:
		RR* coeff = 0;
		long deg, heap_k, heap_m, heaplen;
		RR* chebcoeff = 0;
		Polynomial** poly_heap = 0;

		Polynomial();
		Polynomial(long _deg);
		Polynomial(long _deg, RR* _coeff, string tag);

		~Polynomial();
		void set_polynomial(long _deg, RR* _coeff, string tag);
		void set_zero_polynomial(long _deg);
		void showcoeff();
		void showchebcoeff();
		void copy(Polynomial &poly);
		void power_to_cheb();
		void cheb_to_power();
		RR evaluate(RR value);

		void constmul(RR constant);
		void mulinplace(Polynomial &poly);
		void addinplace(Polynomial &poly);
		void subtinplace(Polynomial &poly);

		void change_variable_scale(RR scale); // f_new(x) = f(x / scale)
		void generate_poly_heap_manual(long k, long m);
		void generate_poly_heap();
		void generate_poly_heap_odd();

		void write_heap_to_file(ostream &out);
		void read_heap_from_file(istream &in);

		void homomorphic_poly_evaluation(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, CompatEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor);
		void homomorphic_poly_evaluation_naive(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, CompatEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor);
		void homomorphic_poly_evaluation_cheb_naive(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, CompatEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor);
	};

	void mul(Polynomial &rtn, Polynomial &a, Polynomial &b);
	void add(Polynomial &rtn, Polynomial &a, Polynomial &b);
	void subt(Polynomial &rtn, Polynomial &a, Polynomial &b);

	void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider);
	void chebyshev(Polynomial &rtn, long deg);
	void second_chebyshev_times_x_for_sine(Polynomial &rtn, long deg);
}