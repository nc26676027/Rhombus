#include "Polynomial.h"
#include "SealCompat.h"
#include <iomanip>

namespace boot
{
	Polynomial::Polynomial() {}

	Polynomial::Polynomial(long _deg) {
		deg = _deg;
		coeff = new RR[deg + 1];
		chebcoeff = new RR[deg + 1];
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = RR(0);
			chebcoeff[i] = RR(0);
		}
	}

	Polynomial::Polynomial(long _deg, RR* _coeff, string tag) {
		deg = _deg;
		coeff = new RR[deg + 1];
		chebcoeff = new RR[deg + 1];
		if (tag == "power") {
			for(int i = 0; i < deg + 1; i++) {
				coeff[i] = _coeff[i];
			}
			power_to_cheb();
		}

		else if (tag == "cheb") {
			for(int i = 0; i < deg + 1; i++) {
				chebcoeff[i] = _coeff[i];
			}
			cheb_to_power();
		}
	}

	Polynomial::~Polynomial() {
		if(coeff) delete[] coeff;
		if(chebcoeff) delete[] chebcoeff;
		if(poly_heap) {
			for (int i = 0; i < heaplen; i++) {
				if(poly_heap[i]) delete poly_heap[i];
			}
		}
	}

	void Polynomial::set_polynomial(long _deg, RR* _coeff, string tag) {
		deg = _deg;
		coeff = new RR[deg + 1];
		chebcoeff = new RR[deg + 1];
		if (tag == "power") {
			for(int i = 0; i < deg + 1; i++) {
				coeff[i] = _coeff[i];
			}
			power_to_cheb();
		}

		else if (tag == "cheb") {
			for(int i = 0; i < deg + 1; i++) {
				chebcoeff[i] = _coeff[i];
			}
			cheb_to_power();
		}
	}

	void Polynomial::set_zero_polynomial(long _deg) {
		deg = _deg;
		coeff = new RR[deg + 1]{};
		chebcoeff = new RR[deg + 1]{};
	}

	void Polynomial::showcoeff() {
		auto orig_output_precision = RR::OutputPrecision();
		RR::SetOutputPrecision(30);
		for(int i = 0; i < deg; i++) {
			cout << "R(" << coeff[i] << ")*x^" << i << " + ";
		}
		cout << "(" << coeff[deg] << ")*x^" << deg << endl;
		cout << endl;
		RR::SetOutputPrecision(orig_output_precision);
	}

	void Polynomial::showchebcoeff() {
		for(int i = 0; i < deg + 1; i++) {
			cout << "chebterm " << i << " : " << chebcoeff[i] << endl;
		}
		cout << endl;
	}
	void Polynomial::copy(Polynomial &poly) {
		deg = poly.deg;
		coeff = new RR[deg + 1];
		chebcoeff = new RR[deg + 1];
		for(int i = 0; i < deg + 1; i++) {
			coeff[i] = poly.coeff[i];
			chebcoeff[i] = poly.chebcoeff[i];
		}
	}

	void Polynomial::power_to_cheb() {
		Polynomial chebbasis, tmp(deg);
		for(int i = 0; i <= deg; i++) {
			tmp.coeff[i] = coeff[i];
		}

		for(int i = deg; i >= 0; i--) {
			chebyshev(chebbasis, i);
			chebcoeff[i] = tmp.coeff[i] / chebbasis.coeff[i];
			for(int j = 0; j <= i; j++) {
				chebbasis.coeff[j] *= chebcoeff[i];
			}
			tmp.subtinplace(chebbasis);
		}
	}

	void Polynomial::cheb_to_power() {
		Polynomial chebbasis, tmp(deg);
		for(int i = 0; i <= deg; i++) {
			chebyshev(chebbasis, i);
			for(int j = 0; j <= i; j++) {
				chebbasis.coeff[j] *= chebcoeff[i];
			}
			tmp.addinplace(chebbasis);
		}

		for(int i = 0; i <= deg; i++) {
			coeff[i] = tmp.coeff[i];
		}
	}

	RR Polynomial::evaluate(RR value) {
		RR rtn = RR(0), term = RR(1);
		for(int i = 0; i <= deg; i++) {
			rtn += coeff[i] * term;
			term *= value;
		}
		
		return rtn;
	}

	void Polynomial::constmul(RR constant) {
		for (int i = 0; i <= deg; i++) {
			coeff[i] *= constant;
			chebcoeff[i] *= constant;
		}
	}

	void Polynomial::mulinplace(Polynomial &poly) {
		Polynomial rtn;
		mul(rtn, (*this), poly);
		copy(rtn);	
	}

	void Polynomial::addinplace(Polynomial &poly) {
		Polynomial rtn;
		add(rtn, (*this), poly);
		copy(rtn);
	}

	void Polynomial::subtinplace(Polynomial &poly) {
		Polynomial rtn;
		subt(rtn, (*this), poly);
		copy(rtn);	
	}

	void Polynomial::change_variable_scale(RR scale) {
		RR div = RR(1);
		for (int i = 0; i <= deg; i++) {
			coeff[i] /= div;
			div *= scale;
		}
		power_to_cheb();
	}

	void Polynomial::generate_poly_heap_manual(long k, long m) {
		heap_k = k;
		heap_m = m;
		heaplen = (1 << (heap_m + 1)) - 1;
		poly_heap = new Polynomial*[heaplen];
		for(long i = 0; i < heaplen; i++) {
			poly_heap[i] = 0;
		}

		poly_heap[0] = new Polynomial();
		poly_heap[0]->copy((*this));

		long chebdeg = heap_k << heap_m;
		long first, last;
		Polynomial chebgiant;
		for(int i = 0; i < heap_m; i++) {
			chebdeg = chebdeg >> 1;
			first = (1 << i) - 1;
			last = (1 << (i + 1)) - 1;
			chebyshev(chebgiant, chebdeg);
			for(int j = first; j < last; j++) {
				if(poly_heap[j]) {
					if(poly_heap[j]->deg < chebdeg) {
						poly_heap[2 * (j + 1)] = new Polynomial();
						poly_heap[2 * (j + 1)]->copy(*poly_heap[j]);
					}

					else {
						poly_heap[2 * (j + 1) - 1] = new Polynomial();
						poly_heap[2 * (j + 1)] = new Polynomial();
						divide_poly(*poly_heap[2 * (j + 1) - 1], *poly_heap[2 * (j + 1)], *poly_heap[j], chebgiant);
						poly_heap[2 * (j + 1) - 1]->power_to_cheb();
						poly_heap[2 * (j + 1)]->power_to_cheb();
					}
				}
			}
		}
	}

	void Polynomial::generate_poly_heap_odd() {
		oddbabycount(heap_k, heap_m, deg);
		generate_poly_heap_manual(heap_k, heap_m);
	}
	void Polynomial::generate_poly_heap() {
		babycount(heap_k, heap_m, deg);
		generate_poly_heap_manual(heap_k, heap_m);
	}

	void Polynomial::write_heap_to_file(ostream &out) {
        RR::SetOutputPrecision(17);
		out << heaplen << endl;

		for (int index = 0; index < heaplen; index++) {
			if(poly_heap[index]) {
				out << index << " " << poly_heap[index]->deg << endl;
				for (int i = 0; i <= poly_heap[index]->deg; i++) {
					out << poly_heap[index]->chebcoeff[i] << endl;
				}
				out << endl;
			}
		}
	}

	void Polynomial::read_heap_from_file(istream &in) {
		long index = 0, in_deg;

		in >> heaplen;

		poly_heap = new Polynomial*[heaplen];

		for(int i = 0; i < heaplen; i++) {
			poly_heap[i] = 0;
		}

		while(index < heaplen - 1) {
			in >> index >> in_deg;
			poly_heap[index] = new Polynomial(in_deg);
			for(int i = 0; i <= in_deg; i++) {
				in >> poly_heap[index]->chebcoeff[i];
			}
			poly_heap[index]->cheb_to_power();
		}

		copy(*poly_heap[0]);
	}

	string to_string(const RR &x) {
		stringstream ss;
		ss << x;
		return ss.str();
	}

	void Polynomial::homomorphic_poly_evaluation(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, CompatEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor) {
		cout << setprecision(16);
		string eval_trace("(x)");
		string left_paret("(");
		
		double eps = 1./cipher.scale();

		if (deg == 1) {
			evaluator.multiply_const(cipher, coeff[1], rtn);
			eval_trace = left_paret + to_string(coeff[1]) + " * " + eval_trace + ")";
            evaluator.rescale_to_next_inplace(rtn);
			evaluator.add_const(rtn, coeff[0], rtn);
			eval_trace = left_paret + to_string(coeff[0]) + " + " + eval_trace + ")";

			cout << eval_trace << endl;
			return;
		} else if (deg == 2) {
			Ciphertext squared;
			evaluator.square(cipher, squared);
			eval_trace = left_paret + eval_trace + "^2)";
            evaluator.relinearize_inplace(squared, relin_keys);
            evaluator.rescale_to_next_inplace(squared);
			evaluator.multiply_const_inplace(squared, coeff[2]);
			eval_trace = left_paret + to_string(coeff[2]) + " * " + eval_trace + ")";
            evaluator.rescale_to_next_inplace(squared);

			if (abs(coeff[1]) >= eps) {
				evaluator.multiply_const(cipher, coeff[1], rtn);
				string temp_trace = left_paret + to_string(coeff[1]) + " * " + "(x))";
                evaluator.rescale_to_next_inplace(rtn);
				evaluator.add_reduced_error(rtn, squared, rtn);
				eval_trace = left_paret + temp_trace + " + " + eval_trace + ")";
			}

			else rtn = squared;

			evaluator.add_const_inplace(rtn, coeff[0]);
			eval_trace = left_paret + to_string(coeff[0]) + " + " + eval_trace + ")";

			cout << eval_trace << endl;
			return;
		} else if (deg == 3) {
			Ciphertext squared, cubic;
			evaluator.square(cipher, squared);
			auto square_trace = left_paret + eval_trace + "^2)";
            evaluator.relinearize_inplace(squared, relin_keys);
            evaluator.rescale_to_next_inplace(squared);
			evaluator.multiply_const(cipher, coeff[3], cubic);
			auto cubic_trace = left_paret + to_string(coeff[3]) + " * " + eval_trace + ")";
            evaluator.rescale_to_next_inplace(cubic);
			evaluator.multiply_inplace_reduced_error(cubic, squared, relin_keys);
			cubic_trace = left_paret + cubic_trace + " * " + square_trace + ")";
            evaluator.rescale_to_next_inplace(cubic);

			if (abs(coeff[1]) >= eps) {
				evaluator.multiply_const(cipher, coeff[1], rtn);
				eval_trace = left_paret + to_string(coeff[1]) + " * " + eval_trace + ")";
                evaluator.rescale_to_next_inplace(rtn);
				evaluator.add_reduced_error(rtn, cubic, rtn);
				eval_trace = left_paret + eval_trace + " + " + cubic_trace + ")";
			} else {
				rtn = cubic;
				eval_trace = cubic_trace;
			}

			if (abs(coeff[2]) >= eps) {
				evaluator.multiply_const_inplace(squared, coeff[2]);
				square_trace = left_paret + to_string(coeff[2]) + " * " + square_trace + ")";
                evaluator.rescale_to_next_inplace(squared);
				evaluator.add_reduced_error(rtn, squared, rtn);
				eval_trace = left_paret + square_trace + " + " + eval_trace + ")";
			}

			evaluator.add_const_inplace(rtn, coeff[0]);
			eval_trace = left_paret + to_string(coeff[0]) + " + " + eval_trace + ")";

			cout << eval_trace << endl;
			return;
		}

		vector<Ciphertext> baby(heap_k, Ciphertext());
		vector<string> baby_trace(heap_k, "");
		vector<bool> babybool(heap_k, false);
		Ciphertext addtmp1, addtmp2;

		Plaintext tmpplain;
		
		baby[1] = cipher;
		babybool[1] = true;
		baby_trace[1] = "(x)";

		for (int i = 2; i < heap_k; i*=2) {
			evaluator.square(baby[i/2], baby[i]);
			baby_trace[i] = left_paret + baby_trace[i / 2] + "^2)";
			evaluator.relinearize_inplace(baby[i], relin_keys);
			evaluator.rescale_to_next_inplace(baby[i]);
			evaluator.double_inplace(baby[i]);
			baby_trace[i] = left_paret + baby_trace[i] + " * 2)";

			evaluator.add_const(baby[i], -1.0, baby[i]);
			baby_trace[i] = left_paret + baby_trace[i] + " - 1)";
			babybool[i] = true;
		}

		long lpow2, res, diff;

		for (int i = 1; i < heap_k; i++) {
			if(!babybool[i]) {
				lpow2 = (1 << (int)floor(log(i) / log(2)));
				res = i - lpow2;
				diff = abs(lpow2 - res);

				evaluator.multiply_reduced_error(baby[lpow2], baby[res], relin_keys, baby[i]);
				baby_trace[i] = left_paret + baby_trace[lpow2] + " * " + baby_trace[res] + ")";
				evaluator.rescale_to_next_inplace(baby[i]);
				evaluator.double_inplace(baby[i]);
				baby_trace[i] = left_paret + baby_trace[i] + " * 2)";
				
				evaluator.sub_reduced_error(baby[i], baby[diff], baby[i]);
				baby_trace[i] = left_paret + baby_trace[i] + " - " + baby_trace[diff] + ")";
				babybool[i] = true;
			}
		}

		vector<Ciphertext> giant(heap_m, Ciphertext());
		vector<string> giant_trace(heap_m, "");
		vector<bool> giantbool(heap_m, false);
		giantbool[0] = true;
		lpow2 = (1 << ((int)ceil(log(heap_k) / log(2)) - 1));
		res = heap_k - lpow2;
		diff = abs(lpow2 - res);

		if (res == 0) {
			giant[0] = baby[lpow2];
			giant_trace[0] = baby_trace[lpow2];
		}

		else if (diff == 0) {
			evaluator.square(baby[lpow2], giant[0]);
			giant_trace[0] = left_paret + baby_trace[lpow2] + "^2)";
			evaluator.relinearize_inplace(giant[0], relin_keys);
			evaluator.rescale_to_next_inplace(giant[0]);
			evaluator.double_inplace(giant[0]);
			giant_trace[0] = left_paret + giant_trace[0] + " * 2)";
			evaluator.add_const(giant[0], -1.0, giant[0]);
			giant_trace[0] = left_paret + giant_trace[0] + " - 1)";
		}

		else {
			evaluator.multiply_reduced_error(baby[lpow2], baby[res], relin_keys, giant[0]);
			giant_trace[0] = left_paret + baby_trace[lpow2] + " * " + baby_trace[res] + ")";
 			evaluator.rescale_to_next_inplace(giant[0]);
			evaluator.double_inplace(giant[0]);
			giant_trace[0] = left_paret + giant_trace[0] + " * 2)";
			evaluator.sub_reduced_error(giant[0], baby[diff], giant[0]);
			giant_trace[0] = left_paret + giant_trace[0] + " - " + baby_trace[diff] + ")";
		}

		for (int i = 1; i < heap_m; i++) {
			giantbool[i] = true;
			
			evaluator.square(giant[i - 1], giant[i]);
			giant_trace[i] = left_paret + giant_trace[i - 1] + "^2)";
			evaluator.relinearize_inplace(giant[i], relin_keys);
			evaluator.rescale_to_next_inplace(giant[i]);
			evaluator.double_inplace(giant[i]);
			giant_trace[i] = left_paret + giant_trace[i] + " * 2)";
			
			evaluator.add_const_inplace(giant[i], -1.0);
			giant_trace[i] = left_paret + giant_trace[i] + " - 1)";
		}

		vector<Ciphertext> cipherheap((1 << (heap_m + 1)) - 1, Ciphertext());
		vector<string> heap_trace((1 << (heap_m + 1)) - 1, "");
		vector<bool> cipherheapbool((1 << (heap_m + 1)) - 1, false);

		long heapfirst = (1 << heap_m) - 1;
		long heaplast = (1 << (heap_m + 1)) - 1;
		
		for (int i = heapfirst; i < heaplast; i++) {
			if(poly_heap[i]) {
				cipherheapbool[i] = true;
				
				evaluator.multiply_const(baby[1], poly_heap[i]->chebcoeff[1], cipherheap[i]);
				heap_trace[i] = left_paret + to_string(poly_heap[i]->chebcoeff[1]) + " * " + baby_trace[1] + ")";
				evaluator.rescale_to_next_inplace(cipherheap[i]);

				if( !(abs(poly_heap[i]->chebcoeff[0]) <= eps) ) {
					evaluator.add_const_inplace(cipherheap[i], poly_heap[i]->chebcoeff[0]);
					heap_trace[i] = left_paret + heap_trace[i] + " + " + to_string(poly_heap[i]->chebcoeff[0]) + ")";
				}
				
				Ciphertext tmp;
				string temp_trace;
				for (int j = 2; j <= poly_heap[i]->deg; j++) {
					if( abs(poly_heap[i]->chebcoeff[j]) <= eps )
						continue;
									
					if (j < heap_k) {
						evaluator.multiply_const(baby[j], poly_heap[i]->chebcoeff[j], tmp);
						temp_trace = left_paret + to_string(poly_heap[i]->chebcoeff[j]) + " * " + baby_trace[j] + ")";
						evaluator.rescale_to_next_inplace(tmp);
					}
					else {
						evaluator.multiply_const(giant[0], poly_heap[i]->chebcoeff[j], tmp);
						temp_trace = left_paret + to_string(poly_heap[i]->chebcoeff[j]) + " * " + giant_trace[0] + ")";
						evaluator.rescale_to_next_inplace(tmp);
					}
			
					evaluator.add_reduced_error(cipherheap[i], tmp, cipherheap[i]);
					heap_trace[i] = left_paret + temp_trace + " + " + heap_trace[i] + ")";
				}
			}
		}
		long depth = heap_m;
		long gindex = 0;

		while (depth != 0) {
			depth--;
			heapfirst = (1 << depth) - 1;
			heaplast = (1 << (depth + 1)) - 1;
			for (int i = heapfirst; i < heaplast; i++) {		
				if(poly_heap[i]) {
					cipherheapbool[i] = true;
					if (!cipherheapbool[2 * (i + 1) - 1]) {
						cipherheap[i] = cipherheap[2 * (i + 1)];
						heap_trace[i] = heap_trace[2 * (i + 1)];
					}
					else {
						evaluator.multiply_reduced_error(cipherheap[2 * (i + 1) - 1], giant[gindex], relin_keys, cipherheap[i]);
						heap_trace[i] = left_paret + heap_trace[2 * (i + 1) - 1] + " * " + giant_trace[gindex] + ")";
						evaluator.rescale_to_next_inplace(cipherheap[i]);
						evaluator.add_reduced_error(cipherheap[i], cipherheap[2 * (i + 1)], cipherheap[i]);
						heap_trace[i] = left_paret + heap_trace[i] + " + " + heap_trace[2 * (i + 1)] + ")";
					}

				}

			}
			gindex++;
		}

		rtn = cipherheap[0];
		// cout << heap_trace[0] << endl;
	}
	void Polynomial::homomorphic_poly_evaluation_naive(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, CompatEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor) {
		vector<pair<size_t, size_t>> dp_table(this->deg + 1, {this->deg, 1}); // (min level, larger component)
		auto dp_at = [&](size_t i) -> pair<size_t, size_t>& { 
			if (i == 0) {
				throw invalid_argument("try to retrieve dp[0]");
			}
			return dp_table.at(i);
		};

		dp_at(1) = make_pair(0, 1);
		for (size_t i = 2; i <= this->deg; i++) {
			auto &min_level = dp_at(i).first;
			auto &larger_comp = dp_at(i).second;
			
			min_level = this->deg;
			for (size_t j = i - 1; j >= (i + 1) / 2; j--) {
				auto rest = i - j;
				if (dp_at(j).first + 1 < min_level && dp_at(rest).first + 1 < min_level) {
					min_level = max(dp_at(j).first, dp_at(rest).first) + 1;
					larger_comp = j;
				}
			}
		}

		vector<Ciphertext> x_powers(this->deg + 1);
		auto pow_at = [&](size_t i) -> Ciphertext& { 
			if (i == 0) {
				throw invalid_argument("try to retrieve x_powers[0]");
			}
			return x_powers.at(i);
		};

		pow_at(1) = cipher;
		for (size_t i = 2; i <= this->deg; i++) {
			auto larger_comp = dp_at(i).second;
			auto rest = i - larger_comp;

			evaluator.multiply_reduced_error(pow_at(larger_comp), pow_at(rest), relin_keys, pow_at(i));
			evaluator.rescale_to_next_inplace(pow_at(i));
		}

		auto &acc(pow_at(this->deg));
		auto parms_to_align = acc.parms_id();
		evaluator.multiply_const_inplace(acc, this->coeff[this->deg]);
		evaluator.rescale_to_next_inplace(acc);
		auto scale_to_align = RR(acc.scale());
		for (size_t i = this->deg - 1; i > 0; i--) {
			while (pow_at(i).parms_id() != parms_to_align) {
				evaluator.mod_switch_to_next_inplace(pow_at(i));
			}
			evaluator.multiply_const_inplace_fit_to_scale(pow_at(i), this->coeff[i], scale_to_align);
			evaluator.rescale_to_next_inplace(pow_at(i));

			evaluator.add_inplace(acc, pow_at(i));
		}
		evaluator.add_const_inplace(acc, this->coeff[0]);

		rtn = acc;

	}

	void Polynomial::homomorphic_poly_evaluation_cheb_naive(SEALContext &context, CKKSEncoder &encoder, Encryptor &encryptor, CompatEvaluator &evaluator, RelinKeys &relin_keys, Ciphertext &rtn, Ciphertext &cipher, Decryptor &decryptor) {
		vector<pair<size_t, size_t>> dp_table(this->deg + 1, {this->deg, 1}); // (min level, larger component)
		auto dp_at = [&](size_t i) -> pair<size_t, size_t>& { 
			return dp_table.at(i);
		};

		dp_at(0) = make_pair(0, 0);
		dp_at(1) = make_pair(0, 1);
		for (size_t i = 2; i <= this->deg; i++) {
			auto &min_level = dp_at(i).first;
			auto &larger_comp = dp_at(i).second;
			
			min_level = this->deg;
			for (size_t j = i - 1; j >= i - j; j--) {
				auto rest = i - j;
				if (dp_at(j).first + 1 < min_level && dp_at(rest).first + 1 < min_level) {
					min_level = max(dp_at(j).first, dp_at(rest).first) + 1;
					larger_comp = j;
				}
			}
		}

		vector<Ciphertext> T(this->deg + 1);
		auto T_at = [&](size_t i) -> Ciphertext& { 
			if (i == 0) {
				throw invalid_argument("try to retrieve T[0] as ciphertext");
			}
			return T.at(i);
		};

		T_at(1) = cipher;
		for (size_t i = 2; i <= this->deg; i++) {
			auto larger_comp = dp_at(i).second;
			auto rest = i - larger_comp;

			evaluator.multiply_reduced_error(T_at(larger_comp), T_at(rest), relin_keys, T_at(i));
			evaluator.double_inplace(T_at(i));
			evaluator.rescale_to_next_inplace(T_at(i));
			if (larger_comp == rest) {
				evaluator.add_const_inplace(T_at(i), -1);
			} else {
				evaluator.sub_inplace_reduced_error(T_at(i), T_at(larger_comp - rest));
			}
		}

		auto &acc(T_at(this->deg));
		auto parms_to_align = acc.parms_id();
		evaluator.multiply_const_inplace(acc, this->chebcoeff[this->deg]);
		evaluator.rescale_to_next_inplace(acc);
		auto scale_to_align = RR(acc.scale());
		for (size_t i = this->deg - 1; i > 0; i--) {
			while (T_at(i).parms_id() != parms_to_align) {
				evaluator.mod_switch_to_next_inplace(T_at(i));
			}
			evaluator.multiply_const_inplace_fit_to_scale(T_at(i), this->chebcoeff[i], scale_to_align);
			evaluator.rescale_to_next_inplace(T_at(i));

			evaluator.add_inplace(acc, T_at(i));
		}
		evaluator.add_const_inplace(acc, this->chebcoeff[0]);

		rtn = acc;

	}

	void mul(Polynomial &rtn, Polynomial &a, Polynomial &b) {
		rtn.set_zero_polynomial(a.deg + b.deg);
		for(long i = 0; i <= rtn.deg; i++) {
			for(long j = 0; j <= i; j++) {
				if(j <= a.deg && i - j <= b.deg) {
					rtn.coeff[i] += a.coeff[j] * b.coeff[i - j];
				}
			}
		}
	}

	void add(Polynomial &rtn, Polynomial &a, Polynomial &b) {
		if(a.deg >= b.deg) {
			rtn.copy(a);
			for(int i = 0; i <= b.deg; i++) {
				rtn.coeff[i] += b.coeff[i];
			}
		}

		else {
			rtn.copy(b);
			for(int i = 0; i < a.deg; i++) {
				rtn.coeff[i] += a.coeff[i];
			}
		}
	}

	void subt(Polynomial &rtn, Polynomial &a, Polynomial &b) {
		if(a.deg >= b.deg) {
			rtn.copy(a);
			for(int i = 0; i <= b.deg; i++) {
				rtn.coeff[i] -= b.coeff[i];
			}
		}

		else {
			rtn.copy(b);
			for(int i = 0; i < b.deg; i++) {
				rtn.coeff[i] *= RR(-1);
			}

			for(int i = 0; i < a.deg; i++) {
				rtn.coeff[i] += a.coeff[i];
			}
		}
	}


	void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider) {
		if(target.deg < divider.deg) {
			quotient.set_zero_polynomial(0);
			quotient.coeff[0] = 0;

			remainder.copy(target);
		}

		else {
			quotient.set_zero_polynomial(target.deg - divider.deg);
			
			Polynomial tmp(target.deg, target.coeff, "power");
			RR currcoeff;
			long currdeg = target.deg;
			for(long i = quotient.deg; i >= 0; i--) {
				currcoeff = tmp.coeff[currdeg] / divider.coeff[divider.deg];
				quotient.coeff[i] = currcoeff;
				tmp.coeff[currdeg] = 0;
				for(int j = 0; j < divider.deg; j++) {
					tmp.coeff[j + i] -= divider.coeff[j] * currcoeff;
				}
				currdeg--;
				
			}
			remainder.set_polynomial(divider.deg - 1, tmp.coeff, "power");
		}
	}

	void chebyshev(Polynomial &rtn, long deg) {
		if(deg == 0) {
			rtn.set_zero_polynomial(0);
			rtn.coeff[0] = RR(1);
		}

		else if(deg == 1) {
			rtn.set_zero_polynomial(1);
			rtn.coeff[1] = RR(1);
		}

		else {
			Polynomial iden2(1);
			iden2.coeff[1] = RR(2);

			Polynomial tmp1(0), tmp2(1), tmp3;

			tmp1.coeff[0] = RR(1);
			tmp2.coeff[1] = RR(1);

			for(int i = 2; i <= deg; i++) {
				mul(tmp3, iden2, tmp2);
				subt(rtn, tmp3, tmp1);
				tmp1.copy(tmp2);
				tmp2.copy(rtn);
			}
		}
	}

	void second_chebyshev_times_x_for_sine(Polynomial &rtn, long deg) {
		if(deg == 1) {
			rtn.set_zero_polynomial(1);
			rtn.coeff[1] = RR(1);
		}

		else if(deg == 3) {
			rtn.set_zero_polynomial(3);
			rtn.coeff[1] = RR(3);
			rtn.coeff[3] = RR(-4);
		}

		else if(deg % 2 == 1 && deg > 0){
			Polynomial mul_fac(2);
			mul_fac.coeff[0] = RR(2);
			mul_fac.coeff[2] = RR(-4);

			Polynomial iden(1);
			iden.coeff[1] = RR(1);

			Polynomial tmp1(0), tmp2(2), tmp3;

			tmp1.coeff[0] = RR(1);
			tmp2.coeff[0] = RR(3);
			tmp2.coeff[2] = RR(-4);

			for(int i = 4; i < deg; i += 2) {
				mul(tmp3, mul_fac, tmp2);
				subt(rtn, tmp3, tmp1);
				tmp1.copy(tmp2);
				tmp2.copy(rtn);
			
			}
			rtn.mulinplace(iden);
		} 
	}
}
