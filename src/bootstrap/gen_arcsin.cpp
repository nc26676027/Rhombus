/**
 * gen_arcsin.cpp
 * Generate arcsin_d127.txt data file required by RBOOT bootstrapping.
 * Uses upgrade_oddbaby + coeff_number to determine the decomposition tree,
 * then generates the coefficients via Taylor expansion of arcsin.
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <NTL/RR.h>
#include "PolyUpdate.h"
#include "comp/program.h"
#include "comp/SEALfunc.h"

using namespace std;
using namespace NTL;
using namespace minicomp;

int main() {
    // Generate the odd-baby decomposition tree for degree 127
    Tree arcsin_tree;
    upgrade_oddbaby(127, arcsin_tree);

    long num_coeff = coeff_number(127, arcsin_tree);
    cout << "arcsin decomposition: " << num_coeff << " coefficients" << endl;

    // Taylor coefficients of arcsin(x) = x + (1/6)x^3 + (3/40)x^5 + ...
    // arcsin(x) = sum_{k=0}^{inf} C(2k,k) / (4^k * (2k+1)) * x^{2k+1}
    RR::SetPrecision(150);
    vector<RR> taylor(128);
    for (long k = 0; k <= 63; k++) {
        RR binom = RR(1);
        for (long i = 0; i < k; i++) {
            binom *= RR(2 * k - i) / RR(i + 1);
        }
        taylor[2 * k + 1] = binom / (power(RR(4), k) * RR(2 * k + 1));
    }

    // Decompose polynomial using the tree structure
    // Allocate decomposition degree array
    long sz = 1 << (arcsin_tree.depth + 1);
    vector<long> decomp_deg(sz, -1);
    vector<long> start_index(sz, -1);
    decomp_deg[1] = 127;
    for (int i = 1; i <= arcsin_tree.depth; i++) {
        for (int j = (1 << i); j < (1 << (i+1)); j++) {
            if (j % 2 == 0) decomp_deg[j] = arcsin_tree.tree[j/2] - 1;
            else decomp_deg[j] = decomp_deg[j/2] - arcsin_tree.tree[j/2];
        }
    }

    long idx = 0;
    for (int i = 1; i < sz; i++) {
        if (arcsin_tree.tree[i] == 0) {
            start_index[i] = idx;
            idx += (decomp_deg[i] + 1);
        }
    }

    // Extract decomposition coefficients
    // For each leaf node, extract the odd-degree Taylor coefficients
    vector<RR> decomp_coeff(num_coeff);
    for (int i = 1; i < sz; i++) {
        if (arcsin_tree.tree[i] == 0) {
            long si = start_index[i];
            long deg = decomp_deg[i];
            // Odd-degree coefficients only (arcsin is an odd function)
            for (long k = 0; k <= deg; k += 2) {
                if (si + k < num_coeff && (k + 1) < 128) {
                    decomp_coeff[si + k] = taylor[k + 1];
                }
            }
        }
    }

    // Write to file
    ofstream ofs("arcsin_d127.txt");
    if (!ofs.good()) {
        cerr << "ERROR: cannot write arcsin_d127.txt" << endl;
        return 1;
    }
    ofs << setprecision(80);
    for (long j = 0; j < num_coeff; j++) {
        ofs << decomp_coeff[j] << "\n";
    }
    ofs.close();

    cout << "✅ Generated arcsin_d127.txt (" << num_coeff << " coefficients)" << endl;
    return 0;
}
