#include "ModularReducer.h"
#include "PolyUpdate.h"
#include "comp/SEALfunc.h"
#include "comp/program.h"
#include <fstream>
#include <iomanip>

ModularReducer::ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg,
    SEALContext &_context,
	CKKSEncoder &_encoder,
    Encryptor &_encryptor,
    CompatEvaluator &_evaluator,
    RelinKeys &_relin_keys,
    Decryptor &_decryptor
    ) : 
	boundary_K(_boundary_K), log_width(_log_width), deg(_deg), num_double_formula(_num_double_formula), inverse_deg(_inverse_deg), context(_context), encoder(_encoder), encryptor(_encryptor), evaluator(_evaluator), relin_keys(_relin_keys), decryptor(_decryptor) {

	// sin_generator = new RemezSin(rmparm, boundary_K, log_width, deg, (1 << num_double_formula));
	// cos_generator = new RemezSin(rmparm, boundary_K, log_width, deg, (1 << num_double_formula));

	// inverse_log_width = -log2(sin(2 * M_PI * pow(2.0, -log_width)));
	// inverse_sin_generator = new RemezArcsin(rmparm, inverse_log_width, inverse_deg);
	// inverse_sin_generator->params.log_scan_step_diff = 12;
	// inverse_sin_generator->params.RR_prec = 1000;
}

void ModularReducer::double_angle_formula(Ciphertext &cipher) {
	evaluator.square_inplace(cipher);
    evaluator.relinearize_inplace(cipher, relin_keys);
    evaluator.rescale_to_next_inplace(cipher);
	evaluator.double_inplace(cipher);
    evaluator.add_const(cipher, -1.0, cipher);
}

void ModularReducer::double_angle_formula_scaled(Ciphertext &cipher, double scale_coeff) {
	evaluator.square_inplace(cipher);
    evaluator.relinearize_inplace(cipher, relin_keys);
    evaluator.rescale_to_next_inplace(cipher);
	evaluator.double_inplace(cipher);
    evaluator.add_const(cipher, -scale_coeff, cipher);
}

void ModularReducer::generate_sin_cos_polynomial() {
	// sin_generator->generate_optimal_poly(sin_polynomial);
    RR coeffs[] = {
        RR(0.92387953231773009), 
        RR(15.027943245428581), 
        RR(-712.36917917753132), 
        RR(-3862.4962173546066), 
        RR(91546.900696525565), 
        RR(297822.73140750612), 
        RR(-4705894.2117281199), 
        RR(-10935240.05946226), 
        RR(129590717.97240609), 
        RR(234215477.02216466), 
        RR(-2220501807.3704581), 
        RR(-3283542596.6943443), 
        RR(25941607670.155979), 
        RR(32459195987.501604), 
        RR(-219808768676.06232), 
        RR(-238362547119.73398), 
        RR(1412383968229.0714), 
        RR(1351413420768.5189), 
        RR(-7117851373948.3598), 
        RR(-6093684508527.8351), 
        RR(28885433488237.621), 
        RR(22374069196578.916), 
        RR(-96413939332866.742), 
        RR(-68186629971363.139), 
        RR(269324742471541.98), 
        RR(175237895712024.71), 
        RR(-638808512782011.86), 
        RR(-384865939301985.97), 
        RR(1302256999237959.8), 
        RR(730509056263588.56), 
        RR(-2304968261174938.7), 
        RR(-1209713598732979.1), 
        RR(3571666231875234.7), 
        RR(1761338358467019.9), 
        RR(-4875744024649225.6), 
        RR(-2268100565747594.1), 
        RR(5886770715832721.1), 
        RR(2592558478326421.3), 
        RR(-6292431181678847.7), 
        RR(-2632767781192529), 
        RR(5939686576589964.4), 
        RR(2369063681350181), 
        RR(-4918259710060236.7), 
        RR(-1876272208761189.2), 
        RR(3532216541468631), 
        RR(1293070069731452.1), 
        RR(-2164304632527603.8), 
        RR(-762689882270018.15), 
        RR(1106207490957092.5), 
        RR(376362391898387.81), 
        RR(-457397159398795.63), 
        RR(-150657282523830.17), 
        RR(146510702076755.46), 
        RR(46835248558682.531), 
        RR(-34007271107894.592), 
        RR(-10574375024775.831), 
        RR(5075533242659.2427), 
        RR(1538211315599.5003), 
        RR(-364844435510.66296), 
        RR(-107962886396.68134)
    };
    sin_polynomial.set_polynomial(deg, coeffs, "power");
	sin_polynomial.generate_poly_heap();
    sin_polynomial.showcoeff();

    cos_polynomial.set_polynomial(deg, coeffs, "power");
	cos_polynomial.generate_poly_heap();
    cos_polynomial.showcoeff();
}

void ModularReducer::generate_inverse_sine_polynomial() {
	// inverse_sin_generator->generate_optimal_poly(inverse_sin_polynomial);
    {
        if (inverse_deg > 127) {
            throw runtime_error("We have only prepared taylor coeffs for terms of degree <= 127.");
        } else if (inverse_deg < 3) { 
            throw std::invalid_argument("We now only allow inverse deg >= 3."); 
        }

        vector<NTL::RR> taylor_coeffs {
            RR(NTL::INIT_VAL_TYPE {}, "2.2947379471529382969744498600366971923620711188040720503135663285591120454439982e-76"),
            RR(NTL::INIT_VAL_TYPE {}, "0.15897595311853101274526061204647097374617576780504320373503002412804013125389352"),
            RR(NTL::INIT_VAL_TYPE {}, "-1.1608808021533176343645051916281395736534684849452555645808255157330666353747941e-73"),
            RR(NTL::INIT_VAL_TYPE {}, "0.056648697109133325293689187145117838558785315198698457540119087978676284245672466"),
            RR(NTL::INIT_VAL_TYPE {}, "9.7846038194851369194745179519276033267815972074371410866890581905249590678060137e-72"),
            RR(NTL::INIT_VAL_TYPE {}, "-1.4962576928053791590280746333430693060007070082898805468523500985795806225772547"),
            RR(NTL::INIT_VAL_TYPE {}, "-3.2702674279438227454484772562766054732640518481367510698173844152325652533644013e-70"),
            RR(NTL::INIT_VAL_TYPE {}, "35.362074308406350474792880925586208455112536866444251515111465396777026303611403"),
            RR(NTL::INIT_VAL_TYPE {}, "5.7551890541776576190524723561565756985691176458704944761227995371007212304187709e-69"),
            RR(NTL::INIT_VAL_TYPE {}, "-471.07953431217396424071402460253047716806523886881267880515901400181061914676474"),
            RR(NTL::INIT_VAL_TYPE {}, "-6.1385995412158862994464629918848808954017674485798557953343092811489608950401273e-68"),
            RR(NTL::INIT_VAL_TYPE {}, "3964.6161629320869022883621817561482183665609895483259339857446167047706330655413"),
            RR(NTL::INIT_VAL_TYPE {}, "4.3063391773465811423959945686875889964018085109242934483393710560231252306784904e-67"),
            RR(NTL::INIT_VAL_TYPE {}, "-22459.541248865152312938065178128097300367927572725676740568338390667504135425985"),
            RR(NTL::INIT_VAL_TYPE {}, "-2.0908474497159231683075942490197489538677345843971745704378214470396584316573739e-66"),
            RR(NTL::INIT_VAL_TYPE {}, "89123.731902018683809289470223194159766482580874743306244672533078890309896168463"),
            RR(NTL::INIT_VAL_TYPE {}, "7.2563980255277427934374583835867474857085908282004973666885313524025023335521144e-66"),
            RR(NTL::INIT_VAL_TYPE {}, "-253768.63233816677971544447178108247093939938513515465922939165725620651828358590"),
            RR(NTL::INIT_VAL_TYPE {}, "-1.8350287759899356541915493275239160064134455598891720786718897728659227461676995e-65"),
            RR(NTL::INIT_VAL_TYPE {}, "524776.79801063595003095752312128769224864972771216682777456495127994560809612862"),
            RR(NTL::INIT_VAL_TYPE {}, "3.4109554196316155049189506446937865662178628136050185783972844603646948740133374e-65"),
            RR(NTL::INIT_VAL_TYPE {}, "-789109.60198359304961497985472868477663171542552926472895053159878609101337416840"),
            RR(NTL::INIT_VAL_TYPE {}, "-4.6549606705774713570078140439899260027092906414990235601470982179195983272398013e-65"),
            RR(NTL::INIT_VAL_TYPE {}, "853745.43290425213388163067627040607005174940281249565369620887280235781821082316"),
            RR(NTL::INIT_VAL_TYPE {}, "4.6067539112626614413484189888682410299940811380246496040148671655469170423923455e-65"),
            RR(NTL::INIT_VAL_TYPE {}, "-647273.70324669429840150943198157914821646349529845866115260376069874133420227858"),
            RR(NTL::INIT_VAL_TYPE {}, "-3.2156205379114234927352893083936135992902282672923966492119385660829250635757931e-65"),
            RR(NTL::INIT_VAL_TYPE {}, "326301.26986669056050070160365511101760051670477775281598819399244409811748101577"),
            RR(NTL::INIT_VAL_TYPE {}, "1.5006065588780182319983124977215458601564535200529088091947918519261862615559678e-65"),
            RR(NTL::INIT_VAL_TYPE {}, "-98218.984331747139107499019920039922801278956938360690938876111674031308704904622"),
            RR(NTL::INIT_VAL_TYPE {}, "-4.2005464365314365056914532097070591744796635056007673884285548727683343842567724e-66"),
            RR(NTL::INIT_VAL_TYPE {}, "13355.853558204680018492786448903280528249651256497222226676393541186558489225246")
        };
        inverse_sin_polynomial_v1.set_polynomial(taylor_coeffs.size() - 1, taylor_coeffs.data(), "power");
        inverse_sin_polynomial_v1.showcoeff();
    }
    inverse_sin_polynomial_v1.generate_poly_heap();    
    double inv_of_scale_for_eval = 1.0 / scale_for_eval;
    inverse_sin_polynomial_v1.constmul(RR(inv_of_scale_for_eval));
    inverse_sin_polynomial_v1.constmul(RR(0.5)); // for relu, comment this for id
    inverse_sin_polynomial_v1.generate_poly_heap_odd();
    inverse_sin_polynomial_v1_original.copy(inverse_sin_polynomial_v1);

    upgrade_oddbaby(127, arcsin_tree);
    string arcsin_coeff_file = "arcsin_d127.txt";
    ifstream tmp_ifs(arcsin_coeff_file);
    if (not tmp_ifs.good()) {
        throw runtime_error(string("No file \"") + arcsin_coeff_file + "\"");
    }

    for (int j = 0; j < coeff_number(127, arcsin_tree); j++) {
        RR temp;
        tmp_ifs >> temp;

        // temp *= RR(0.5); // for abs
        temp *= RR(0.25); // for relu
        temp *= RR(inv_of_scale_for_eval);
        arcsin_decomp_coeff.emplace_back(temp);
        // cout << temp << " ";
    }
    // cout << endl;
    arcsin_decomp_coeff_original = arcsin_decomp_coeff;

    // abs_lift = 0.125 * inv_of_scale_for_eval; // no 0.5 for abs
    abs_lift = 0.125 * inv_of_scale_for_eval * 0.5; // one 0.5 for relu
}

double scale_for_boost_relu_range = 2.0;
double dyn_scale_for_turn_back_q = std::pow(2.0, 46) / 69834700000000.0 / scale_for_boost_relu_range; // for cnn first boot (magic num, to check again)
bool after_cnn_first_boot = false;

void ModularReducer::scaling_for_turn_back_q() {

    inverse_sin_polynomial_v1.copy(inverse_sin_polynomial_v1_original);
    inverse_sin_polynomial_v1.constmul(to_RR(dyn_scale_for_turn_back_q));
    inverse_sin_polynomial_v1.generate_poly_heap_odd();

    arcsin_decomp_coeff = arcsin_decomp_coeff_original;
    for (int j = 0; j < arcsin_decomp_coeff.size(); j++) {
        arcsin_decomp_coeff[j] *= RR(dyn_scale_for_turn_back_q);
    }

    double inv_of_scale_for_eval = 1.0 / scale_for_eval;
    // abs_lift = 0.125 * inv_of_scale_for_eval * dyn_scale_for_turn_back_q; // no 0.5 for abs
    abs_lift = 0.125 * inv_of_scale_for_eval * 0.5 * dyn_scale_for_turn_back_q; // one 0.5 for relu
}

void ModularReducer::write_polynomials() {
    ofstream sin_out("sine.txt"), cos_out("cosine.txt"), inverse_out1("inverse_sine_v1.txt"), inverse_out2("inverse_sine_v2.txt");
    sin_polynomial.write_heap_to_file(sin_out);
    cos_polynomial.write_heap_to_file(cos_out);
    inverse_sin_polynomial_v1.write_heap_to_file(inverse_out1);
    // inverse_sin_polynomial_v2.write_heap_to_file(inverse_out2);
    sin_out.close();
    cos_out.close();
    inverse_out1.close();
    // inverse_out2.close();
}

void ModularReducer::modular_reduction(Ciphertext &rtn, Ciphertext &cipher) {
	Ciphertext sin_tmp1, sin_tmp2, cos_tmp1, cos_tmp2;
	sin_tmp1 = cipher;
    cos_tmp1 = cipher;
    Ciphertext sin_rtn, cos_rtn;

    cout << "[DIAG] before sin poly, cipher #q=" << cos_tmp1.coeff_modulus_size() << " scale=" << log2(cos_tmp1.scale()) << endl;
    try {
        sin_polynomial.homomorphic_poly_evaluation(context, encoder, encryptor, evaluator, relin_keys, cos_tmp2, cos_tmp1, decryptor);
    } catch (const exception& e) {
        cout << "[DIAG] FAIL in sin_polynomial: " << e.what() << endl;
        throw;
    }
    cout << "[DIAG] after sin poly, #q = " << cos_tmp2.coeff_modulus_size() << " scale=" << log2(cos_tmp2.scale()) << endl;
    {
        for (int i = 0; i < num_double_formula + 1; i++) {
            cout << "[DIAG] double_angle_formula iteration " << i << ", before: #q=" << cos_tmp2.coeff_modulus_size() << endl;
            double_angle_formula(cos_tmp2);
            cout << "[DIAG] double_angle_formula iteration " << i << ", after:  #q=" << cos_tmp2.coeff_modulus_size() << endl;
            if (i == num_double_formula - 1) {
                sin_tmp2 = cos_tmp2;
            }
        }
        cout << "[DIAG] after double angle for cos, #q = " << cos_tmp2.coeff_modulus_size() << endl;
    }    
    std::thread thread1([&]() 
    {
        PublicKey dummy_pk;
        SecretKey dummy_sk;
        cout << "[DIAG] thread1: before eval_polynomial_integrate, cos_tmp2 #q=" << cos_tmp2.coeff_modulus_size() << endl;
        eval_polynomial_integrate(encryptor, evaluator, decryptor, encoder, dummy_pk, dummy_sk, relin_keys, cos_rtn, cos_tmp2, 127, arcsin_decomp_coeff, arcsin_tree);
        cout << "[DIAG] after eval arcsin, #q = " << cos_rtn.coeff_modulus_size() << endl;
        evaluator.add_const_inplace(cos_rtn, abs_lift);
        cout << "[DIAG] thread1 done, cos_rtn #q=" << cos_rtn.coeff_modulus_size() << endl;
    }
    );

    std::thread thread2([&]() 
    {
        cout << "[DIAG] thread2: before inverse_sin poly, sin_tmp2 #q=" << sin_tmp2.coeff_modulus_size() << endl;
        inverse_sin_polynomial_v1.homomorphic_poly_evaluation(context, encoder, encryptor, evaluator, relin_keys, sin_rtn, sin_tmp2, decryptor);
        cout << "[DIAG] thread2: after inverse_sin poly, sin_rtn #q=" << sin_rtn.coeff_modulus_size() << endl;
        evaluator.mod_switch_to_next_inplace(sin_rtn);
        cout << "[DIAG] thread2: after mod_switch, sin_rtn #q=" << sin_rtn.coeff_modulus_size() << endl;
    }
    );
    
    thread1.join(); // rtn = cos_rtn; return;
    thread2.join(); // evaluator.mod_switch_to_next_inplace(sin_rtn); rtn = sin_rtn; return;

    cout << "[mod_red] cos_rtn.scale=" << log2(cos_rtn.scale()) << " sin_rtn.scale=" << log2(sin_rtn.scale()) << endl;
    int64_t ql = (*context.get_context_data(sin_rtn.parms_id())->parms().coeff_modulus().rbegin()).value();
    double scale_adj = cos_rtn.scale() * ql / sin_rtn.scale() / sin_rtn.scale();
    cout << "[mod_red] ql=" << log2((double)ql) << " scale_adj=" << log2(abs(scale_adj)) << " adj_sign=" << (scale_adj >= 0 ? "+" : "-") << endl;
    evaluator.multiply_const_inplace(sin_rtn, scale_adj);
    evaluator.rescale_to_next_inplace(sin_rtn);
    sin_rtn.scale() = cos_rtn.scale();
    evaluator.add(sin_rtn, cos_rtn, rtn);

}
