#define ISOSPEC_TESTS_SKIP_MAIN true
#include "from_formula_layered.cpp"
#include "from_formula_ordered.cpp"
#include "from_formula_threshold.cpp"
#include "from_formula_threshold_simple.cpp"
#include "element_zero.cpp"
#include "unity-build.cpp"
#include "empty_iso.cpp"
#include <vector>

#if !defined(ISOSPEC_TESTS_MEMSAN)
#define TEST(formula, prob, function) \
std::cout << "Testing " << formula << " prob: " << prob << " function: " << #function << "..." << std::flush; \
tmp = function(formula, prob, false); \
std::cout << " " << tmp << " confs." << std::endl; \
total += tmp;
#else
#define TEST(formula, prob, function) \
tmp = function(formula, prob, false); \
total += tmp;
#endif

int main()
{
        bool zero_ok = false;
        try{
            test_zero(0.1, false);
        }
        catch(std::invalid_argument&)
        {
            zero_ok = true;
        }
        assert(zero_ok);
        test_empty_and_print();
        #if !defined(ISOSPEC_SKIP_SLOW_TESTS)
	char test_formulas[] = "P1 P2 H1 H2 O1 O2 H2O1 C0 P0 C10000P10 F10C10000P10 P10F10O100 C100O0P100 C100 P100 C1 H10C10O10N10S5 Se1 Se10 Sn1 Sn4 Sn4C1 C2H6O1 C1000 C1H1O2N2Se1Sn1P1 P1C1Sn1 Se5 Sn5 Se2Sn2C2O2N2S2B2He2U2Na2Cl2";
        #else
        char test_formulas[] = "P1 P2 H1 H2 O1 O2 H2O1 C0 P0 C10000P10 F10C10000P10 P10F10O100 C100O0P100 C100 P100 C1 H3C3O3N3S3 Se1 Se3 Sn1 Sn3C1 C2H6O1 C1000 C1H1O2N2Se1Sn1P1 P1C1Sn1 Se5";
        #endif
	size_t tf_len = strlen(test_formulas);
	std::vector<const char*> formulas;
	formulas.push_back(test_formulas);
	std::vector<float> probs = {0.0, 0.1, 0.5, 0.01, 0.9, 0.99, 0.01, 0.0001, 0.999, 0.362, 0.852348};

	for(size_t ii=0; ii<tf_len; ii++)
	{
		if(test_formulas[ii] == ' ')
		{
			test_formulas[ii] = '\0';
			formulas.push_back(&test_formulas[ii+1]);
		}
	}

	size_t total = 0;
	size_t tmp;
	for(auto it_formula = formulas.begin(); it_formula != formulas.end(); it_formula++)
		for(auto it_prob = probs.begin(); it_prob != probs.end(); it_prob++)
		{
			TEST(*it_formula, *it_prob, test_threshold_simple);
			TEST(*it_formula, *it_prob, test_threshold);
			TEST(*it_formula, *it_prob, test_layered_tabulator);
			TEST(*it_formula, *it_prob, test_ordered);
		}

        #if !defined(ISOSPEC_TESTS_MEMSAN)
	std::cout << "Total confs considered: " << total << std::endl;
        #endif
}
