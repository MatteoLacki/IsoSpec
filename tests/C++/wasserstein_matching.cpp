#include <iostream>
#include "unity-build.cpp"
#include <cassert>

using namespace IsoSpec;

int main()
{
    double masses1[] = {1.0, 2.0};
    double probs1[]  = {0.5, 0.5};
    double masses2[] = {0.9, 2.01};
    double probs2[]  = {0.4, 0.6};

    FixedEnvelope F1(masses1, probs1, 2);
    FixedEnvelope F2(masses2, probs2, 2);

    auto [m1, m2, flow] = F1.WassersteinMatch(F2, 0.05);
    std::cout << m1 << " " << m2 << " " << flow << std::endl;

    assert(m1 == 0.5 && m2 == 0.5 && flow == 0.5);

    F1.release_everything();
    F2.release_everything();
}

