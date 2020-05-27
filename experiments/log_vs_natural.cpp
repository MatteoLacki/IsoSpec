#include <cstddef>
#include <cmath>
#include <iostream>

//extern double a, b, c, d, e, f;
int main()
{
    double a = 1.3;
    double b = 1.2;
    double c = 3.4;
    double d = 2.3;
    double e = 0.1;
    double f = 2.3;

    for(size_t ii = 0; ii < 1000000000; ii++)
    {
        double res = exp(a);//+b+c+d+e);
        //double res = a+b+c+d+e;
        //double res = a*b*c*d*e;
        f += res;
    }
    std::cout << f << std::endl;
}
