#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define M_PI 3.14159265359

#include <algorithm>
#include <complex>
#include <random>
#include <iostream>
#include <ostream>


static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0, 1);



int main() {
    int N = 1000000;
    double sigma = 1;
    double s = 0;
    for (int i = 0; i < N; i++) {
        double u1 = uniform(engine);
        double u2 = uniform(engine);
        double xi = sigma * cos(2*M_PI*u1)*sqrt(-2*log(u2));
        if ((xi > -M_PI/2) & (xi < M_PI)) {
            double p = 1/(sigma*sqrt(2*M_PI)) * exp(-xi*xi/(2*sigma*sigma));
            s += std::pow(cos(xi), 10) / p / N;
            }
        }
    std::cout << s << std::endl;

    double s2 = 0;
    for (int i = 0; i < N; i++) {
        double u1 = uniform(engine);
        double u2 = uniform(engine);
        double u3 = uniform(engine);
        double u4 = uniform(engine);
        double x1 = sigma * cos(2*M_PI*u1)*sqrt(-2*log(u2));
        double x2 = sigma * sin(2*M_PI*u1)*sqrt(-2*log(u2));
        double x3 = sigma * cos(2*M_PI*u3)*sqrt(-2*log(u4));
        double x4 = sigma * sin(2*M_PI*u3)*sqrt(-2*log(u4));
        if (((x1 > -M_PI/2) && (x1 < M_PI/2)) && ((x2 > -M_PI/2) && (x2 < M_PI/2)) && ((x3 > -M_PI/2) && (x3 < M_PI/2)) && ((x4 > -M_PI/2) && (x4 < M_PI/2))) {
            double p1 = 1/(sigma*sqrt(2*M_PI)) * exp(-x1*x1/(2*sigma*sigma));
            double p2 = 1/(sigma*sqrt(2*M_PI)) * exp(-x2*x2/(2*sigma*sigma));
            double p3 = 1/(sigma*sqrt(2*M_PI)) * exp(-x3*x3/(2*sigma*sigma));
            double p4 = 1/(sigma*sqrt(2*M_PI)) * exp(-x4*x4/(2*sigma*sigma));
            s2 += std::pow(cos(x1 + x2 + x3 + x4), 2) / p1 / p2 / p3 / p4 / N;
            }
        }
    std::cout << s2 << std::endl;

    return 0;
}
