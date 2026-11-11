// ============================================================================
// Derived from the Dawson Integral (dawsn.c) in the Cephes Mathematical Library.
// Copyright 1984, 1987, 2000 by Stephen L. Moshier.
// 
// Permission to use, copy, modify, and distribute this software for any
// purpose without fee is hereby granted, provided that this entire notice 
// is included in all copies of any software which is or includes a copy 
// or modification of this software.
//
// Optimized into a branchless NaN-free version CMOV configuration using IEEE 754 hardware operations.
// ============================================================================
#include "pch.h"
#include "framework.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include<complex>
// TODO: This is an example of a library function

// ============================================================================
// CEPHES CORE MATHEMATICAL POLYNOMIAL ENGINE
// ============================================================================



// Ensure these Cephes helper definitions are placed above your function
double cephes_polevl(double x, const double coef[], int N) {
    double ans = coef[0];
    for (int i = 1; i <= N; i++) {
        ans = ans * x + coef[i];
    }
    return ans;
}
double cephes_p1evl(double x, const double coef[], int N) {
    double ans = x + coef[0];
    for (int i = 1; i < N; i++) {
        ans = ans * x + coef[i];
    }
    return ans;
}
double dawsn(double x) {
    double PI = 4 * std::atan(1);
    double e1_x = std::pow(std::exp(800.0 * (double)(x - 3.25)), 2.0);
    double  e2_x = std::pow(std::exp(800 * (x + 3.25)), 2);
    double  e3_x = std::pow(std::exp(800 * (x - 6.25)), 2);///
    double  e4_x = std::pow(std::exp(800 * (x + 6.25)), 2);///
    double  R5_x = (std::exp(-pow(e1_x, 2))) * (-1 * std::exp(-pow(e2_x, 2)) + 1);//////
    double  R6_x = ((-1 * std::exp(-pow(e1_x, 2)) + 1) + (std::exp(-pow(e2_x, 2)))) * (std::exp(-pow(e3_x, 2))) + ((-1 * std::exp(-pow(e1_x, 2)) + 1) + (-1 * std::exp(-pow(e2_x, 2)))) * (std::exp(-pow(e4_x, 2)));
    double  R7_x = (-1 * std::exp(-pow(e3_x, 2)) + 1) + (std::exp(-pow(e4_x, 2)));
    double  R8_x = (std::exp(-pow(e3_x, 2))) * ((-1 * std::exp(-pow(e4_x, 2)) + 1));
    double  R9_x = (2 * std::atan(std::pow(10, 100) * x)) / PI;
    const double AN[] = {
        1.13681498971755972054e-11,
        8.49262267667473811108e-10,
        1.94434204175553054283e-08,
        9.53151741254484363489e-07,
        3.07828309874913200438e-06,
        3.52513368520288738649e-04,
       -8.50149846724410912031e-04,
        4.22618223005546594270e-02,
       -9.17480371773452345351e-02,
        9.99999999999999994612e-01
    };

    const double AD[] = {
        2.40372073066762605484e-11,
        1.48864681368493396752e-09,
        5.21265281010541664570e-08,
        1.27258478273186970203e-06,
        2.32490249820789513991e-05,
        3.25524741826057911661e-04,
        3.48805814657162590916e-03,
        2.79448531198828973716e-02,
        1.58874241960120565368e-01,
        5.74918629489320327824e-01,
        1.00000000000000000539e0
    }; 

    
    double t_A = std::pow(std::abs(x), 2 * R5_x);
    double F3_X = std::pow(std::abs(x), R5_x) * R5_x * (cephes_polevl(t_A, AN, 9) / cephes_polevl(t_A, AD, 10));
    

    const double BN[] = {
        5.08955156417900903354e-01,
        -2.44754418142697847934e-01,
        9.41512335303534411857e-02,
        -2.18711255142039025206e-02,      
        3.66207612329569181322e-03,
        -4.23209114460388756528e-04,
        3.59641304793896631888e-05,
        -2.14640351719968974225e-06,
        9.10010780076391431042e-08,
        -2.40274520828250956942e-09,
        3.59233385440928410398e-11
    };

    const double BD[] = {
        -6.31839869873368190192e-01,
        2.36706788228248691528e-01,
        -5.31806367003223277662e-02,
        8.48041718586295374409e-03,
        -9.47996768486665330168e-04,
        7.81025592944552338085e-05,
        -4.55875153252442634831e-06,
        1.89100358111421846170e-07,
        -4.91324691331920606875e-09,
        7.18466403235734541950e-11
    };
    double t_BC = 1.0 / (std::pow(std::abs(x) + R5_x, 2));
    double F5_X = (1.0 / (2.0 * std::abs(x) + R5_x)) + (t_BC / (2.0 * std::abs(x) + R5_x)) *
        (cephes_polevl(t_BC, BN, 10) / cephes_p1evl(t_BC, BD, 10));

    const double CN[] = {
        -5.90592860534773254987e-01,
        6.29235242724368800674e-01,
        -1.72858975380388136411e-01,
        1.64837047825189632310e-02,
        -4.86827613020462700845e-04
    };

    const double CD[] = {
        -2.69820057197544900361e+00,
        1.73270799045947845857e+00,
        -3.93708582281939493482e-01,
        3.44278924041233391079e-02,
        -9.73655226040941223894e-04
    };

    double t1_BC = 1.0 / (std::pow(std::abs(x) + R8_x, 2));
    double F7_X = (1.0 / (2.0 * std::abs(x) + R8_x)) + (t_BC / (2.0 * std::abs(x) + R8_x)) *
        (cephes_polevl(t1_BC, CN, 4) / cephes_p1evl(t1_BC, CD, 5));


    double result = R9_x * (F3_X + R6_x * F5_X + R7_x * F7_X);
    return result;
}

void fndawsn()
{
	
    std::cout << "Real Part values at break points:\n";
    std::cout << "Image Part(-6.25):-0.162132188122\n";
    std::cout << "Image Part(-3.25):-0.325141829121\n";
    std::cout << "Image Part(3.25):0.325141829121\n";
    std::cout << "Image Part(6.25):0.162132188122\n";
}

std::complex<double> fried_conte_Z(double x) {
    const double SQRT_PI = 1.77245385090551602730; 
    double real_part = -2.0 * dawsn(x);
    double imag_part = SQRT_PI * std::exp(-x * x);
    return std::complex<double>(real_part, imag_part);
}
