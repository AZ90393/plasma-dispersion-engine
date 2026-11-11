// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// ============================================================================
// Derived from the Dawson Integral (dawsn.c) in the Cephes Mathematical Library.
// Copyright 1984, 1987, 2000 by Stephen L. Moshier.
// 
// Permission to use, copy, modify, and distribute this software for any
// purpose without fee is hereby granted, provided that this entire notice 
// is included in all copies of any software which is or includes a copy 
// or modification of this software.
//
// Optimized into a branchless NaN-free version  CMOV configuration using IEEE 754 hardware operations.
// ============================================================================
#include <iostream>
#include"dawsn.h"
#include<iomanip>
#include<complex>
#include <chrono>
int main()
{
    fndawsn();
    // test velocity vector(e.g,wave phase velocity)
    double normalized_velocity = 1.5;
    std::complex<double>Z_val = fried_conte_Z(normalized_velocity);
    //compute plasma dospersion directly using your modified code
    std::cout << std::fixed << std::setprecision(16);
    std::cout << "------------------------------------------\n";
    std::cout << "input velocity(x)  :"<< normalized_velocity << "\n";
    std::cout << "Fried-Conte Z(x)  : " << Z_val << "\n";
    std::cout << "-> Real Part  : " << Z_val.real() << "(via dawsn rational engine)\n";
    std::cout << "-> Image Part : " << Z_val.imag() << "(Landau Damping component)\n";
    return 0;
}


