#pragma once
#include<complex>
void fndawsn();
double dawsn(double x);
std::complex<double>fried_conte_Z(double x);

// dawsn.cpp : Defines the functions for the static library.
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
// Optimized into a branchless CMOV configuration using IEEE 754 hardware operations.
// ============================================================================