# Cross-Language Plasma Dispersion Engine

**Copyright (c) 2026 RMIO Enterprise, LLC**  
*Licensed under the MIT License.*

## Overview
A high-performance, cross-language (C++ and Python) numerical library for evaluating the **Fried-Conte plasma dispersion relation**. This engine models electromagnetic wave transitions as phase velocities accelerate across kinetic thermal bounds into the superluminal fluid regime ($v_{ph} > c$).

## Repository Architecture
* **Branchless C++ Core (`dawsn.cpp`, `dawsn.h`, `project.cpp`):** A standalone static library engineered to compile into hardware Conditional Moves (`CMOV`),  NaN-free and eliminating CPU branch-prediction stalls during heavy simulation loops.
* **Interactive Domain Visualizer (`branchless.py`):** A Python prototyping module featuring an interactive slider to dynamically inspect real-part Dawson decay and imaginary Landau damping across wave speeds.
* **Commercial & Academic Integrity (`LICENSE.txt`):** Secured under RMIO Enterprise, LLC copyright with explicit attribution to the classic Cephes Mathematical Library rational approximations by Stephen L. Moshier.

## Usage
Download and extract `plasma-dispersion-engine.zip` from the file list above to inspect the raw C++ source files, run the interactive Python visualizer, or integrate the static library into high-performance computing (HPC) pipelines.
