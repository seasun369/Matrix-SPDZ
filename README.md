# Implementations for Dishonest Majority MPC over Matrix Rings

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Build Status](https://github.com/YOUR_USERNAME/YOUR_REPOSITORY/actions/workflows/cmake.yml/badge.svg)](https://github.com/YOUR_USERNAME/YOUR_REPOSITORY/actions/workflows/cmake.yml) ## Overview

This repository provides the C++ implementation accompanying our research paper:

**"[Full Title of Your Paper: Dishonest Majority Multiparty Computation over Matrix Rings]"** *Authors: [Author 1, Author 2, ...]* *Published in/Presented at: [Conference/Journal Name or ePrint Archive Link, Year]*

The code focuses on secure multiparty computation (MPC) protocols designed to operate directly over **matrix rings** (e.g., $M_n(R)$ for some base ring $R$, specify if applicable, e.g., $M_n(\mathbb{Z}_{p^k})$) in the **dishonest majority setting** (where up to $n-1$ parties out of $n$ can be corrupt).

This library is built upon the robust and efficient Oblivious Transfer (OT) library, [**libOTe**](https://github.com/osu-crypto/libOTe).

## Features

The primary contributions implemented in this repository include:

* **Vector Oblivious Polynomial Evaluation (VOPE) / Vector OLE over Matrix Rings:** Efficient protocols for VOPE/VOLE tailored for matrix ring elements. *[Optional: Briefly clarify what VOPE means in your context if ambiguous, e.g., related to OLE or specific polynomial evaluation]*.
* **Matrix Multiplication Triples Generation:** Secure protocols to generate multiplication triples $(A, B, C)$ such that $AB = C$, where $A, B, C$ are secret-shared matrices over the target ring. These are fundamental building blocks for evaluating matrix multiplications within MPC.
* **Supporting Infrastructure:** Utilities and framework components necessary for performing MPC over the specified matrix rings using `libOTe`.

## Dependencies

* [**libOTe**](https://github.com/osu-crypto/libOTe): Included as a submodule or required as an external dependency. (Our build system assumes it's a submodule).
* **CMake:** Version 3.15 or higher recommended.
* **C++ Compiler:** A modern C++ compiler supporting C++17 (e.g., GCC 9+ or Clang 9+).
* **Boost:** Required by `libOTe`. (Usually Boost.System, Boost.Thread, etc.)
* *[Optional: Add any other specific dependencies, e.g., GMP, NTL]*

## Building

We use CMake for building the project. `libOTe` is included as a Git submodule.

```bash
# 1. Clone the repository recursively to fetch submodules (like libOTe)
git clone --recursive [https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git](https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git)
cd YOUR_REPOSITORY_NAME

# 2. Configure the build using CMake
#    Create a build directory
mkdir build && cd build
#    Run CMake (adjust options if needed, e.g., -DCMAKE_BUILD_TYPE=Release)
cmake ..

# 3. Compile the project
#    Use '-j' for parallel compilation (e.g., -j4 for 4 cores)
make -j

# The executables/libraries will be located in the 'build/' directory
# (e.g., build/bin/, build/lib/, or directly in build/ depending on CMakeLists.txt)
