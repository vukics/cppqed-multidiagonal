# cppqed-multidiagonal

[![CI](https://github.com/vukics/cppqed-multidiagonal/actions/workflows/ci.yml/badge.svg)](https://github.com/vukics/cppqed-multidiagonal/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-github.io-blue)](https://vukics.github.io/cppqed-multidiagonal/)
[![License](https://img.shields.io/badge/license-BSL--1.0-green)](LICENSE.txt)

A C++20 library implementing sparse multidiagonal operator algebra for quantum systems of arbitrary arity, extracted from the [C++QED](https://github.com/vukics/cppqed) framework.

## What is a multidiagonal operator?

A multidiagonal operator over a Hilbert space of dimension $N$ is a sparse matrix storing an arbitrary set of diagonals at arbitrary offsets. In the context of quantum optics and cavity QED, Hamiltonians and jump operators naturally take this form. For example, the harmonic oscillator ladder operators $a$, $a^\dagger$ are single off-diagonals; number operators $a^\dagger a$ are main diagonals; nonlinear terms like $(a^\dagger)^2 a^2$ are further off-diagonals.

For composite systems of arity $M$, `MultiDiagonal<RANK>` stores tensor-product-structured operators of the form

$$H/i = \bigotimes_{m=0}^{M-1} \sum_{i_m \in \mathbb{K}_m} \sum_{n_m=0}^{N_m-1-i_m} \alpha^{i_m}_{m,n} \left|n_m+i_m\right\rangle\left\langle n_m\right|$$

with $\mathbb{K}_m$ an **arbitrary** set of offsets per subsystem axis.

## Why not `Tridiagonal`?

The predecessor class `Tridiagonal<RANK>` in C++QED v2 stored exactly three diagonals per axis — at positions $-K$, $0$, $+K$ for a single shared offset $K$ — and was **not closed under composition**: the product of two tridiagonal operators is in general pentadiagonal, and with different offsets is not tridiagonal at all. `MultiDiagonal` lifts all these restrictions — any number of diagonals at independent arbitrary offsets, and the operator algebra is closed under composition (`operator|`), direct product (`operator*`), addition, and Hermitian conjugation.

## Features

- **Arbitrary rank**: `MultiDiagonal<RANK>` for any compile-time arity $M$
- **Closed algebra**: composition, direct product, addition, scalar multiplication, Hermitian conjugation
- **Efficient application**: applies as $H/i$ on a state vector via strided multi-index iteration, with compile-time loop unrolling via Boost.Hana
- **Move-only semantics**: no accidental copies of large diagonal arrays
- **JSON serialization**: via Boost.JSON `tag_invoke`, useful for testing and debugging
- **Time dependence**: interaction-picture time dependence ($\alpha_{m,n} e^{\delta_{m,n} t}$) is handled by `InteractionPictureDiagonal` in C++QED proper, keeping this library unconditionally time-independent

## Requirements

- C++20 compiler (GCC ≥ 12, Clang ≥ 15)
- [Boost](https://www.boost.org/) ≥ 1.80 (Hana, JSON, Serialization)
- [Eigen3](https://eigen.tuxfamily.org/) ≥ 3.4
- CMake ≥ 3.24

## Building

```bash
cmake -B build
cmake --build build
```

## Integration into a CMake project

### As a submodule (recommended for C++QED development)

```bash
git submodule add https://github.com/vukics/cppqed-multidiagonal extern/cppqed-multidiagonal
```

In your `CMakeLists.txt`:

```cmake
add_subdirectory(extern/cppqed-multidiagonal)
target_link_libraries(MyTarget PUBLIC cppqed::multidiagonal)
```

### Via `find_package` (after installation)

```bash
cmake --install build
```

```cmake
find_package(cppqed-multidiagonal REQUIRED)
target_link_libraries(MyTarget PUBLIC cppqed::multidiagonal)
```

## Repository structure

```
include/
  Traits.h            — FWD, Boost.Hana/JSON setup, passByValue_v, ReferenceMF
  MultiArray.h        — generic rank-N array with strided view
  MultiArrayComplex.h — dcomp specializations, directProduct, Eigen interop
  MultiDiagonal.h     — multidiagonal operator algebra
src/
  MultiDiagonal.cc    — rank-1 identity and composition (operator|)
cmake/
  cppqed-multidiagonalConfig.cmake.in
```

## License

Distributed under the [Boost Software License, Version 1.0](LICENSE.txt).

## Authors

András Vukics — extracted from C++QED, originally developed 2006–2026.
