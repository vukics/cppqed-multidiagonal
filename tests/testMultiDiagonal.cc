// Part of the cppqed-multidiagonal test suite. BSL-1.0.
#include "MultiDiagonal.h"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace cppqedutils;
using quantumoperator::MultiDiagonal;

namespace {

constexpr double tol = 1e-12;

bool approxEqual(const CMatrix& A, const CMatrix& B)
{
  return A.rows()==B.rows() && A.cols()==B.cols() && (A-B).norm() < tol;
}

CMatrix kron(const CMatrix& A, const CMatrix& B)   // A = slow factor, B = fast factor
{
  CMatrix res{A.rows()*B.rows(), A.cols()*B.cols()};
  for (Eigen::Index i=0;i<A.rows();++i) for (Eigen::Index j=0;j<A.cols();++j)
    res.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
  return res;
}

/// The oracle: probe with unit vectors; column ordering = the flat storage ordering (axis 0 fastest).
template <size_t RANK>
CMatrix densify(const MultiDiagonal<RANK>& md)
{
  const auto dims = calculateAndCheckDimensions(md);
  const size_t D = multiarray::calculateExtent(dims);
  CMatrix res{CMatrix::Zero(D,D)};
  MultiArray<dcomp,RANK> psi{dims, multiarray::zeroInit<dcomp>}, dpsidt{dims, multiarray::zeroInit<dcomp>};
  for (size_t k=0; k<D; ++k) {
    psi.dataStorage()[k]=1.;
    for (dcomp& v : dpsidt.dataStorage()) v=0.;
    md(0., psi, dpsidt);
    res.col(k) = vectorize(dpsidt);
    psi.dataStorage()[k]=0.;
  }
  return res;
}

MultiDiagonal<1> upperDiag(size_t offset, std::vector<dcomp> v)
{
  MultiDiagonal<1> res;
  const size_t s = v.size();
  res.diagonals[1].emplace(
    MultiDiagonal<1>::Offsets{offset},
    MultiDiagonal<1>::Diagonal{ multiarray::fromStorage(std::move(v)) } );
  return res;
}

MultiDiagonal<1> destroy(size_t dim)                // ⟨n|a|n+1⟩ = √(n+1)
{
  std::vector<dcomp> v(dim-1);
  for (size_t k=0;k<dim-1;++k) v[k]=std::sqrt(double(k+1));
  return upperDiag(1, std::move(v));
}

CMatrix denseDestroy(size_t dim)
{
  CMatrix a{CMatrix::Zero(dim,dim)};
  for (size_t n=0;n+1<dim;++n) a(n,n+1)=std::sqrt(double(n+1));
  return a;
}

} // anonymous


TEST_CASE("MultiDiagonal is move-only", "[multidiagonal]")
{
  STATIC_REQUIRE(!std::is_copy_constructible_v<MultiDiagonal<1>>);
  STATIC_REQUIRE( std::is_move_constructible_v<MultiDiagonal<1>>);
}


TEST_CASE("identity densifies to the identity matrix", "[multidiagonal]")
{
  CHECK(approxEqual(densify(quantumoperator::multidiagonal::identity(5)), CMatrix::Identity(5,5)));
}


TEST_CASE("ladder operator via the densify oracle", "[multidiagonal]")
{
  CHECK(approxEqual(densify(destroy(6)), denseDestroy(6)));
}


TEST_CASE("application accumulates and ignores t", "[multidiagonal]")
{
  auto id = quantumoperator::multidiagonal::identity(3);
  MultiArray<dcomp,1> psi{{3}, multiarray::zeroInit<dcomp>}, d{{3}, multiarray::zeroInit<dcomp>};
  psi(1)=2.+1i;
  id(0., psi, d); id(1., psi, d);
  CHECK(d(1) == 2.*psi(1));
}


#ifndef NDEBUG
TEST_CASE("application checks state dimensions in debug", "[multidiagonal]")
{
  auto a = destroy(4);
  MultiArray<dcomp,1> psi{{5}, multiarray::zeroInit<dcomp>}, d{{5}, multiarray::zeroInit<dcomp>};
  CHECK_THROWS_AS(a(0., psi, d), std::runtime_error);
}
#endif


TEST_CASE("Hermitian conjugation matches the dense adjoint; involution", "[multidiagonal]")
{
  auto a = destroy(6);
  CHECK(approxEqual(densify(hermitianConjugateOf(a)), denseDestroy(6).adjoint()));
  CHECK(approxEqual(densify(hermitianConjugateOf(hermitianConjugateOf(a))), denseDestroy(6)));
}


TEST_CASE("addition, subtraction, scalar algebra, real/imag parts", "[multidiagonal]")
{
  auto a = destroy(5);
  const CMatrix da = denseDestroy(5);

  CHECK(approxEqual(densify(a + hermitianConjugateOf(a)), da + da.adjoint()));   // inserts a new diagonal
  CHECK(approxEqual(densify(a + a), 2.*da));                                     // accumulates on an existing one
  CHECK(densify(a - a).norm() < tol);
  CHECK(approxEqual(densify(2.*a/2.), da));
  CHECK(approxEqual(densify(twoTimesRealPartOf(a)), da + da.adjoint()));
  CHECK(approxEqual(densify(twoTimesImagPartOf(a)), da - da.adjoint()));
}


TEST_CASE("copy is deep", "[multidiagonal]")
{
  auto a = destroy(4);
  auto b = copy(a);
  b *= 3.;
  CHECK(approxEqual(densify(a), denseDestroy(4)));
}


TEST_CASE("calculateAndCheckDimensions", "[multidiagonal]")
{
  CHECK((calculateAndCheckDimensions(MultiDiagonal<1>{}) == Extents<1>{}));
  CHECK((calculateAndCheckDimensions(destroy(7)) == Extents<1>{7}));

  MultiDiagonal<1> bad;
  {
    bad.diagonals[1].emplace(MultiDiagonal<1>::Offsets{1},
      MultiDiagonal<1>::Diagonal{ multiarray::fromStorage( std::vector<dcomp>(3,1.) ) } );
  }
  {
    bad.diagonals[0].emplace(MultiDiagonal<1>::Offsets{1},
      MultiDiagonal<1>::Diagonal{ multiarray::fromStorage( std::vector<dcomp>(4,1.) ) } );
  }
  CHECK_THROWS_AS(calculateAndCheckDimensions(bad), std::runtime_error);
}


// ── composition: branches that are correct today are STRICT ──────────────────

TEST_CASE("composition: number operator a†|a (lower·upper)", "[multidiagonal][composition]")
{
  auto a = destroy(6);
  CHECK(approxEqual(densify(hermitianConjugateOf(a) | a), denseDestroy(6).adjoint()*denseDestroy(6)));
}


TEST_CASE("composition: like directions (upper·upper, lower·lower)", "[multidiagonal][composition]")
{
  auto a = destroy(6); auto ad = hermitianConjugateOf(a);
  const CMatrix da = denseDestroy(6);
  CHECK(approxEqual(densify(a | a),   (da*da).eval()));
  CHECK(approxEqual(densify(ad | ad), (da.adjoint()*da.adjoint()).eval()));
}


TEST_CASE("composition: identity is two-sided neutral on all branches", "[multidiagonal][composition]")
{
  auto a = destroy(5); auto ad = hermitianConjugateOf(a);
  auto id = quantumoperator::multidiagonal::identity(5);
  const CMatrix da = denseDestroy(5);
  CHECK(approxEqual(densify(id | a),  da));
  CHECK(approxEqual(densify(a | id),  da));
  CHECK(approxEqual(densify(id | ad), da.adjoint()));
  CHECK(approxEqual(densify(ad | id), da.adjoint()));
}


TEST_CASE("composition associativity on (a†,a,a)", "[multidiagonal][composition]")
{
  auto a = destroy(6); auto ad = hermitianConjugateOf(a);
  CHECK(approxEqual(densify((ad|a)|a), densify(ad|(a|a))));
}


// ── composition: the upper·lower branch — CONTRACT for the signed-offset refactor (finding #1)

TEST_CASE("composition: a|a† keeps the full Hilbert dimension — finding #1", "[multidiagonal][composition][contract][!mayfail]")
{
  auto a = destroy(4);
  auto aad = a | hermitianConjugateOf(a);
  CHECK((calculateAndCheckDimensions(aad) == Extents<1>{4}));   // currently 3: main diagonal loses its trailing zero
  CHECK(approxEqual(densify(aad), (denseDestroy(4)*denseDestroy(4).adjoint()).eval()));
}


TEST_CASE("commutator [a,a†] == 1 — finding #1", "[multidiagonal][composition][contract][!mayfail]")
{
  auto a = destroy(4); auto ad = hermitianConjugateOf(a);
  auto comm = (a|ad) - (ad|a);                                  // currently a dimension-mismatch in debug
  CHECK(approxEqual(densify(comm), CMatrix::Identity(4,4)));
}


TEST_CASE("(A|B)† == B†|A† on the n̂ witness", "[multidiagonal][composition]")
{
  auto a = destroy(5);
  auto n = hermitianConjugateOf(a) | a;                          // Hermitian
  CHECK(approxEqual(densify(hermitianConjugateOf(n)), densify(hermitianConjugateOf(a) | a)));
}


// ── direct product ───────────────────────────────────────────────────────────

TEST_CASE("direct product densifies to a Kronecker product (axis 0 fastest)", "[multidiagonal][directproduct]")
{
  auto a  = destroy(3);
  auto id = quantumoperator::multidiagonal::identity(2);
  CHECK(approxEqual(densify(a * id), kron(CMatrix::Identity(2,2), denseDestroy(3))));  // a on axis 0 (fast)
  CHECK(approxEqual(densify(id * a), kron(denseDestroy(3), CMatrix::Identity(2,2))));  // a on axis 1 (slow)
}


TEST_CASE("Jaynes–Cummings coupling σ⁺⊗a + σ⁻⊗a† is Hermitian and correct", "[multidiagonal][directproduct]")
{
  const size_t dim=4;
  auto sminus = upperDiag(1, {1.});                              // ⟨0|σ⁻|1⟩ = 1
  auto a      = destroy(dim);
  auto jc = hermitianConjugateOf(sminus)*a + sminus*hermitianConjugateOf(a);

  CMatrix ds{CMatrix::Zero(2,2)}; ds(0,1)=1.;                    // dense σ⁻
  const CMatrix da = denseDestroy(dim);
  const CMatrix ref = kron(da, ds.adjoint()) + kron(da.adjoint(), ds);  // qubit fast (axis 0), mode slow (axis 1)

  auto dense = densify(jc);
  CHECK(approxEqual(dense, dense.adjoint().eval()));
  CHECK(approxEqual(dense, ref));
}


TEST_CASE("JSON serialization of MultiDiagonal — finding #2", "[multidiagonal][json][contract][!mayfail]")
{
  auto a = destroy(3);
  json::value jv;
  CHECK_NOTHROW(jv = json::value_from(a));                       // currently as_object() on a null value throws
}