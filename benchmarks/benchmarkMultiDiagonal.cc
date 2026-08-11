// Part of the cppqed-multidiagonal benchmark suite. BSL-1.0.
// Emits CSV to stdout: benchmark,impl,N,d,time_per_op_ns
// Cross-validates every configuration against Eigen before timing (exit 1 on mismatch).
#include "MultiDiagonal.h"

#include <Eigen/SparseCore>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <span>
#include <string_view>
#include <vector>

using namespace cppqedutils;
using quantumoperator::MultiDiagonal;

namespace {

// ── harness ──────────────────────────────────────────────────────────────────

#if defined(__GNUC__) || defined(__clang__)
template <typename T> void doNotOptimize(T const& v) { asm volatile("" : : "r,m"(v) : "memory"); }
#else
volatile const void* dnoSink_;
template <typename T> void doNotOptimize(T const& v) { dnoSink_ = &v; }
#endif

using Clock = std::chrono::steady_clock;

template <typename F>
double batchNs(F&& f, size_t iters)
{
  const auto t0 = Clock::now();
  for (size_t i=0; i<iters; ++i) f();
  return std::chrono::duration<double,std::nano>(Clock::now()-t0).count();
}

/// doubling calibration to ≥ targetBatchNs, then median of `reps` batches
template <typename F>
double medianPerOpNs(F&& f, double targetBatchNs, size_t reps)
{
  size_t iters=1;
  for (double batch=batchNs(f,iters); batch<targetBatchNs; ) { iters*=2; batch=batchNs(f,iters); }
  std::vector<double> t(reps);
  for (auto& v : t) v = batchNs(f,iters)/double(iters);
  std::nth_element(t.begin(), t.begin()+t.size()/2, t.end());
  return t[t.size()/2];
}

void emit(const char* bench, const char* impl, size_t N, size_t d, double ns)
{
  std::printf("%s,%s,%zu,%zu,%.6g\n", bench, impl, N, d, ns);
  std::fflush(stdout);
}

// ── operator construction (identical coefficients on all sides) ─────────────

dcomp coeff(ptrdiff_t o, size_t k) { return { std::cos(0.1*double(k)+double(o)), std::sin(0.05*double(k)-double(o)) }; }

MultiDiagonal<1> single(ptrdiff_t o, std::vector<dcomp> v)
{
  MultiDiagonal<1> md;
  md.diagonals.emplace(MultiDiagonal<1>::Offsets{o}, multiarray::fromStorage(std::move(v)));
  return md;
}

MultiDiagonal<1> makeMD(size_t N, std::span<const ptrdiff_t> offsets)
{
  MultiDiagonal<1> md;
  for (ptrdiff_t o : offsets) {
    std::vector<dcomp> v(N-size_t(std::abs(o)));
    for (size_t k=0; k<v.size(); ++k) v[k]=coeff(o,k);
    md.diagonals.emplace(MultiDiagonal<1>::Offsets{o}, multiarray::fromStorage(std::move(v)));
  }
  return md;
}

template <int Options>
Eigen::SparseMatrix<dcomp,Options> makeSparse(size_t N, std::span<const ptrdiff_t> offsets)
{
  std::vector<Eigen::Triplet<dcomp>> tr;
  for (ptrdiff_t o : offsets) for (size_t k=0; k<N-size_t(std::abs(o)); ++k)
    tr.emplace_back(int(k+size_t(std::max(-o,ptrdiff_t{0}))),
                    int(k+size_t(std::max( o,ptrdiff_t{0}))),
                    coeff(o,k));
  Eigen::SparseMatrix<dcomp,Options> S(N,N);
  S.setFromTriplets(tr.begin(),tr.end());
  S.makeCompressed();
  return S;
}

MultiArray<dcomp,1> randomState(size_t N, unsigned seed)
{
  std::mt19937 gen{seed};
  std::uniform_real_distribution<double> d{-1.,1.};
  std::vector<dcomp> v(N);
  for (auto& x : v) x={d(gen),d(gen)};
  return multiarray::fromStorage(std::move(v));
}

void validate(std::string_view what, size_t N, const CVector& got, const CVector& ref)
{
  const double rel = (got-ref).norm()/ref.norm();
  if (!(rel < 1e-10)) {
    std::fprintf(stderr, "VALIDATION FAILED: %s N=%zu relative deviation %g\n", std::string(what).c_str(), N, rel);
    std::exit(1);
  }
}

std::vector<size_t> sweep(bool smoke)
{
  return smoke ? std::vector<size_t>{16,256}
               : std::vector<size_t>{16,64,256,1024,4096,16384,65536,262144,1048576};
}

// ── benchmark 1: application, rank 1 ────────────────────────────────────────

void applicationBench(bool smoke)
{
  static constexpr std::array<ptrdiff_t,5> offsets{-2,-1,0,1,2};
  const double target = smoke ? 2e6 : 2e7;  const int reps = smoke ? 2 : 5;
  const size_t denseMax = 4096;

  for (size_t N : sweep(smoke)) {
    auto md  = makeMD(N,offsets);
    auto S   = makeSparse<Eigen::RowMajor>(N,offsets);
    auto psi = randomState(N,12345u);
    MultiArray<dcomp,1> dpsidt{{N},multiarray::zeroInit<dcomp>};
    const CVector x = vectorize(psi);           // copy of the same data for Eigen
    CVector y{CVector::Zero(Eigen::Index(N))};

    md(0.,psi,dpsidt);
    y.noalias() = S*x;
    validate("apply1",N,CVector(vectorize(dpsidt)),y);

    emit("apply1","multidiagonal",N,offsets.size(),
         medianPerOpNs([&]{ md(0.,psi,dpsidt); doNotOptimize(dpsidt.dataView.data()); },target,reps));
    emit("apply1","eigen_sparse",N,offsets.size(),
         medianPerOpNs([&]{ y.noalias() += S*x; doNotOptimize(y.data()); },target,reps));

    if (N<=denseMax) {
      CMatrix D{CMatrix::Zero(Eigen::Index(N),Eigen::Index(N))};
      for (ptrdiff_t o : offsets) for (size_t k=0; k<N-size_t(std::abs(o)); ++k)
        D(Eigen::Index(k+size_t(std::max(-o,ptrdiff_t{0}))),Eigen::Index(k+size_t(std::max(o,ptrdiff_t{0}))))=coeff(o,k);
      emit("apply1","eigen_dense",N,offsets.size(),
           medianPerOpNs([&]{ y.noalias() += D*x; doNotOptimize(y.data()); },target,reps));
    }
  }
}

// ── benchmark 2: composition, rank 1 ────────────────────────────────────────

void compositionBench(bool smoke)
{
  static constexpr std::array<ptrdiff_t,3> offsets{-1,0,1};
  const double target = smoke ? 2e6 : 2e7;  const int reps = smoke ? 2 : 5;
  const size_t denseMax = 1024;

  for (size_t N : sweep(smoke)) {
    auto a = makeMD(N,offsets), b = makeMD(N,offsets);
    auto SA = makeSparse<Eigen::ColMajor>(N,offsets); auto SB = SA;
    auto psi = randomState(N,23456u);
    const CVector x = vectorize(psi);

    { // validate (A|B)ψ against A·(B·ψ)
      auto c = a|b;
      MultiArray<dcomp,1> dpsidt{{N},multiarray::zeroInit<dcomp>};
      c(0.,psi,dpsidt);
      const CVector ref = SA*(SB*x).eval();
      validate("compose1",N,CVector(vectorize(dpsidt)),ref);
    }

    emit("compose1","multidiagonal",N,offsets.size(),
         medianPerOpNs([&]{ auto c = a|b; doNotOptimize(c.diagonals.size()); },target,reps));
    emit("compose1","eigen_sparse",N,offsets.size(),
         medianPerOpNs([&]{ Eigen::SparseMatrix<dcomp> C = SA*SB; doNotOptimize(C.nonZeros()); },target,reps));

    if (N<=denseMax) {
      CMatrix DA{CMatrix::Zero(Eigen::Index(N),Eigen::Index(N))};
      for (ptrdiff_t o : offsets) for (size_t k=0; k<N-size_t(std::abs(o)); ++k)
        DA(Eigen::Index(k+size_t(std::max(-o,ptrdiff_t{0}))),Eigen::Index(k+size_t(std::max(o,ptrdiff_t{0}))))=coeff(o,k);
      CMatrix C(N,N);   // preallocated — generous to dense
      emit("compose1","eigen_dense",N,offsets.size(),
           medianPerOpNs([&]{ C.noalias() = DA*DA; doNotOptimize(C.data()); },target,reps));
    }
  }
}

// ── benchmark 3: application, rank-2 Jaynes–Cummings ────────────────────────

MultiDiagonal<1> ladder(size_t N)
{
  std::vector<dcomp> v(N-1);
  for (size_t k=0; k<N-1; ++k) v[k]=std::sqrt(double(k+1));
  return single(+1,std::move(v));
}

MultiDiagonal<1> number(size_t N)
{
  std::vector<dcomp> v(N);
  for (size_t k=0; k<N; ++k) v[k]=double(k);
  return single(0,std::move(v));
}

void jaynesCummingsBench(bool smoke)
{
  const double delta=0.7, g=1.1;  const dcomp eta{0.3,0.2};
  const double target = smoke ? 2e6 : 2e7;  const int reps = smoke ? 2 : 5;

  for (size_t N : sweep(smoke)) {
    // multidiagonal, rank 2: qubit = axis 0 (fast), mode = axis 1 (slow)
    auto jc = g*( single(-1,{1.}) * ladder(N) )                          // σ⁺ ⊗ a
            + g*( single(+1,{1.}) * hermitianConjugateOf(ladder(N)) )    // σ⁻ ⊗ a†
            + quantumoperator::multidiagonal::identity(2)
              * ( delta*number(N) + eta*hermitianConjugateOf(ladder(N)) + conj(eta)*ladder(N) );

    // flattened sparse on dimension 2N, flat index = i_qubit + 2*i_mode
    struct Ent { int r,c; dcomp v; };
    const std::vector<Ent> q_sp{{1,0,1.}}, q_sm{{0,1,1.}}, q_id{{0,0,1.},{1,1,1.}};
    std::vector<Ent> m_a, m_ad, m_n;
    for (size_t k=0; k<N-1; ++k) { m_a.push_back({int(k),int(k+1),std::sqrt(double(k+1))});
                                   m_ad.push_back({int(k+1),int(k),std::sqrt(double(k+1))}); }
    for (size_t k=0; k<N; ++k) m_n.push_back({int(k),int(k),double(k)});

    std::vector<Eigen::Triplet<dcomp>> tr;
    auto addKron = [&] (const std::vector<Ent>& q, const std::vector<Ent>& m, dcomp s) {
      for (const auto& qe : q) for (const auto& me : m)
        tr.emplace_back(qe.r+2*me.r, qe.c+2*me.c, s*qe.v*me.v);
    };
    addKron(q_sp,m_a ,g); addKron(q_sm,m_ad,g);
    addKron(q_id,m_n ,delta); addKron(q_id,m_ad,eta); addKron(q_id,m_a,conj(eta));

    Eigen::SparseMatrix<dcomp,Eigen::RowMajor> S(Eigen::Index(2*N),Eigen::Index(2*N));
    S.setFromTriplets(tr.begin(),tr.end());
    S.makeCompressed();

    std::mt19937 gen{34567u};
    std::uniform_real_distribution<double> dist{-1.,1.};
    MultiArray<dcomp,2> psi{{2,N},[&](size_t e){ std::vector<dcomp> v(e); for (auto& x : v) x={dist(gen),dist(gen)}; return v; }};
    MultiArray<dcomp,2> dpsidt{{2,N},multiarray::zeroInit<dcomp>};
    const CVector x = vectorize(psi);
    CVector y{CVector::Zero(Eigen::Index(2*N))};

    jc(0.,psi,dpsidt);
    y.noalias() = S*x;
    validate("applyJC",N,CVector(vectorize(dpsidt)),y);

    emit("applyJC","multidiagonal",N,jc.diagonals.size(),
         medianPerOpNs([&]{ jc(0.,psi,dpsidt); doNotOptimize(dpsidt.dataView.data()); },target,reps));
    emit("applyJC","eigen_sparse",N,jc.diagonals.size(),
         medianPerOpNs([&]{ y.noalias() += S*x; doNotOptimize(y.data()); },target,reps));
  }
}

} // anonymous


int main(int argc, char** argv)
{
  bool smoke=false;
  for (int i=1; i<argc; ++i) if (std::string_view{argv[i]}=="--smoke") smoke=true;

#ifndef NDEBUG
  std::fprintf(stderr,"WARNING: non-Release build — bounds checks active, timings meaningless.\n");
#endif

  std::printf("benchmark,impl,N,d,time_per_op_ns\n");
  applicationBench(smoke);
  compositionBench(smoke);
  jaynesCummingsBench(smoke);
}