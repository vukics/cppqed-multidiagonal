option(MULTIDIAGONAL_BENCHMARK_NATIVE "Compile benchmarks with -march=native" ON)

add_executable(cppqed-multidiagonal-benchmarks benchmarkMultiDiagonal.cc)

target_link_libraries(cppqed-multidiagonal-benchmarks PRIVATE cppqed::multidiagonal)

target_compile_options(cppqed-multidiagonal-benchmarks PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
  $<$<AND:$<BOOL:${MULTIDIAGONAL_BENCHMARK_NATIVE}>,$<CXX_COMPILER_ID:GNU,Clang>>:-march=native>)