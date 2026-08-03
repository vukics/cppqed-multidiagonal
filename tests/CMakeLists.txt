include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.7.1)        # pin; bump to the current v3.x at will
FetchContent_MakeAvailable(Catch2)

add_executable(cppqed-multidiagonal-tests
  testMultiArray.cc
#  testSliceIterator.cc
#  testMultiArrayComplex.cc
#  testMultiDiagonal.cc
)

target_link_libraries(cppqed-multidiagonal-tests PRIVATE
  cppqed::multidiagonal          # brings Boost/Eigen transitively (PUBLIC on the lib)
  Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(cppqed-multidiagonal-tests)