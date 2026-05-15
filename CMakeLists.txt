cmake_minimum_required(VERSION 3.24)

project(cppqed-multidiagonal
  VERSION 0.1.0
  DESCRIPTION "Multidiagonal sparse operator algebra for C++QED"
  HOMEPAGE_URL "https://github.com/vukics/cppqed-multidiagonal"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Dependencies ────────────────────────────────────────────────────────────
find_package(Boost REQUIRED COMPONENTS json serialization)
find_package(Eigen3 REQUIRED NO_MODULE)

# Boost.Hana is header-only and part of the Boost distribution;
# it is pulled in transitively via Boost::headers.

# ── Library target ──────────────────────────────────────────────────────────
add_library(cppqed-multidiagonal STATIC
  src/MultiDiagonal.cc)

add_library(cppqed::multidiagonal ALIAS cppqed-multidiagonal)

target_include_directories(cppqed-multidiagonal
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(cppqed-multidiagonal
  PUBLIC
    Boost::headers      # Hana, JSON
    Boost::json
    Boost::serialization  # transitional — will be replaced by HDF5
    Eigen3::Eigen)

target_compile_options(cppqed-multidiagonal
  PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>)

# ── Install ─────────────────────────────────────────────────────────────────
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS cppqed-multidiagonal
  EXPORT cppqed-multidiagonalTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT cppqed-multidiagonalTargets
  FILE cppqed-multidiagonalTargets.cmake
  NAMESPACE cppqed::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cppqed-multidiagonal)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/cppqed-multidiagonalConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  cmake/cppqed-multidiagonalConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cppqed-multidiagonalConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cppqed-multidiagonal)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/cppqed-multidiagonalConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/cppqed-multidiagonalConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cppqed-multidiagonal)
