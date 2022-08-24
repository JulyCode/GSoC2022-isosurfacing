# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.14)
project( Isosurfacing )

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(CGAL REQUIRED)

create_single_source_cgal_program( "Isosurfacing_3/main.cpp" )
create_single_source_cgal_program( "Isosurfacing_3/benchmarks/benchmark.cpp" )

option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

if(ACTIVATE_CONCURRENCY)
  find_package(TBB)
  include(CGAL_TBB_support)
  if(TARGET CGAL::TBB_support)
    target_link_libraries(main PUBLIC CGAL::TBB_support)
    target_link_libraries(benchmark PUBLIC CGAL::TBB_support)
  endif()
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries( main PRIVATE OpenMP::OpenMP_CXX )
  target_link_libraries( benchmark PRIVATE OpenMP::OpenMP_CXX )
endif()

# enable all warnings
if(MSVC)
  target_compile_options(main PRIVATE /W4)
else()
  target_compile_options(main PRIVATE -Wall -Wextra -pedantic)
endif()

add_subdirectory(Isosurfacing_3/examples/Isosurfacing_3)
