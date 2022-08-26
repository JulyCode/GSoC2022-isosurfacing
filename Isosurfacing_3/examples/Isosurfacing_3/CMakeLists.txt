# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

#cmake_minimum_required(VERSION 3.1...3.14)
#project( Isosurfacing_examples )

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(CGAL REQUIRED)

option(ACTIVATE_CONCURRENCY "Enable concurrency" ON)

if(ACTIVATE_CONCURRENCY)
  find_package(TBB)
  include(CGAL_TBB_support)

  find_package(OpenMP)
endif()

file(GLOB EXECUTABLES *.cpp)

foreach(exec ${EXECUTABLES})
  get_filename_component(exec_name ${exec} NAME_WE)

  create_single_source_cgal_program( ${exec} )

  if(ACTIVATE_CONCURRENCY)
    if(TARGET CGAL::TBB_support)
      target_link_libraries(${exec_name} PUBLIC CGAL::TBB_support)
    endif()
    if(OpenMP_CXX_FOUND)
      target_link_libraries(${exec_name} PRIVATE OpenMP::OpenMP_CXX )
    endif()
  endif()

  # enable all warnings
  if(MSVC)
    target_compile_options(${exec_name} PRIVATE /W3)
  else()
    target_compile_options(${exec_name} PRIVATE -Wall -Wextra -pedantic)
  endif()

endforeach()
