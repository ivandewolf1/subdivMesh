# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.6 )
project( gradedSubdiv VERSION 0.0.1 )
add_library( subtri SHARED
  subtri.C
  subtri.h
  )
include_directories( ../../lwtv/src ../src)
add_executable(tester tester.C)
