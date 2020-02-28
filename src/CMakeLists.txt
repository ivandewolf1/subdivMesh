# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.6 )
project( gradedSubdiv VERSION 0.0.1 )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(LINK_TYPE SHARED CACHE STING "Choose if subti is linked dynamically (Default) or statically -DLINK_TYPE=STATIC")
message(STATUS "link: ${LINK_TYPE}")

add_library( subtri ${_LINK_TYPE}
  subtri.C
  subtri.h
  )
include_directories( ../../lwtv/src ../src)
add_executable(tester tester.C)
