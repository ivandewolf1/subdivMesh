# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.6 )
list( APPEND CMAKE_PREFIX_PATH "$ENV{HFS}/toolkit/cmake" )
project( gradedSubdiv VERSION 0.0.1 )
# Locate Houdini's libraries and header files.
# Registers an imported library target named 'Houdini'.
find_package( Houdini REQUIRED )
# Add a library with two source files.
set( houSrc houdini/ )
set( triSrc subTri/ )
set( subtriBasicSop SOP_SubtriBasic )
add_library( ${subtriBasicSop} SHARED
  ${houSrc}SOP_SubtriBasic.C
  ${houSrc}SOP_SubtriBasic.h
  )
set( CMAKE_CXX_FLAGS -FPIC )
add_library( subtri SHARED
  ${triSrc}subtri.C
  ${triSrc}subtri.h
  )
include_directories( lwtv )
# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${subtriBasicSop} PUBLIC subtri Houdini )
target_include_directories( ${subtriBasicSop} PUBLIC ${triSrc} )
# Configure several common target properties, such as its output directory.
houdini_configure_target( ${subtriBasicSop} )
