# CMakeLists.txt to build the NBodylib library
#
# ICRAR - International Centre for Radio Astronomy Research
# (c) UWA - The University of Western Australia, 2018
# Copyright by UWA (in the framework of the ICRAR)
# All rights reserved
#
# Contributed by Rodrigo Tobar
#
# This file is part of VELOCIraptor.

cmake_minimum_required(VERSION 3.10)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(_subdirs Analysis Cosmology InitCond KDTree Math NBody)

# We have the version here
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/VERSION NBODY_VERSION)
string(REPLACE "." ";" VERSION_LIST ${NBODY_VERSION})
list(GET VERSION_LIST 0 NBODY_VERSION_MAJOR)
list(GET VERSION_LIST 1 NBODY_VERSION_MINOR)

# We have only C++ sources, but until cmake 3.4 the check_symbol_exists
# macro didn't work properly if only defining C++ as the project language
if( ${CMAKE_MAJOR_VERSION} EQUAL 3 AND ${CMAKE_MINOR_VERSION} LESS 4 )
	set(NBODY_LANGS C CXX)
else()
	set(NBODY_LANGS CXX)
endif()

project(NBodyLib VERSION ${NBODY_VERSION} LANGUAGES ${NBODY_LANGS})
set(CMAKE_CXX_STANDARD 17)

#include cmake stuff
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
# include some general cmake utility functions
include (NBodylibCMakeUtils)
# complation summary
include(NBodylibCompilationMessages)

include(CTest)
enable_testing()


# Are we part of an outside build?
# If so, we want to let the caller know what to include, define, etc at the end
set(_export OFF)
get_directory_property(_hasParent PARENT_DIRECTORY)
if (_hasParent)
	set(_export ON)
endif()

# Options users can give on the command line via -D
macro(nbody_option optname optdesc status )
	if (NOT DEFINED NBODY_${optname})
		option(NBODY_${optname} "${optdesc}" "${status}")
	endif()
endmacro()

# Optional dependencies
nbody_option(OPENMP "Attempt to include OpenMP support in NBodyLib" ON)
nbody_option(ENABLE_OPENMP_GPU "Attempt to include OpenMP GPU support in NBodyLib" OFF)
nbody_option(OPENACC "Attempt to include OpenACC support in NBodyLib" OFF)

# Precision options
nbody_option(SINGLE_PRECISION  "Use single point precision to store all properties and perform all calculations" OFF)
nbody_option(LONG_INT          "Use long ints to represent all integers. Needed if dealing with more than MAXINT number of particles" OFF)

# Particle class details options
nbody_option(NO_MASS                    "Do not store the mass as all particles are the same mass" OFF)
nbody_option(SINGLE_PARTICLE_PRECISION  "Do not store the mass as all particles are the same mass" OFF)
nbody_option(UNSIGNED_PARTICLE_PIDS     "Use unsigned particle PIDs" OFF)
nbody_option(UNSIGNED_PARTICLE_IDS      "Use unsigned particle IDs" OFF)
nbody_option(USE_GAS                    "Particle class has gas" OFF)
nbody_option(USE_STARS                  "Particle class has stars" OFF)
nbody_option(USE_BARYONS                "Gas and Stars" OFF)
nbody_option(USE_HYDRO                  "All particle types (gas, stars, black holes, etc)" OFF)
nbody_option(USE_EXTRA_DM_PROPERTIES    "DM particle have extra properties" OFF)
nbody_option(USE_SWIFT_INTERFACE		"SWIFT Interface" OFF)
nbody_option(USE_EXTRA_INPUT_INFO		"Extra Input Info" OFF)
nbody_option(USE_EXTRA_FOF_INFO		    "Extra FOF Info" OFF)
nbody_option(USE_LARGE_KDTREE		    "Require large mem KDTree as there are more than > max 32-bit integer entries" OFF)

# Process high-level options
if (NBODY_USE_BARYONS)
	set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
endif()
if (NBODY_USE_HYDRO)
	set(NBODY_USE_GAS ON)
	set(NBODY_USE_STARS ON)
	set(NBODY_USE_BH ON)
endif()


# Flags used for compiling nbodylib.
# They are optionally exported at the end if there is an outer scope
set(NBODY_INCLUDE_DIRS "")
set(NBODY_DEFINES "")
set(NBODY_LIBS "")
set(NBODY_CXX_FLAGS "")
set(NBODY_LINK_FLAGS "")
set(NBODY_DOC_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/doc)

# Define macros depending on the user's choices
macro(nbody_option_defines varname define)
	if (NBODY_${varname})
		list(APPEND NBODY_DEFINES ${define})
		set(NBODY_HAS_${varname} Yes)
	else()
		set(NBODY_HAS_${varname} No)
	endif()
endmacro()

nbody_option_defines(SINGLE_PRECISION          SINGLEPRECISION)
nbody_option_defines(LONG_INT                  LONGINT)
nbody_option_defines(NO_MASS                   NOMASS)
nbody_option_defines(SINGLE_PARTICLE_PRECISION LOWPRECISIONPOS)
nbody_option_defines(UNSIGNED_PARTICLE_PIDS    PARTICLEUPIDS)
nbody_option_defines(UNSIGNED_PARTICLE_IDS     PARTICLEUIDS)
nbody_option_defines(LARGEKDTREE               LARGETREE)
nbody_option_defines(USE_GAS                   GASON)
nbody_option_defines(USE_STARS                 STARON)
nbody_option_defines(USE_BH                    BHON)
nbody_option_defines(USE_EXTRA_DM_PROPERTIES   EXTRADMON)
nbody_option_defines(USE_SWIFT_INTERFACE       SWIFTINTERFACE)
nbody_option_defines(USE_EXTRA_INPUT_INFO      EXTRAINPUTINFO)
nbody_option_defines(USE_EXTRA_FOF_INFO        EXTRAFOFINFO)
nbody_option_defines(USE_LARGE_KDTREE          LARGETREE)

#
# How we find GSL and set it up
#
macro(find_gsl)
	find_package(GSL REQUIRED)
	if( ${GSL_VERSION} GREATER 2.2 OR ${GSL_VERSION} EQUAL 2.2)
		list(APPEND NBODY_DEFINES HAVE_GSL22)
	endif()
	list(APPEND NBODY_INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
	list(APPEND NBODY_LIBS ${GSL_LIBRARIES})
endmacro()

# How we find OpenMP and set it up
#
macro(find_openmp)
	find_package(OpenMP COMPONENTS CXX)
	if (OPENMP_FOUND)
		list(APPEND NBODY_CXX_FLAGS ${OpenMP_CXX_FLAGS})
		list(APPEND NBODY_LINK_FLAGS ${OpenMP_CXX_FLAGS})
		list(APPEND NBODY_DEFINES USEOPENMP USEOMP)
		set(NBODY_HAS_OPENMP Yes)
		if (NBODY_ENABLE_OPENMP_GPU)
			if (${OpenMP_CXX_VERSION} GREATER 4.0 OR {OpenMP_CXX_VERSION} EQUAL 4.0)
				set(NBODY_HAS_OPENMPGPU Yes)
				#list(APPEND NBODY_CXX_FLAGS -foffload=-lm)
				#list(APPEND NBODY_LINK_FLAGS -foffload=-lm)
				list(APPEND NBODY_DEFINES USEOPENMPTARGET)
			endif()
		endif()
	endif()
endmacro()

macro(find_openacc)
	find_package(OpenACC COMPONENTS CXX)
	if (OPENMP_FOUND AND NBODY_ENABLE_OPENACC)
		list(APPEND NBODY_CXX_FLAGS ${OpenACC_CXXX_FLAGS})
		list(APPEND NBODY_LINK_FLAGS ${OpenACC_CXX_FLAGS})
		list(APPEND NBODY_DEFINES USEOPENACC)
		set(NBODY_HAS_OPENACC Yes)
	endif()
endmacro()

#
# Now go and find our dependencies, depending on whether
# we actually need them or not
#
find_gsl()

set(NBODY_HAS_OPENMP No)
set(NBODY_HAS_OPENMPGPU No)
set(NBODY_HAS_OPENACC No)
if (NBODY_OPENMP)
	find_openmp()
endif()
if (NBODY_OPENACC)
	find_openacc()
endif()


# Add subdirectories to list of include directories
foreach(_subdir IN ITEMS ${_subdirs})
	list(APPEND NBODY_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/src/${_subdir}")
endforeach()

# compilation state
nbody_compilation_summary()


#
# Macro used by subdirectories to build each a static library
#
macro(add_nbody_lib libname sources)
	add_library(${libname} OBJECT ${sources})
	target_compile_definitions(${libname} PRIVATE ${NBODY_DEFINES})
	target_include_directories(${libname} PRIVATE ${NBODY_INCLUDE_DIRS})
	if (NBODY_CXX_FLAGS)
		set_target_properties(${libname} PROPERTIES COMPILE_FLAGS ${NBODY_CXX_FLAGS})
	endif()
	sanitizer_options(${libname})
endmacro()

# Include all subdirectories,
# where the magic happens
foreach(_subdir IN ITEMS ${_subdirs})
	add_subdirectory(src/${_subdir})
endforeach()

# The nbodylib OBJECT library, which can then be used in the same fashion
# to bring all objects together into a final executable/library
add_library(nbodylib_iface INTERFACE)
target_sources(nbodylib_iface INTERFACE
    $<TARGET_OBJECTS:Analysis> $<TARGET_OBJECTS:Cosmology>
    $<TARGET_OBJECTS:InitCond> $<TARGET_OBJECTS:KDTree>
    $<TARGET_OBJECTS:Math> $<TARGET_OBJECTS:NBody>
)
add_library(NBodylib::nbodylib_iface ALIAS nbodylib_iface)

# The final static libnbodylib.a library, in case it's needed on its own
add_library(nbodylib STATIC)
target_link_libraries(nbodylib nbodylib_iface)
add_library(NBodylib::nbodylib ALIAS nbodylib)
sanitizer_options(nbodylib)

if (BUILD_TESTING)
	add_subdirectory(src/tests)
endif()


# Export the include directories and definitions, if necessary
# If building on our own, add the "doc" target
if (_export)
	set(NBODYLIB_VERSION "${NBODY_VERSION_MAJOR}.${NBODY_VERSION_MINOR}" PARENT_SCOPE)
	set(NBODYLIB_VERSION_MAJOR "${NBODY_VERSION_MAJOR}" PARENT_SCOPE)
	set(NBODYLIB_VERSION_MINOR "${NBODY_VERSION_MINOR}" PARENT_SCOPE)
	set(NBODYLIB_INCLUDE_DIRS "${NBODY_INCLUDE_DIRS}" PARENT_SCOPE)
	set(NBODYLIB_DEFINES "${NBODY_DEFINES}" PARENT_SCOPE)
	set(NBODYLIB_HAS_OPENMP "${NBODY_HAS_OPENMP}" PARENT_SCOPE)
	set(NBODYLIB_HAS_OPENMPGPU "${NBODY_HAS_OPENMPGPU}" PARENT_SCOPE)
	set(NBODYLIB_HAS_OPENACC "${NBODY_HAS_OPENACC}" PARENT_SCOPE)
	set(NBODYLIB_CXX_FLAGS "${NBODY_CXX_FLAGS}" PARENT_SCOPE)
	set(NBODYLIB_LINK_FLAGS "${NBODY_LINK_FLAGS}" PARENT_SCOPE)
	set(NBODYLIB_LIBS "${NBODY_LIBS}" PARENT_SCOPE)
	set(NBODYLIB_DOC_DIRS "${NBODY_DOC_DIRS}" PARENT_SCOPE)
else()
	try_add_doc_target("${NBODY_DOC_DIRS}")
endif()
