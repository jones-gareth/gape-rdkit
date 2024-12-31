project(Gape)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_VS_PLATFORM_TOOLSET v140)
if (NOT RDBASE)
    set(RDBASE $ENV{RDBASE})
endif()
if (NOT RDLIB) 
    set(RDLIB $ENV{RDLIB})
endif()


MESSAGE(STATUS "RDBASE is  ${RDBASE}")
MESSAGE(STATUS "RDLIB  is  ${RDLIB}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# boost
# if (NOT DEFINED BOOST_ROOT)
# set(Boost_NO_SYSTEM_PATHS ON)
# set(BOOST_ROOT /home/packages/boost)
# endif()

# Windows needs dynamic libs for boost
# see https://github.com/boostorg/program_options/issues/58 for MSVC and boost::program_options issues
set(Boost_USE_STATIC_LIBS  OFF)
find_package(Boost REQUIRED COMPONENTS iostreams filesystem serialization system program_options)

# Eigen

#[[
# use this to include a local eigen
if (NOT DEFINED EIGEN_ROOT)
    set(EIGEN_ROOT /home/packages/eigen-3.3.8)
endif ()
set(EIGEN_INCLUDE_DIR "${EIGEN_ROOT}")
#]]

include(ExternalProject)
set(EIGEN_BUILD_DIR   ${CMAKE_BINARY_DIR}/eigen)
set(EIGEN_INSTALL_DIR ${CMAKE_SOURCE_DIR}/eigen_include)
set(EIGEN_INCLUDE_DIR ${EIGEN_INSTALL_DIR})
ExternalProject_Add(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
        PREFIX ${EIGEN_BUILD_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND
        ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/Eigen ${EIGEN_INCLUDE_DIR}/Eigen
        && ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/unsupported ${EIGEN_INCLUDE_DIR}/unsupported
        TEST_COMMAND ""
)

# include directories
include_directories(
        "src"
        "/usr/local/include"
        ${Boost_INCLUDE_DIRS}
        ${RDBASE}/Code
        ${EIGEN_INCLUDE_DIR}
        ${RDBASE}/External/rapidjson-1.1.0/include
        ${RDBASE}/External/catch/catch/single_include)
link_directories ( ${RDLIB} )

# Set the directories under src which contain cpp files to be compiled
set(srcDirs util ga mol test gape)
set(cSrcDirs '')

# exes to be default compiled
set(allExes gape)

if (MSVC)
    set(LIB_PREFIX "")
    set(OS_LIBS "")
    # see https://github.com/boostorg/program_options/issues/58 for MSVC and boost::program_options issues
    add_definitions( -DBOOST_ALL_NO_LIB )
    # force all boost libraries to dynamic link (we already disabled
    # autolinking, so I don't know why we need this, but we do!)
    add_definitions( -DBOOST_ALL_DYN_LINK )
elseif (UNIX)
    set(LIB_PREFIX "RDKit")
    set(OS_LIBS, "m")
endif()

set(RDKit_LIBS ${LIB_PREFIX}ChemReactions ${LIB_PREFIX}FileParsers ${LIB_PREFIX}SmilesParse ${LIB_PREFIX}Depictor
        ${LIB_PREFIX}RDGeometryLib ${LIB_PREFIX}RDGeneral ${LIB_PREFIX}SubstructMatch ${LIB_PREFIX}Subgraphs
        ${LIB_PREFIX}MolDraw2D ${LIB_PREFIX}GraphMol ${LIB_PREFIX}DistGeometry ${LIB_PREFIX}DistGeomHelpers
        ${LIB_PREFIX}MolAlign ${LIB_PREFIX}Optimizer ${LIB_PREFIX}ForceField ${LIB_PREFIX}ForceFieldHelpers
        ${LIB_PREFIX}Alignment ${LIB_PREFIX}ForceField  ${LIB_PREFIX}MolTransforms ${LIB_PREFIX}EigenSolvers
        ${LIB_PREFIX}GenericGroups ${LIB_PREFIX}coordgen ${LIB_PREFIX}RingDecomposerLib ${LIB_PREFIX}Fingerprints
        ${LIB_PREFIX}DataStructs)

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

# libraries to link against
MESSAGE(STATUS, "Boost libraries are ${Boost_LIBRARIES}")
set( LIBS ${Boost_LIBRARIES}  ${RDKit_THREAD_LIBS}  ${OS_LIBS})


# program files (end in *.cxx)- can be in source directories, or prog directory
set(progDirs prog test)

# find all sources
set(allSources "")
foreach (_dir ${srcDirs})
    set(_pattern "src/${_dir}/*.cpp")
    file(GLOB _srcs ${_pattern})
    #MESSAGE( DBG "_dir is ${_dir} _pattern = ${_pattern} _srcs is ${_srcs}")
    list(APPEND allSources ${_srcs})
endforeach ()
MESSAGE(STATUS "allSources is ${allSources} ")

# find all C sources
set(allCSources ""
        src/gape/FeatureOverlay.cpp
        src/gape/FeatureOverlay.h)
foreach (_dir ${cSrcDirs})
    set(_pattern "src/${_dir}/*.c")
    file(GLOB _srcs ${_pattern})
    #MESSAGE( DBG "_dir is ${_dir} _pattern = ${_pattern} _srcs is ${_srcs}")
    list(APPEND allCSources ${_srcs})
endforeach ()
MESSAGE(STATUS "allCSources is ${allCSources} ")

#[[
# find all include directories
set(allIncludes "")
foreach (_dir ${srcDirs})
    #MESSAGE( DBG "_dir is ${_dir}")
    list(APPEND allIncludes "src" "src/${_dir}")
endforeach ()
MESSAGE(STATUS "allIncludes is  ${allIncludes}")
include_directories(${allIncludes})
#]]

# create library
add_library(GapeLib ${allSources} ${allCSources})
add_dependencies(GapeLib eigen)

# add all programs
add_custom_target(build_all)
MESSAGE(STATUS "progDirs is ${progDirs}")
foreach (_dir ${progDirs})
    set(_pattern "src/${_dir}/*.cxx")
    file(GLOB APP_SOURCES ${_pattern})
    foreach (exesourcefile ${APP_SOURCES})
        # cut off .cxx.
        get_filename_component(exename ${exesourcefile} NAME_WE)
        add_executable(${exename} ${exesourcefile})
        MESSAGE(STATUS " added program " ${exename})
        target_link_libraries(${exename} GapeLib ${RDKit_LIBS} ${LIBS} )
    endforeach (exesourcefile ${APP_SOURCES})
endforeach ()

# Add test programs (in test source directory)
enable_testing()
add_custom_target(build_and_test ${CMAKE_CTEST_COMMAND} -V)

set(_pattern "src/test/*.cxx")
file(GLOB APP_SOURCES ${_pattern})
foreach (exesourcefile ${APP_SOURCES})
    get_filename_component(exename ${exesourcefile} NAME_WE)
    add_test(${exename} ${exename})
    MESSAGE(STATUS " added test " ${exename})
    add_dependencies(build_and_test ${exename})
endforeach (exesourcefile ${APP_SOURCES})

