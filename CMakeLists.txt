project(Gape)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)
set(RDBASE $ENV{RDBASE})

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif ()

# Google test on Windows:
# md build
# cd build
# cmake .. -DCMAKE_INSTALL_PREFIX=c:/gareth/src/googletest/install -DCMAKE_BUILD_TYPE=Release
# cmake --build .
# cd lib
# mv Debug Release
# cd ..
# cmake --install .
# then (from RDKit Anaconda build prompt)
# cd C:\gareth\src\gape-rdkit\out\build\x64-debug
# cmake ../../.. -DGTEST_ROOT=/gareth/src/googletest/install
#
#
# need gtest (Download googletest and 'cmake . ; cmake --build .')

if (NOT DEFINED GTEST_ROOT)
    set(GTEST_ROOT /home/packages/googletest)
    set(GTEST_INCLUDE_DIR /home/packages/googletest/googletest/include)
endif ()
find_package(GTest REQUIRED)

# boost
# if (NOT DEFINED BOOST_ROOT)
# set(Boost_NO_SYSTEM_PATHS ON)
# set(BOOST_ROOT /home/packages/boost)
# endif()
find_package(Boost REQUIRED COMPONENTS iostreams filesystem system)

# Eigen
if (NOT DEFINED EIGEN_ROOT)
    set(EIGEN_ROOT /home/packages/eigen-3.3.8)
endif ()
set(EIGEN_INCLUDE_DIR "${EIGEN_ROOT}")

# Set the directories under src which contain cpp files to be compiled
set(srcDirs util ga mol test gape)

set(cSrcDirs '')

# exes to be default compiled
set(allExes TaffForcefieldTest TaffOperations)

# libraries to link against 
set(coreLibs
        ${Boost_IOSTREAMS_LIBRARY}
        ${Boost_IOSTREAMS_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${GTEST_BOTH_LIBRARIES}
        pthread
        m
        dl
        ${RDBASE}/lib/libRDKitRDGeneral.so
        ${RDBASE}/lib/libRDKitGraphMol.so
        ${RDBASE}/lib/libRDKitForceField.so
        ${RDBASE}/lib/libRDKitForceFieldHelpers.so
        ${RDBASE}/lib/libRDKitDistGeomHelpers.so
        ${RDBASE}/lib/libRDKitFileParsers.so
        )

# program files (end in *.cxx)- can be in source directories, or prog directory 
set(progDirs ${srcDirs} prog)

# library paths- add boost
set(LIBS ${LIBS} ${Boost_LIBRARY_DIRS})

# include directories
include_directories(
        "/usr/local/include"
        ${Boost_INCLUDE_DIRS}
        ${GTEST_INCLUDE_DIRS}
        ${RDBASE}/Code
        ${EIGEN_INCLUDE_DIR}
        ${RDBASE}/External/rapidjson-1.1.0/include)
link_directories(${LIBS})

# cflags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC -Wall -g -ggdb -std=c++17 -fmessage-length=0 ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -Wall -std=c++17 -fmessage-length=0 ")

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
set(allCSources "" src/gape/RotatableBond.cpp src/gape/RotatableBond.h src/gape/Gape.cpp src/gape/Gape.h)
foreach (_dir ${cSrcDirs})
    set(_pattern "src/${_dir}/*.c")
    file(GLOB _srcs ${_pattern})
    #MESSAGE( DBG "_dir is ${_dir} _pattern = ${_pattern} _srcs is ${_srcs}")
    list(APPEND allCSources ${_srcs})
endforeach ()
MESSAGE(STATUS "allCSources is ${allCSources} ")

# find all include directories
set(allIncludes "")
foreach (_dir ${srcDirs})
    #MESSAGE( DBG "_dir is ${_dir}")
    list(APPEND allIncludes "src" "src/${_dir}")
endforeach ()
MESSAGE(STATUS "allIncludes is  ${allIncludes}")
include_directories(${allIncludes})

# create library
add_library(GarethLib ${allSources} ${allCSources})

# add all programs
add_custom_target(build_all)
MESSAGE(STATUS "progDirs is ${progDirs}")
foreach (_dir ${progDirs})
    set(_pattern "src/${_dir}/*.cxx")
    file(GLOB APP_SOURCES ${_pattern})
    foreach (exesourcefile ${APP_SOURCES})
        # cut off .cxx.
        get_filename_component(exename ${exesourcefile} NAME_WE)
        #string( REPLACE ".cxx" "" exename ${exesourcefile} )
        list(FIND allExes "${exename}" _index)
        if (${_index} GREATER -1)
            add_executable(${exename} ${exesourcefile})
            MESSAGE(STATUS " added program " ${exename} " to all")
        else ()
            add_executable(${exename} EXCLUDE_FROM_ALL ${exesourcefile})
            MESSAGE(STATUS " added program " ${exename})
        endif ()
        add_dependencies(build_all ${exename})
        target_link_libraries(${exename} GarethLib ${coreLibs})
    endforeach (exesourcefile ${APP_SOURCES})
endforeach ()

# Add gtest programs (in test source directory)
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


# add all programs- manaual specification

# List of programs - each program has a main method and is located in src/prog
#set(programs gape testChrom testGa loadMol sdfParser gravity testEnum conformerMatcher)
#foreach(program ${programs})
#  add_executable(${program} "src/prog/${program}.cxx")
#  target_link_libraries(${program} ${coreLibs} GarethLib)
#  MESSAGE( STATUS " added program " ${program} )
#endforeach(program)

