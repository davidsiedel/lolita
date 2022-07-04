cmake_minimum_required(VERSION 3.19)
project(lolita
        VERSION 0.0.1
        DESCRIPTION "C++ implementation of the Hybrid High Order method"
        HOMEPAGE_URL "https://github.com/davidsiedel/feta"
        LANGUAGES CXX)

# Additional modules path for cmake
#set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}/cmake)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CheckLanguage)

set(CMAKE_CXX_STANDARD 23)
set(CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=reference_point++2a -fsanitize=address")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

######################################################################
## Find and use OpemMP
#OPTION (USE_OpenMP "Use OpenMP" ON)
#IF(USE_OpenMP)
#    FIND_PACKAGE(OpenMP)
#    IF(OPENMP_FOUND)
#        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    ENDIF()
#ENDIF()
#
#if ( WIN32 )
#    message(FATAL_ERROR "Sorry, Windows is not supported.")
#endif()

######################################################################
## Helpers for compiler identification
#if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#    set(COMPILER_IS_CLANG TRUE)
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#    set(COMPILER_IS_GNU TRUE)
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
#    set(COMPILER_IS_INTEL)
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
#    set(COMPILER_IS_PGI)
#endif ()
#
#if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
#    cmake_host_system_information(RESULT OSVER QUERY OS_RELEASE)
#    if(OSVER VERSION_LESS "10.14")
#        message(FATAL_ERROR "DiSk++ depends on libraries that use C++17 \
#                 features available starting with Mac OS X 10.14. You should \
#                 either upgrade your system or use a compiler different than \
#                 AppleClang.")
#    endif()
#endif()

######################################################################
## Package: Threading support
#set(THREADS_PREFER_PTHREAD_FLAG ON)
#find_package(Threads REQUIRED)
#if (COMPILER_IS_CLANG OR COMPILER_IS_INTEL OR COMPILER_IS_GNU)
#    set(LINK_LIBS ${LINK_LIBS} Threads::Threads)
#endif()

######################################################################
## Package: Intel MKL support
find_package(MKL REQUIRED)
if (MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIRS})
    set(INCL_DIRS ${INCL_DIRS} ${MKL_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${MKL_LIBRARIES})
    set(HAVE_INTEL_MKL TRUE)
    add_definitions(-DHAVE_INTEL_MKL)
endif()

######################################################################
## Package: Eigen3
find_package(Eigen3 REQUIRED)
if (Eigen3_FOUND)
    message(STATUS "Eigen3 found")
    set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)
    set(INCL_DIRS ${INCL_DIRS} ${EIGEN3_INCLUDE_DIR})
    include_directories(${EIGEN3_INCLUDE_DIR})
    set(HAVE_EIGEN3 TRUE)
    add_definitions(-DHAVE_EIGEN3)
else()
    message(STATUS "Eigen3 not found")
    set(HAVE_EIGEN3 FALSE)
endif()

######################################################################
## Package: Boost
#find_package(Boost REQUIRED COMPONENTS python38 numpy38)
#if (Boost_FOUND)
#    message(STATUS "Boost found")
#    set(LINK_LIBS ${LINK_LIBS} ${Boost_LIBRARIES})
#    set(INCL_DIRS ${INCL_DIRS} ${Boost_INCLUDE_DIRS})
#    set(HAVE_BOOST TRUE)
#    add_definitions(-DHAVE_BOOST)
#else()
#    message(STATUS "Boost not found")
#    set(HAVE_BOOST FALSE)
#endif()

######################################################################
## Package: PythonLibs
#find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
#if (Python3_FOUND)
#    message(STATUS "Python3 and Numpy found")
#    set(LINK_LIBS ${LINK_LIBS} ${Python3_LIBRARIES})
#    set(LINK_LIBS ${LINK_LIBS} Python3::NumPy)
#    set(INCL_DIRS ${INCL_DIRS} ${Python3_INCLUDE_DIRS})
#    set(HAVE_PYTHON3 TRUE)
#    add_definitions(-DHAVE_HAVE_PYTHON3)
#else()
#    message(STATUS "Python3 not found")
#    set(HAVE_PYTHON3 FALSE)
#endif()


######################################################################
## Package: MGIS
#find_package(MFrontGenericInterface REQUIRED)
#if (MFrontGenericInterface_FOUND)
#    message(STATUS "MGIS found")
#    #    include_directories(${MGIS_INCLUDE_DIRS})
#    set(LINK_LIBS ${LINK_LIBS} mgis::MFrontGenericInterface)
#    set(HAVE_MGIS TRUE)
#    add_definitions(-DHAVE_MGIS)
#else()
#    message(STATUS "MGIS not found")
#    set(HAVE_MGIS FALSE)
#endif()
find_package(MGIS REQUIRED)
if (MGIS_FOUND)
    message(STATUS "MGIS found")
    include_directories(${MGIS_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${MGIS_LIBRARIES})
    set(HAVE_MGIS TRUE)
    add_definitions(-DHAVE_MGIS)
else()
    message(STATUS "MGIS not found")
endif()

######################################################################
## Package: fmt
#find_package(fmt REQUIRED)
#if (fmt_FOUND)
#    message(STATUS "FMT found")
#    #    include_directories(${MGIS_INCLUDE_DIRS})
#    set(LINK_LIBS ${LINK_LIBS} fmt::fmt)
#    set(HAVE_FMT TRUE)
#    add_definitions(-DHAVE_FMT)
#else()
#    message(STATUS "FMT not found")
#    set(HAVE_FMT FALSE)
#endif()

######################################################################
## Optimization: Prefer 512 bit vectors on certain architectures.
## See https://reviews.llvm.org/D67259 for more info
#option(OPT_PREFER_512bit "Prefer 512 bit vectors with AVX512 (Clang > 10 & GCC)" OFF)
#if (OPT_PREFER_512bit)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mprefer-vector-width=512")
#endif()

######################################################################
## Optimization: Disallow denormals and activate -ffast-math (violates IEEE754)
#option(OPT_AGGRESSIVE_FP "Enable DAZ, FTZ and -ffast-math" ON)
#if (OPT_AGGRESSIVE_FP)
#    add_definitions(-DDISALLOW_DENORMALS)
#    if (COMPILER_IS_CLANG OR COMPILER_IS_GNU)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
#    endif()
#endif()

######################################################################
## Optimization: Enable vectorizer output
#option(OPT_VECTORIZER_REMARKS "Enable vectorizer remarks" OFF)
#if (OPT_VECTORIZER_REMARKS)
#    if (COMPILER_IS_CLANG)
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass-missed=loop-vectorize")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass-analysis=loop-vectorize")
#    endif()
#
#    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report-phase=vec -qopt-report=2")
#    endif()
#
#    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized")
#    endif()
#
#    if (CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minfo")
#    endif()
#endif()

######################################################################
## Compiler settings: Settings for different build types
#set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address -fpermissive")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -g -fpermissive -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASEASSERT "-O3 -march=native -g -fpermissive")
#
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

# add_subdirectory(src)


######################################################################
## Google tests
if (TEST_MODE STREQUAL "test")
    add_subdirectory(tests)
endif()

add_subdirectory(applications)