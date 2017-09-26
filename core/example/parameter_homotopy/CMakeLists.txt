cmake_minimum_required (VERSION 3.4)

project (parameter_homotopy_example)


IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE debug)
ENDIF()

message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g")



include_directories (include)

set(MY_HEADERS
	include/parameter_homotopy.hpp
	)

set(MY_SOURCES
	src/main.cpp
	)

include(GenerateExportHeader)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

find_library(B2_LIBRARIES 
	NAMES "bertini2"
)

find_library(GMP_LIBRARIES 
	NAMES "gmp"
)

find_library(MPFR_LIBRARIES 
	NAMES "mpfr"
)

#Prep for compiling against boost
find_package(Boost REQUIRED
			COMPONENTS system log)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})



find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${B2_INCLUDE_DIRS})



add_executable(parameter_homotopy_example ${MY_SOURCES})


target_link_libraries (parameter_homotopy_example ${B2_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} Eigen3::Eigen ${Boost_LIBRARIES})


#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltcmalloc")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lprofiler")
