cmake_minimum_required(VERSION 3.22)
project(pybertini)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(BUILD_SHARED_LIBS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_PHOENIX_STL_TUPLE_H_")


option(ENABLE_UNIT_TESTING "Turn on building and running of unit tests while compiling" OFF)
option(BMP_EXPRESSION_TEMPLATES "Enable Boost.Multiprecision expression templates for mpfr_float, mpfr_complex, mpz_int, mpq_rational" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type set to Release by default, since it was not specified.")
endif(NOT CMAKE_BUILD_TYPE)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

include(CMakePrintHelpers) # so we can print cmake variable values

if(BMP_EXPRESSION_TEMPLATES)
    add_compile_definitions(BMP_EXPRESSION_TEMPLATES=1)
    message(STATUS "BMP expression templates: ON")
else()
    message(STATUS "BMP expression templates: OFF")
endif()

if (WIN32)
    add_compile_options(/bigobj)
    add_compile_options(/EHsc)
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
else ()
    add_compile_options(
        # put things for debug compile here
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Debug>:-Wall>

        # and things for release here
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-g>
        # using O2 or O3 cuts runtime by over half.
	)
endif ()

#sometimes you need to see a lot, and this lets you see ... so much
#include(CMakePrintHelpers)


# note about building using conda-provided tools:
# In order to find conda, run 'conda activate' and then use 'cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX' when cmaking

# All source files to be compiled
# We can either explicitly list all files or use glob, we chose to explicitly list files and not to glob
#file(GLOB SOURCES src/*.cpp)



#Builds a C++ library and the python bindings around it

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MPC REQUIRED)






find_package(eigenpy 3.3 REQUIRED CONFIG)
# Eigenpy finds Boostpython, having our call to find boost libraries before the eigenpy call erases our call, so our call for boost libraries is after this line





# Find python and Boost - both are required dependencies
# so that we can `python_add_library`
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

find_package(Eigen3 REQUIRED NO_MODULE)

find_package(MPI COMPONENTS CXX)
if(MPI_CXX_FOUND)
    message(STATUS "MPI found: ${MPI_CXX_COMPILER}")
endif()

if (PROJECT_IS_TOP_LEVEL)
    find_package(bertini2 REQUIRED)
endif()

set(_PB_BOOST_COMPONENTS
    serialization
    wserialization
    unit_test_framework
    filesystem
    chrono
    regex
    timer
    log
    thread
    log_setup
    python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}
)

find_package(Boost 1.82 REQUIRED COMPONENTS ${_PB_BOOST_COMPONENTS})

# Boost 1.89.0 removed the boost_system library (it's header-only now).
if(Boost_VERSION_STRING VERSION_LESS "1.89.0")
    find_package(Boost REQUIRED COMPONENTS ${_PB_BOOST_COMPONENTS} system)
endif()








#################  done finding things.  now we list our sources


set(PYBERTINI_HEADERS
    include/bertini_python.hpp
    include/eigenpy_interaction.hpp
    include/function_tree_export.hpp
    include/mpfr_export.hpp
    include/random_export.hpp
    include/node_export.hpp
    include/symbol_export.hpp
    include/operator_export.hpp
    include/root_export.hpp
    include/system_export.hpp
    include/tracker_export.hpp
    include/endgame_export.hpp
    include/parser_export.hpp
    include/generic_observer.hpp
    include/generic_observable.hpp
    include/tracker_observers.hpp
    include/endgame_observers.hpp
    include/detail.hpp
    include/logging.hpp
    include/info.hpp
    include/parallel_export.hpp
)

set(PYBERTINI_SOURCES
    src/eigenpy_interaction.cpp
    src/logging.cpp
    src/detail.cpp
    src/containers.cpp
    # tracker: coordinator + config split into separate TU for build-time parallelism
    src/tracker_export.cpp
    src/tracker_config_export.cpp
    # endgame: coordinator + per-tracker-family split
    src/endgame_export.cpp
    src/endgame_config_export.cpp
    src/endgame_double_export.cpp
    src/endgame_mp_export.cpp
    src/endgame_amp_export.cpp
    src/random_export.cpp
    src/mpfr_export.cpp
    src/node_export.cpp
    src/symbol_export.cpp
    src/operator_export.cpp
    src/root_export.cpp
    src/system_export.cpp
    src/parser_export.cpp
    src/generic_observable.cpp
    src/generic_observer.cpp
    src/tracker_observers.cpp
    src/endgame_observers.cpp
    # zero_dim: coordinator + configs + per-tracker-family split
    src/zero_dim_export.cpp
    src/zero_dim_configs_export.cpp
    src/zero_dim_double_export.cpp
    src/zero_dim_mp_export.cpp
    src/zero_dim_amp_export.cpp
    # NID: coordinator + datatypes + per-tracker-family split
    src/numerical_irreducible_decomposition_export.cpp
    src/nid_datatypes_export.cpp
    src/nid_double_export.cpp
    src/nid_mp_export.cpp
    src/nid_amp_export.cpp
    src/bertini_python.cpp
    src/info.cpp
    src/parallel_export.cpp
)



######   set up some stuff about the project
add_project_dependency(Eigen3 REQUIRED PKG_CONFIG_REQUIRES "eigen3 >= 3.4.0")


python_add_library(_pybertini MODULE WITH_SOABI)

target_include_directories(_pybertini PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_include_directories(_pybertini
  SYSTEM PUBLIC
    ${GMP_INCLUDES}
    ${MPFR_INCLUDES}
    ${MPC_INCLUDES}
)

target_sources(_pybertini PRIVATE ${PYBERTINI_SOURCES} ${PYBERTINI_HEADERS})

if(NOT PROJECT_IS_TOP_LEVEL)
    target_link_libraries(_pybertini PUBLIC bertini2)
else()
    target_link_libraries(_pybertini PUBLIC ${Bertini2_LIBRARIES})
endif()

# which libraries are needed
target_link_libraries(_pybertini PUBLIC ${GMP_LIBRARIES})
target_link_libraries(_pybertini PUBLIC ${MPFR_LIBRARIES})
target_link_libraries(_pybertini PUBLIC ${MPC_LIBRARIES})


# which directories contain headers we need.
# NOTE: the SYSTEM keyword must precede the scope keyword (SYSTEM PUBLIC, not
# PUBLIC SYSTEM) -- otherwise CMake reads "SYSTEM" as a relative include path and
# these external headers are NOT treated as system, leaking their warnings into
# our build output. (Eigen comes via the imported Eigen3::Eigen target below,
# which is already treated as SYSTEM.)
target_include_directories(_pybertini SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(_pybertini SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS})
target_include_directories(_pybertini SYSTEM PUBLIC ${Python3_NumPy_INCLUDE_DIRS})
# eigenpy is linked as an imported target below, but its install prefix (e.g.
# /usr/local/include) is not reliably treated as system, so its headers (e.g.
# numpy-allocator.hpp) leak -Wshorten-64-to-32 etc. Mark its include dirs SYSTEM
# explicitly so this third-party dependency's warnings stay out of our output.
target_include_directories(_pybertini SYSTEM PUBLIC
    $<TARGET_PROPERTY:eigenpy::eigenpy,INTERFACE_INCLUDE_DIRECTORIES>)
target_include_directories(_pybertini PUBLIC ${Bertini2_INCLUDES})

if(WIN32)
  target_compile_options(_pybertini
                         PRIVATE $<$<CXX_COMPILER_ID:MSVC>:-bigobj -MP>)
  target_compile_definitions(_pybertini PUBLIC "HAVE_SNPRINTF")
endif()

# which libraries do we need to link to
target_link_libraries(_pybertini PUBLIC Eigen3::Eigen)
target_link_libraries(_pybertini PUBLIC eigenpy::eigenpy)
target_link_libraries(_pybertini PUBLIC ${Boost_LIBRARIES})
target_link_boost_python(_pybertini PUBLIC)

if(MPI_CXX_FOUND)
    target_link_libraries(_pybertini PUBLIC MPI::MPI_CXX)
    target_compile_definitions(_pybertini PUBLIC BERTINI2_HAVE_MPI)
endif()

# eigenpy's own templates (e.g. numpy-allocator.hpp's npy_intp->int) emit
# conversion/shorten warnings when instantiated from our binding code. Marking
# eigenpy's include dir SYSTEM is not enough: clang still reports warnings raised
# during template instantiation requested from non-system code. --system-header-
# prefix=eigenpy/ tells clang to treat any eigenpy/ header as a system header,
# suppressing exactly that third-party noise -- without a blanket -Wno-conversion
# that would also hide genuine conversion bugs in OUR binding code. Clang-only
# flag (GCC/MSVC don't accept it).
target_compile_options(_pybertini PRIVATE
    $<$<CXX_COMPILER_ID:AppleClang,Clang>:--system-header-prefix=eigenpy/>)

# lld-link (Windows) won't pull data-only object files (e.g. explicit_predictors.obj)
# from the static bertini2 archive unless every object is forced in. CMake 3.24+
# LINK_LIBRARY:WHOLE_ARCHIVE places the archive in the linker input list with the
# appropriate /WHOLEARCHIVE: flag (vs. target_link_options which puts it in options only).
# Strip the normal bertini2 link entry first to avoid the archive appearing twice.
# /FORCE:MULTIPLE handles Boost.Serialization duplicate template instantiations.
if(MSVC AND NOT PROJECT_IS_TOP_LEVEL)
  get_target_property(_ll _pybertini LINK_LIBRARIES)
  if(_ll)
    list(REMOVE_ITEM _ll bertini2)
    set_target_properties(_pybertini PROPERTIES LINK_LIBRARIES "${_ll}")
  endif()
  target_link_options(_pybertini PRIVATE /FORCE:MULTIPLE)
  target_link_libraries(_pybertini PRIVATE "$<LINK_LIBRARY:WHOLE_ARCHIVE,bertini2>")
endif()

# Precompile the common root header shared by all 24 python_bindings source files.
# python_common.hpp pulls in Boost.Python (~50 headers), Boost.Operators, Eigen,
# and Bertini's multiprecision headers. Without PCH every TU re-parses that stack
# from scratch, causing peak-RAM spikes during parallel builds and inflating log output.
target_precompile_headers(_pybertini PRIVATE "include/python_common.hpp")

cmake_print_variables(SKBUILD)
if(${SKBUILD})
    # see https://stackoverflow.com/questions/1242904/finding-python-site-packages-directory-with-cmake
    install(TARGETS _pybertini LIBRARY DESTINATION .)
else()
    message(STATUS "Will install python bindings to: ${Python_SITELIB}/bertini")
    install(TARGETS _pybertini LIBRARY DESTINATION ${Python_SITELIB}/bertini)
endif()

