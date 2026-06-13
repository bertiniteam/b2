cmake_minimum_required(VERSION 3.22)
project(bertini2-core)

include(CMakePrintHelpers) # so we can print cmake variable values
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# set(BUILD_SHARED_LIBS ON)
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# In order to find conda, run 'conda activate' and then use 'cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX' when cmaking

# All source files to be compiled
# We can either explicitly list all files or use glob, we chose to explicitly list files and not to glob
#file(GLOB SOURCES src/*.cpp)

option(ENABLE_UNIT_TESTING "Turn on building and running of unit tests while compiling" OFF)
option(BMP_EXPRESSION_TEMPLATES "Enable Boost.Multiprecision expression templates for mpfr_float, mpfr_complex, mpz_int, mpq_rational" ON)

find_package(MPI)
option(ENABLE_MPI "Enable MPI parallelism (requires OpenMPI or MPICH; no Boost.MPI needed)" ${MPI_FOUND})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type set to Release by default, since it was not specified.")
endif(NOT CMAKE_BUILD_TYPE)


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "ENABLE_UNIT_TESTING: ${ENABLE_UNIT_TESTING}")

if(BMP_EXPRESSION_TEMPLATES)
    add_compile_definitions(BMP_EXPRESSION_TEMPLATES=1)
    message(STATUS "BMP expression templates: ON")
else()
    message(STATUS "BMP expression templates: OFF")
endif()

if (WIN32)
    add_compile_options(/bigobj)
    add_compile_options(/EHsc)
    set(CMAKE_BUILD_TYPE Release)
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


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MPC REQUIRED)

# include_directories(${GMP_INCLUDES})
# include_directories(${MPC_INCLUDES})

# set_boost_default_options()
# export_boost_default_options()

set(_BOOST_REQUIRED_COMPONENTS
    serialization
    unit_test_framework
    filesystem
    graph
    chrono
    regex
    timer
    log
    log_setup
    thread
)

if(ENABLE_MPI AND NOT MPI_FOUND)
    message(FATAL_ERROR "ENABLE_MPI=ON but no MPI installation was found.")
endif()
if(ENABLE_MPI)
    message(STATUS "MPI enabled: ${MPI_CXX_COMPILER}")
endif()

find_package(Boost REQUIRED COMPONENTS ${_BOOST_REQUIRED_COMPONENTS})

# Boost 1.89.0 removed the boost_system library (it's header-only now).
# Re-find with ALL components so Boost_LIBRARIES is not overwritten.
if(Boost_VERSION_STRING VERSION_LESS "1.89.0")
    find_package(Boost REQUIRED COMPONENTS ${_BOOST_REQUIRED_COMPONENTS} system)
endif()

find_package(Eigen3 REQUIRED NO_MODULE)

# This was from a makemodule.am, take the lists and define them in seperate blocks with the variables they already have, then add the variable name to library_headrs
set(endgames_headers
	include/bertini2/endgames/amp_endgame.hpp
	include/bertini2/endgames/base_endgame.hpp
	include/bertini2/endgames/cauchy.hpp
	include/bertini2/endgames/config.hpp
	include/bertini2/endgames/events.hpp
	include/bertini2/endgames/fixed_prec_endgame.hpp
	include/bertini2/endgames/interpolation.hpp
	include/bertini2/endgames/observers.hpp
	include/bertini2/endgames/powerseries.hpp
	include/bertini2/endgames/prec_base.hpp
)

set(basics_rootinclude_headers
    include/bertini2/have_bertini.hpp
    include/bertini2/mpfr_extensions.hpp
    include/bertini2/mpfr_complex.hpp
    include/bertini2/forbid_mixed_arithmetic.hpp
    include/bertini2/double_extensions.hpp
    include/bertini2/random.hpp
    include/bertini2/num_traits.hpp
    include/bertini2/classic.hpp
    include/bertini2/eigen_extensions.hpp
    include/bertini2/eigen_serialization_addon.hpp
    include/bertini2/logging.hpp
    include/bertini2/endgames.hpp
    include/bertini2/bertini.hpp
)

set(common_headers
    include/bertini2/common/config.hpp
    include/bertini2/common/stream_enum.hpp
)

set(detail_headers
    include/bertini2/detail/configured.hpp
    include/bertini2/detail/events.hpp
    include/bertini2/detail/visitable.hpp
    include/bertini2/detail/visitor.hpp
    include/bertini2/detail/observer.hpp
    include/bertini2/detail/observable.hpp
    include/bertini2/detail/is_template_parameter.hpp
    include/bertini2/detail/enable_permuted_arguments.hpp
    include/bertini2/detail/typelist.hpp
)

set(function_tree_headers
    include/bertini2/function_tree.hpp
    include/bertini2/function_tree/node.hpp
    include/bertini2/function_tree/gather.hpp
    include/bertini2/function_tree/forward_declares.hpp
    include/bertini2/function_tree/simplify.hpp
    include/bertini2/function_tree/operators/operator.hpp
    include/bertini2/function_tree/symbols/symbol.hpp
    include/bertini2/function_tree/symbols/variable.hpp
    include/bertini2/function_tree/symbols/differential.hpp
    include/bertini2/function_tree/symbols/special_number.hpp
    include/bertini2/function_tree/symbols/number.hpp
    include/bertini2/function_tree/symbols/linear_product.hpp
    include/bertini2/function_tree/roots/function.hpp
    include/bertini2/function_tree/roots/jacobian.hpp
    include/bertini2/function_tree/operators/arithmetic.hpp
    include/bertini2/function_tree/operators/trig.hpp
)

set(function_tree_headers_rootinclude_HEADERS
	include/bertini2/function_tree.hpp
)

set(settings_headers
    include/bertini2/settings/configIni_parse.hpp
)

set(functiontreeinclude_HEADERS
	include/bertini2/function_tree/node.hpp
	include/bertini2/function_tree/gather.hpp
	include/bertini2/function_tree/factory.hpp
	include/bertini2/function_tree/forward_declares.hpp
	include/bertini2/function_tree/simplify.hpp
)

set(functiontree_operators_HEADERS
	include/bertini2/function_tree/operators/operator.hpp
	include/bertini2/function_tree/operators/arithmetic.hpp
	include/bertini2/function_tree/operators/trig.hpp
)

set(functiontree_symbols_HEADERS
	include/bertini2/function_tree/symbols/symbol.hpp
	include/bertini2/function_tree/symbols/variable.hpp
	include/bertini2/function_tree/symbols/differential.hpp
	include/bertini2/function_tree/symbols/special_number.hpp
	include/bertini2/function_tree/symbols/number.hpp
    include/bertini2/function_tree/symbols/linear_product.hpp
)

set(functiontree_roots_HEADERS
	include/bertini2/function_tree/roots/function.hpp
	include/bertini2/function_tree/roots/jacobian.hpp
)

set(io_headers
    include/bertini2/io/file_utilities.hpp
    include/bertini2/io/generators.hpp
    include/bertini2/io/parsing.hpp
    include/bertini2/io/splash.hpp
)

set(io_parsing_headers
    include/bertini2/io/parsing/classic_utilities.hpp
    include/bertini2/io/parsing/function_parsers.hpp
    include/bertini2/io/parsing/function_rules.hpp
    include/bertini2/io/parsing/number_parsers.hpp
    include/bertini2/io/parsing/number_rules.hpp
    include/bertini2/io/parsing/qi_files.hpp
    include/bertini2/io/parsing/settings_parsers.hpp
    include/bertini2/io/parsing/settings_rules.hpp
    include/bertini2/io/parsing/system_parsers.hpp
    include/bertini2/io/parsing/system_rules.hpp
)

set(io_parsing_settings_headers
    include/bertini2/io/parsing/settings_parsers/algorithm.hpp
    include/bertini2/io/parsing/settings_parsers/base.hpp
    include/bertini2/io/parsing/settings_parsers/endgames.hpp
    include/bertini2/io/parsing/settings_parsers/tracking.hpp
)

set(nag_algorithms_headers
	include/bertini2/nag_algorithms/midpath_check.hpp
	include/bertini2/nag_algorithms/numerical_irreducible_decomposition.hpp
	include/bertini2/nag_algorithms/output.hpp
	include/bertini2/nag_algorithms/sharpen.hpp
	include/bertini2/nag_algorithms/trace.hpp
	include/bertini2/nag_algorithms/zero_dim_solve.hpp
)

set(nag_algorithms_common_headers
    include/bertini2/nag_algorithms/common/algorithm_base.hpp
    include/bertini2/nag_algorithms/common/config.hpp
    include/bertini2/nag_algorithms/common/policies.hpp
)

set(nag_datatypes_headers
	include/bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp
	include/bertini2/nag_datatypes/witness_set.hpp
)

set(nag_datatypes_common_headers
	include/bertini2/nag_datatypes/common/policies.hpp
)

set(parallel_headers
    include/bertini2/parallel/initialize_finalize.hpp
    include/bertini2/parallel/path_result.hpp
    include/bertini2/parallel/manager.hpp
    include/bertini2/parallel/worker.hpp
)

set(parallel_rootinclude_HEADERS
    include/bertini2/parallel.hpp
)

set(pool_headers
    include/bertini2/pool/pool.hpp
    include/bertini2/pool/system.hpp
)

set(system_headers
	include/bertini2/system/patch.hpp
	include/bertini2/system/precon.hpp
	include/bertini2/system/slice.hpp
	include/bertini2/system/start_base.hpp
	include/bertini2/system/start_systems.hpp
	include/bertini2/system/straight_line_program.hpp
	include/bertini2/system/system.hpp
)

set(system_start_headers
	include/bertini2/system/start/total_degree.hpp
	include/bertini2/system/start/mhom.hpp
	include/bertini2/system/start/user.hpp
	include/bertini2/system/start/utility.hpp
)

set(system_rootinclude_HEADERS
	include/bertini2/system.hpp
)

set(trackers_HEADERS
	include/bertini2/trackers/adaptive_precision_utilities.hpp
	include/bertini2/trackers/amp_criteria.hpp
	include/bertini2/trackers/amp_tracker.hpp
	include/bertini2/trackers/base_predictor.hpp
	include/bertini2/trackers/base_tracker.hpp
	include/bertini2/trackers/config.hpp
	include/bertini2/trackers/events.hpp
	include/bertini2/trackers/explicit_predictors.hpp
	include/bertini2/trackers/fixed_precision_tracker.hpp
	include/bertini2/trackers/fixed_precision_utilities.hpp
	include/bertini2/trackers/newton_correct.hpp
	include/bertini2/trackers/newton_corrector.hpp
	include/bertini2/trackers/observers.hpp
	include/bertini2/trackers/ode_predictors.hpp
	include/bertini2/trackers/predict.hpp
	include/bertini2/trackers/step.hpp
	include/bertini2/trackers/tracker.hpp
)

set(tracking_rootinclude_HEADERS
	include/bertini2/tracking.hpp
)

set(BERTINI2_LIBRARY_HEADERS
    ${endgames_headers}
    ${basics_headers}
    ${common_headers}
    ${detail_headers}
    ${function_tree_headers}
    ${io_headers}
    ${io_parsing_headers}
    ${io_parsing_settings_headers}
    ${nag_algorithms_headers}
    ${nag_algorithms_common_headers}
    ${nag_datatypes_headers}
    ${nag_datatypes_common_headers}
    ${parallel_headers}
    ${parallel_rootinclude_HEADERS}
    ${pool_headers}
    ${system_headers}
    ${system_start_headers}
    ${system_rootinclude_HEADERS}
    ${systeminclude_HEADERS}
    ${startinclude_HEADERS}
    ${trackers_HEADERS}
    ${tracking_header_files}
    ${tracking_rootinclude_HEADERS}
    ${trackersinclude_HEADERS}
)

set(basics_sources
    src/basics/random.cpp
    src/basics/have_bertini.cpp
)

set(function_tree_sources
    src/function_tree/node.cpp
    src/function_tree/gather.cpp
    src/function_tree/simplify.cpp
    src/function_tree/operators/arithmetic.cpp
    src/function_tree/operators/trig.cpp
    src/function_tree/linear_product.cpp
    src/function_tree/operators/operator.cpp
    src/function_tree/symbols/special_number.cpp
    src/function_tree/symbols/differential.cpp
    src/function_tree/symbols/symbol.cpp
    src/function_tree/symbols/variable.cpp
    src/function_tree/symbols/number.cpp
    src/function_tree/roots/jacobian.cpp
    src/function_tree/roots/function.cpp
)

set(parallel_sources
    src/parallel/parallel.cpp
    src/parallel/initialize_finalize.cpp
)

set(system_source_files
    src/system/precon.cpp
    src/system/slice.cpp
    src/system/start_base.cpp
    src/system/system.cpp
    src/system/straight_line_program.cpp
    src/system/start/total_degree.cpp
    src/system/start/mhom.cpp
    src/system/start/user.cpp
)

set(tracking_source_files
    src/tracking/explicit_predictors.cpp
)

# explicit template instantiation TUs: the one home for the heavy
# endgame/zerodim instantiations (see ADR-0014); consumers get extern
# declarations from the umbrella headers and just link.
set(eti_source_files
    src/eti/endgames_eti.cpp
    src/eti/zero_dim_eti.cpp
    src/eti/zero_dim_blackbox_eti.cpp
)

set(BERTINI2_LIBRARY_SOURCES
    ${basics_sources}
    ${function_tree_sources}
    ${parallel_sources}
    ${system_source_files}
    ${tracking_source_files}
    ${eti_source_files}
)

set(BERTINI2_EXE_SOURCES
    src/blackbox/bertini.cpp
    src/blackbox/main_mode_switch.cpp
    src/blackbox/argc_argv.cpp
    src/blackbox/algorithm_builder.cpp
)

set(BERTINI2_EXE_HEADERS
    include/bertini2/blackbox/main_mode_switch.hpp
    include/bertini2/blackbox/argc_argv.hpp
    include/bertini2/blackbox/config.hpp
    include/bertini2/blackbox/algorithm_builder.hpp
    include/bertini2/blackbox/global_configs.hpp
    include/bertini2/blackbox/switches_zerodim.hpp
    include/bertini2/blackbox/user_homotopy.hpp
)


# this is correct if run from core, but not if from a higher level...
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# see https://stackoverflow.com/questions/25676277/cmake-target-include-directories-prints-an-error-when-i-try-to-add-the-source
# for more on these generator expressions, $<INSTALL_INTERFACE, etc.
add_library(bertini2 ${BERTINI2_LIBRARY_SOURCES}
                            ${BERTINI2_LIBRARY_HEADERS})
add_library(bertini2::bertini2 ALIAS bertini2)

# The (static) library is linked into the shared _pybertini Python module, so
# its objects must be position-independent. Previously this requirement was
# accidentally satisfied because every consumer recompiled the library sources
# itself (with -fPIC); now that the library is compiled once and linked, mark
# it PIC explicitly.
set_target_properties(bertini2 PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(bertini2
  SYSTEM PUBLIC
    ${GMP_INCLUDES}
    ${MPFR_INCLUDES}
    ${MPC_INCLUDES}
)

target_include_directories(bertini2 PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>  # <prefix>/include
    )

if(WIN32)
  target_compile_options(bertini2
                         PRIVATE $<$<CXX_COMPILER_ID:MSVC>:-bigobj -MP>)
  target_compile_definitions(bertini2 PUBLIC "HAVE_SNPRINTF")
endif()

target_link_libraries(bertini2 PUBLIC ${GMP_LIBRARIES})
target_link_libraries(bertini2 PUBLIC ${MPFR_LIBRARIES})
target_link_libraries(bertini2 PUBLIC ${MPC_LIBRARIES})
target_link_libraries(bertini2 PUBLIC Eigen3::Eigen)
target_link_libraries(bertini2 PUBLIC ${Boost_LIBRARIES})

if(ENABLE_MPI)
    target_compile_definitions(bertini2 PUBLIC BERTINI2_HAVE_MPI)
    target_link_libraries(bertini2 PUBLIC MPI::MPI_CXX)
    message(STATUS "bertini2 built with plain C MPI support (no Boost.MPI)")
endif()


if(SUFFIX_SO_VERSION)
  set_target_properties(bertini2 PROPERTIES SOVERSION ${PROJECT_VERSION})
endif(SUFFIX_SO_VERSION)

add_executable(bertini2_exe)
target_sources(bertini2_exe PUBLIC ${BERTINI2_EXE_SOURCES} ${BERTINI2_EXE_HEADERS})
target_link_libraries(bertini2_exe ${Boost_LIBRARIES} bertini2)
set_target_properties(bertini2_exe PROPERTIES OUTPUT_NAME bertini2)
# On MSVC/lld-link, linking an executable with OUTPUT_NAME matching the static library
# (bertini2) causes lld-link to write bertini2.lib as an import library, overwriting
# the static archive. /NOIMPLIB suppresses the import library entirely.
if(MSVC)
  target_link_options(bertini2_exe PRIVATE /NOIMPLIB)
endif()

# Precompile the four blackbox headers that collectively cover all exe source files.
# algorithm_builder.hpp is the heaviest (pulls in full tracking + NAG algorithm
# machinery); the others are lighter but still include significant Bertini subsystems.
target_precompile_headers(bertini2_exe PRIVATE
    <bertini2/eigen_extensions.hpp>
    <bertini2/blackbox/argc_argv.hpp>
    <bertini2/blackbox/main_mode_switch.hpp>
    <bertini2/blackbox/algorithm_builder.hpp>
    <bertini2/parallel.hpp>
)

# todo: this should be made a devmode thing
#target_compile_options(bertini2 PRIVATE -Wall -Wextra)

install(
  TARGETS bertini2
  EXPORT ${TARGETS_EXPORT_NAME}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  INCLUDES
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_header_group(${PROJECT_NAME}_HEADERS)
add_source_group(${PROJECT_NAME}_SOURCES)

install(
    TARGETS bertini2
    ARCHIVE DESTINATION "lib"
    LIBRARY DESTINATION "lib"
    COMPONENT library
)

install(
    FILES ${BERTINI2_EXE_HEADERS}
    DESTINATION "include/bertini2/blackbox"
)

install(
    FILES ${function_tree_headers_rootinclude_HEADERS} ${parallel_rootinclude_HEADERS} ${system_rootinclude_HEADERS} ${tracking_rootinclude_HEADERS} ${parallel_rootinclude_HEADERS} ${system_rootinclude_HEADERS} ${tracking_rootinclude_HEADERS}
    DESTINATION "include/bertini2"
)

install(
    FILES ${CMAKE_BINARY_DIR}/include/bertini2/config.hpp
    DESTINATION "include/bertini2"
)

install(
    FILES ${endgames_headers}
    DESTINATION "include/bertini2/endgames"
)

install(
    FILES ${trackers_headers}
    DESTINATION "include/bertini2/trackers"
)

install(
    FILES ${basics_rootinclude_headers}
    DESTINATION "include/bertini2"
)

install(
    FILES ${common_headers}
    DESTINATION "include/bertini2/common"
)

install(
    FILES ${detail_headers}
    DESTINATION "include/bertini2/detail"
)

install(
    FILES ${function_tree_headers}
    DESTINATION "include/bertini2/function_tree"
)

install(
    FILES ${function_tree_headers_rootinclude_HEADERS}#
    DESTINATION "include/bertini2"
)

install(
    FILES ${functiontree_operators_HEADERS}
    DESTINATION "include/bertini2/function_tree/operators"
)

install(
    FILES ${functiontree_symbols_HEADERS}
    DESTINATION "include/bertini2/function_tree/symbols"
)

install(
    FILES ${functiontree_roots_HEADERS}
    DESTINATION "include/bertini2/function_tree/roots/"
)

install(
    FILES ${io_headers}
    DESTINATION "include/bertini2/io"
)

install(
    FILES ${io_parsing_headers}
    DESTINATION "include/bertini2/io/parsing"
)

install(
    FILES ${io_parsing_settings_headers}
    DESTINATION "include/bertini2/io/parsing/settings_parsers"
)

install(
    FILES ${nag_algorithms_headers}
    DESTINATION "include/bertini2/nag_algorithms"
)

install(
    FILES ${nag_algorithms_common_headers}
    DESTINATION "include/bertini2/nag_algorithms/common"
)

install(
    FILES ${nag_datatypes_headers}
    DESTINATION "include/bertini2/nag_datatypes"
)

install(
    FILES ${nag_datatypes_common_headers}
    DESTINATION "include/bertini2/nag_datatypes/common"
)

install(
    FILES ${parallel_headers}
    DESTINATION "include/bertini2/parallel"
)

install(
    FILES ${parallel_rootinclude_HEADERS}
    DESTINATION "include/bertini2"
)

install(
    FILES ${pool_headers}
    DESTINATION "include/bertini2/pool"
)

install(
    FILES ${system_headers}
    DESTINATION "include/bertini2/system"
)

install(
    FILES ${system_start_headers}
    DESTINATION "include/bertini2/system/start"
)

install(
    FILES ${system_rootinclude_HEADERS}
    DESTINATION "include/bertini2"
)

install(
    FILES ${tracking_rootinclude_HEADERS}
    DESTINATION "include/bertini2"
)

install(
    FILES ${trackers_HEADERS}
    DESTINATION "include/bertini2/trackers"
)

install(
    FILES ${settings_headers}
    DESTINATION "include/bertini2/settings"
)

install(
    FILES ${tracking_header_files}
    DESTINATION "include/bertini2/tracking"
)

###################
#
#
#       UNIT TESTING
#
#
#############





if (ENABLE_UNIT_TESTING)

	set(B2_CLASS_TEST_SOURCES
		test/classes/boost_multiprecision_test.cpp
		test/classes/fundamentals_test.cpp
		test/classes/precision_propagation_test.cpp
		test/classes/eigen_test.cpp
		test/classes/complex_test.cpp
		test/classes/function_tree_test.cpp
		test/classes/function_tree_transform.cpp
		test/classes/system_test.cpp
		test/classes/slp_test.cpp
		test/classes/differentiate_test.cpp
		test/classes/differentiate_wrt_var.cpp
		test/classes/homogenization_test.cpp
		test/classes/start_system_test.cpp
		test/classes/node_serialization_test.cpp
		test/classes/patch_test.cpp
		test/classes/slice_test.cpp
		test/classes/m_hom_start_system.cpp
		test/classes/products_of_linears_block_test.cpp
		test/classes/class_test.cpp
	)


	add_executable(test_classes ${B2_CLASS_TEST_SOURCES})
	target_include_directories(test_classes PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_classes ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_classes COMMAND ${CMAKE_BINARY_DIR}/core/test_classes)


	set(B2_BLACKBOX_TEST
		test/blackbox/blackbox.cpp
		test/blackbox/zerodim.cpp
		test/blackbox/parsing.cpp
		test/blackbox/user_homotopy.cpp
	)

	add_executable(test_blackbox ${B2_BLACKBOX_TEST})
	target_include_directories(test_blackbox PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
	target_link_libraries(test_blackbox ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_blackbox COMMAND ${CMAKE_BINARY_DIR}/core/test_blackbox)


	set(B2_CLASSIC_TEST
		test/classic/classic_parsing_test.cpp
		test/classic/classic_test.cpp
	)

	add_executable(test_classic ${B2_CLASSIC_TEST})
	target_include_directories(test_classic PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_classic ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_classic COMMAND ${CMAKE_BINARY_DIR}/core/test_classic)

	set(B2_ENDGAMES_TEST
		test/endgames/endgames_test.cpp
		test/endgames/generic_interpolation.hpp
		test/endgames/interpolation.cpp
		test/endgames/generic_pseg_test.hpp
		test/endgames/amp_powerseries_test.cpp
		test/endgames/fixed_double_powerseries_test.cpp
		test/endgames/fixed_multiple_powerseries_test.cpp
		test/endgames/generic_cauchy_test.hpp
		test/endgames/amp_cauchy_test.cpp
		test/endgames/fixed_double_cauchy_test.cpp
		test/endgames/fixed_multiple_cauchy_test.cpp
	)

	add_executable(test_endgames ${B2_ENDGAMES_TEST})
	target_include_directories(test_endgames PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_endgames ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_endgames COMMAND ${CMAKE_BINARY_DIR}/core/test_endgames)

	set(B2_GENERATING_TEST
		test/generating/mpfr_float.cpp
		test/generating/mpfr_complex.cpp
		test/generating/double.cpp
		test/generating/std_complex.cpp
		test/generating/generating_test.cpp
	)

	add_executable(test_generating ${B2_GENERATING_TEST})
	target_include_directories(test_generating PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_generating ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_generating COMMAND ${CMAKE_BINARY_DIR}/core/test_generating)

	set(B2_NAG_ALGORITHMS_TEST
		test/nag_algorithms/nag_algorithms_test.cpp
		test/nag_algorithms/zero_dim.cpp
		test/nag_algorithms/numerical_irreducible_decomposition.cpp
		test/nag_algorithms/trace.cpp
	)

	add_executable(test_nag_algorithms ${B2_NAG_ALGORITHMS_TEST})
	target_include_directories(test_nag_algorithms PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_nag_algorithms ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_nag_algorithms COMMAND ${CMAKE_BINARY_DIR}/core/test_nag_algorithms)

	set(B2_NAG_DATATYPES_TEST
		test/nag_datatypes/witness_set.cpp
		test/nag_datatypes/nag_datatypes_test.cpp
		test/nag_datatypes/numerical_irreducible_decomposition.cpp
	)

	add_executable(test_nag_datatypes ${B2_NAG_DATATYPES_TEST})
	target_include_directories(test_nag_datatypes PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_nag_datatypes ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_nag_datatypes COMMAND ${CMAKE_BINARY_DIR}/core/test_nag_datatypes)

	set(B2_POOLS_TEST
		test/pools/pool_test.cpp
	)

	add_executable(test_pool ${B2_POOLS_TEST})
	target_include_directories(test_pool PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_pool ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_pool COMMAND ${CMAKE_BINARY_DIR}/core/test_pool)

	set(B2_SETTINGS_TEST
		test/settings/settings_test.cpp
	)

	add_executable(test_settings ${B2_SETTINGS_TEST})
	target_include_directories(test_settings PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_settings ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_settings COMMAND ${CMAKE_BINARY_DIR}/core/test_settings)

	set(B2_TRACKING_BASICS_TEST
		test/tracking_basics/newton_correct_test.cpp
		test/tracking_basics/euler_test.cpp
		test/tracking_basics/heun_test.cpp
		test/tracking_basics/higher_predictor_test.cpp
		test/tracking_basics/tracking_basics_test.cpp
		test/tracking_basics/fixed_precision_tracker_test.cpp
		test/tracking_basics/amp_criteria_test.cpp
		test/tracking_basics/amp_tracker_test.cpp
		test/tracking_basics/path_observers.cpp
	)

	add_executable(test_tracking_basics ${B2_TRACKING_BASICS_TEST})
	target_include_directories(test_tracking_basics PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
	target_link_libraries(test_tracking_basics ${Boost_LIBRARIES} bertini2)
	add_test(NAME test_tracking_basics COMMAND ${CMAKE_BINARY_DIR}/core/test_tracking_basics)

	enable_testing()

endif (ENABLE_UNIT_TESTING)

# lld-link (Windows) won't pull data-only object files (e.g. explicit_predictors.obj,
# which holds all ExplicitRKPredictor Butcher-tableau statics) unless every object in
# the archive is forced in. CMake 3.24+ LINK_LIBRARY:WHOLE_ARCHIVE places the archive
# in the linker's input list (not just options) with the appropriate /WHOLEARCHIVE: flag.
#
# We strip the normal "bertini2" link entry first to avoid having the archive appear
# twice. $<LINK_LIBRARY:WHOLE_ARCHIVE,bertini2> re-adds it and propagates transitive
# interface deps (gmp, mpfr, boost, MPI, etc.) automatically.
# /FORCE:MULTIPLE: Boost.Serialization void_cast_register and oserializer templates are
# instantiated in multiple TUs within bertini2 itself; /WHOLEARCHIVE exposes all copies
# and /FORCE:MULTIPLE picks the first (same as weak/COMDAT on Linux/macOS).
if(MSVC AND ENABLE_UNIT_TESTING)
  foreach(_t test_classes test_blackbox test_classic test_endgames test_generating
              test_nag_algorithms test_nag_datatypes test_pool test_settings
              test_tracking_basics)
    if(TARGET ${_t})
      get_target_property(_ll ${_t} LINK_LIBRARIES)
      if(_ll)
        list(REMOVE_ITEM _ll bertini2)
        set_target_properties(${_t} PROPERTIES LINK_LIBRARIES "${_ll}")
      endif()
      target_link_options(${_t} PRIVATE /FORCE:MULTIPLE)
      target_link_libraries(${_t} "$<LINK_LIBRARY:WHOLE_ARCHIVE,bertini2>")
    endif()
  endforeach()
endif()
