# Try to find the FLINT library (Fast Library for Number Theory).
# See https://flintlib.org/
#
# FLINT 3 includes Arb, the arbitrary-precision ball arithmetic b2 uses for
# interval evaluation, so the minimum is 3.0.0 unless the caller asks for newer:
#   find_package(FLINT 3.0.0)
#
# Once done this will define
#
#  FLINT_FOUND - system has FLINT with an acceptable version
#  FLINT_INCLUDES - the directory holding flint/flint.h
#  FLINT_LIBRARIES - the FLINT library
#  FLINT_VERSION - FLINT version
#
# Modeled on FindMPFR.cmake (BSD-licensed, from KDE/Eigen).

find_path(FLINT_INCLUDES
  NAMES
  flint/flint.h
  PATHS
  $ENV{FLINTDIR}/include
  ${INCLUDE_INSTALL_DIR}
)

if(NOT FLINT_FIND_VERSION)
  set(FLINT_FIND_VERSION "3.0.0")
endif()

if(FLINT_INCLUDES)

  file(READ "${FLINT_INCLUDES}/flint/flint.h" _flint_version_header)

  string(REGEX MATCH "define[ \t]+__FLINT_VERSION[ \t]+([0-9]+)" _flint_major_version_match "${_flint_version_header}")
  set(FLINT_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+__FLINT_VERSION_MINOR[ \t]+([0-9]+)" _flint_minor_version_match "${_flint_version_header}")
  set(FLINT_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+__FLINT_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _flint_patchlevel_version_match "${_flint_version_header}")
  set(FLINT_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(FLINT_VERSION ${FLINT_MAJOR_VERSION}.${FLINT_MINOR_VERSION}.${FLINT_PATCHLEVEL_VERSION})

  if(${FLINT_VERSION} VERSION_LESS ${FLINT_FIND_VERSION})
    set(FLINT_VERSION_OK FALSE)
    message(STATUS "FLINT version ${FLINT_VERSION} found in ${FLINT_INCLUDES}, "
                   "but at least version ${FLINT_FIND_VERSION} is required")
  else()
    set(FLINT_VERSION_OK TRUE)
  endif()

endif(FLINT_INCLUDES)

# The library's file name differs by platform and packager: libflint.so / .dylib on unix,
# flint.lib or libflint.dll.a from conda-forge on Windows.  Name every spelling.
find_library(FLINT_LIBRARIES
  NAMES flint libflint flint.lib libflint.lib libflint.dll.a
  PATHS $ENV{FLINTDIR}/lib ${LIB_INSTALL_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLINT DEFAULT_MSG
                                  FLINT_INCLUDES FLINT_LIBRARIES FLINT_VERSION_OK)
mark_as_advanced(FLINT_INCLUDES FLINT_LIBRARIES)
