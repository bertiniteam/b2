# Borrowed from https://github.com/libigl/eigen/tree/master
# By Jack Hagen (December 2023)
# Try to find the Multiple Precision Complex (MPC)

if (MPC_INCLUDES AND MPC_LIBRARIES)
  set(MPC_FIND_QUIETLY TRUE)
endif (MPC_INCLUDES AND MPC_LIBRARIES)

find_path(MPC_INCLUDES
  NAMES
  mpc.h
  PATHS
  $ENV{MPC_INC}
  ${INCLUDE_INSTALL_DIR}
)

if(WIN32)
    find_library(MPC_LIBRARIES libmpc.dll.a PATHS $ENV{MPC_LIB} ${LIB_INSTALL_DIR})
else()
    find_library(MPC_LIBRARIES mpc PATHS $ENV{MPC_LIB} ${LIB_INSTALL_DIR})
endif()

include(FindPackageHandleStandardArgs)

# Makes sure that mpc_include and mpc_libraries are valid
# https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
find_package_handle_standard_args(MPC DEFAULT_MSG
                                  MPC_INCLUDES MPC_LIBRARIES)
mark_as_advanced(MPC_INCLUDES MPC_LIBRARIES)
