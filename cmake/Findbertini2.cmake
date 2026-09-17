if(WIN32)
    find_path(Bertini2_INCLUDE_DIR
            NAMES bertini2
            PATHS "$ENV{BERTINI2_DIR}/Include")
    find_library(Bertini2_LIBRARY
            NAMES bertini2
            PATHS "$ENV{BERTINI2_DIR}/Lib")
elseif(UNIX)
    find_path(Bertini2_INCLUDE_DIR
            NAMES bertini2/bertini.hpp
            PATHS "$ENV{BERTINI2_DIR}/include" ${INCLUDE_INSTALL_DIR})
    find_library(Bertini2_LIBRARY
            NAMES bertini2
            PATHS "$ENV{BERTINI2_DIR}/lib" ${LIB_INSTALL_DIR})
endif()
set(Bertini2_INCLUDES ${Bertini2_INCLUDE_DIR})
set(Bertini2_LIBRARIES ${Bertini2_LIBRARY})
