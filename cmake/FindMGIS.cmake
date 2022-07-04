
include(FindPackageHandleStandardArgs)

find_library(MGIS_MFRONT_LIBRARY
        NAMES   libMFrontGenericInterface.so libMFrontGenericInterface.dylib
        PATHS   ${MGIS_ROOT}/lib
)

find_path(MGIS_INCLUDE_DIRS
        #    NAMES   Integrate.hxx
        NAMES	MGIS/Behaviour/Integrate.hxx
        PATHS   ${MGIS_ROOT}/include
)

if (MGIS_MFRONT_LIBRARY AND MGIS_INCLUDE_DIRS)
    set(MGIS_FOUND TRUE)
    set(MGIS_LIBRARIES ${MGIS_MFRONT_LIBRARY} )
endif ()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MGIS DEFAULT_MSG MGIS_LIBRARIES MGIS_INCLUDE_DIRS)

