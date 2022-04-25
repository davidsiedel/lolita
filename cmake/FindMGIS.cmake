
include(FindPackageHandleStandardArgs)

find_library(MGIS_MFRONT_LIBRARY
        NAMES   libMFrontGenericInterface.so libMFrontGenericInterface.dylib
        PATHS	/home/dsiedel/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/mgis-master-vdgcnh7wvfthpcmgsghwsxilhgh5bkn3/lib
        )

#if(MGIS_MFRONT_LIBRARY)
#    message("MGIS_MFRONT_LIBRARY found")
#else()
#    message("MGIS_MFRONT_LIBRARY not found")
#endif()


find_path(MGIS_INCLUDE_DIRS
        #    NAMES   Integrate.hxx
        NAMES	MGIS/Behaviour/Integrate.hxx
        #    PATHS   /home/dsiedel/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/mgis-master-rl6bdjtij4rspaz636wx3uf46wghnnib/include/MGIS/Behaviour /home/dsiedel/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/mgis-master-rl6bdjtij4rspaz636wx3uf46wghnnib/include/MGIS /home/dsiedel/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/mgis-master-rl6bdjtij4rspaz636wx3uf46wghnnib/include
        PATHS   /home/dsiedel/spack/opt/spack/linux-ubuntu20.04-skylake/gcc-9.3.0/mgis-master-vdgcnh7wvfthpcmgsghwsxilhgh5bkn3/include
        )

#if(MGIS_INCLUDE_DIRS)
#    message("MGIS_INCLUDE_DIRS found")
#else()
#    message("MGIS_INCLUDE_DIRS not found")
#endif()



if (MGIS_MFRONT_LIBRARY AND MGIS_INCLUDE_DIRS)
    set(MGIS_FOUND TRUE)
    set(MGIS_LIBRARIES ${MGIS_MFRONT_LIBRARY} )
endif ()

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MGIS DEFAULT_MSG MGIS_LIBRARIES MGIS_INCLUDE_DIRS)

