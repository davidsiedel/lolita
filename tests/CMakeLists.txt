## Download and unpack googletest at configure time
#configure_file(cmake/CMakeLists.txt.in googletest-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -Notations "${CMAKE_GENERATOR}" .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#        RESULT_VARIABLE result
#        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
#if(result)
#    message(FATAL_ERROR "Build step for googletest failed: ${result}")
#endif()
## Prevent overriding the parent project'local_point compiler/linker
## settings on Windows
#set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#
## Add googletest directly to our build. This defines
## the gtest and gtest_main targets.
#add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
#        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
#        EXCLUDE_FROM_ALL)
#-----------------------------------------------------------------------------------------------------------------------
enable_testing()
find_package(GTest REQUIRED)

#add_executable(foo foo.cc)
#target_link_libraries(foo GTest::gtest GTest::gtest_main)
#
#add_test(AllTestsInFoo foo)
#-----------------------------------------------------------------------------------------------------------------------


# Now simply link against gtest or gtest_main as needed. Eg
add_executable(${PROJECT_NAME}_test
        main.cxx
        test_lolita.cxx
        test_mgis.cxx)
add_test(NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME}_test)
#target_link_libraries(${PROJECT_NAME}_test
#        PUBLIC ${PROJECT_NAME} gtest)
target_link_libraries(${PROJECT_NAME}_test
        PUBLIC GTest::gtest
        PUBLIC GTest::gtest_main
#        PUBLIC ${PROJECT_NAME} gtest
        PUBLIC ${LINK_LIBS}
        )
target_include_directories(${PROJECT_NAME}_test
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        #        PRIVATE ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}
        PUBLIC ${INCL_DIRS}
        #        PUBLIC Eigen3::Eigen
        #        ${Boost_LIBRARIES}
        #        ${PYTHON_LIBRARIES}
        #        PUBLIC Python3::NumPy
        )
