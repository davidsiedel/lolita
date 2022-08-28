enable_testing()
find_package(GTest REQUIRED)

add_executable(${PROJECT_NAME}_test
        main.cxx
        test_0.cxx
        # test_deb.cxx
        # test_triangle.cxx
        # test_quadrangle.cxx
        test_misc.cxx
        # test_integ.cxx
        # test_nnn.cxx
)
add_test(NAME ${PROJECT_NAME}_test COMMAND ${PROJECT_NAME}_test)
target_link_libraries(${PROJECT_NAME}_test
        PUBLIC GTest::gtest
        PUBLIC GTest::gtest_main
        PUBLIC ${LINK_LIBS}
)
target_include_directories(${PROJECT_NAME}_test
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        PUBLIC ${INCL_DIRS}
)
