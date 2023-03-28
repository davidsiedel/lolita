add_executable(${PROJECT_NAME}_test2
    test.cxx
)

target_link_libraries(${PROJECT_NAME}_test2
    PUBLIC ${LINK_LIBS}
)

target_include_directories(${PROJECT_NAME}_test2
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include_>
    $<INSTALL_INTERFACE:include_>
    PUBLIC ${INCL_DIRS}
)