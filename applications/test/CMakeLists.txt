# add_subdirectory(app)
add_executable(${PROJECT_NAME}_test1
    app.cxx
)

target_link_libraries(${PROJECT_NAME}_test1
    PUBLIC ${LINK_LIBS}
)

target_include_directories(${PROJECT_NAME}_test1
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PUBLIC ${INCL_DIRS}
)