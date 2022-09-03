add_executable(${PROJECT_NAME}_linear_elastic_tensile_square
    app.cxx
)

target_link_libraries(${PROJECT_NAME}_linear_elastic_tensile_square
    PUBLIC ${LINK_LIBS}
)

target_include_directories(${PROJECT_NAME}_linear_elastic_tensile_square
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PUBLIC ${INCL_DIRS}
)