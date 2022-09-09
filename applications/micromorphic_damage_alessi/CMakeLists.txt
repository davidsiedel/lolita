# add_subdirectory(app)
add_executable(${PROJECT_NAME}_micromorphic_damage_alessi
    app.cxx
)

target_link_libraries(${PROJECT_NAME}_micromorphic_damage_alessi
    PUBLIC ${LINK_LIBS}
)

target_include_directories(${PROJECT_NAME}_micromorphic_damage_alessi
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PUBLIC ${INCL_DIRS}
)