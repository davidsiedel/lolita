add_executable(${PROJECT_NAME}_app
        app.cxx
)

target_link_libraries(${PROJECT_NAME}_app
        PUBLIC ${LINK_LIBS}
)

target_include_directories(${PROJECT_NAME}_app
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        PUBLIC ${INCL_DIRS}
)
