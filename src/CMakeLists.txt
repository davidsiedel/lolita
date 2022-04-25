#add_library(${PROJECT_NAME} SHARED
#
#        )
#
#target_include_directories(${PROJECT_NAME}
#        PUBLIC
#        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
#        $<INSTALL_INTERFACE:include>
#        #        PRIVATE ${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}
#        PUBLIC ${INCL_DIRS}
#        )
#target_link_libraries(${PROJECT_NAME}
#        PUBLIC ${LINK_LIBS}
#        #        PUBLIC Eigen3::Eigen
#        #        ${Boost_LIBRARIES}
#        #        ${PYTHON_LIBRARIES}
#        #        PUBLIC Python3::NumPy
#        )
#set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")