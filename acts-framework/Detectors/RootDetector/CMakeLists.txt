add_library(ACTFWRootDetector INTERFACE)
target_include_directories(ACTFWRootDetector INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWRootDetector INTERFACE ActsCore ActsTGeoPlugin)
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
