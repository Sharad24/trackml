add_library(ACTFWPropagation INTERFACE)
target_include_directories(ACTFWPropagation INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTFWPropagation INTERFACE ActsCore ACTFramework)

install(TARGETS ACTFWPropagation LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})