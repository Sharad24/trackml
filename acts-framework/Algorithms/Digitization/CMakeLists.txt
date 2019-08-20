file(GLOB_RECURSE src_files "src/*.*pp")

add_library(ACTFWDigitization SHARED ${src_files})
target_include_directories(ACTFWDigitization PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWDigitization PRIVATE ACTFramework ActsCore ActsDigitizationPlugin IdentificationPlugin)

install(TARGETS ACTFWDigitization LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
