file(GLOB_RECURSE src_files "src/*.*pp")

add_library(ACTFWGenericDetector SHARED ${src_files})
target_include_directories(ACTFWGenericDetector PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWGenericDetector PUBLIC ActsCore IdentificationPlugin ActsDigitizationPlugin)
target_link_libraries(ACTFWGenericDetector PUBLIC ACTFramework)

install(TARGETS ACTFWGenericDetector LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
