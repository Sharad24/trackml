include(${Geant4_USE_FILE})

file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTFWMaterialMapping SHARED ${src_files})
target_include_directories(ACTFWMaterialMapping PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWMaterialMapping PUBLIC ActsCore ActsMaterialMappingPlugin)
target_link_libraries(ACTFWMaterialMapping PUBLIC ACTFramework ACTFWGeant4Plugin)

install(TARGETS ACTFWMaterialMapping LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
