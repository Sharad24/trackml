file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTFWCsvPlugin SHARED ${src_files})
target_include_directories(ACTFWCsvPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWCsvPlugin PRIVATE ACTFramework ActsCore ActsDigitizationPlugin IdentificationPlugin)
target_link_libraries(ACTFWCsvPlugin PRIVATE Threads::Threads)

install(TARGETS ACTFWCsvPlugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
