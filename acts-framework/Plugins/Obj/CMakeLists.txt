file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTFWObjPlugin SHARED ${src_files})
target_include_directories(ACTFWObjPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWObjPlugin PUBLIC ActsCore ACTFramework)
target_link_libraries(ACTFWObjPlugin PUBLIC Threads::Threads)

install(TARGETS ACTFWObjPlugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
