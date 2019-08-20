file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTFramework SHARED ${src_files})
target_include_directories (ACTFramework PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTFramework PUBLIC ${ROOT_INCLUDE_DIRS})
target_include_directories(ACTFramework PUBLIC ${TBB_INCLUDE_DIRS})
target_link_libraries(ACTFramework PUBLIC ActsCore)
target_link_libraries(ACTFramework PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(ACTFramework PUBLIC ${TBB_LIBRARIES})

install(TARGETS ACTFramework LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
