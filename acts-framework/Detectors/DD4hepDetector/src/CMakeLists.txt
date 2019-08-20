file(GLOB_RECURSE src_files "src/*.cpp")

add_library (ACTFWDD4hepDetector SHARED ${src_files})
target_include_directories(ACTFWDD4hepDetector PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTFWDD4hepDetector PUBLIC ${DD4hep_INCLUDE_DIRS})
target_include_directories(ACTFWDD4hepDetector PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ACTFWDD4hepDetector PRIVATE ACTFramework)
target_link_libraries(ACTFWDD4hepDetector PUBLIC ActsCore Acts::DD4hepPlugin)
target_link_libraries(ACTFWDD4hepDetector PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(ACTFWDD4hepDetector PUBLIC ${DD4hep_LIBRARIES})
 
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(ACTFWDD4hepDetector PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
dd4hep_set_version(ACTFWDD4hepDetector MAJOR 1 MINOR 0 PATCH 0)
dd4hep_generate_rootmap(ACTFWDD4hepDetector)

install(TARGETS ACTFWDD4hepDetector LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
