file (GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_executable(ACTFWWhiteBoardExample ${src_files})
target_include_directories(ACTFWWhiteBoardExample PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWWhiteBoardExample PUBLIC ActsCore)
target_link_libraries(ACTFWWhiteBoardExample PRIVATE ACTFramework ACTFWExamplesCommon)
target_link_libraries(ACTFWWhiteBoardExample PRIVATE ${Boost_LIBRARIES})

install(TARGETS ACTFWWhiteBoardExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
