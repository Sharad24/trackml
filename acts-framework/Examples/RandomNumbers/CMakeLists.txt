file (GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_executable(ACTFWRandomNumberExample ${src_files})
target_include_directories(ACTFWRandomNumberExample PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWRandomNumberExample PUBLIC ActsCore)
target_link_libraries(ACTFWRandomNumberExample PRIVATE ACTFramework ACTFWExamplesCommon)
target_link_libraries(ACTFWRandomNumberExample PRIVATE ${Boost_LIBRARIES})


install(TARGETS ACTFWRandomNumberExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
