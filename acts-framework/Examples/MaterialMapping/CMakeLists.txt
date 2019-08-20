include(${Geant4_USE_FILE})

file (GLOB_RECURSE src_files "src/*.cpp")

add_executable(ACTFWGeantinoRecordingExample src/GeantinoRecordingExample.cpp)
add_executable(ACTFWMaterialMappingExample src/MaterialMappingExample.cpp)

foreach(_example ACTFWGeantinoRecordingExample ACTFWMaterialMappingExample)
  target_include_directories(${_example} PRIVATE ${DD4hep_INCLUDE_DIRS})
  target_link_libraries(${_example} PRIVATE ACTFWDD4hepG4Plugin)
  target_link_libraries(${_example} PRIVATE ACTFWExtrapolation)
  target_link_libraries(${_example} PRIVATE ACTFWMaterialMapping)
  target_link_libraries(${_example} PRIVATE ACTFWRootPlugin)
  target_link_libraries(${_example} PRIVATE ${DD4hep_LIBRARIES})
  target_link_libraries(${_example} PRIVATE ${Geant4_LIBRARIES})
  target_link_libraries(${_example} PRIVATE ${Boost_LIBRARIES})
  install(TARGETS ${_example} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()
