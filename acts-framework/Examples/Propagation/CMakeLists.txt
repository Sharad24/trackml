file (GLOB_RECURSE src_files "src/*.*pp")

add_executable(ACTFWGenericPropagationExample src/GenericPropagationExample.cpp)
target_include_directories(ACTFWGenericPropagationExample PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWGenericPropagationExample PRIVATE ActsCore)
target_link_libraries(ACTFWGenericPropagationExample PRIVATE ACTFramework ACTFWExamplesCommon ACTFWPropagation ACTFWGenericDetector)
target_link_libraries(ACTFWGenericPropagationExample PRIVATE ACTFWRootPlugin ACTFWObjPlugin ACTFWBFieldPlugin)
target_link_libraries(ACTFWGenericPropagationExample PRIVATE ${Boost_LIBRARIES})

install(TARGETS ACTFWGenericPropagationExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# TGEO based detector
if (USE_TGEO)
  add_executable(ACTFWRootPropagationExample src/RootPropagationExample.cpp)
  target_include_directories(ACTFWRootPropagationExample PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(ACTFWRootPropagationExample PRIVATE ActsCore)
  target_link_libraries(ACTFWRootPropagationExample PRIVATE ACTFramework ACTFWExamplesCommon ACTFWPropagation ACTFWRootDetector)
  target_link_libraries(ACTFWRootPropagationExample PRIVATE ACTFWRootPlugin ACTFWObjPlugin ACTFWBFieldPlugin)
  target_link_libraries(ACTFWRootPropagationExample PRIVATE ${Boost_LIBRARIES})

  install(TARGETS ACTFWRootPropagationExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# DD4hep detector
if (USE_DD4HEP)
  add_executable(ACTFWDD4hepPropagationExample src/DD4hepPropagationExample.cpp)
  target_include_directories(ACTFWDD4hepPropagationExample PRIVATE ${DD4hep_INCLUDE_DIRS})
  target_include_directories(ACTFWDD4hepPropagationExample PRIVATE ${Boost_INCLUDE_DIRS})  
  target_link_libraries(ACTFWDD4hepPropagationExample PRIVATE ActsCore)
  target_link_libraries(ACTFWDD4hepPropagationExample PRIVATE ACTFramework ACTFWExamplesCommon ACTFWPropagation ACTFWDD4hepDetector)
  target_link_libraries(ACTFWDD4hepPropagationExample PRIVATE ACTFWRootPlugin ACTFWObjPlugin ACTFWBFieldPlugin)
  target_link_libraries(ACTFWDD4hepPropagationExample PRIVATE ${DD4hep_LIBRARIES})
  target_link_libraries(ACTFWDD4hepPropagationExample PRIVATE ${Boost_LIBRARIES})

  install(TARGETS ACTFWDD4hepPropagationExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
