# generic detector
add_executable(ACTFWGenericGeometryExample src/GenericGeometryExample)
target_include_directories(ACTFWGenericGeometryExample PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ACTFWGenericGeometryExample PRIVATE ActsCore)
target_link_libraries(ACTFWGenericGeometryExample PRIVATE ACTFramework ACTFWExamplesCommon)
target_link_libraries(ACTFWGenericGeometryExample PRIVATE ACTFWObjPlugin ACTFWCsvPlugin ACTFWGenericDetector)
target_link_libraries(ACTFWGenericGeometryExample PRIVATE ${Boost_LIBRARIES})

install(TARGETS ACTFWGenericGeometryExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# TGEO based detector
if (USE_TGEO)
  add_executable(ACTFWRootGeometryExample src/RootGeometryExample)
  target_include_directories(ACTFWRootGeometryExample PUBLIC ${Boost_INCLUDE_DIRS})
  target_link_libraries(ACTFWRootGeometryExample PRIVATE ActsCore)
  target_link_libraries(ACTFWRootGeometryExample PRIVATE ACTFramework ACTFWExamplesCommon)
  target_link_libraries(ACTFWRootGeometryExample PRIVATE ACTFWObjPlugin ACTFWCsvPlugin ACTFWRootDetector)
  target_link_libraries(ACTFWRootGeometryExample PRIVATE ${Boost_LIBRARIES})

  install(TARGETS ACTFWRootGeometryExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

# DD4hep detector
if (USE_DD4HEP)
  add_executable(ACTFWDD4hepGeometryExample src/DD4hepGeometryExample.cpp)
  target_include_directories(ACTFWDD4hepGeometryExample PRIVATE ${DD4hep_INCLUDE_DIRS})
  target_include_directories(ACTFWDD4hepGeometryExample PRIVATE ${Boost_INCLUDE_DIRS})  
  target_link_libraries(ACTFWDD4hepGeometryExample PRIVATE ActsCore)
  target_link_libraries(ACTFWDD4hepGeometryExample PRIVATE ACTFramework ACTFWExamplesCommon)
  target_link_libraries(ACTFWDD4hepGeometryExample PRIVATE ACTFWCsvPlugin ACTFWObjPlugin ACTFWDD4hepDetector)
  target_link_libraries(ACTFWDD4hepGeometryExample PRIVATE ${DD4hep_LIBRARIES})
  target_link_libraries(ACTFWDD4hepGeometryExample PRIVATE ${Boost_LIBRARIES})

  install(TARGETS ACTFWDD4hepGeometryExample RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
