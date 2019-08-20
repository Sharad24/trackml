add_executable(ActsParticleGun ParticleGun.cpp)
target_include_directories(ActsParticleGun PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(ActsParticleGun
  PRIVATE ActsFrameworkGenerators ACTFWExamplesCommon ACTFWRootPlugin
    ACTFWCsvPlugin ${Boost_LIBRARIES})

install(TARGETS ActsParticleGun RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(USE_PYTHIA8)
  add_executable(ActsPythia8 Pythia8.cpp)
  target_include_directories(ActsPythia8 PRIVATE ${Boost_INCLUDE_DIRS})
  target_link_libraries(ActsPythia8
    PRIVATE ActsFrameworkGenerators ActsFrameworkPythia8 ACTFWExamplesCommon
      ACTFWRootPlugin ACTFWCsvPlugin ${Boost_LIBRARIES})

  install(TARGETS ActsPythia8 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
