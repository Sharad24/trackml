add_library(ActsFrameworkPythia8 SHARED
  ACTFW/Generators/Pythia8ProcessGenerator.cpp
  ACTFW/Options/Pythia8Options.cpp)

target_include_directories(ActsFrameworkPythia8
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_include_directories(ActsFrameworkPythia8 PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(ActsFrameworkPythia8 PUBLIC ActsCore ActsFrameworkGenerators ACTFramework Pythia8)
target_link_libraries(ActsFrameworkPythia8 PUBLIC ${Boost_LIBRARIES})
