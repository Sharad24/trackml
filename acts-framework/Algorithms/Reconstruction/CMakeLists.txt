add_library(ACTFWReconstruction SHARED
    ACTFW/Reconstruction/Empty.cpp)
target_include_directories(
    ACTFWReconstruction
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(ACTFWReconstruction PRIVATE ACTFramework ActsCore)

install(TARGETS ACTFWReconstruction LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_subdirectory(ACTFW/Top-quarks)
#add_library(PolarModule SHARED ACTFW/Top-quarks/PolarModule.cpp)
#target_include_directories(
#    PolarModule
#    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
# target_include_directories(ACTFramework PUBLIC ${TBB_INCLUDE_DIRS})
# target_link_libraries(ACTFramework PUBLIC ${TBB_LIBRARIES})
#target_link_libraries(PolarModule PRIVATE ACTFramework ActsCore ${TBB_LIBRARIES})


#add_library(ACTFWTopquarks SHARED
#    ACTFW/Top-quarks/main.cpp ACTFW/Top-quarks/Point.hpp ACTFW/Top-quarks/Triple.hpp ACTFW/Top-quarks/Layer.hpp)
#target_include_directories(
#    ACTFWTopquarks
#    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#target_link_libraries(ACTFWTopquarks PRIVATE ACTFramework ActsCore PolarModule ${TBB_LIBRARIES})

#install(TARGETS ACTFWTopquarks LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

#add_library(DAGbasedNNTracker SHARED ACTFW/DAGbasedNNTracker/extern.cxx ACTFW/DAGbasedNNTracker/parallel.h ACTFW/DAGbasedNNTracker/Graph.h ACTFW/DAGbasedNNTracker/Parameters.h ACTFW/DAGbasedNNTracker/Point.h ACTFW/DAGbasedNNTracker/XMLP.h ACTFW/DAGbasedNNTracker/Tracker.h ACTFW/DAGbasedNNTracker/Reconstruction.h ACTFW/DAGbasedNNTracker/PolarModule.h)
#target_include_directories(DAGbasedNNTracker PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})

#add_library(PolarModule.o ACTFW/DAGbasedNNTracker/PolarModule.cxx)
#target_include_directories(PolarModule.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(Reconstruction.o ACTFW/DAGbasedNNTracker/Reconstruction.cxx)
#target_include_directories(Reconstruction.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#target_link_libraries(Reconstruction.o PolarModule.o)
#add_library(XMLP.o ACTFW/DAGbasedNNTracker/XMLP.cxx)
#target_include_directories(XMLP.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(Tracker.o ACTFW/DAGbasedNNTracker/Tracker.cxx)
#target_include_directories(Tracker.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})

#target_link_libraries(DAGbasedNNTracker Reconstruction.o Tracker.o XMLP.o Reconstruction.o)
#add_library(DAGbasedNNTracker SHARED ACTFW/DAGbasedNNTracker/extern.cxx ACTFW/DAGbasedNNTracker/parallel.h ACTFW/DAGbasedNNTracker/Graph.h ACTFW/DAGbasedNNTracker/Parameters.h ACTFW/DAGbasedNNTracker/Point.h ACTFW/DAGbasedNNTracker/XMLP.h ACTFW/DAGbasedNNTracker/Tracker.h ACTFW/DAGbasedNNTracker/Reconstruction.h ACTFW/DAGbasedNNTracker/PolarModule.h)
#target_include_directories(DAGbasedNNTracker PUBLIC ACTFW/DAGbasedNNTracker)
#install(TARGETS DAGbasedNNTracker LIBRARY DESTINATION lib/)

#add_library(ACTFWDAGbasedNNTracker SHARED
#    ACTFW/DAGbasedNNTracker/main.cpp ACTFW/DAGbasedNNTracker/Tracker.h)
#target_include_directories(
#    ACTFWDAGbasedNNTracker
#    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#target_link_libraries(ACTFWDAGbasedNNTracker PRIVATE ACTFramework ActsCore DAGbasedNNTracker ${TBB_LIBRARIES})

#install(TARGETS ACTFWDAGbasedNNTracker LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

#add_executable(eval ACTFW/DAGbasedNNTracker/eval.cxx)
#target_include_directories(eval PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#target_link_libraries(eval DAGbasedNNTracker)
add_subdirectory(ACTFW/DAGbasedNNTracker)
add_subdirectory(ACTFW/MikadoTracker)
add_subdirectory(ACTFW/pytorch)

#add_subdirectory(ACTFW/MikadoTracker/analysis)
#add_library(MikadoTracker SHARED ACTFW/MikadoTracker/reconstruction.cxx ACTFW/MikadoTracker/AccuracyEvaluator.h ACTFW/MikadoTracker/Cuts.h ACTFW/MikadoTracker/DataStructures.h ACTFW/MikadoTracker/Engine.h ACTFW/MikadoTracker/EventReader.h ACTFW/MikadoTracker/Geo.h ACTFW/MikadoTracker/Learning.h ACTFW/MikadoTracker/SearchLayer.h ACTFW/MikadoTracker/Tracker.h ACTFW/MikadoTracker/TrackModelPhysical.h ACTFW/MikadoTracker/TrackSelector.h ACTFW/MikadoTracker/util.h)
#target_include_directories(MikadoTracker PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})

#add_library(TrackModelPhysical.o ACTFW/MikadoTracker/TrackModelPhysical.cxx)
#target_include_directories(TrackModelPhysical.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(SearchLayer.o ACTFW/MikadoTracker/SearchLayer.cxx)
#target_include_directories(SearchLayer.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(MTracker.o ACTFW/MikadoTracker/Tracker.cxx)
#target_include_directories(MTracker.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(MReconstruction.o ACTFW/MikadoTracker/Reconstruction.cxx)
#target_include_directories(MReconstruction.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(Geo.o ACTFW/MikadoTracker/Geo.cxx)
#target_include_directories(Geo.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(Cuts.o ACTFW/MikadoTracker/Cuts.cxx)
#target_include_directories(Cuts.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(DataStructures.o ACTFW/MikadoTracker/DataStructures.cxx)
#target_include_directories(DataStructures.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(Learning.o ACTFW/MikadoTracker/Learning.cxx)
#target_include_directories(Learning.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(AccuracyEvaluator.o ACTFW/MikadoTracker/AccuracyEvaluator.cxx)
#target_include_directories(AccuracyEvaluator.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(Engine.o ACTFW/MikadoTracker/Engine.cxx)
#target_include_directories(Engine.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(TrackSelector.o ACTFW/MikadoTracker/TrackSelector.cxx)
#target_include_directories(TrackSelector.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#add_library(EventReader.o ACTFW/MikadoTracker/EventReader.cxx)
#target_include_directories(EventReader.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})

#target_link_libraries(MikadoTracker TrackModelPhysical.o SearchLayer.o MTracker.o MReconstruction.o Geo.o Cuts.o PolynomFit.o DataStructures.o Learning.o AccuracyEvaluator.o Engine.o TrackSelector.o EventReader.o TrackFitTest.o AnalyseGeometry.o AnalyseField.o)
#target_include_directories(MikadoTracker PUBLIC . analysis)
#install(TARGETS MikadoTracker LIBRARY DESTINATION lib/)

#add_library(ACTFWMikadoTracker SHARED
#    ACTFW/MikadoTracker/main.cpp ACTFW/MikadoTracker/AccuracyEvaluator.h ACTFW/MikadoTracker/Cuts.h ACTFW/MikadoTracker/DataStructures.h ACTFW/MikadoTracker/Engine.h ACTFW/MikadoTracker/EventReader.h ACTFW/MikadoTracker/Geo.h ACTFW/MikadoTracker/Learning.h ACTFW/MikadoTracker/SearchLayer.h ACTFW/MikadoTracker/Tracker.h ACTFW/MikadoTracker/TrackModelPhysical.h ACTFW/MikadoTracker/TrackSelector.h ACTFW/MikadoTracker/util.h)
#target_include_directories(
#    ACTFWMikadoTracker
#    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
#target_link_libraries(ACTFWMikadoTracker PRIVATE ACTFramework ActsCore MikadoTracker ${TBB_LIBRARIES})

#install(TARGETS ACTFWMikadoTracker LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

#add_executable(reco ACTFW/MikadoTracker/reconstruction.cxx)
#target_link_libraries(reco MikadoTracker)
