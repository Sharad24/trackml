#cmake_minimum_required(VERSION 3.5)
#project(trackml LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
add_subdirectory(analysis)

add_library(MikadoTracker SHARED reconstruction.cxx AccuracyEvaluator.h Cuts.h DataStructures.h Engine.h EventReader.h Geo.h Learning.h SearchLayer.h Tracker.h TrackModelPhysical.h TrackSelector.h util.h)

add_library(TrackModelPhysical.o TrackModelPhysical.cxx)
add_library(SearchLayer.o SearchLayer.cxx)
add_library(MTracker.o Tracker.cxx)
add_library(MReconstruction.o Reconstruction.cxx)
add_library(Geo.o Geo.cxx)
add_library(Cuts.o Cuts.cxx)
add_library(DataStructures.o DataStructures.cxx)
add_library(Learning.o Learning.cxx)
add_library(AccuracyEvaluator.o AccuracyEvaluator.cxx)
add_library(Engine.o Engine.cxx)
add_library(TrackSelector.o TrackSelector.cxx)
add_library(EventReader.o EventReader.cxx)

target_link_libraries(MikadoTracker TrackModelPhysical.o SearchLayer.o MTracker.o MReconstruction.o Geo.o Cuts.o PolynomFit.o DataStructures.o Learning.o AccuracyEvaluator.o Engine.o TrackSelector.o EventReader.o TrackFitTest.o AnalyseGeometry.o AnalyseField.o)
target_include_directories(MikadoTracker PUBLIC . analysis)
install(TARGETS MikadoTracker LIBRARY DESTINATION lib/)

#add_executable(reco reconstruction.cxx)
#target_link_libraries(reco MikadoTracker)

add_library(ACTFWMikadoTracker SHARED
    main.cpp AccuracyEvaluator.h Cuts.h DataStructures.h Engine.h EventReader.h Geo.h Learning.h SearchLayer.h Tracker.h TrackModelPhysical.h TrackSelector.h util.h)
target_include_directories(
    ACTFWMikadoTracker
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
target_link_libraries(ACTFWMikadoTracker PRIVATE ACTFramework ActsCore MikadoTracker ${TBB_LIBRARIES})

install(TARGETS ACTFWMikadoTracker LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
