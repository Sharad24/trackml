#cmake_minimum_required(VERSION 3.5)
#project(trackml LANGUAGES CXX)

#set(CMAKE_CXX_STANDARD 17)

add_library(DAGbasedNNTracker SHARED extern.cxx parallel.h Graph.h Parameters.h Point.h XMLP.h Tracker.h Reconstruction.h PolarModule.h)


add_library(PolarModule.o PolarModule.cxx)
target_include_directories(PolarModule.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
add_library(Reconstruction.o Reconstruction.cxx)
target_include_directories(PolarModule.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
target_link_libraries(Reconstruction.o PolarModule.o)
add_library(XMLP.o XMLP.cxx)
target_include_directories(PolarModule.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
add_library(Tracker.o Tracker.cxx)
target_include_directories(PolarModule.o PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})

target_link_libraries(DAGbasedNNTracker Reconstruction.o Tracker.o XMLP.o Reconstruction.o)
#add_library(DAGbasedNNTracker SHARED extern.cxx parallel.h Graph.h Parameters.h Point.h XMLP.h Tracker.h Reconstruction.h PolarModule.h)
target_include_directories(DAGbasedNNTracker PUBLIC .)
install(TARGETS DAGbasedNNTracker LIBRARY DESTINATION lib/)

#add_executable(eval eval.cxx)
#target_link_libraries(eval DAGbasedNNTracker)

add_library(ACTFWDAGbasedNNTracker SHARED
    main.cpp Tracker.h)
target_include_directories(
    ACTFWDAGbasedNNTracker
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
target_link_libraries(ACTFWDAGbasedNNTracker PRIVATE ACTFramework ActsCore DAGbasedNNTracker ${TBB_LIBRARIES})

install(TARGETS ACTFWDAGbasedNNTracker LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
