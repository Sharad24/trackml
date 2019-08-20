include(${Geant4_USE_FILE})

file(GLOB_RECURSE src_files "src/*.cpp" "include/*.*pp")

add_library(ACTFWDD4hepG4Plugin SHARED ${src_files})
target_include_directories(ACTFWDD4hepG4Plugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTFWDD4hepG4Plugin PUBLIC ${DD4hep_INCLUDE_DIRS})
target_include_directories(ACTFWDD4hepG4Plugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ACTFWDD4hepG4Plugin PUBLIC ACTFramework)
target_link_libraries(ACTFWDD4hepG4Plugin PUBLIC ${DD4hep_LIBRARIES} ${DD4hep_DDG4_LIBRARY})
target_link_libraries(ACTFWDD4hepG4Plugin PUBLIC ${Geant4_LIBRARIES})
target_link_libraries(ACTFWDD4hepG4Plugin PUBLIC ${ROOT_LIBRARIES})

install(TARGETS ACTFWDD4hepG4Plugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
