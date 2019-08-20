set(srcs include/ACTFW/Plugins/Root/RootExCellWriter.hpp
         include/ACTFW/Plugins/Root/RootExCellWriter.ipp
         include/ACTFW/Plugins/Root/RootPlanarClusterWriter.hpp
         #include/ACTFW/Plugins/Root/RootIndexedMaterialReader.hpp
         #include/ACTFW/Plugins/Root/RootIndexedMaterialWriter.hpp
         #include/ACTFW/Plugins/Root/RootMaterialTrackReader.hpp
         #include/ACTFW/Plugins/Root/RootMaterialTrackWriter.hpp
         include/ACTFW/Plugins/Root/RootPropagationStepsWriter.hpp
         include/ACTFW/Plugins/Root/RootParticleWriter.hpp
         include/ACTFW/Plugins/Root/RootSimHitWriter.hpp
         #src/RootIndexedMaterialReader.cpp
         #src/RootIndexedMaterialWriter.cpp
         #src/RootMaterialTrackReader.cpp
         #src/RootMaterialTrackWriter.cpp
         src/RootPlanarClusterWriter.cpp
         src/RootParticleWriter.cpp
         src/RootSimHitWriter.cpp
         src/RootPropagationStepsWriter.cpp)
# generate root dictionary for material mapping - needed to write this class to a root file
#if (DEFINED Acts_FOUND)
#root_generate_dictionary(
#  G__Material
#  ${Acts_INCLUDE_DIR}/Acts/Material/Material.hpp
#  ${Acts_INCLUDE_DIR}/Acts/Plugins/MaterialMapping/MaterialStep.hpp
#  ${Acts_INCLUDE_DIR}/Acts/Plugins/MaterialMapping/MaterialTrack.hpp
#  LINKDEF include/ACTFW/Plugins/Root/Material_LinkDef.h
#  OPTIONS -I${Acts_INCLUDE_DIR})
#else()
#root_generate_dictionary(
#  G__Material
#  ${PROJECT_SOURCE_DIR}/external/acts-core/Core/include/Acts/Material/Material.hpp
#  ${PROJECT_SOURCE_DIR}/external/acts-core/Plugins/MaterialPlugin/include/Acts/Plugins/MaterialMapping/MaterialStep.hpp
#  ${PROJECT_SOURCE_DIR}/external/acts-core/Plugins/MaterialPlugin/include/Acts/Plugins/MaterialMapping/MaterialTrack.hpp
#  LINKDEF include/ACTFW/Plugins/Root/Material_LinkDef.h
#  OPTIONS -I${PROJECT_SOURCE_DIR}/external/acts-core/Core/include
#          -I${PROJECT_SOURCE_DIR}/external/acts-core/Plugins/MaterialPlugin/include)
#endif()

add_library(ACTFWRootPlugin SHARED ${srcs}) # G__Material.cxx)
target_include_directories(ACTFWRootPlugin PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_include_directories(ACTFWRootPlugin PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ACTFWRootPlugin PUBLIC ActsCore ActsDigitizationPlugin IdentificationPlugin) # ActsMaterialPlugin)
target_link_libraries(ACTFWRootPlugin PUBLIC ACTFramework ACTFWPropagation)
target_link_libraries(ACTFWRootPlugin PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(ACTFWRootPlugin PUBLIC Threads::Threads)

install(TARGETS ACTFWRootPlugin LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ACTFW DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
