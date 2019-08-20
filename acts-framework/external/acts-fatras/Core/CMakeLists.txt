add_library(FatrasCore SHARED src/RandomNumberDistributions.cpp)

# prefix output library file name with Acts
set_target_properties(FatrasCore PROPERTIES OUTPUT_NAME "FatrasCore")

# main library
target_include_directories(FatrasCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include>)
target_link_libraries(FatrasCore PUBLIC ActsCore)
set_target_properties(FatrasCore PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED TRUE)

# install library and include files to default directories
install(
  TARGETS FatrasCore
  EXPORT FatrasTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(
  DIRECTORY include/Fatras
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Fatras)
target_include_directories(FatrasCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
