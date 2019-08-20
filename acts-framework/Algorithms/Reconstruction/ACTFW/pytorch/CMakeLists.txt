find_package(Torch REQUIRED)

add_library(ACTFWpytorch SHARED main.cpp)
target_include_directories(ACTFWpytorch PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS} ${TORCH_LIBRARIES} "${TORCH_INCLUDE_DIRS}")
target_link_libraries(ACTFWpytorch PRIVATE ACTFramework ActsCore ${TBB_LIBRARIES} "${TORCH_LIBRARIES}")

#add_executable(exec main.cpp)
#target_link_libraries(exec PUBLIC "${TORCH_LIBRARIES}" ACTFramework ActsCore ${TBB_LIBRARIES})
