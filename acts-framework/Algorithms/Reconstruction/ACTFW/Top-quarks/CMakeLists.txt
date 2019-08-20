add_library(PolarModule SHARED PolarModule.cpp)
target_include_directories(
    PolarModule
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
# target_include_directories(ACTFramework PUBLIC ${TBB_INCLUDE_DIRS})
# target_link_libraries(ACTFramework PUBLIC ${TBB_LIBRARIES})
target_link_libraries(PolarModule PRIVATE ACTFramework ActsCore ${TBB_LIBRARIES})


add_library(ACTFWTopquarks SHARED
    main.cpp Point.hpp Triple.hpp Layer.hpp)
target_include_directories(
    ACTFWTopquarks
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include> ${TBB_INCLUDE_DIRS})
target_link_libraries(ACTFWTopquarks PRIVATE ACTFramework ActsCore PolarModule ${TBB_LIBRARIES})

install(TARGETS ACTFWTopquarks LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
