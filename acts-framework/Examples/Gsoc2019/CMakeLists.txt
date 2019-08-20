add_executable(
    ACTFWGsoc2019Reconstruction Reconstruction.cpp)
target_include_directories(ACTFWGsoc2019Reconstruction PUBLIC ../Algorithms/Reconstruction "{TORCH_LIBRARIES}" )
target_link_libraries(
    ACTFWGsoc2019Reconstruction
    PRIVATE ACTFramework ACTFWExamplesCommon
        ActsFrameworkGenerators
        ACTFWGenericDetector
        ACTFWFatras
        ACTFWDigitization
        FatrasCore
        ActsCore
        ACTFWDAGbasedNNTracker
        ACTFWTopquarks
        ACTFWMikadoTracker
        ACTFWpytorch
        Boost::program_options
        ${TBB_LIBRARIES}
        ${TORCH_LIBRARIES})
install(
    TARGETS ACTFWGsoc2019Reconstruction
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
