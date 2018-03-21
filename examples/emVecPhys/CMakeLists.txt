find_package(benchmark REQUIRED)

include_directories(
        ${CMAKE_SOURCE_DIR}/core/numa/inc
        ${CMAKE_SOURCE_DIR}/core/concurrency/inc
        ${CMAKE_SOURCE_DIR}/core/interfaces/inc
        ${CMAKE_SOURCE_DIR}/core/base/inc

        ${CMAKE_SOURCE_DIR}/run/userapp/inc
        ${CMAKE_SOURCE_DIR}/run/scheduler/inc
        ${CMAKE_SOURCE_DIR}/magneticfield/inc
        ${CMAKE_SOURCE_DIR}/physics/generators/inc
        ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
        ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
        ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
        ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
        ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/inc
        ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
        ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
)



add_executable(SauterGavrilaBench SauterGavrilaBench.cxx)
target_link_libraries(SauterGavrilaBench Geant_v RealPhysics benchmark pthread)

add_executable(SauterGavrila_SampleShell_benchmark SauterGavrila_SampleShell_benchmark.cxx)
target_link_libraries(SauterGavrila_SampleShell_benchmark Geant_v RealPhysics benchmark pthread)


