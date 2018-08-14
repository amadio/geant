
find_package(benchmark)

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

        Common/
)

set(COMMON_HEADERS Common/Hist.h Common/Common.h)
set(__extralibs ${PROFILE_LIBRARIES} benchmark::benchmark pthread)


if(${benchmark_FOUND})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks/vecphys)

    add_executable(RngBench Misc/VecRngBench.cxx ${COMMON_HEADERS})
    target_link_libraries(RngBench Geant_v ${__extralibs})

    add_executable(MathBench Misc/VecMathBench.cxx ${COMMON_HEADERS})
    target_link_libraries(MathBench Geant_v ${__extralibs})

    add_executable(KleinNishBench KleinNishina/KleinNishinaBench.cxx Common/CommonBenchmark.h ${COMMON_HEADERS})
    target_link_libraries(KleinNishBench Geant_v RealPhysics ${__extralibs})

    add_executable(VectorPairBench VectorPairTests/VectorPairBench.cxx Common/CommonBenchmark.h ${COMMON_HEADERS})
    target_link_libraries(VectorPairBench Geant_v RealPhysics ${__extralibs})

    add_executable(VectorBremsBench VectorBremsTests/VectorBremsBench.cxx Common/CommonBenchmark.h ${COMMON_HEADERS})
    target_link_libraries(VectorBremsBench Geant_v RealPhysics ${__extralibs})

    add_executable(VectorPos2GammaBench PositronTo2Gamma/PositronTo2GammaBench.cxx Common/CommonBenchmark.h ${COMMON_HEADERS})
    target_link_libraries(VectorPos2GammaBench Geant_v RealPhysics ${__extralibs})

    add_executable(MollerBhabhaBench MollerBhabha/MollerBhabhaBench.cxx Common/CommonBenchmark.h ${COMMON_HEADERS})
    target_link_libraries(MollerBhabhaBench Geant_v RealPhysics ${__extralibs})

    add_executable(SauterGavrilaBench SauterGavrila/SGBench.cxx Common/CommonBenchmark.h ${COMMON_HEADERS})
    target_link_libraries(SauterGavrilaBench Geant_v RealPhysics ${__extralibs})

endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/vecphys)

add_executable(KleinNishValid KleinNishina/KleinNishinaValidate.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(KleinNishValid Geant_v ${PROFILE_LIBRARIES} RealPhysics)

add_executable(PairValid VectorPairTests/VectorPairValidate.cxx VectorPairTests/VectorPairTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PairValid Geant_v ${PROFILE_LIBRARIES} RealPhysics)

add_executable(BremsValid VectorBremsTests/VectorBremsValidate.cxx VectorBremsTests/VectorBremsTestCommon.h ${COMMON_HEADERS})
target_link_libraries(BremsValid Geant_v ${PROFILE_LIBRARIES} RealPhysics)

add_executable(Pos2GammaValid PositronTo2Gamma/PositronTo2GammaValidate.cxx PositronTo2Gamma/PositronTo2GammaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(Pos2GammaValid Geant_v ${PROFILE_LIBRARIES} RealPhysics)

add_executable(MollerBhabhaValid MollerBhabha/MollerBhabhaValidate.cxx MollerBhabha/MollerBhabhaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(MollerBhabhaValid Geant_v ${PROFILE_LIBRARIES} RealPhysics)

add_executable(SauterGavrilaValid SauterGavrila/SGValidate.cxx SauterGavrila/SGTestCommon.h ${COMMON_HEADERS})
target_link_libraries(SauterGavrilaValid Geant_v ${PROFILE_LIBRARIES} RealPhysics)


