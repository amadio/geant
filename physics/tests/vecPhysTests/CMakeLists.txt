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
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm5/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm3/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullCMS/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullLHCb/inc
#        ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
)

set(COMMON_HEADERS Common/Hist.h Common/Common.h)

if(${benchmark_FOUND})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/benchmarks/vecphys)

    add_executable(RngBench Misc/VecRngBench.cxx ${COMMON_HEADERS})
    target_link_libraries(RngBench Geant_v benchmark pthread)

    add_executable(MathBench Misc/VecMathBench.cxx ${COMMON_HEADERS})
    target_link_libraries(MathBench Geant_v benchmark pthread)

    add_executable(KleinNishBench KleinNishina/KleinNishinaBench.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(KleinNishBench Geant_v RealPhysics benchmark pthread)

    add_executable(KleinNishAliasBench KleinNishina/KleinNishinaAliasBench.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(KleinNishAliasBench Geant_v RealPhysics benchmark pthread)

    add_executable(KleinNishRejBench KleinNishina/KleinNishinaRejBench.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(KleinNishRejBench Geant_v RealPhysics benchmark pthread)

    add_executable(BetheHeitlerBench BetheHeitler/BetheHeitlerBench.cxx BetheHeitler/BetheHeitlerTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(BetheHeitlerBench Geant_v RealPhysics benchmark pthread)

    add_executable(RelativisticPairBench RelativisticPair/RelativisticPairBench.cxx RelativisticPair/RelativisticPairTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(RelativisticPairBench Geant_v RealPhysics benchmark pthread)

    add_executable(PositronTo2GammaBench PositronTo2Gamma/PositronTo2GammaBench.cxx PositronTo2Gamma/PositronTo2GammaTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(PositronTo2GammaBench Geant_v RealPhysics benchmark pthread)

    add_executable(MollerBhabhaBench MollerBhabha/MollerBhabhaBench.cxx MollerBhabha/MollerBhabhaTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(MollerBhabhaBench Geant_v RealPhysics benchmark pthread)

    add_executable(SeltzerBergerBench SeltzerBerger/SeltzerBergerBench.cxx SeltzerBerger/SeltzerBergerTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(SeltzerBergerBench Geant_v RealPhysics benchmark pthread)

    add_executable(RelativisticBremsBench RelativisticBrems/RelativisticBremsBench.cxx RelativisticBrems/RelativisticBremsTestCommon.h ${COMMON_HEADERS})
    target_link_libraries(RelativisticBremsBench Geant_v RealPhysics benchmark pthread)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/tests/vecphys)

add_executable(KleinNishValid KleinNishina/KleinNishinaValidate.cxx ${COMMON_HEADERS})
target_link_libraries(KleinNishValid Geant_v RealPhysics)

add_executable(KleinNishAliasValid KleinNishina/KleinNishinaAliasValid.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(KleinNishAliasValid Geant_v RealPhysics)

add_executable(BetheHeitlerValid BetheHeitler/BetheHeitlerValidate.cxx BetheHeitler/BetheHeitlerTestCommon.h ${COMMON_HEADERS})
target_link_libraries(BetheHeitlerValid Geant_v RealPhysics)

add_executable(RelativisticPairValid RelativisticPair/RelativisticPairValidate.cxx RelativisticPair/RelativisticPairTestCommon.h ${COMMON_HEADERS})
target_link_libraries(RelativisticPairValid Geant_v RealPhysics)

add_executable(PositronTo2GammaValid PositronTo2Gamma/PositronTo2GammaValidate.cxx ${COMMON_HEADERS})
target_link_libraries(PositronTo2GammaValid Geant_v RealPhysics)

add_executable(MollerBhabhaValid MollerBhabha/MollerBhabhaValidate.cxx MollerBhabha/MollerBhabhaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(MollerBhabhaValid Geant_v RealPhysics)

add_executable(SeltzerBergerValid SeltzerBerger/SeltzerBergerValidate.cxx SeltzerBerger/SeltzerBergerTestCommon.h ${COMMON_HEADERS})
target_link_libraries(SeltzerBergerValid Geant_v RealPhysics)

add_executable(RelativisticBremsValid RelativisticBrems/RelativisticBremsValidate.cxx RelativisticBrems/RelativisticBremsTestCommon.h ${COMMON_HEADERS})
target_link_libraries(RelativisticBremsValid Geant_v RealPhysics)

