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

        Common/
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm5/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm3/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullCMS/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullLHCb/inc
#        ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
)

set(COMMON_HEADERS Common/Hist.h Common/Common.h)

add_executable(PhysVecRngBench Misc/VecRngBench.cxx ${COMMON_HEADERS})
target_link_libraries(PhysVecRngBench Geant_v benchmark pthread)

add_executable(PhysVecMathBench Misc/VecMathBench.cxx ${COMMON_HEADERS})
target_link_libraries(PhysVecMathBench Geant_v benchmark pthread)

add_executable(PhysVecKleinNishBench KleinNishina/KleinNishinaBench.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecKleinNishBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecKleinNishValid KleinNishina/KleinNishinaValidate.cxx ${COMMON_HEADERS})
target_link_libraries(PhysVecKleinNishValid Geant_v RealPhysics)

add_executable(PhysVecKleinNishAliasBench KleinNishina/KleinNishinaAliasBench.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecKleinNishAliasBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecKleinNishAliasValid KleinNishina/KleinNishinaAliasValid.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecKleinNishAliasValid Geant_v RealPhysics)

add_executable(PhysVecKleinNishRejBench KleinNishina/KleinNishinaRejBench.cxx KleinNishina/KleinNishinaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecKleinNishRejBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecBetheHeitlerValid BetheHeitler/BetheHeitlerValidate.cxx BetheHeitler/BetheHeitlerTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecBetheHeitlerValid Geant_v RealPhysics)

add_executable(PhysVecBetheHeitlerBench BetheHeitler/BetheHeitlerBench.cxx BetheHeitler/BetheHeitlerTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecBetheHeitlerBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecRelativisticPairValid RelativisticPair/RelativisticPairValidate.cxx RelativisticPair/RelativisticPairTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecRelativisticPairValid Geant_v RealPhysics)

add_executable(PhysVecRelativisticPairBench RelativisticPair/RelativisticPairBench.cxx RelativisticPair/RelativisticPairTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecRelativisticPairBench Geant_v RealPhysics benchmark pthread)


add_executable(PhysVecPositronTo2GammaBench PositronTo2Gamma/PositronTo2GammaBench.cxx PositronTo2Gamma/PositronTo2GammaTestCommon.h ${COMMON_HEADERS})
target_link_libraries(PhysVecPositronTo2GammaBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecPositronTo2GammaValid PositronTo2Gamma/PositronTo2GammaValidate.cxx ${COMMON_HEADERS})
target_link_libraries(PhysVecPositronTo2GammaValid Geant_v RealPhysics)
