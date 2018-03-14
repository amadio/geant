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
        ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
        ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm5/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/TestEm3/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullCMS/inc
#        ${CMAKE_SOURCE_DIR}/examples/physics/GeantV/FullLHCb/inc
#        ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
)

add_executable(PhysVecRngBench VecRng/VecRngBench.cxx)
target_link_libraries(PhysVecRngBench Geant_v benchmark pthread)

add_executable(PhysVecKleinNishBench KleinNishina/KleinNishinaBench.cxx KleinNishina/KleinNishinaTestCommon.h)
target_link_libraries(PhysVecKleinNishBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecKleinNishValid KleinNishina/KleinNishinaValidate.cxx)
target_link_libraries(PhysVecKleinNishValid Geant_v RealPhysics)

add_executable(PhysVecKleinNishAliasBench KleinNishina/KleinNishinaAliasBench.cxx KleinNishina/KleinNishinaTestCommon.h)
target_link_libraries(PhysVecKleinNishAliasBench Geant_v RealPhysics benchmark pthread)

add_executable(PhysVecKleinNishAliasValid KleinNishina/KleinNishinaAliasValid.cxx KleinNishina/KleinNishinaTestCommon.h)
target_link_libraries(PhysVecKleinNishAliasValid Geant_v RealPhysics)
