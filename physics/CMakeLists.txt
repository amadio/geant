#----------------------------------------------------------------------------------------------
# Add source, header files & include directories
#----------------------------------------------------------------------------------------------

if(WITH_GEANT4)
# DO NOT INCLUDE EXPLICITLY THE FILE BELOW WHICH OVERWRITES CXX_FLAGS(DEBUG)
# Rather append Geant4_CXX_FLAGS(_DEBUG) to existing ones
#  include(${Geant4_USE_FILE})
endif()

file(GLOB sources
  ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/src/*.cxx
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/hadron/models/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/hadron/processes/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/src/*.cc
  ${CMAKE_SOURCE_DIR}/physics/hadron/cross_sections/src/*.cc
)

file(GLOB headers
  ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/PhysicalConstants.h
  ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/SystemOfUnits.h
  ${CMAKE_SOURCE_DIR}/core/base/inc/Geant/VectorPhysicsTypes.h
  ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/hadron/models/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/hadron/processes/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc/Geant/*.h
  ${CMAKE_SOURCE_DIR}/physics/hadron/cross_sections/inc/Geant/*.h
)

include_directories(
   ${CMAKE_SOURCE_DIR}/physics/kernel/material/inc
   ${CMAKE_SOURCE_DIR}/physics/kernel/cuts/inc
   ${CMAKE_SOURCE_DIR}/core/base/inc
   ${CMAKE_SOURCE_DIR}/physics/kernel/particles/inc
   ${CMAKE_SOURCE_DIR}/physics/kernel/management/inc
   ${CMAKE_SOURCE_DIR}/physics/electromagnetic/models/inc
   ${CMAKE_SOURCE_DIR}/physics/hadron/models/inc
   ${CMAKE_SOURCE_DIR}/physics/electromagnetic/processes/inc
   ${CMAKE_SOURCE_DIR}/physics/hadron/processes/inc
   ${CMAKE_SOURCE_DIR}/physics/kernel/utils/inc
   ${CMAKE_SOURCE_DIR}/physics/hadron/cross_sections/inc
   ${CMAKE_SOURCE_DIR}/core/concurrency/inc
   ${CMAKE_SOURCE_DIR}/core/interfaces/inc
   ${CMAKE_SOURCE_DIR}/core/numa/inc
)

#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
set(_extraLibs "")

if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(realphysicsDict ${headers} MODULE RealPhysics LINKDEF rootdict/realphysicsLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  set(sources ${sources} realphysicsDict.cxx)
  set(_extraLibs "${ROOT_LIBRARIES}")
endif()

#----------------------------------------------------------------------------------------------
# Add vec-phys sources, include directories and lib if vec-phys build was required explicitly
#----------------------------------------------------------------------------------------------
if(USE_VECPHYS)
  set(sources ${sources}
      ${CMAKE_SOURCE_DIR}/realphysics/physics/vecem/src/GUGammaComptonProcess.cc
      ${CMAKE_SOURCE_DIR}/realphysics/physics/vecem/src/GUGammaPhotoElectricProcess.cc
      ${CMAKE_SOURCE_DIR}/realphysics/physics/vecem/src/GUGammaConversionProcess.cc
      ${CMAKE_SOURCE_DIR}/realphysics/physics/vecem/src/GUBetheHeitlerConversionModel.cc
      ${CMAKE_SOURCE_DIR}/realphysics/physics/vecem/src/GUKleinNishinaComptonModel.cc
      ${CMAKE_SOURCE_DIR}/realphysics/physics/vecem/src/GUSauterGavrilaModel.cc
  )
  include_directories(${CMAKE_SOURCE_DIR}/vecphys/inc physics/vecem/inc)
  set(_extraLibs ${_extraLibs} "vecphys")
endif()

#----------------------------------------------------------------------------------------------
# Add target library and links
#----------------------------------------------------------------------------------------------
add_library(RealPhysics SHARED ${sources})
target_link_libraries(RealPhysics -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Material Geant_v ${_extraLibs}
                      ${VECGEOM_LIBRARIES} ${Geant4_LIBRARIES} ${VecCoreLib_LIBRARIES})
add_custom_target(realphyslib DEPENDS RealPhysics)

#----------------------------------------------------------------------------------------------
# Install the files
#
install(FILES ${headers} DESTINATION inc/Geant)

#----------------------------------------------------------------------------------------------
if(APPLE)
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRealPhysics${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS RealPhysics DESTINATION lib)
endif()
if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRealPhysics.rootmap
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libRealPhysics_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add source files & include directories to doxygen
#----------------------------------------------------------------------------------------------
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
  set(DOXYFILE_CITE_BIB_FILES "${DOXYFILE_CITE_BIB_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/physics/models/EM/doc/biblio.bib" PARENT_SCOPE)
endif()
