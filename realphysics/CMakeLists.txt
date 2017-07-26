#----------------------------------------------------------------------------------------------
# Add source, header files & include directories
#----------------------------------------------------------------------------------------------
file(GLOB sources
  ${CMAKE_SOURCE_DIR}/realphysics/dummies/src/*.cc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/management/src/*.cxx
  ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/src/*.cc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/src/*.cc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/utils/src/*.cc
  ${CMAKE_SOURCE_DIR}/realphysics/physics/cross_sections/HAD/src/*.cc
)

file(GLOB headers
  ${CMAKE_SOURCE_DIR}/base/inc/PhysicalConstants.h
  ${CMAKE_SOURCE_DIR}/base/inc/SystemOfUnits.h
  ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc/*.h
  ${CMAKE_SOURCE_DIR}/realphysics/physics/management/inc/*.h
  ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/inc/*.h
  ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/inc/*.h
  ${CMAKE_SOURCE_DIR}/realphysics/physics/utils/inc/*.h
  ${CMAKE_SOURCE_DIR}/realphysics/physics/cross_sections/HAD/inc/*.h
)

include_directories(
   ${CMAKE_SOURCE_DIR}/material/inc
   ${CMAKE_SOURCE_DIR}/base/inc
   ${CMAKE_SOURCE_DIR}/realphysics/dummies/inc
   ${CMAKE_SOURCE_DIR}/realphysics/physics/management/inc
   ${CMAKE_SOURCE_DIR}/realphysics/physics/models/EM/inc
   ${CMAKE_SOURCE_DIR}/realphysics/physics/processes/EM/inc
   ${CMAKE_SOURCE_DIR}/realphysics/physics/utils/inc
   ${CMAKE_SOURCE_DIR}/realphysics/physics/cross_sections/HAD/inc
   ${CMAKE_SOURCE_DIR}/numa/inc
)

#----------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------
set(_extraLibs "")

if(USE_ROOT)
  ROOT_GENERATE_DICTIONARY(realphysicsDict ${headers} MODULE RealPhysics LINKDEF realphysicsLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
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
                      ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})

add_headers("${headers}")
add_custom_target(realphyslib DEPENDS RealPhysics)

#----------------------------------------------------------------------------------------------
# Install the files
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
