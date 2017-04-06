#----------------------------------------------------------------------------------------------
# Add source files & include directories
#
file(GLOB sources src/T*.cxx src/GunGenerator.cxx src/TestProcess.cxx)
file(GLOB headers inc/T*.h inc/GunGenerator.h inc/TestProcess.h)

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx src/HepMCTruth.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h inc/HepMCTruth.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc)
include_directories(${CMAKE_SOURCE_DIR}/physics/particle/inc)
include_directories(${CMAKE_SOURCE_DIR}/numa/inc)

#----------------------------------------------------------------------------------------------
# Add the library
#
if (USE_ROOT)
   if (USE_UMESIMD)
      ROOT_GENERATE_DICTIONARY(xsecDict ${headers} MODULE Xsec LINKDEF inc/xsecLinkDef.h OPTIONS "-DVECGEOM_UMESIMD=1")
      add_library(Xsec SHARED ${sources} xsecDict.cxx)
      target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v NewPhysics 
        ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${HEPMC_LIBRARIES})
   else()
      ROOT_GENERATE_DICTIONARY(xsecDict ${headers} MODULE Xsec LINKDEF inc/xsecLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
      add_library(Xsec SHARED ${sources} xsecDict.cxx)
      target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v NewPhysics 
         ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})
   endif()
  
else()
  add_library(Xsec SHARED ${sources})
  
  target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} Geant_v NewPhysics
    ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})
endif()

set(XSEC_LINK_LIBS "Xsec ${HEPMC_LIBRARIES}" PARENT_SCOPE)

add_dependencies(Xsec Geantlib_v)
add_headers("${headers}")
add_custom_target(xseclib DEPENDS Xsec)

#----------------------------------------------------------------------------------------------
# Add the test executable
#
add_executable(testserialRead testserialRead.cc)
target_link_libraries(testserialRead -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} GeantNuma Geant_v Xsec Vmagfield NewPhysics)
add_dependencies(testserialRead Geant_v Xsec)

#----------------------------------------------------------------------------------------------
# CUDA section
#
if (CUDA)
  set(SRC_CPP_RELATIVE
    xsec/src/TEFstate.cxx
    xsec/src/TPartIndex.cxx
    xsec/src/TPXsec.cxx
    xsec/src/TEXsec.cxx
    xsec/src/TPDecay.cxx
    xsec/src/TFinState.cxx
    xsec/src/TPFstate.cxx
    #    xsec/src/TMXsec.cxx
    #    xsec/src/TTabPhysMgr.cxx
    #    xsec/src/TTabPhysProcess.cxx
    )

  foreach(SRC ${SRC_CPP_RELATIVE})
    set(SRC_CPP ${SRC_CPP} ${CMAKE_SOURCE_DIR}/${SRC})
  endforeach()
  

  foreach(SRC_FILE ${SRC_CPP})
    get_filename_component(SRC_FILENAME ${SRC_FILE} NAME_WE)

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      COMMAND ${CMAKE_COMMAND} -E copy ${SRC_FILE}
      ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      DEPENDS ${SRC_FILE}
      )
    
    set(
      SRC_CUDA ${SRC_CUDA}
      ${CMAKE_CURRENT_BINARY_DIR}/cuda_src/${SRC_FILENAME}.cu
      )
  endforeach()

  cuda_add_library(xseccuda
    ${SRC_CUDA}
    SHARED
    OPTIONS ${CUDA_ARCH}
    )  

  # Note: the repeat below is due the author of  cuda_add_library_depend
  # not knowing how to automatically go from the target to the real
  # file in add_custom_command
  cuda_add_library_depend(xseccuda ../lib/libNewPhysicsCuda_static.a ${VECGEOM_CUDA_STATIC_LIBRARY} )
  target_link_libraries(xseccuda NewPhysicsCuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL})
  
  cuda_add_library(xseccudatest
    testserialReadCuda.cu
    SHARED
    OPTIONS ${CUDA_ARCH}
    )

  # Note: the repeat below is due the author of  cuda_add_library_depend
  # not knowing how to automatically go from the target to the real
  # file in add_custom_command
  cuda_add_library_depend(xseccudatest xseccuda_static ../lib/libxseccuda_static.a 
    ../lib/libNewPhysicsCuda_static.a ../lib/libGeant_v_static.a ${VECGEOM_CUDA_STATIC_LIBRARY})
  target_link_libraries(xseccudatest -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} xseccuda Geant_v)
  
  add_executable(testserialReadCuda testserialReadCuda.cpp)
  target_link_libraries(testserialReadCuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} xseccudatest xseccuda Xsec Geant_v)

endif()

#----------------------------------------------------------------------------------------------
# Install the files
#
if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Xsec DESTINATION lib)
endif()

if(USE_ROOT)
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec.rootmap 
    ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec_rdict.pcm DESTINATION lib)
endif()

#----------------------------------------------------------------------------------------------
# Add this directory to Doxygen
#
if(DOXYGEN_FOUND)
  set(DOXYFILE_SOURCE_DIRS "${DOXYFILE_SOURCE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
endif()

