# Setup the project


#----------------------------------------------------------------------------
# Find HepMC (required package)
#
find_package(HepMC QUIET)
if(HepMC_FOUND)
  message(STATUS "HepMC includes: ${HEPMC_INCLUDE_DIR}")
  include_directories(AFTER SYSTEM ${HEPMC_INCLUDE_DIR})
endif()

# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/T*.cxx src/GunGenerator.cxx src/TestProcess.cxx)
file(GLOB headers inc/T*.h inc/GunGenerator.h inc/TestProcess.h)

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc)

ROOT_GENERATE_DICTIONARY(xsecDict ${headers} MODULE Xsec LINKDEF inc/xsecLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(Xsec SHARED ${sources} xsecDict.cxx)

set(HEPMC_LIBS "")
if(HepMC_FOUND)
  set(HEPMC_LIBS "${HEPMC_LIBRARIES}")
endif()

target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
     ${ROOT_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES})

set(XSEC_LINK_LIBS "-lXsec" ${HEPMC_LIBRARIES} PARENT_SCOPE)

add_dependencies(Xsec Geantlib_v)

add_headers("${headers}")

add_custom_target(xseclib DEPENDS Xsec)

#add_executable(testserialRead testserialRead.cc)

#target_link_libraries(testserialRead -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lXsec
#        ${HEPMC_LIBRARIES} ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
#add_dependencies(testserialRead Geant_v Xsec)
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

    ADD_CUSTOM_COMMAND(
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

cuda_add_library_depend(xseccuda ${VECGEOM_CUDA_STATIC_LIBRARY} )
target_link_libraries(xseccuda ${VECGEOM_LIBRARIES} ${VECGEOM_LIBRARIES_EXTERNAL})

add_executable(testserialRead testserialRead.cc)
target_link_libraries(testserialRead xseccuda Xsec Geant_v)

cuda_add_library(xseccudatest
    testserialReadCuda.cu
    SHARED
    OPTIONS ${CUDA_ARCH}
)

  # Note: the repeat below is due the author of  cuda_add_library_depend
  # not knowing how to automatically go from the target to the real
  # file in add_custom_command
cuda_add_library_depend(xseccudatest xseccuda_static ../lib/libxseccuda_static.a ../lib/libGeant_v_static.a ${VECGEOM_CUDA_STATIC_LIBRARY})
target_link_libraries(xseccudatest xseccuda)

add_executable(testserialReadCuda testserialReadCuda.cpp)
target_link_libraries(testserialReadCuda xseccudatest xseccuda Xsec Geant_v)

endif()



if(APPLE) 
# for some reason the installation process alters the libraries on Mac, and the program fails
  install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec${CMAKE_SHARED_LIBRARY_SUFFIX} DESTINATION lib)
else()
  install(TARGETS Xsec DESTINATION lib)
endif()
if (NOT MIC)
install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec.rootmap 
              ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libXsec_rdict.pcm DESTINATION lib)
endif()
