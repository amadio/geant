#---------------------------------------------------------------------------
# Setup the project


#----------------------------------------------------------------------------
# Find HepMC (required package)
#
find_package(HepMC QUIET)
if(HepMC_FOUND)
  message(STATUS "HepMC includes: ${HEPMC_INCLUDE_DIR}")
  include_directories(AFTER SYSTEM ${HEPMC_INCLUDE_DIR})
endif()

#set(HepMC_FOUND NO)

#----------------------------------------------------------------------------
# Setup include directory for this project
#


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources src/T*.cxx src/GunGenerator.cxx src/TestProcess.cxx)
file(GLOB headers inc/T*.h inc/GunGenerator.h inc/TestProcess.h)

#set(XSEC_HEADERS ${headers} PARENT_SCOPE)
#set(XSEC_SOURCES ${sources} PARENT_SCOPE)

#set(headers ${headers} ${VC_INCLUDE_DIR}/base/Vector.h)

if(HepMC_FOUND)
  add_definitions(-DHEPMC)
  file(GLOB _mcg_cxx src/HepMCGenerator.cxx)
  file(GLOB _mcg_h   inc/HepMCGenerator.h)
  set(sources ${sources} ${_mcg_cxx})
  set(headers ${headers} ${_mcg_h})
endif()

include_directories(inc)


if (NOT MIC)
  ROOT_GENERATE_DICTIONARY(xsecDict ${headers} MODULE libXsec LINKDEF inc/xsecLinkDef.h OPTIONS "-DVC_COMMON_SIMDIZE_H_=1")
  #----------------------------------------------------------------------------
  # Add the executable, and link it to the Geant4 libraries
  #
  add_library(Xsec SHARED ${sources} xsecDict.cxx)

  set(ROOT_LIBS "")
  if(USE_ROOT OR NOT USE_VECGEOM_NAVIGATOR)
    set(ROOT_LIBS "${ROOT_LIBRARIES} -lEG -lGui -lGeom")
  endif()
else()
  add_library(Xsec SHARED ${sources})
endif()

set(HEPMC_LIBS "")
if(HepMC_FOUND)
  set(HEPMC_LIBS "${HEPMC_LIBRARIES}")
endif()

if (NOT MIC)
target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v  
     ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})
else()
target_link_libraries(Xsec -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v 
     ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})
endif()

set(XSEC_LINK_LIBS "-lXsec" ${HEPMC_LIBRARIES} PARENT_SCOPE)

add_dependencies(Xsec Geantlib_v)

add_headers("${headers}")

add_custom_target(xseclib DEPENDS Xsec)

#add_executable(testserialRead testserialRead.cc)

#target_link_libraries(testserialRead -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v -lXsec
#        ${HEPMC_LIBRARIES} ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES})
#add_dependencies(testserialRead Geant_v Xsec)
if (CUDA)


 set(SRC_CPP_RELATIVE
#      xsec/testserialReadCuda.cu
     xsec/src/TPartIndex.cxx
#     xsec/src/TEFstate.cxx
#     xsec/src/TFinState.cxx
#     xsec/src/TEXsec.cxx
#     xsec/src/TMXsec.cxx
#     xsec/src/TPDecay.cxx
#     xsec/src/TPFstate.cxx
 #    xsec/src/TPXsec.cxx
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
    OPTIONS ${CUDA_ARCH}
 #   ${CUDA_NVCC_FLAGS}

)  
target_link_libraries(xseccuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -lGeant_v -lXsec ${ROOT_LIBS} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBS})

#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#message (STATUS "dir='${dirs}'")
#endforeach()

add_custom_target(testserial ALL DEPENDS xseccuda)
#target_include_directories(xseccuda PUBLIC ${dirs})
ADD_CUSTOM_COMMAND(
    TARGET testserial
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "Separate Compilation for testserialReadCuda"
    COMMAND ${CUDA_HOST_COMPILER} -c --std=c++11 ${VECGEOM_DEFINITIONS} ${GEANTV_COMPILE_DEFINITIONS} ${CMAKE_SOURCE_DIR}/xsec/testserialReadCuda.cpp -o testserialReadCuda.cpp.o -I${VECGEOM_INSTALL_DIR}/include -I${CUDA_TOOLKIT_ROOT_DIR}/include -lxseccuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND ${CMAKE_COMMAND} -E echo "Separate Compilation for testserialReadCuda nvcc ${GEANTV_COMPILE_DEFINITIONS}"
    COMMAND ${CUDA_NVCC_EXECUTABLE} ${CUDA_NVCC_FLAGS} ${VECGEOM_DEFINITIONS} -DUSE_VECGEOM_NAVIGATOR -dc -arch=sm_30 ${CMAKE_SOURCE_DIR}/xsec/testserialReadCuda.cu -o testserial.o -I${VECGEOM_INSTALL_DIR}/include -I${CMAKE_SOURCE_DIR}/include -I${CMAKE_SOURCE_DIR}/xsec/inc -I${CMAKE_SOURCE_DIR}/base/inc -lxseccuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND ${CMAKE_COMMAND} -E echo "Separate Compilation for testserialReadCuda linking"
    COMMAND ${CUDA_NVCC_EXECUTABLE} ${CUDA_NVCC_FLAGS} -arch=sm_30 ${VECGEOM_INSTALL_DIR}/lib/libvecgeom.a testserialReadCuda.cpp.o testserial.o ${VECGEOM_INSTALL_DIR}/librawvecgeomcuda.a -o runtest -lxseccuda -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    )

message (STATUS "VECGEOM LIBS ${VECGEOM_LIBRARIES}")


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
