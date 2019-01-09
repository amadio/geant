# generate and install benchmark scripts
set(_examples TestEm3
              FullCMS)
set(_mcs      G4 GV)

foreach(_example ${_examples})
foreach(_mc ${_mcs})
  configure_file(${_example}/${_example}_${_mc}.run.in
                 ${CMAKE_BINARY_DIR}/bin/benchmarks/${_example}/${_example}_${_mc}.run
                 IMMEDIATE @ONLY)
  install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/benchmarks/${_example}/${_example}_${_mc}.run DESTINATION bin/benchmarks/${_example})
endforeach()
endforeach()

file(GLOB config_testem3 TestEm3/*.sh)
file(GLOB config_fullcms FullCMS/*.sh)
install(PROGRAMS ${config_testem3} DESTINATION bin/benchmarks/TestEm3)
install(FILES README.md DESTINATION bin/benchmarks/TestEm3)
install(PROGRAMS ${config_fullcms} DESTINATION bin/benchmarks/FullCMS)
install(FILES README.md DESTINATION bin/benchmarks/FullCMS)
