# generate and install benchmark scripts
set(_macros TestEm3/TestEm3_G4 
            TestEm3/TestEm3_GV 
            FullCMS/FullCMS_G4 
            FullCMS/FullCMS_GV)

foreach(_macro ${_macros})
  configure_file(${_macro}.run.in
                 ${CMAKE_BINARY_DIR}/bin/benchmarks/${_macro}.run
                 IMMEDIATE @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/bin/benchmarks/${_macro}.run DESTINATION bin/benchmarks)
endforeach()

file(GLOB config_testem3 TestEm3/*.sh)
file(GLOB config_fullcms FullCMS/*.sh)
install(PROGRAMS ${config_testem3} DESTINATION bin/benchmarks/TestEm3)
install(FILES README.md DESTINATION bin/benchmarks/TestEm3)
install(PROGRAMS ${config_fullcms} DESTINATION bin/benchmarks/FullCMS)
install(FILES README.md DESTINATION bin/benchmarks/FullCMS)
