# generate and install benchmark scripts
configure_file(
    "TestEm3/TestEm3.run.in"
    "${CMAKE_INSTALL_PREFIX}/bin/benchmarks/TestEm3/TestEm3.run"
    IMMEDIATE @ONLY)

configure_file(
    "TestEm3/TestEm3_GV.run.in"
    "${CMAKE_INSTALL_PREFIX}/bin/benchmarks/TestEm3/TestEm3_GV.run"
    IMMEDIATE @ONLY)

configure_file(
    "FullCMS/FullCMS.run.in"
    "${CMAKE_INSTALL_PREFIX}/bin/benchmarks/FullCMS/FullCMS.run"
    IMMEDIATE @ONLY)

configure_file(
    "FullCMS/FullCMS_GV.run.in"
    "${CMAKE_INSTALL_PREFIX}/bin/benchmarks/FullCMS/FullCMS_GV.run"
    IMMEDIATE @ONLY)

file(GLOB config_testem3 TestEm3/*.sh)
file(GLOB config_fullcms FullCMS/*.sh)
install(FILES ${config_testem3} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/benchmarks/TestEm3)
install(FILES ${config_fullcms} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/benchmarks/FullCMS)
