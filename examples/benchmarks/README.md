GeantV benchmarks
=================

The scripts in the subfolders TestEm3 and FullCMS are installed in INSTALLATION_PATH/bin/benchmarks. They allow comparing the performance 
of the simulation for the corresponding examples from INSTALLATION_PATH/bin/examples between GeantV with multiple options and Geant4.

To run a given GeantV configuration:
   
  bash [example]_GV.run [NREP] [F] [P] [G] [M]

where:
  example = TestEm3 or FullCMS in the corresponding folders
  NREP    = Number of repetitions. The run time average and standard deviation are computed and printed out
  F       = Magnetic field basketizing option applied to all models
  P       = Physics basketizing option
  G       = Geometry basketizing option
  M       = Multiple scattering basketizing option

The F/P/G/M possible values for the basketizing options:
  0    = basketizing off
  1    = basketizing on calling vectorized code
  2    = basketizing with scalar dispatch in a loop

To run the corresponding Geant4 configuration with NREP repetitions:

  bash [example]_G4.run [NREP]

To run benchmarks in all configurations with NREP repetitions:

  bash bench_all.sh [NREP]

The scripts FullCMS_GV.run and TestEm3_GV.run can be used to do profiling based on gperftools. To do this GeantV has to be
compiled with USE_PERFTOOLS option ON. The above GeaqntV scripts should be run like:

  [example]_GV.run 1

The method RunSimulation gets automatically profiled, resulting in the file prof_F_P_G_M.out

Inspect the profile using the pprof tool from gperftools:
  pprof [--focus=method_name] [–nodecount=maxvis_nodes] [--nodefraction=0.001] [--edgefraction=0.001] [--ps] [profile_file]

This will display profiling info for method_name, showing maximum maxvis_nodes, dropping nodes and edges with <0.1% hits, displaying a graph
using ghostview.
