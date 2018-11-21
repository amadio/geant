GeantV benchmarks
=================

The scripts in the subfolders TestEm3 and FullCMS are installed in INSTALLATION_PATH/bin/benchmarks. They allow comparing the performance 
of the simulation for the corresponding examples from INSTALLATION_PATH/bin/examples between GeantV with multiple options and Geant4.

To run benchmarks in all configurations: bench_all.sh
To do statistical benchmarking for a given configuration of GeantV: FullCMS_GV.run bench_GV.sh <N>
  where N is the number of repetitions. The configuration in the script bench_GV.sh should be modified to match the desired options.

The same scripts  FullCMS_GV.run and TestEm3_GV.run can be used to do profiling based on gperftools. To do this GeantV has to be
compiled with USE_PERFTOOLS option ON. The above script should be run like:

FullCMS_GV.run bench_GV.sh 1

The method RunSimulation gets automatically profiled, resulting in the file $GEANT_PERFTOOLS_FILE (see bench_GV.sh)

Inspect the profile using the pprof tool from gperftools:
  pprof [--focus=method_name] [–nodecount=maxvis_nodes] [--nodefraction=0.001] [--edgefraction=0.001] [--ps]

This will display profiling info for method_name, showing maximum maxvis_nodes, dropping nodes and edges with <0.1% hits, displaying a graph
using ghostview.

