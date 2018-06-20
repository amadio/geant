Usage of gperftools profiler with GeantV on linux
=================================================

1. Download and install.

For extensive instructions, see: https://github.com/gperftools/gperftools/blob/master/INSTALL
Two options are possible:

a. Download the development package from your preferred distibution (e.g. libgoogle-perftools-dev on Ubuntu)
or:
b. Install manually from the github repository:
Make sure libunwind is installed on your system. This is critical to determine the call-chain of a program
git clone https://github.com/gperftools/gperftools.git
cd gperftools
./autogen.sh
./configure --enable-frame-pointers --with-tcmalloc-pagesize=32
sudo make -j install (in case you want libs in /usr/local/lib)

2. Compile GeantV with gperftools support

Just add -DGPERFTOOLS=ON and compile, watching that the tool detection went fine. For project dependencies it is 
recommended to compile with -DCMAKE_CXX_FLAGS="-fno-omit-frame-pointer -gdwarf"

3. Profiling is automatically enabled for the method geant::RunManager::RunSimulation. In future some API will be added to enable
user-defined profiling blocks. The name of the output file of the profiler is by default gperfprof.out, but can be changed:
- either use RunManager::SetProfilingFile() or:
- define the environment variable: GEANT_PERFTOOLS_FILE
After running the GeantV executable MyExec, the profiling file will be automatically produced. The difference in 
run time compared to the non-profiling case is very small.

4. Inspect the profile

Install gv (ghostview) for producing the output in a nice graphics form:
pprof --gv MyExec gperfprof.out
There are other options and profiling tweaks described here: https://github.com/gperftools/gperftools/blob/master/README
