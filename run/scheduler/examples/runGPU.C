#define COPROCESSOR_REQUEST true

#include "run.C"

void runGPU(int nthreads=4,
            bool performance=true,
            const char *geomfile="ExN03.root",
            const char *xsec="xsec_FTFP_BERT.root",
            const char *fstate="fstate_FTFP_BERT.root"
            )
{
   run(nthreads,performance,geomfile,xsec,fstate,true);
}
