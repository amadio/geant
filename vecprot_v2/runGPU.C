#define COPROCESSOR_REQUEST true

#include "run.C"

void runGPU(Int_t nthreads=4,
            const char *geomfile="ExN03.root",
            const char *xsec="xsec_FTFP_BERT.root",
            const char *fstate="fstate_FTFP_BERT.root"
            )
{
   run(nthreads,geomfile,xsec,fstate,true);
}
