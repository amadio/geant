#define COPROCESSOR_REQUEST true

#include "runCMS.C"

void runCMSGPU(int nthreads=4,
         bool performance=true,
	 const char *geomfile="../cmstrack/cms2015.root",
	 const char *xsec="xsec_FTFP_BERT_G496p02_1mev.root",
	 const char *fstate="fstate_FTFP_BERT_G496p02_1mev.root"
            )
{
   runCMS(nthreads,performance,geomfile,xsec,fstate,true);
}
