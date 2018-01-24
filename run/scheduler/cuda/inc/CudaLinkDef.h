#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#ifdef GEANTCUDA_REPLACE
#include "../../vecprot_v2/inc/LinkDef.h"
// #include "../../xsec/inc/xsmcLinkDef.h"
#endif

#pragma link C++ class CoprocessorBroker+;
#pragma link C++ class CoprocessorBroker::TaskData;

#endif
