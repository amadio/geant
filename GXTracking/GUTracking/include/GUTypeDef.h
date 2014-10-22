#ifndef GUTypeDef_hh
#define GUTypeDef_hh 1

//  Temporary definitions for CUDA, OpenCL 'qualifications'
//

// - Geant4 type definition
typedef bool G4bool;
typedef int G4int;
typedef double G4double;

#ifdef USE_OPENCL
    #define GLOBALFUNC __kernel
    #define GLOBALTYPE __global
    #define SHAREDTYPE __local
    #define CONSTTYPE  __constant
    #define NULL       ((void*)0)
    #define GNULL ((GLOBALTYPE void*)0)
    #define FQUALIFIER 
#elif USE_MIC
    #define GLOBALFUNC __attribute__ ((target(mic)))
    #define FQUALIFIER __attribute__ ((target(mic)))
    #define GLOBALTYPE 
    #define SHAREDTYPE 
    #define CONSTTYPE  __attribute__ ((target(mic)))
    #define GNULL 0
#else
    #define GLOBALFUNC __global__
    #define GLOBALTYPE 
    #define SHAREDTYPE __shared__
    #define CONSTTYPE  
    #define GNULL 0
    #define FQUALIFIER __host__ __device__
#endif

// Variable Type Qualifer 
#define GEOMETRYLOC GLOBALTYPE
#define GEOMETRYNULL GNULL

#ifndef __CUDA_ARCH__
  #define VARTYPE CONSTTYPE const
#else
  #define VARTYPE __constant__
#endif


#endif
