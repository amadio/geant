#ifndef GXSAMPLINGTEXTURE_H
#define GXSAMPLINGTEXTURE_H 1

#include <cuda.h>
#include <cuda_runtime.h>

#ifdef __CUDACC__

texture<int2, 1> texPDFX;
texture<int2, 1> texPDFY;

static __inline__ __device__ G4double FetchToDouble(texture<int2, 1> T, int i)
{
  // for double precision float as a texture, use int2 and cast it to double 
  int2 v = tex1Dfetch(T,i);
  return __hiloint2double(v.y, v.x);
}

class GXSamplingTexture {
public:

  GXSamplingTexture() {}

  ~GXSamplingTexture() {}

  void Load(unsigned int size, G4double *PDFX_d, G4double *PDFY_d)
  {
    cudaBindTexture(0, texPDFX, PDFX_d, size);
    cudaBindTexture(0, texPDFY, PDFY_d, size);
  }
};

#endif
#endif
