#ifndef GXMAGNETICFIELDTEXTURE_H
#define GXMAGNETICFIELDTEXTURE_H 1

#include <cuda.h>
#include <cuda_runtime.h>
#include "GXFieldMap.h"

#ifdef __CUDACC__
texture<float, cudaTextureType2D, cudaReadModeElementType> texZ;
texture<float, cudaTextureType2D, cudaReadModeElementType> texR;

class MagneticFieldTexture {
public:
  cudaArray *texZarray;
  cudaArray *texRarray;

  MagneticFieldTexture() : texZarray(0),texRarray(0) {}

  ~MagneticFieldTexture() {
    cudaUnbindTexture(texZ);
    cudaUnbindTexture(texR);
    cudaFreeArray(texZarray);
    cudaFreeArray(texRarray);
  }

  void Load(unsigned int Zsize, unsigned int Rsize, float * fieldMapZ, float *fieldMapR)
  {
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

    cudaMallocArray(&texZarray, &channelDesc, Zsize, Rsize);
    cudaMemcpyToArray(texZarray, 0, 0, fieldMapZ, Zsize*Rsize*sizeof(float), cudaMemcpyHostToDevice);

    // set texture parameters
    texZ.addressMode[0] = cudaAddressModeWrap;
    texZ.addressMode[1] = cudaAddressModeWrap;
    texZ.filterMode = cudaFilterModeLinear;
    texZ.normalized = false;

    cudaBindTextureToArray(texZ, texZarray, channelDesc);

    int Rsize2 = Rsize * 2;
    int Zsize2 = (Zsize + 1) / 2;
    cudaMallocArray(&texRarray, &channelDesc, Rsize2, Zsize2);
    cudaMemcpyToArray(texRarray, 0, 0, fieldMapR, Rsize2*Zsize2*sizeof(float), cudaMemcpyHostToDevice);

    // set texture parameters
    texR.addressMode[0] = cudaAddressModeWrap;
    texR.addressMode[1] = cudaAddressModeWrap;
    texR.filterMode = cudaFilterModeLinear;
    texR.normalized = false;

    cudaBindTextureToArray(texR, texRarray, channelDesc);
  }
};
#endif
#endif
