#include <stdio.h>

void deviceProperty(int idevice);

int deviceQuery() {

  printf("CUDA Device Query\n");  

  int count = 0;
  cudaGetDeviceCount(&count);

  if (count == 0) {
    printf("There is no CUDA capable devices\n");
  }
  else {
    printf("There is %d CUDA capable devices\n",count);
    for (int id = 0 ; id < count ; id++) deviceProperty(id);
  }

  return 0;
}

void deviceProperty(int id) {

  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop,id);
  
  printf("Device %d Name        : %s\n",id,prop.name);
  printf("Capability Maj/Minor : %d %d\n",prop.major, prop.minor);
  printf("Clock Rate           : %d\n",prop.clockRate);
  printf("Global Memory        : %ld\n",prop.totalGlobalMem);
  printf("Constant Memory      : %ld\n",prop.totalConstMem);
  printf("Memory Pitch         : %ld\n",prop.memPitch);
  printf("Texture Alignment    : %ld\n",prop.textureAlignment);
  printf("MultiProcess count   : %d\n",prop.multiProcessorCount);
  printf("Shared Memory/Block  : %ld\n",prop.sharedMemPerBlock);
  printf("Registers/Block      : %ld\n",prop.regsPerBlock);
  printf("Threads in Warp      : %d\n",prop.warpSize);
  printf("Max Thread per Block : %d\n",prop.maxThreadsPerBlock);
  printf("Max Thread Dimension : (%d,%d,%d)\n",prop.maxThreadsDim[0],
  	                         prop.maxThreadsDim[1],prop.maxThreadsDim[2]);
  printf("Max Grid Size        : (%d,%d,%d)\n",prop.maxGridSize[0],
  	                         prop.maxGridSize[1],prop.maxGridSize[2]);
}