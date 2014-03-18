#include "soa_kernel.h"
#include "stdio.h"
//struct of array (4-dimension)

__global__ void soa4_kernel(StructSoA4 soa, size_t numTracks)
{
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

  while (tid < numTracks) {
    (soa.x)[tid]  += 2.0;
    (soa.y)[tid]  += 3.0;
    (soa.z)[tid]  += 4.0;
    (soa.t)[tid]  += 5.0;
    //    printf("GPU (soa.t)[%d] = %f\n",tid,(soa.t)[tid]);
    tid += blockDim.x * gridDim.x;
  }
}

void soa4_gpu(StructSoA4 soa, size_t numTracks, int NB, int NT)
{
  soa4_kernel<<< NB, NT >>>(soa, numTracks);
}

void soa4_cpu(StructSoA4 soa, size_t numTracks)
{
  for (size_t tid = 0; tid < numTracks; tid++) {
    (soa.x)[tid]  += 2.0;
    (soa.y)[tid]  += 3.0;
    (soa.z)[tid]  += 4.0;
    (soa.t)[tid]  += 5.0;
    //   printf("CPU (soa.t)[%d] = %f\n",tid,(soa.t)[tid]);
  }
}

//struct of array (8-dimension)

__global__ void soa8_kernel(StructSoA8 soa, size_t numTracks)
{
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

  while (tid < numTracks) {
    (soa.x)[tid] += 2.0;
    (soa.y)[tid] += 3.0;
    (soa.z)[tid] += 4.0;
    (soa.t)[tid] += 5.0;
    (soa.u)[tid] += 2.0;
    (soa.v)[tid] += 3.0;
    (soa.w)[tid] += 4.0;
    (soa.s)[tid] += 5.0;
    tid += blockDim.x * gridDim.x;
  }
}

void soa8_gpu(StructSoA8 soa, size_t numTracks, int NB, int NT)
{
  soa8_kernel<<< NB, NT >>>(soa, numTracks);
}

void soa8_cpu(StructSoA8 soa, size_t numTracks)
{
  for (size_t tid = 0; tid < numTracks; tid++) {
    (soa.x)[tid] += 2.0;
    (soa.y)[tid] += 3.0;
    (soa.z)[tid] += 4.0;
    (soa.t)[tid] += 5.0;
    (soa.u)[tid] += 2.0;
    (soa.v)[tid] += 3.0;
    (soa.w)[tid] += 4.0;
    (soa.s)[tid] += 5.0;
  }
}

//array of struct (4-dimension)

__global__ void aos4_kernel(StructAoS4 *aos, size_t numTracks)
{
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

  while (tid < numTracks) {
    aos[tid].x  += 2.0;
    aos[tid].y  += 3.0;
    aos[tid].z  += 4.0;
    aos[tid].t  += 5.0;
    tid += blockDim.x * gridDim.x;
  }
}

void aos4_gpu(StructAoS4 *aos, size_t numTracks, int NB, int NT)
{
  aos4_kernel<<< NB, NT >>>(aos, numTracks);
}

void aos4_cpu(StructAoS4 *aos, size_t numTracks)
{
  for (size_t tid = 0; tid < numTracks; tid++) {
    aos[tid].x  += 2.0;
    aos[tid].y  += 3.0;
    aos[tid].z  += 4.0;
    aos[tid].t  += 5.0;
  }
}

//array of struct (8-dimension)

__global__ void aos8_kernel(StructAoS8 *aos, size_t numTracks)
{
  unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;

  while (tid < numTracks) {
    aos[tid].x += 2.0;
    aos[tid].y += 3.0;
    aos[tid].z += 4.0;
    aos[tid].t += 5.0;
    aos[tid].u += 2.0;
    aos[tid].v += 3.0;
    aos[tid].w += 4.0;
    aos[tid].s += 5.0;
    tid += blockDim.x * gridDim.x;
  }
}

void aos8_gpu(StructAoS8 *aos, size_t numTracks, int NB, int NT)
{
  aos8_kernel<<< NB, NT >>>(aos, numTracks);
}

void aos8_cpu(StructAoS8 *aos, size_t numTracks)
{
  for (size_t tid = 0; tid < numTracks; tid++) {
    aos[tid].x += 2.0;
    aos[tid].y += 3.0;
    aos[tid].z += 4.0;
    aos[tid].t += 5.0;
    aos[tid].u += 2.0;
    aos[tid].v += 3.0;
    aos[tid].w += 4.0;
    aos[tid].s += 5.0;
  }
}
