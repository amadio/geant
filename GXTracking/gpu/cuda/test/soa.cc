#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <ctime>
#include "soa_kernel.h"
#include "stdio.h"

using namespace std;

double Rndm(); 
void InitRandomState(unsigned seed=0);
float ElapsedTime(cudaEvent_t& start, cudaEvent_t& stop); 
void ProcessOneEvent(int ievt, size_t NTracks, 
		     unsigned int NBlocks, unsigned int NThreads);

int main(int argc, char* argv[]) {

  //default configuration
  unsigned int NBlocks  =  32;
  unsigned int NThreads = 128;
  size_t NTracks  = 32*128*16;

  if(argc >= 2) NBlocks  = atoi(argv[1]);
  if(argc >= 3) NThreads = atoi(argv[2]);
  if(argc >= 4) NTracks  = atoi(argv[3]);

  printf("Processing %d Tracks with (NBlocks,NThreads)= (%d,%d)\n",
	 NTracks,NBlocks,NThreads);

  //event-loop
  int NEvents = 10;

  for (int ievt = 0 ; ievt < NEvents ; ++ievt) {
    ProcessOneEvent(ievt, NTracks, NBlocks, NThreads);
  }

  return 0;
}

void ProcessOneEvent(int ievt, size_t NTracks,
		     unsigned int NBlocks, unsigned int NThreads) 
{
  cudaEvent_t start, stop;
  cudaEventCreate (&start);
  cudaEventCreate (&stop);

  // 1. Allocate host array

  //1.1 array of struct (AoS)

  StructAoS4 *aos4_c = (StructAoS4 *) malloc (NTracks*sizeof(StructAoS4));
  StructAoS4 *aos4_h = (StructAoS4 *) malloc (NTracks*sizeof(StructAoS4));
  StructAoS8 *aos8_c = (StructAoS8 *) malloc (NTracks*sizeof(StructAoS8));
  StructAoS8 *aos8_h = (StructAoS8 *) malloc (NTracks*sizeof(StructAoS8));

  //1.2 struct of array (SoA)

  StructSoA4 soa4_h;
  StructSoA4 soa4_c;
  StructSoA8 soa8_h;
  StructSoA8 soa8_c;

  G4double *h_x = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_y = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_z = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_t = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_u = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_v = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_w = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *h_s = (G4double *) malloc (NTracks*sizeof(G4double));

  G4double *c_x = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_y = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_z = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_t = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_u = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_v = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_w = (G4double *) malloc (NTracks*sizeof(G4double));
  G4double *c_s = (G4double *) malloc (NTracks*sizeof(G4double));

  //1.3 populate array elements

  for(int i = 0 ; i < NTracks ; ++i) {
    aos8_h[i].x = aos8_c[i].x = h_x[i] = c_x[i] = 1.0+10.0*Rndm();
    aos8_h[i].y = aos8_c[i].y = h_y[i] = c_y[i] = 2.0+10.0*Rndm();
    aos8_h[i].z = aos8_c[i].z = h_z[i] = c_z[i] = 3.0+10.0*Rndm();
    aos8_h[i].t = aos8_c[i].t = h_t[i] = c_t[i] = 4.0+10.0*Rndm();
    aos8_h[i].u = aos8_c[i].u = h_u[i] = c_u[i] = 1.0+10.0*Rndm();
    aos8_h[i].v = aos8_c[i].v = h_v[i] = c_v[i] = 2.0+10.0*Rndm();
    aos8_h[i].w = aos8_c[i].w = h_w[i] = c_w[i] = 3.0+10.0*Rndm();
    aos8_h[i].s = aos8_c[i].s = h_s[i] = c_s[i] = 4.0+10.0*Rndm();

    aos4_h[i].x = aos4_c[i].x = c_x[i];
    aos4_h[i].y = aos4_c[i].y = c_y[i];
    aos4_h[i].z = aos4_c[i].z = c_z[i];
    aos4_h[i].t = aos4_c[i].t = c_t[i];
  }

  soa8_c.x = c_x;
  soa8_c.y = c_y;
  soa8_c.z = c_z;
  soa8_c.t = c_t;
  soa8_c.u = c_u;
  soa8_c.v = c_v;
  soa8_c.w = c_w;
  soa8_c.s = c_s;

  soa4_c.x = c_x;
  soa4_c.y = c_y;
  soa4_c.z = c_z;
  soa4_c.t = c_t;

  // 2. Allocate device array

  //2.1 array of struct (AoS)

  StructAoS4 *aos4_d;
  StructAoS8 *aos8_d;
  cudaMalloc((void**)&aos4_d, NTracks*sizeof(StructAoS4));
  cudaMalloc((void**)&aos8_d, NTracks*sizeof(StructAoS8));
  cudaMemcpy(aos4_d,aos4_h, NTracks*sizeof(StructAoS4), cudaMemcpyHostToDevice);
  cudaMemcpy(aos8_d,aos8_h, NTracks*sizeof(StructAoS8), cudaMemcpyHostToDevice);

  //2.2 struct of array (SoA)
  G4double *d_x;
  G4double *d_y;
  G4double *d_z;
  G4double *d_t;
  G4double *d_u;
  G4double *d_v;
  G4double *d_w;
  G4double *d_s;

  cudaMalloc((void**) &(d_x), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_y), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_z), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_t), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_u), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_v), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_w), sizeof(G4double)*NTracks);
  cudaMalloc((void**) &(d_s), sizeof(G4double)*NTracks);

  //Copy array contents from host to device.
  cudaMemcpy(d_x, h_x, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_y, h_y, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_z, h_z, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_t, h_t, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_u, h_u, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_v, h_v, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_w, h_w, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);
  cudaMemcpy(d_s, h_s, sizeof(G4double)*NTracks, cudaMemcpyHostToDevice);

  //2.3 Point to device pointer in host struct.
  soa8_h.x = d_x;
  soa8_h.y = d_y;
  soa8_h.z = d_z;
  soa8_h.t = d_t;
  soa8_h.u = d_u;
  soa8_h.v = d_v;
  soa8_h.w = d_w;
  soa8_h.s = d_s;

  soa4_h.x = d_x;
  soa4_h.y = d_y;
  soa4_h.z = d_z;
  soa4_h.t = d_t;

  // 3. execute GPU kernels 

  //elapsed time for AoS4
  cudaEventRecord (start,0);
  aos4_gpu(aos4_d, NTracks, NBlocks, NThreads);
  float elapsedGPUAoS4 = ElapsedTime(start,stop);

  //D2H
  cudaMemcpy(aos4_h,aos4_d, sizeof(StructAoS4)*NTracks, cudaMemcpyDeviceToHost);

  //elapsed time for SoA4
  cudaEventRecord (start,0);
  soa4_gpu(soa4_h, NTracks, NBlocks, NThreads);
  float elapsedGPUSoA4 = ElapsedTime(start,stop);

  //copy pointer from device to host.
  cudaMemcpy(h_x, d_x, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_y, d_y, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_z, d_z, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_t, d_t, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);

  //elapsed time for AoS8
  cudaEventRecord (start,0);
  aos8_gpu(aos8_d, NTracks, NBlocks, NThreads);
  float elapsedGPUAoS8 = ElapsedTime(start,stop);

  //D2H
  cudaMemcpy(aos8_h,aos8_d, sizeof(StructAoS8)*NTracks, cudaMemcpyDeviceToHost);

  //elapsed time for SoA8
  cudaEventRecord (start,0);
  soa8_gpu(soa8_h, NTracks, NBlocks, NThreads);
  float elapsedGPUSoA8 = ElapsedTime(start,stop);

  //copy pointer from device to host.
  cudaMemcpy(h_x, d_x, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_y, d_y, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_z, d_z, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_t, d_t, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_u, d_u, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_v, d_v, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_w, d_w, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);
  cudaMemcpy(h_s, d_s, sizeof(G4double)*NTracks, cudaMemcpyDeviceToHost);

  // 4. execute CPU programs

  //elapsed time for CPU AoS4
  cudaEventRecord (start,0);
  aos4_cpu(aos4_c, NTracks);
  float elapsedCPUAoS4 = ElapsedTime(start,stop);

  //elapsed time for CPU SoA4
  cudaEventRecord (start,0);
  soa4_cpu(soa4_c, NTracks);
  float elapsedCPUSoA4 = ElapsedTime(start,stop);

  //elapsed time for CPU AoS8
  cudaEventRecord (start,0);
  aos8_cpu(aos8_c, NTracks);
  float elapsedCPUAoS8 = ElapsedTime(start,stop);

  //elapsed time for CPU SoA8
  cudaEventRecord (start,0);
  soa8_cpu(soa8_c, NTracks);
  float elapsedCPUSoA8 = ElapsedTime(start,stop);

  // 5. print performance

  printf("Event> %d elapsedTimeAoS4 (GPU CPU Ratio) = ( %f %f %f ) ms\n",
	 ievt, elapsedGPUAoS4,elapsedCPUAoS4,elapsedCPUAoS4/elapsedGPUAoS4);
  printf("Event> %d elapsedTimeSoA4 (GPU CPU Ratio) = ( %f %f %f ) ms\n",
	 ievt, elapsedGPUSoA4,elapsedCPUSoA4,elapsedCPUSoA4/elapsedGPUSoA4);

  printf("Event> %d elapsedTimeAoS8 (GPU CPU Ratio) = ( %f %f %f ) ms\n",
	 ievt, elapsedGPUAoS8,elapsedCPUAoS8,elapsedCPUAoS8/elapsedGPUAoS8);
  printf("Event> %d elapsedTimeSoA8 (GPU CPU Ratio) = ( %f %f %f ) ms\n",
	 ievt, elapsedGPUSoA8,elapsedCPUSoA8,elapsedCPUSoA8/elapsedGPUSoA8);

  // 6. validation
  //  for(int i = 0 ; i < NTracks ; ++i) {
  //    printf("t(GPU:SoA,GPU:AoS,CPU:SoA,CPU:AoS) = (%f,%f,%f,%f)\n",
  //               h_t[i],aos4_h[i].t,(soa4_c.t)[i],aos4_c[i].t);
  //  }

  //clean up

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  free(aos4_h);
  free(aos4_c);
  free(aos8_h);
  free(aos8_c);
  cudaFree(aos4_d);
  cudaFree(aos8_d);

  free(h_x);
  free(h_y);
  free(h_z);
  free(h_t);
  free(h_u);
  free(h_v);
  free(h_w);
  free(h_s);

  free(c_x);
  free(c_y);
  free(c_z);
  free(c_t);
  free(c_u);
  free(c_v);
  free(c_w);
  free(c_s);

  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_z);
  cudaFree(d_t);
  cudaFree(d_u);
  cudaFree(d_v);
  cudaFree(d_w);
  cudaFree(d_s);

}

double Rndm() 
{ //random (0,1]
  return ( static_cast<double> (rand()) )/(RAND_MAX);
}

void InitRandomState(unsigned seed)
{ //set the initial seed for random number generator
  if( seed==0 ) {
    time_t curtime;    
    time(&curtime);    
    srand(static_cast<unsigned> (curtime));
  }
  else {
    srand(seed);
  }
}

float ElapsedTime(cudaEvent_t& start, cudaEvent_t& stop) 
{
  float elapsedTime = 0.0;
  cudaEventRecord (stop,0);
  cudaEventSynchronize (stop);
  cudaEventElapsedTime (&elapsedTime,start,stop);
  return elapsedTime; 
}
