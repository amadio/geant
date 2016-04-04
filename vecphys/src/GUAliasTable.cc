#include "GUAliasTable.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST_DEVICE
GUAliasTable::GUAliasTable(int ngrid)
{
  fNGrid = ngrid;
  Allocate(ngrid);
}

VECCORE_CUDA_HOST_DEVICE
GUAliasTable::GUAliasTable(const GUAliasTable& table)
{
  Deallocate();
  Allocate(table.fNGrid);
  CopyData(table);
}

VECCORE_CUDA_HOST_DEVICE
GUAliasTable& GUAliasTable::operator=(const GUAliasTable& table)
{
  if (this != &table) {
    Deallocate();
    Allocate(table.fNGrid);
    CopyData(table);
  }
  return *this;
}

VECCORE_CUDA_HOST_DEVICE
GUAliasTable::~GUAliasTable()
{
  Deallocate();
}

VECCORE_CUDA_HOST_DEVICE
void GUAliasTable::Allocate(int ngrid) {
  //try
  fpdf   = new Precision [fNGrid];
  fProbQ = new Precision [fNGrid];
  fAlias = new int    [fNGrid];

  for (int i = 0 ; i < fNGrid ; ++i) {
    fpdf[i]   = -1;
    fProbQ[i] = -1;
    fAlias[i] = -1;
  }
}

VECCORE_CUDA_HOST_DEVICE
void GUAliasTable::Deallocate() {
  delete [] fpdf;
  delete [] fProbQ;
  delete [] fAlias;
}

VECCORE_CUDA_HOST_DEVICE
void GUAliasTable::CopyData(const GUAliasTable& table) {
  int i;
  fNGrid = table.fNGrid;
  for(i = 0 ; i < fNGrid  ; ++i) {
    fpdf[i] =   table.fpdf[i];
    fProbQ[i] = table.fProbQ[i];
    fAlias[i] = table.fAlias[i];
  }
}

VECCORE_CUDA_HOST_DEVICE
int GUAliasTable::SizeOfTable() {
  return sizeof(int) + SizeOfGrid()*(2.*sizeof(Precision)+sizeof(int));
}

VECCORE_CUDA_HOST_DEVICE
void GUAliasTable::PrintInfo() {
  printf("Size(NGrid,Table) = (%d,%d)\n",SizeOfGrid(),SizeOfTable());
}

#ifdef VECCORE_NVCC
void GUAliasTable::Relocate(void *devPtr)
{
  //Implement/use a general way to (byte-wise) copy a object to GPU
  Precision *d_fpdf;
  Precision *d_fProbQ;
  int       *d_fAlias;

  cudaMalloc((void**) &(d_fpdf)  , sizeof(Precision)*fNGrid);
  cudaMalloc((void**) &(d_fProbQ), sizeof(Precision)*fNGrid);
  cudaMalloc((void**) &(d_fAlias), sizeof(int)*fNGrid);

  //Copy array contents from host to device.
  cudaMemcpy(d_fpdf,  fpdf,   sizeof(Precision)*fNGrid, cudaMemcpyHostToDevice);
  cudaMemcpy(d_fProbQ,fProbQ, sizeof(Precision)*fNGrid, cudaMemcpyHostToDevice);
  cudaMemcpy(d_fAlias,fAlias, sizeof(int)*fNGrid, cudaMemcpyHostToDevice);

  //point to device pointer in host struct.
  fpdf =   d_fpdf;
  fAlias = d_fAlias;
  fProbQ = d_fProbQ;

  cudaMemcpy(devPtr, this, SizeOfTable(), cudaMemcpyHostToDevice);
}

#endif

} // end namespace impl
} // end namespace vecphys
