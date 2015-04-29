#include "GUAliasTableManager.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
GUAliasTableManager::GUAliasTableManager(int numberOfElement) 
  : fNElement(0)
{
  for (int i= 0; i < maximumZ ; ++i) fIndex[i] = -1;   
  fAliasTables = (GUAliasTable**)malloc(numberOfElement*sizeof(GUAliasTable*)); 
}

VECPHYS_CUDA_HEADER_HOST
GUAliasTableManager::~GUAliasTableManager() 
{
  free(fAliasTables);
}

VECPHYS_CUDA_HEADER_HOST
void GUAliasTableManager::SetTableIndex(int Z) { 
  assert(Z > -1 && Z < maximumZ );
  if(fIndex[Z] == -1 ) {
    fIndex[Z] = fNElement;
    fNElement++;
  }
}

VECPHYS_CUDA_HEADER_BOTH
GUAliasTable* GUAliasTableManager::GetAliasTable(int Z) { 
  return fAliasTables[GetTableIndex(Z)];
}

VECPHYS_CUDA_HEADER_BOTH
int GUAliasTableManager::GetTableIndex(int Z) { 
  assert(Z > -1 && Z < maximumZ ) ;
  return fIndex[Z];
}

VECPHYS_CUDA_HEADER_BOTH
int GUAliasTableManager::GetNumberOfElements() {
  return fNElement;
}

VECPHYS_CUDA_HEADER_HOST
void GUAliasTableManager::AddAliasTable(int Z, GUAliasTable* table) { 
  assert( Z > -1 && Z < maximumZ ) ;
  fAliasTables[fNElement] = table;
  SetTableIndex(Z);  
}

VECPHYS_CUDA_HEADER_HOST
int GUAliasTableManager::SizeOfManager() { 
  return sizeof(GUAliasTable*)*fNElement + (maximumZ+1)*sizeof(int);
}

#ifdef VECPHYS_NVCC
void GUAliasTableManager::Relocate(void *devPtr) 
{
  //device pointers in device memory
  GUAliasTable **fAliasTables_d;     
  cudaMalloc((void**)&fAliasTables_d, sizeof(GUAliasTable*)*fNElement);

  //device pointers in host memory
  GUAliasTable* tables_d[fNElement]; 

  // relocate pointers of this to the corresponding device pointers 
  for(int i = 0 ; i < fNElement ; i++) {
    cudaMalloc((void**)&tables_d[i], fAliasTables[i]->SizeOfTable());
    fAliasTables[i]->Relocate(tables_d[i]);
  }

  // copy the pointer to alias table pointers from the host to the device
  cudaMemcpy(fAliasTables_d, tables_d, sizeof(GUAliasTable*)*fNElement, 
	     cudaMemcpyHostToDevice);
  fAliasTables = fAliasTables_d;

  // copy the manager from host to device.
  cudaMemcpy(devPtr, this, SizeOfManager(), cudaMemcpyHostToDevice);
}
#endif

} // end namespace impl
} // end namespace vecphys
