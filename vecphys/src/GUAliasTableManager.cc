#include "GUAliasTableManager.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST
GUAliasTableManager::GUAliasTableManager(int nelements, int /*ngrid*/)
  : fNElement(0)
{
  for (int i = 0; i < maximumZ ; ++i) fIndex[i] = -1;
  //better to have another member for the total number of elements = nelements
  fAliasTables = new GUAliasTable*[nelements];

}

VECCORE_CUDA_HOST
GUAliasTableManager::~GUAliasTableManager()
{
  for(int i = 0 ; i < fNElement ; ++i) delete fAliasTables[i];
  delete [] fAliasTables;
}

VECCORE_CUDA_HOST
void GUAliasTableManager::SetTableIndex(int Z) {
  assert(Z > -1 && Z < maximumZ );
  if(fIndex[Z] == -1 ) {
    fIndex[Z] = fNElement;
    fNElement++;
  }
 }

VECCORE_CUDA_HOST_DEVICE
GUAliasTable* GUAliasTableManager::GetAliasTable(int Z) {
  return fAliasTables[GetTableIndex(Z)];
}

VECCORE_CUDA_HOST_DEVICE
int GUAliasTableManager::GetTableIndex(int Z) {
  assert(Z > -1 && Z < maximumZ ) ;
  return fIndex[Z];
}

VECCORE_CUDA_HOST_DEVICE
int GUAliasTableManager::GetNumberOfElements() {
  return fNElement;
}

VECCORE_CUDA_HOST
void GUAliasTableManager::AddAliasTable(int Z, GUAliasTable* table) {
  assert( Z > -1 && Z < maximumZ ) ;

  fAliasTables[fNElement] = table;
  SetTableIndex(Z);
}

VECCORE_CUDA_HOST
int GUAliasTableManager::SizeOfManager() {
  return sizeof(GUAliasTable*)*fNElement + (maximumZ+1)*sizeof(int);
}

#ifdef VECCORE_NVCC
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
