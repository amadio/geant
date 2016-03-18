#include "MaterialHandler.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

MaterialHandler* MaterialHandler::fInstance = 0;

VECPHYS_CUDA_HEADER_HOST
MaterialHandler* MaterialHandler::Instance()
{
  if (fInstance == 0) fInstance = new MaterialHandler();
  return fInstance;
}

VECPHYS_CUDA_HEADER_HOST
MaterialHandler::MaterialHandler() 
{
  //test mode: 0 all elements in the element table, 1 single element
  fElementMode = 0; 

  //initialize the element array
  fNumberOfElements = 0;
  for(int i = 0 ; i < maximumZ ; ++i)  fElementArray[i] = 0;

  //build the element array
  BuildElementTable();
}

VECPHYS_CUDA_HEADER_HOST
MaterialHandler::~MaterialHandler() 
{
  ;
}

VECPHYS_CUDA_HEADER_HOST
void MaterialHandler::BuildElementTable() 
{
  //This should interface with the global material manager of GeantV so that 
  //the element arrary is properly filled with all elements of detector
  //materials. Temporarily, build a table based on John's arrary
  
  constexpr int NumFx = 16;
  int Element[NumFx] = { 82, 74, 8, 7, 6, 13, 18, 22, 26, 27, 30, 48, 54, 64, 79, 91};  
                      // Pb   W  O  N  C  Al  Ar, Ti  Fe  Cu  Zn  Cd  Xe  Gd  Au  Pa

  for(int ie = 0 ; ie < NumFx ; ++ie) AddElement(Element[ie]);   

}

VECPHYS_CUDA_HEADER_HOST
void MaterialHandler::AddElement(int element) {
  //check validity of the element
  if(element > 0 && element < maximumZ) {
    //seach whether this element already exists
    bool found = false;
    for(int i = 0 ; i < fNumberOfElements ; ++i) {
      if (fElementArray[i] == element) {
        found = true;
        break;
      }    
    } 

    //add a new element to the array
    if(!found) { 
      fElementArray[fNumberOfElements] = element;
      fNumberOfElements++;
    }
  }
}

VECPHYS_CUDA_HEADER_HOST
void MaterialHandler::PrepareTargetElements(int *targetElements, int ntracks, int elementMode)
{
  //only two modes for now based on John's original method
  static int noCalls=0 ;
  noCalls++;
  
  bool report = (noCalls == 1 );
  
  if(elementMode == 0 ) { // all elements in the material table
    if( report )
      printf(" Generating Target Elements from table of %d elements - mode # =  %d\n",
             fNumberOfElements,elementMode);
    int indEl;
    for(int i = 0 ; i < ntracks ; ++i) {
      indEl = ( i % fNumberOfElements ) ;
      targetElements[i] = fElementArray[ indEl ]; 
    }
  }
  else if( elementMode == 1 ) { // using a single element
    if( report ) 
      printf(" Using *Constant* Target Element Z = %d - mode # = %d\n",
	     fElementArray[0],elementMode);
    
    for(int i = 0 ; i < ntracks ; ++i) {
      targetElements[i] = fElementArray[0];
    }
  }
  else {
    printf(" Illeagal - mode # = %d\n",elementMode);
    assert(0);
  }
}

} // end namespace impl
} // end namespace vecphys
