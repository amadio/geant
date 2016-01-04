#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <ctime>

#include "base/Global.h"

#include <Vc/Vc>
#include "backend/vc/Backend.h"
#include "backend/Backend.h"
//#include "GatherAlias.h"
#include "ReadVectorData.h"

#include "CMSmagField.h"

using namespace std;
/*
class ReadVectorData
{
public: 
  vector<double> fRadius, fPhi, fZ, fBr, fBz, fBphi;
  ReadVectorData(string inputMap){
    dataFile = inputMap;
    PleaseReadData();
  };
  
  ~ReadVectorData(){}

private:
  string dataFile;
  void PleaseReadData(){
    string line;
    string s1,s2,s3,s4,s5,s0;
    double d1,d2,d3,d4,d5,d0;
    ifstream pFile(dataFile);
    if (pFile.is_open()){
      while(getline(pFile,line)){
        stringstream ss(line);
        ss>> s0>> d1>> s1>> d0>> s2>> d2>> s3>> d3>> s4>> d4>> s5>> d5;
        fRadius.push_back(d1);
        fPhi.push_back(d0);
        fZ.push_back(d2);
        fBz.push_back(d3);
        fBr.push_back(d4);
        fBphi.push_back(d5);
      }
      pFile.close();
    }
    else{
      cout<<"Unable to open file"<<endl;
    }
  }

};



template<>
void GatherAlias(typename kVc::Index_t    index, typename kVc::Double_t  &probNA){
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  //vecgeom::kVectorSize == vecphys::kVc::kSize
  ReadVectorData dataMap("/home/ananya/Work/MagFieldRoutine/cms2015.txt");
  for(int i = 0; i < vecgeom::kVectorSize ; ++i) 
  {
    int ind = index[i];

    if(ind < 0) {
      // printf("Warning: negative index! - in GUPhotoElectronSauterGavrila\n");
      ind = 0;
    }

    //assert( z > 0  && z <= fMaxZelement );
    //    assert( ind >= 0 && ind < fAliasTable[z]->SizeOfGrid() );
  
  //Say index has 5 indices a,b,c,d,e i.e. 5 particles
  //ind=a.  Now for first particle, need to fill probNA with Br[a.. .. e]
  //so say
  probNA[i] = dataMap.fBr[ind];
  //... but access to fBr needed which we will have obviously
  


    // int tableIndex = fAliasTableManager->GetTableIndex(z);
    // probNA[i]=   (fAliasTableManager->GetAliasTable(tableIndex))->fProbQ[ ind ];
    // aliasInd[i]= (fAliasTableManager->GetAliasTable(tableIndex))->fAlias[ ind ];
  }
};

*/
int main(){

  ReadVectorData r1("/home/ananya/Work/MagFieldRoutine/cms2015.txt");

  typedef vecgeom::kVc::precision_v Double_t;
  typedef vecgeom::kVc::Index_t Index_t;
  typedef vecgeom::kVc::Int_t Int_t;
  Index_t i1; 
  Double_t probNAtest;
  srand(time(NULL));
  for (int i = 0; i < vecgeom::kVectorSize; ++i)
  {
    i1[i] = i + ((double) rand()/(RAND_MAX))*100;
  }

  cout<<"i1 randomizes is: "<<i1<<endl;
  cout<<"probNAtest before gatheralias is: "<<probNAtest<<endl;
    int count  = 1e+8;
    for (int i = 0; i < count; ++i)
    {
       r1.GatherAlias<vecgeom::kVc>(i1, probNAtest);
    }
  
  cout<<"probNAtest after gatheralias is: "<<probNAtest<<endl;
  // typedef vecgeom::Vector3D<vecgeom::Precision> Vec_t;
  // typedef vecgeom::Vector3D<vecgeom::kVc::precision_v> VEC_t;


  // int temp_arr[150];
  // double *arr = new double[150];
  // for (int i = 0; i < 150; ++i)
  // {
  //   temp_arr[i] = i;
  // }

  //vecgeom::gather(arr, temp_arr);

  vecgeom::kVc::precision_v v1; 
  double arr2[] = {12,14,16,18,20};
  int ind2[] = {0,1,2,3,4,5};
  cout<<"v1 before gather: "<<v1<<endl;
  v1.gather(arr2, ind2);
  cout<<"v1 after gather: "<<v1<<endl;

  Int_t i2;
  i2[0] = 2;
  i2[1] = 3;
  i2[2] = 1;
  i2[3] = 0;
  cout<<"i2 looks like: "<<i2<<endl;
  v1.gather(arr2, i2);
  cout<<"v1 after new gather: "<<v1<<endl;

  Int_t i3 = (Int_t) i1;
  v1.gather(arr2,i3);
  cout<<"i3 looks like: "<<i3<<endl;
  cout<<"v1 after 3rd gather is: "<<v1<<endl;

  vecgeom::kVc::precision_v v3; 
  vecgeom::kVc::precision_v arr3; 
  arr3[0] = 5.0;
  arr3[1] = 51.0;
  arr3[2] = 115.0;
  arr3[3] = 15.0;

  // typedef vecgeom::Vector3D<vecgeom::kVc::precision_v> ThreeVector;
  // ThreeVector v4;
  // v4.gather(arr2, i3);
  // cout<<"v3 now looks like: "<<v3<<endl;


  return 0;
}

/*

// Scalar method - to be used ONLY for 'scalar-type' backends
//                  i.e. currently: scalar & CUDA
template<class Backend>
inline
void GUAliasSampler::
GatherAlias(typename Backend::Index_t    index,
            typename Backend::Index_t    zElement,
            typename Backend::Double_t  &probNA,
            typename Backend::Double_t  &aliasInd
           ) const
{
#ifdef CHECK
  if( zElement <= 0  || zElement > fMaxZelement )
  {
    printf(" Illegal zElement = %d\n",zElement);
  }
#endif
  //  assert( (zElement > 0)  && (zElement <= fMaxZelement) );

  int     intIndex= (int) index;
   
#ifdef CHECK
  //  int     tableSize= fAliasTable[zElement]->SizeOfGrid();
  //  INDEX_CHECK( intIndex, 0, tableSize, "Index", "TableSize" );
#endif
  //  assert( (intIndex >= 0) && (intIndex < tableSize) );

  int tableIndex = fAliasTableManager->GetTableIndex(zElement);
  probNA=   (fAliasTableManager->GetAliasTable(tableIndex))->fProbQ[ intIndex ];
  aliasInd= (fAliasTableManager->GetAliasTable(tableIndex))->fAlias[ intIndex ];
}

// Specialisation for all vector-type backends - Vc for now
#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::
GatherAlias<kVc>(typename kVc::Index_t    index, 
                 typename kVc::Index_t    zElement,
                 typename kVc::Double_t  &probNA,  
                 typename kVc::Double_t  &aliasInd

                ) const 
{
  //gather for alias table lookups - (backend type has no ptr arithmetic)
  for(int i = 0; i < kVc::kSize ; ++i) 
  {
    int z= zElement[i];
    int ind = index[i];

    if(ind < 0) {
      // printf("Warning: negative index! - in GUPhotoElectronSauterGavrila\n");
      ind = 0;
    }

    assert( z > 0  && z <= fMaxZelement );
    //    assert( ind >= 0 && ind < fAliasTable[z]->SizeOfGrid() );

    int tableIndex = fAliasTableManager->GetTableIndex(z);
    probNA[i]=   (fAliasTableManager->GetAliasTable(tableIndex))->fProbQ[ ind ];
    aliasInd[i]= (fAliasTableManager->GetAliasTable(tableIndex))->fAlias[ ind ];
  }
}
#endif

*/
