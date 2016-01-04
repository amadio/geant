
#ifndef _TOYCLASS_H_
#define _TOYCLASS_H_

#include <iostream>
#include "AlignedBase.h"
#include <Vc/Vc>
#include <stdlib.h>
#include "base/Vector3D.h" //To utlimately include kVc, kScalar backends etc


using namespace std;

class ToyClass : public AlignedBase
{
  typedef typename vecgeom::kVc::bool_v      Bool_v;
  typedef typename vecgeom::kVc::precision_v Double_v;
  typedef typename vecgeom::kVc::int_v       Int_v;
public:
  ToyClass(){};

  //used to set total number of tracks user might give
  //Can be modified later. Makes things easier for now.
  ToyClass(int n);

  ~ToyClass();

  //Main method here
  int ToyMethod(double nTracks[], double finalResults[]);

  //Called by ToyMethod; kinda like OneGoodStep
  void SimpleMethod(Double_v &fPreProcVc, 
                    Double_v &outputSimpleMethod);

  //Called in constructor. Set input no. of tracks
  void SetInputTracks(int n);

  //Initialize the Vc vector in beginning if not initialized. Assume 
  //uninitialized here. Later on replacement when work done for it 
  void InitializeVcVector(double nTracks[]);

  void PreProcess(double nTracks[], 
                  double fPreprocData[]);

  //Print statements for debugging
  void PrintCurrentData(Double_v fPreProcVc, 
                        Double_v outputSimpleMethod,
                        Bool_v Done                );

private:
  int kSizeOfVcVector = 4; //can use vecgeom::kVectorSize
  int fInputTotalTracks;
  double hstep = 0.1;

  //Stores indices, used to take care of sth sth.. random drops and pickups
  // Int_v fIndex;
  int *fIndex;

  //Stores scalar preprocessed data 
  //probably needs to be vectorized but leave that for later as of now 
  //this is kinda mix of scalar and vector version
  //but let it be for the time being 
  double *fPreprocData;

  Int_v    fNoGoodSteps;
  Double_v fYVcVector;
  Double_v fPreProcVc;
  Double_v fState;

};


ToyClass::ToyClass(int n)
{
  SetInputTracks(n);
  fPreprocData = new double[n];
  fIndex       = new int[kSizeOfVcVector];
}

ToyClass::~ToyClass()
{
  delete fPreprocData;
  cout<<"----ToyClass destructor being called"<<endl;
}

void ToyClass::SetInputTracks(int n)
{
  fInputTotalTracks = n;
}

void ToyClass::InitializeVcVector(double nTracks[])
{
  for (int i = 0; i < kSizeOfVcVector; ++i)
  {
    fPreProcVc[i] = nTracks[i];
    fIndex    [i] = i;    
    fState    [i] = hstep; 
  }
}

void ToyClass::PreProcess(double nTracks[], double fPreprocData[])
{  
  double factor = 0.4;

  for (int i = 0; i < fInputTotalTracks; ++i)
  {
    fPreprocData[i] = nTracks[i] * factor;
  }
}

void ToyClass::PrintCurrentData(typename vecgeom::kVc::precision_v fPreProcVc, 
                                typename vecgeom::kVc::precision_v outputSimpleMethod, 
                                typename vecgeom::kVc::bool_v      Done                )
{
  cout<<"\n----Currently:"<<endl;
  cout<<"----fPreProcVc is         : "<<fPreProcVc<<endl;
  cout<<"----outputSimpleMethod is: "<<outputSimpleMethod<<endl;
  cout<<"----Done is              : "<<Done<<endl;
}

void ToyClass::SimpleMethod(typename vecgeom::kVc::precision_v &fPreProcVc, 
                            typename vecgeom::kVc::precision_v &outputSimpleMethod)
{
  // srand(time(NULL));
  for (int i = 0; i < kSizeOfVcVector; ++i)
  {
    outputSimpleMethod[i] =  (float) rand()/(RAND_MAX) ;
  }
  // outputSimpleMethod = fPreProcVc/2.; //some random processing
  fPreProcVc /= 3.;

}


int ToyClass::ToyMethod(double nTracks[], double finalResults[])
{
  PreProcess(nTracks, fPreprocData);

  InitializeVcVector(fPreprocData);

  int trackNextInput = 4; 

  //now start working from preprocessed data

  //What simple work to do?
  //Say if number < 0.5
  // add 1 to the number and store it in finalResults
  //else keep subtracting .1 and check
  //push in new tracks into the empty slot

  typedef typename vecgeom::kVc::bool_v Bool_v;
  typedef typename vecgeom::kVc::precision_v Double_v;
  Bool_v isDone(0.), isGoodStep(0.);
  Double_v outputSimpleMethod;
  // PrintCurrentData(fPreProcVc, outputSimpleMethod, isDone);
  // SimpleMethod(fPreProcVc, outputSimpleMethod);
  // isDone = outputSimpleMethod < fState;
  // PrintCurrentData(fPreProcVc, outputSimpleMethod, isDone);

  // do{
  // trackNextInput = 16;
  while(!(vecgeom::IsFull(isDone)) || (trackNextInput < fInputTotalTracks))
  {
    // cout<<"\n----trackNextInput is: "<<trackNextInput<<endl;
    // if( !(vecgeom::IsFull(isDone)) ) // && trackNextInput < 16 )    
    // {
      SimpleMethod(fPreProcVc, outputSimpleMethod);
      isDone = outputSimpleMethod < fState;
      PrintCurrentData(fPreProcVc, outputSimpleMethod, isDone);
    // }

    //return if all are done
    // if(vecgeom::IsFull(isDone)) return;

    //what if all are not done?
    //if none is done, then simply continue
    //But if say one is done, then we need to insert new tracks in it's position
    //This executes only if at least one is done 
    //Do it only if anything is left in nTracks, otherwise no extra tracks left
    // cout<<  !(vecgeom::IsFull(isDone)) <<" "<< vecgeom::IsEmpty(isDone) <<endl;
    // if( !(vecgeom::IsFull(isDone)) && !(vecgeom::IsEmpty(isDone)))
    if(!(vecgeom::IsEmpty(isDone)))
    {
      cout<<"here1"<<endl;
      for (int i = 0; i < kSizeOfVcVector; ++i)
      {
        cout<<"here2"<<endl;
        if(isDone[i]==1 && fIndex[i] != -1) 
        {
          cout<<"here3"<<endl;
          finalResults[fIndex[i]] = outputSimpleMethod[i]; //store the output

          if(trackNextInput < fInputTotalTracks)
          {
            cout<<"\n----trackNextInput is: "<<trackNextInput<<endl;
            // cout<<"here4"<<endl;
            //insert new track
            //And store output of this one in finalResults (in scalar array for now)
            //Probably need to define an index to take care of these random pickups and drops
            //Size of index = 16/n
            fPreProcVc [i]          = nTracks[trackNextInput]; //sending in next one
            fIndex     [i]          = trackNextInput;
            isDone[i]     = 0;
            trackNextInput++;
          }
          else
          {
            cout<<"here5"<<endl;
            fIndex[i]               = -1;
            //need to stop processing for other things or do something else of the kind
          }
        }
      }
    }

/*    if( !(vecgeom::IsFull(isDone)) ) // && trackNextInput < 16 )    
    {
      SimpleMethod(fPreProcVc, outputSimpleMethod);
      isDone = outputSimpleMethod < 0.1;
      PrintCurrentData(fPreProcVc, outputSimpleMethod, isDone);
    }*/

/*    SimpleMethod(fPreProcVc, outputSimpleMethod);
    isDone = outputSimpleMethod < 0.1;
    PrintCurrentData(fPreProcVc, outputSimpleMethod, isDone);*/

  }
  // } while(!(vecgeom::IsFull(isDone)) || (trackNextInput < fInputTotalTracks));

  cout<<"----Exited while loop"<<endl;

  //All are done. Store data if not stored already
  //now being taken care of in the do while loop since data is stored after isDone
/*  for (int i = 0; i < kSizeOfVcVector; ++i)
  {
    cout<<"fIndex is: "<<fIndex[i]<<endl;
    if(fIndex[i] != -1)
    {
      cout<<"here?"<<endl;
      finalResults[fIndex[i]] = outputSimpleMethod[i];
    }
  }
*/

  cout<<"----here we are at the end"<<endl;
  cout<<"\n----Output obtained is: "<<endl;
  for (int i = 0; i < fInputTotalTracks; ++i)
  {
    cout<<i<<" "<<finalResults[i]<<" "<<endl;
  }

  return 1;
}




#endif