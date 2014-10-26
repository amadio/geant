#ifndef VTfileio_h
#define VTfileio_h

#include <string>

class TFile;
class TTree;
#include "TObjArray.h"
#include "TObjString.h"

#include "G4String.hh"
 
class VTfileio {
 public:
   static VTfileio *I() {if(!fgInstance) fgInstance=new VTfileio; return fgInstance;}
   ~VTfileio();
   void NewTree(const char* name);
   TFile* OutFile() const {return fOutFile;}
   void OutFile(TFile* f) {fOutFile=f;}
   TTree* CurTree() const {return fCurTree;}
   void CurTree(TTree* f) {fCurTree=f;}
   int VolumeIndex(const char* volname) const;
   int ProcessIndex(const char* procname) const;
   TObjArray *GetVolumeDictionary() {return fVolumeDictionary;}
   TObjArray *GetProcessDictionary() {return fProcessDictionary;}
   void WriteDictionaries();
   void Fill(double x, double y, double z, double px, double py, double pz, Short_t pid,
	     UShort_t lvid, double safety, double snext, double step, UChar_t surfid, 
	     UChar_t process, UChar_t begend, UInt_t trid, UInt_t trpid);
   Bool_t IsNewEvent() {if(fNewEvent) {fNewEvent=kFALSE; return kTRUE;} 
      else return kFALSE;}
   void SetNewEvent() {fNewEvent=kTRUE;}
   void SetPrimaries(int prim) {fNprimaries=prim;}
   int  GetPrimaries() const {return fNprimaries;}
 private:
   VTfileio();
   static VTfileio *fgInstance;
   TFile *fOutFile; // output file
   TTree* fCurTree; // current tree
   double fX;            // x position
   double fY;            // y position
   double fZ;            // z position
   double fPx;           // x momentum
   double fPy;           // y momentum
   double fPz;           // z momentum
   Short_t fPID;         // PDG particle id
   UShort_t  fLVid;      // logical volume id
   double fSafety;       // safety
   double fSnext;        // snext
   double fStep;         // step
   UChar_t fSurfid;      // surface id
   UChar_t fProcess;     // Process
   UChar_t fBegEnd;      // Beginning or end of track
   UInt_t  fTrid;        // Track ID
   UInt_t  fTrPid;       // Track Parend ID
   //
   Bool_t  fNewEvent;    // if new event
   Int_t   fNprimaries;  // Number of primaries
   TObjArray *fVolumeDictionary; // dictionary of volumes
   TObjArray *fProcessDictionary; // dictionary of processes
};

#endif
