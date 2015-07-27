#ifndef VTfileio_h
#define VTfileio_h

#include <string>

class TFile;
class TTree;
#include "THashList.h"
#include "TNamed.h"

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

   inline int VolumeIndex(const char* volname) const {
      TObject *obj=fVolumeDictionary->FindObject(volname);
      if(obj) return obj->GetUniqueID();
      else return -1;
   }

   inline int ShapeIndex(const char* shapename) const {
      TObject *obj=fShapeDictionary->FindObject(shapename);
      if(obj) return obj->GetUniqueID();
      else return -1;
   }

   inline int ProcessIndex(const char* procname) {
      TObject *obj=fProcessDictionary->FindObject(procname);
      if(obj) return obj->GetUniqueID();
      int nproc = fProcessDictionary->GetEntries();
      TNamed *nam = new TNamed(procname,procname);
      nam->SetUniqueID(nproc);
      fProcessDictionary->Add(nam);
      return nproc;
   }

   int ProcessIndex(const char* procname) const;
   THashList *GetVolumeDictionary() {return fVolumeDictionary;}
   void AddVolume(const char *volname);
   void AddShape(const char *shapename);
   THashList *GetProcessDictionary() {return fProcessDictionary;}
   void WriteDictionaries();
   void Fill(double x, double y, double z, double px, double py, double pz, Short_t pid,
	     UShort_t lvid, UShort_t shapeid, double safety, double snext, double step, UChar_t surfid, 
	     UChar_t process, UChar_t begend, unsigned int trid, unsigned int trpid, double cputime,
	     double cpustep);
   Bool_t IsNewEvent() {if(fNewEvent) {fNewEvent=kFALSE; return kTRUE;} 
      else return kFALSE;}
   void SetNewEvent() {fNewEvent=kTRUE;}
   void SetPrimaries(int prim) {fNprimaries=prim;}
   int  GetPrimaries() const {return fNprimaries;}
 private:
   VTfileio(const VTfileio &); // Not implemented
   VTfileio& operator=(const VTfileio &); // Not implemented

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
   UShort_t fLVid;       // logical volume id
   UShort_t fShapeid;    // shape id 
   double fSafety;       // safety
   double fSnext;        // snext
   double fStep;         // step
   UChar_t fSurfid;      // surface id
   UChar_t fProcess;     // Process
   UChar_t fBegEnd;      // Beginning or end of track
   unsigned int  fTrid;        // Track ID
   unsigned int  fTrPid;       // Track Parend ID
   double fCPUtime;    // CPU time used since start of track
   double fCPUstep;    // CPU time used for current step
 
   //
   Bool_t  fNewEvent;    // if new event
   int   fNprimaries;  // Number of primaries
   THashList *fVolumeDictionary;  // dictionary of volumes
   THashList *fShapeDictionary;   // dictionary of shapes
   THashList *fProcessDictionary; // dictionary of processes
};

#endif
