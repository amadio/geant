#include "VTfileio.h"

#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"

#include "G4String.hh"

VTfileio* VTfileio::fgInstance = 0;

VTfileio::VTfileio(): fOutFile(0), fCurTree(0), 
		      fVolumeDictionary(new TObjArray),
		      fProcessDictionary(new TObjArray)
{fProcessDictionary->Add(new TObjString("Transportation"));}

VTfileio::~VTfileio() {
   if(fCurTree) fCurTree->Write(); 
   delete fCurTree; 
   WriteDictionaries();
   fOutFile->Close();
   fVolumeDictionary->Delete();
   fProcessDictionary->Delete();
   delete fVolumeDictionary;
   delete fProcessDictionary;
}

int VTfileio::VolumeIndex(const char* volname) const {
   for(int iv=0; iv<fVolumeDictionary->GetEntries(); ++iv)
      if(((TObjString*) fVolumeDictionary->At(iv))->GetString() == volname) return iv;
   return -1;
}

void VTfileio::WriteDictionaries() {
   fVolumeDictionary->Write("LogicalVolumes",TObject::kSingleKey);
   fProcessDictionary->Write("ProcessDictionary",TObject::kSingleKey);
}

int VTfileio::ProcessIndex(const char* procname) const {
   int nproc = fProcessDictionary->GetEntries();
   for(int ip=0; ip<nproc; ++ip)
      if(((TObjString*) fProcessDictionary->At(ip))->GetString() == procname) return ip;
   fProcessDictionary->Add(new TObjString(procname));
   return nproc;
}

void VTfileio::NewTree(const char* name) 
    {if(fCurTree) fCurTree->Write(); delete fCurTree; 
       fCurTree = new TTree(name,"Event Tree");
       fCurTree->Branch("x",&fX,"x/D");
       fCurTree->Branch("y",&fY,"y/D");
       fCurTree->Branch("z",&fZ,"z/D");
       fCurTree->Branch("px",&fPx,"px/D");
       fCurTree->Branch("py",&fPy,"py/D");
       fCurTree->Branch("pz",&fPz,"px/D");
       fCurTree->Branch("pid",&fPID,"pid/S");
       fCurTree->Branch("lvid",&fLVid,"lvid/s");
       fCurTree->Branch("safety",&fSafety,"safety/D");
       fCurTree->Branch("snext",&fSnext,"snext/D");
       fCurTree->Branch("step",&fStep,"step/D");
       fCurTree->Branch("surfid",&fSurfid,"surfid/b");
       fCurTree->Branch("process",&fProcess,"process/b");
       fCurTree->Branch("begend",&fBegEnd,"begend/b");
       fCurTree->Branch("trid",&fTrid,"trid/i");
       fCurTree->Branch("trpid",&fTrPid,"trpid/i");
    }

void VTfileio::Fill(double x, double y, double z, double px, double py, double pz, Short_t pid,
		    UShort_t lvid, double safety, double snext, double step, UChar_t surfid, 
		    UChar_t process, UChar_t begend, UInt_t trid, UInt_t trpid) {
   fX = x;
   fY = y;
   fZ = z;
   fPx = px;
   fPy = py;
   fPz = pz;
   fPID = pid;
   fLVid = lvid;
   fSafety = safety;
   fSnext = snext;
   fStep = step;
   fSurfid = surfid;
   fProcess = process;
   fBegEnd = begend;
   fTrid = trid;
   fTrPid = trpid;
   fCurTree->Fill();
}
