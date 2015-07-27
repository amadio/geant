#include "VTfileio.h"

#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"

#include "G4String.hh"

VTfileio* VTfileio::fgInstance = 0;

VTfileio::VTfileio(): fOutFile(0), 
		      fCurTree(0),
		      fX(0),
		      fY(0),
		      fZ(0),
		      fPx(0),
		      fPy(0),
		      fPz(0),
		      fPID(0),
		      fLVid(0),
		      fShapeid(0),
		      fSafety(0),
		      fSnext(0),
		      fStep(0),
		      fSurfid(0),
		      fProcess(0),
		      fBegEnd(0),
		      fTrid(0),
		      fTrPid(0),
		      fCPUtime(0),
		      fCPUstep(0),
		      fNewEvent(0),
		      fNprimaries(0),
		      fVolumeDictionary(new THashList),
		      fShapeDictionary(new THashList),
		      fProcessDictionary(new THashList)
{
   TNamed *nam = new TNamed("Transportation","Transportation");
   nam->SetUniqueID(0);
   fProcessDictionary->Add(nam);
}

VTfileio::~VTfileio() {
   if(fCurTree) fCurTree->Write(); 
   WriteDictionaries();
   fOutFile->Close();
   //   delete fCurTree; 
   fVolumeDictionary->Delete();
   fProcessDictionary->Delete();
   fShapeDictionary->Delete();
   delete fVolumeDictionary;
   delete fProcessDictionary;
   delete fShapeDictionary;
}

void VTfileio::AddVolume(const char* volname) {
   TNamed *obj = new TNamed(volname,volname);
   obj->SetUniqueID(fVolumeDictionary->GetEntries());
   fVolumeDictionary->Add(obj);
}

void VTfileio::AddShape(const char* shapename) {
   TNamed *obj = (TNamed *) fShapeDictionary->FindObject(shapename);
   if(obj) return;
   obj = new TNamed(shapename,shapename);
   obj->SetUniqueID(fShapeDictionary->GetEntries());
   fShapeDictionary->Add(obj);
}

void VTfileio::WriteDictionaries() {
   fVolumeDictionary->Write("LogicalVolumes",TObject::kSingleKey);
   fProcessDictionary->Write("ProcessDictionary",TObject::kSingleKey);
   fShapeDictionary->Write("ShapeDictionary",TObject::kSingleKey);
}

void VTfileio::NewTree(const char* name) {
   if(fCurTree) fCurTree->Write(); delete fCurTree; 
   fCurTree = new TTree(name,"Event Tree");
   fCurTree->Branch("x",&fX,"x/D");
   fCurTree->Branch("y",&fY,"y/D");
   fCurTree->Branch("z",&fZ,"z/D");
   fCurTree->Branch("px",&fPx,"px/D");
   fCurTree->Branch("py",&fPy,"py/D");
   fCurTree->Branch("pz",&fPz,"px/D");
   fCurTree->Branch("pid",&fPID,"pid/S");
   fCurTree->Branch("lvid",&fLVid,"lvid/s");
   fCurTree->Branch("shapeid",&fShapeid,"shapeid/s");
   fCurTree->Branch("safety",&fSafety,"safety/D");
   fCurTree->Branch("snext",&fSnext,"snext/D");
   fCurTree->Branch("step",&fStep,"step/D");
   fCurTree->Branch("surfid",&fSurfid,"surfid/b");
   fCurTree->Branch("process",&fProcess,"process/b");
   fCurTree->Branch("begend",&fBegEnd,"begend/b");
   fCurTree->Branch("trid",&fTrid,"trid/i");
   fCurTree->Branch("trpid",&fTrPid,"trpid/i");
   fCurTree->Branch("cputime",&fCPUtime,"cputime/D");
   fCurTree->Branch("cpustep",&fCPUstep,"cpustep/D");
}

void VTfileio::Fill(double x, double y, double z, double px, double py, double pz, Short_t pid,
		    UShort_t lvid, UShort_t shapeid, double safety, double snext, double step, UChar_t surfid, 
		    UChar_t process, UChar_t begend, UInt_t trid, UInt_t trpid, double cputime,
		    double cpustep) {
   fX = x;
   fY = y;
   fZ = z;
   fPx = px;
   fPy = py;
   fPz = pz;
   fPID = pid;
   fLVid = lvid;
   fShapeid = shapeid;
   fSafety = safety;
   fSnext = snext;
   fStep = step;
   fSurfid = surfid;
   fProcess = process;
   fBegEnd = begend;
   fTrid = trid;
   fTrPid = trpid;
   fCPUtime = cputime;
   fCPUstep = cpustep;
   fCurTree->Fill();
}
