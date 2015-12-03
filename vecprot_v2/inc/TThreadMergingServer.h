#ifndef GEANT_TTHREADMERGINGSERVER
#define GEANT_TTHREADMERGINGSERVER

#include "TBufferFile.h"
#include "TBenchmark.h"
#include "TSocket.h"
#include "TH2.h"
#include "TTree.h"
#include "TMemFile.h"
#include "TRandom.h"
#include "TError.h"
#include "TFileMerger.h"

#include "TServerSocket.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TMonitor.h"

#include "TFileCacheWrite.h"
#include "TSystem.h"
#include "THashTable.h"

#include "TMath.h"
#include "TTimeStamp.h"

#include "TKey.h"

#include "dcqueue.h"

namespace Geant {
  
struct ClientInfo 
{
   TFile      *fFile;      // This object does *not* own the file, it will be owned by the owner of the ClientInfo.
   TString    fLocalName;
   UInt_t     fContactsCount;
   TTimeStamp fLastContact;
   Double_t   fTimeSincePrevContact;

   ClientInfo() : fFile(0), fLocalName(), fContactsCount(0), fTimeSincePrevContact(0) {}
   ClientInfo(const char *filename, UInt_t clientId) : fFile(0), fContactsCount(0), fTimeSincePrevContact(0) {
      fLocalName.Form("%s-%d-%d",filename,clientId,gSystem->GetPid());
   }

  void Set(TFile *file);
};

struct ThreadFileMerger : public TObject
{
   typedef std::vector<ClientInfo> ClientColl_t;

   TString       fFilename;
   TBits         fClientsContact;       //
   UInt_t        fNClientsContact;      //
   ClientColl_t  fClients;
   TTimeStamp    fLastMerge;
   TFileMerger   fMerger;

   ThreadFileMerger(const char *filename, Bool_t writeCache = kFALSE) : fFilename(filename), fNClientsContact(0), fMerger(kFALSE,kTRUE)
   {
      // Default constructor.

      fMerger.SetPrintLevel(0);
      fMerger.OutputFile(filename,"RECREATE");
      if (writeCache) new TFileCacheWrite(fMerger.GetOutputFile(),32*1024*1024);
   }

   ~ThreadFileMerger()
   {
      // Destructor.

      for(unsigned int f = 0 ; f < fClients.size(); ++f) {
         fprintf(stderr,"Client %d reported %u times\n",f,fClients[f].fContactsCount);
      }
      for( ClientColl_t::iterator iter = fClients.begin();
          iter != fClients.end();
          ++iter)
      {
         delete iter->fFile;
      }
   }

   ULong_t  Hash() const
   {
      // Return hash value for this object.
      return fFilename.Hash();
   }

   const char *GetName() const
   {
      // Return the name of the object which is the name of the output file.
      return fFilename;
   }

   Bool_t InitialMerge(TFile *input);

   Bool_t Merge();

   Bool_t NeedFinalMerge();
 
   Bool_t NeedMerge(Float_t clientThreshold);

   void RegisterClient(UInt_t clientId, TFile *file);

   ClassDef(ThreadFileMerger,0);
};

struct TThreadMergingServer : public TObject
{
  dcqueue<TBufferFile*>* fOutput; // this pointer is NOT owned by the server
  bool finish;
  
 TThreadMergingServer(dcqueue<TBufferFile*>* queue):finish(false)
    {   
      fOutput = queue;
    }
  
  ~TThreadMergingServer(){}
  void Listen();
  void Finish(){finish=true;}
};

} // namespace Geant

#endif
