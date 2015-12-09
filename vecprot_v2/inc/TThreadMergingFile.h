// @(#)root/net:$Id$
// Author: Philippe Canal October 2011.

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
#ifdef USE_ROOT

#ifndef ROOT_TThreadMergingFile
#define ROOT_TThreadMergingFile


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TThreadMergingFile                                                 //
//                                                                      //
// Specialization of TMemFile to connect to a parallel file merger.     //
// Upon a call to UploadAndReset, the content already written to the    //
// file is upload to the server and the object implementing the function//
// ResetAfterMerge (like TTree) are reset.                              //
// The parallel file merger will then collate the information coming    //
// from this client and any other client in to the file described by    //
// the filename of this object.                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMemFile
#include "TMemFile.h"
#endif
#ifndef ROOT_TBufferFile
#include "TBufferFile.h"
#endif
#ifndef ROOT_TUrl
#include "TUrl.h"
#endif

#include "dcqueue.h"

class TArrayC;

namespace Geant {

class TThreadMergingFile : public TMemFile
{
private:
  // queue of TBuffer objects
  dcqueue<TBufferFile*>* fQueue;
  TArrayC *fClassSent;      // Record which StreamerInfo we already sent.
  
 public:
  TThreadMergingFile(const char *filename, dcqueue<TBufferFile*>* qeue, Option_t *option = "", const char *ftitle = "", Int_t compress = 1);
  ~TThreadMergingFile();
  
  virtual void   Close(Option_t *option="");
  Bool_t UploadAndReset();
  virtual Int_t  Write(const char *name=0, Int_t opt=0, Int_t bufsiz=0);
  virtual Int_t  Write(const char *name=0, Int_t opt=0, Int_t bufsiz=0) const;
  virtual void   WriteStreamerInfo();
  
  ClassDef(TThreadMergingFile,1);  // TFile specialization that will semi-automatically upload its content to a merging server.
};

} // namespace Geant

#endif // ROOT_TThreadMergingFile
#endif
