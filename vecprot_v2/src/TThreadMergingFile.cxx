// @(#)root/net:$Id$
// Author: Philippe Canal October 2011.

/*************************************************************************
 * Copyright (C) 1995-2011, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

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

#include "TThreadMergingFile.h"
#include "TSocket.h"
#include "TArrayC.h"

using namespace Geant;

////////////////////////////////////////////////////////////////////////////////
/// Constructor.
/// We do no yet open any connection to the server.  This will be done at the
/// time the first upload will be requested.

TThreadMergingFile::TThreadMergingFile(const char *filename, dcqueue<TBufferFile*>* queue, Option_t *option /* = "" */,
                                           const char *ftitle /* = "" */, Int_t compress /* = 1 */) :
  TMemFile(filename,option,ftitle,compress),fQueue(queue),fClassSent(0)
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TThreadMergingFile::~TThreadMergingFile()
{
   // We need to call Close, right here so that it is executed _before_
   // the data member of TThreadMergingFile are destructed.
   Close();
   delete fClassSent;
}

////////////////////////////////////////////////////////////////////////////////

void TThreadMergingFile::Close(Option_t *option)
{
   TMemFile::Close(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Upload the current file data to the merging server.
/// Reset the file and return true in case of success.

Bool_t TThreadMergingFile::UploadAndReset()
{
  
  TBufferFile* fBuffer = new TBufferFile(TBuffer::kWrite);
  
  fBuffer->WriteTString(GetName());
  fBuffer->WriteLong64(GetEND());
  CopyTo(*fBuffer);

  //  Printf("UploadAndReset: bufferLength %x GetBuffer %p GetEnd %x", fBuffer->Length(), fBuffer->Buffer(), GetEND());
  
  // adding to common queue to be merged afterwards

  fQueue->push(fBuffer);

   // Record the StreamerInfo we sent over.
   Int_t isize = fClassIndex->GetSize();
   if (!fClassSent) {
      fClassSent = new TArrayC(isize);
   } else {
      if (isize > fClassSent->GetSize()) {
         fClassSent->Set(isize);
      }
   }
   for(Int_t c = 0; c < isize; ++c) {
      if (fClassIndex->fArray[c]) {
         fClassSent->fArray[c] = 1;
      }
   }
   ResetAfterMerge(0);
   
   
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Write memory objects to this file and upload them to the parallel merge server.
/// Then reset all the resetable object (those with a ResetAfterMerge routine,
/// like TTree).
///
/// Loop on all objects in memory (including subdirectories).
/// A new key is created in the KEYS linked list for each object.
/// The list of keys is then saved on the file (via WriteKeys)
/// as a single data record.
/// For values of opt see TObject::Write().
/// The directory header info is rewritten on the directory header record.
/// The linked list of FREE segments is written.
/// The file header is written (bytes 1->fBEGIN).

Int_t TThreadMergingFile::Write(const char *, Int_t opt, Int_t bufsiz)
{
   Int_t nbytes = TMemFile::Write(0,opt,bufsiz);

   if (nbytes) {
      UploadAndReset();
   }
   return nbytes;
}

////////////////////////////////////////////////////////////////////////////////
/// One can not save a const TDirectory object.

Int_t TThreadMergingFile::Write(const char *n, Int_t opt, Int_t bufsize) const
{
   Error("Write const","A const TFile object should not be saved. We try to proceed anyway.");
   return const_cast<TThreadMergingFile*>(this)->Write(n, opt, bufsize);
}

////////////////////////////////////////////////////////////////////////////////
/// Write the list of TStreamerInfo as a single object in this file
/// The class Streamer description for all classes written to this file
/// is saved. See class TStreamerInfo.

void TThreadMergingFile::WriteStreamerInfo()
{
   if (!fWritable) return;
   if (!fClassIndex) return;
   //no need to update the index if no new classes added to the file
   if (fClassIndex->fArray[0] == 0) return;

   // clear fClassIndex for anything we already sent.
   if (fClassSent) {
      Int_t isize = fClassIndex->GetSize();
      Int_t ssize = fClassSent->GetSize();
      for(Int_t c = 0; c < isize && c < ssize; ++c) {
         if (fClassSent->fArray[c]) {
            fClassIndex->fArray[c] = 0;
         }
      }
   }

   TMemFile::WriteStreamerInfo();
}
