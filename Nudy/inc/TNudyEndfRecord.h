#ifndef __TNudyEndfRecord__
#define __TNudyEndfRecord__

#include "TObject.h"
#include "Riostream.h"

class TNudyEndfRecord : public TObject{
 public:
  TNudyEndfRecord();
  virtual ~TNudyEndfRecord(){}
  virtual void DumpENDF(Int_t /*mat*/,Int_t /*mf*/,Int_t /*mt*/,Int_t& /*ns*/,Int_t) {}
 private:

  ClassDef(TNudyEndfRecord,1)
};
#endif
