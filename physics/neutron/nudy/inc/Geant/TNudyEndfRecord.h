#ifndef __TNudyEndfRecord__
#define __TNudyEndfRecord__

#include "TObject.h"
#include "Riostream.h"

namespace Nudy {
class TNudyEndfRecord : public TObject {
public:
  TNudyEndfRecord();
  virtual ~TNudyEndfRecord() {}
  virtual void DumpENDF(int /*mat*/, int /*mf*/, int /*mt*/, int & /*ns*/, int) {}

private:
  ClassDef(TNudyEndfRecord, 1)
};

} // namespace
#endif
