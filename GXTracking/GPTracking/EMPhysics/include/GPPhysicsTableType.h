#ifndef GPPhysicsTableType_H
#define GPPhysicsTableType_H 1

enum GPPhysicsTableType {
  kNullTableType = -1, 
  kLambda_eBrem,        // Lambda.eBrem.e-.asc
  kLambda_eIoni,        // Lambda.eIoni.e-.asc
  kRange_eIoni,         // Range.eIoni.e-.asc
  kDEDX_eIoni,          // DEDX.eIoni.e-.asc
  kInverseRange_eIoni,  // InverseRange.eIoni.e-.asc
  kLambda_msc,          // Lambda.msc.e-.asc
  kLambda_compt,        // Lambda.compt.gamma.asc
  kLambda_conv,         // Lambda.conv.gamma.asc
  kLambdaPrim_phot,     // LambdaPrim.phot.gamma.asc
  kNumberPhysicsTable
};

static const char* GPPhysicsTableName[kNumberPhysicsTable] = {
  "Lambda.eBrem.e-.asc",
  "Lambda.eIoni.e-.asc",
  "Range.eIoni.e-.asc",
  "DEDX.eIoni.e-.asc",
  "InverseRange.eIoni.e-.asc",
  "Lambda.msc.e-.asc",
  "Lambda.compt.gamma.asc",
  "Lambda.conv.gamma.asc",
  "LambdaPrim.phot.gamma.asc"
};

#endif
