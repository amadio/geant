#ifndef TrialStats_Def
#define TrialStats_Def

#include "Geant/FormattedReporter.h"

//  Keep simple statistics for good, bad and total steps.
//  Also steps in which a lane was inactive ('keep' steps.)
//
//  J. Apostolakis  2019.05.14

template <typename Real_v>
class TrialStats {
  public:
   TrialStats() = default;
   ~TrialStats() = default;
   TrialStats( const TrialStats &) = default;

   void Update(vecCore::Mask_v<Real_v> active,
               vecCore::Mask_v<Real_v> goodStep );   // Add to stats 

   // Access methods
   vecCore::Index_v<Real_v> GetNoSteps()     const { return fNumTotalSteps; }
   vecCore::Index_v<Real_v> GetNoGoodSteps() const { return fNumGoodSteps;  }
   vecCore::Index_v<Real_v> GetNoBadSteps()  const { return fNumBadSteps;  }
   vecCore::Index_v<Real_v> GetNoKeepSteps() const { return fNumKeepSteps; } 

   unsigned int GetNoSteps(     unsigned int i ) const { return vecCore::Get( fNumTotalSteps, i ); }
   unsigned int GetNoGoodSteps( unsigned int i ) const { return vecCore::Get( fNumGoodSteps, i );  }
   unsigned int GetNoBadSteps(  unsigned int i ) const { return vecCore::Get( fNumBadSteps, i );  }
   unsigned int GetNoKeepSteps( unsigned int i ) const { return vecCore::Get( fNumKeepSteps, i ); } 
   
   //  Reset the count of selected lanes - after updating totals for lane
   void fSumAndReset(vecCore::Mask_v<Real_v> resetLane );
   void fSumAndResetLane(unsigned int i);
   // void fSumAndReset();
   
   //  Reset the count of selected lanes - no update of totals
   void ResetLane(unsigned int i);
   void Reset(vecCore::Mask_v<Real_v> resetLane );
   void Reset();

   //  Reset both counts and running totals
   void FullReset();

   // Access sums
   vecCore::Index_v<Real_v> GetSumSteps()     const { return fSumTotalSteps; }
   vecCore::Index_v<Real_v> GetSumGoodSteps() const { return fSumGoodSteps;  }
   vecCore::Index_v<Real_v> GetSumBadSteps()  const { return fSumBadSteps;  }
   vecCore::Index_v<Real_v> GetSumKeepSteps() const { return fSumKeepSteps; } 
   
   void PrintSums() const;
   //  Prints only the Sum statistics
   void PrintStats( bool full = false ) const;
   //  Prints the current statistics - and optionally the Sums
   // void PrintAllStats() const;
   
   void PrintLaneStats( unsigned int lane ) const;
   
  public:
   vecCore::Index_v<Real_v> fNumTotalSteps= vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fNumGoodSteps = vecCore::Index_v<Real_v>(0);   
   vecCore::Index_v<Real_v> fNumBadSteps  = vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fNumKeepSteps = vecCore::Index_v<Real_v>(0);

   // Old stats: numSmallSteps(0), numInitialSmallSteps(0);
   //
   // Extension idea: count 'initial' bad steps - ones at start of integration
   //  2019.05.15
   // vecCore::Index_v<Real_v> fNumInitialBadSteps = vecCore::Index_v<Real_v>(0);   
   
   // Additional capability:  sums
   vecCore::Index_v<Real_v> fSumTotalSteps= vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fSumBadSteps  = vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fSumGoodSteps = vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fSumKeepSteps = vecCore::Index_v<Real_v>(0);
};

template <typename Real_v>
void TrialStats<Real_v>::Update(vecCore::Mask_v<Real_v> active,
                                vecCore::Mask_v<Real_v> goodStep
      )
{
   vecCore::MaskedAssign( fNumTotalSteps, active,              fNumTotalSteps + 1 );
   vecCore::MaskedAssign( fNumGoodSteps,  goodStep,            fNumGoodSteps  + 1 );
   vecCore::MaskedAssign( fNumBadSteps,   active && !goodStep, fNumBadSteps + 1 );
   vecCore::MaskedAssign( fNumKeepSteps, !active,              fNumKeepSteps + 1 );
}

template <typename Real_v>
void TrialStats<Real_v>::ResetLane(unsigned int i)
{
   assert( i <= vecCore::VectorSize<Real_v>() );
   vecCore::Set( fNumTotalSteps, i, 0 );
   vecCore::Set( fNumBadSteps,   i, 0 );
   vecCore::Set( fNumGoodSteps,  i, 0 );
   vecCore::Set( fNumKeepSteps,  i, 0 );   
}

template <typename Real_v>
void TrialStats<Real_v>::Reset(vecCore::Mask_v<Real_v> changeMask )
{
   vecCore::MaskedAssign( fNumTotalSteps, changeMask,  0 );
   vecCore::MaskedAssign( fNumBadSteps,   changeMask,  0 );
   vecCore::MaskedAssign( fNumGoodSteps,  changeMask,  0 );
   vecCore::MaskedAssign( fNumKeepSteps,  changeMask,  0 );
}

template <typename Real_v>
void TrialStats<Real_v>::fSumAndReset(vecCore::Mask_v<Real_v> changeMask )
{
   vecCore::MaskedAssign( fSumTotalSteps,  changeMask, fSumTotalSteps + fNumTotalSteps );
   vecCore::MaskedAssign( fSumBadSteps,    changeMask, fSumBadSteps   + fNumBadSteps );
   vecCore::MaskedAssign( fSumGoodSteps,   changeMask, fSumGoodSteps  + fNumGoodSteps );
   vecCore::MaskedAssign( fSumKeepSteps,   changeMask, fSumKeepSteps  + fNumKeepSteps );
   
   Reset( changeMask );
}

template <typename Real_v>
void TrialStats<Real_v>::fSumAndResetLane(unsigned int i)
{
   assert( i <= vecCore::VectorSize<Real_v>() );
   vecCore::Mask_v<Real_v> changeMask= vecCore::Mask_v<Real_v>(false);

   Set( changeMask, i, true );
   
   fSumAndReset( changeMask );
}

template <typename Real_v>
void TrialStats<Real_v>::Reset()
{
   fNumTotalSteps= vecCore::Index_v<Real_v>(0);
   fNumBadSteps  = vecCore::Index_v<Real_v>(0);
      
   fNumGoodSteps = vecCore::Index_v<Real_v>(0);
   fNumKeepSteps = vecCore::Index_v<Real_v>(0);
}


template <typename Real_v>
void TrialStats<Real_v>::FullReset()
{
   Reset();
      
   fSumTotalSteps= vecCore::Index_v<Real_v>(0);
   fSumBadSteps  = vecCore::Index_v<Real_v>(0);
      
   fSumGoodSteps = vecCore::Index_v<Real_v>(0);
   fSumKeepSteps = vecCore::Index_v<Real_v>(0);
}


template <typename Real_v>
void TrialStats<Real_v>::PrintSums() const
{
   FormattedReporter::ReportRowOfInts<Real_v>(  "fSum Total Steps",  fSumTotalSteps );
   FormattedReporter::ReportRowOfInts<Real_v>(  "fSum Good  Steps",  fSumGoodSteps  );
   FormattedReporter::ReportRowOfInts<Real_v>(  "fSum Bad   Steps",  fSumBadSteps  );
   FormattedReporter::ReportRowOfInts<Real_v>(  "fSum Keep  Steps",  fSumKeepSteps );
}

template <typename Real_v>
void TrialStats<Real_v>::PrintStats( bool full ) const
{
   FormattedReporter::ReportRowOfInts<Real_v>(  "no Total Steps",  fNumTotalSteps );
   FormattedReporter::ReportRowOfInts<Real_v>(  "no Good  Steps",  fNumGoodSteps  );
   FormattedReporter::ReportRowOfInts<Real_v>(  "no Bad   Steps",  fNumBadSteps  );
   FormattedReporter::ReportRowOfInts<Real_v>(  "no Keep  Steps",  fNumKeepSteps );

   if( full )
   {
      PrintSums();
   }
}

template <typename Real_v>
void TrialStats<Real_v>::PrintLaneStats( unsigned int lane ) const
{
   std::cout << "Information for lane " << lane << "  Num Steps:  "
             << "  Total = " << GetNoSteps( lane ) << " " 
             << "  Good  = " << GetNoGoodSteps(lane)    << " " 
             << "  Bad   = " << GetNoBadSteps(lane)     << " " 
             << "  Keep  = " << GetNoKeepSteps(lane)    << std::endl;
}



#endif
