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
   TrialStats();
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
   void SumAndReset(vecCore::Mask_v<Real_v> resetLane );
   void SumAndResetLane(unsigned int i);
   void SumAndResetLaneOld(unsigned int i);
   // void fSumAndReset();
   
   //  Reset the count of selected lanes - no update of totals
   void ResetLane(unsigned int i);
   void ResetLanes(vecCore::Mask_v<Real_v> resetLane );
   void Reset();

   //  Reset both counts and running totals
   void FullReset();

   // Access sums
   vecCore::Index_v<Real_v> GetSumSteps()     const { return fSumTotalSteps; }
   vecCore::Index_v<Real_v> GetSumGoodSteps() const { return fSumGoodSteps;  }
   vecCore::Index_v<Real_v> GetSumBadSteps()  const { return fSumBadSteps;  }
   vecCore::Index_v<Real_v> GetSumKeepSteps() const { return fSumKeepSteps; } 
   
   void PrintStats( bool full = false ) const;
   //  Prints the current statistics - and optionally the Sums

   void PrintSums() const;
   //  Prints only the Sum statistics - per Lane
   void PrintSummary() const;
   //  Print Grand Totals
   
   // void PrintAllStats() const;
   void PrintState()    const;  // All .. buffers, sums, summary
   void PrintBuffers()  const;   // Just the current buffers
   
   void PrintLaneStats( unsigned int lane ) const;

   void PrintNumActive() const;
   void ResetNumActive();

   // Check that the expected 'conservation law' is obeyed - else report! 
   void CheckSums() const; 
   
 private:
   // STATE
   // 1 - current counters
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
   // 2 - Sum counters
   vecCore::Index_v<Real_v> fSumTotalSteps= vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fSumBadSteps  = vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fSumGoodSteps = vecCore::Index_v<Real_v>(0);
   vecCore::Index_v<Real_v> fSumKeepSteps = vecCore::Index_v<Real_v>(0);

   static constexpr int VecSize = vecCore::VectorSize<Real_v>();
   unsigned long int fNumActive[ VecSize + 1 ];
   unsigned long int fNumUpdateCalls = 0;
   // static vecCore::Index_v<Real_v> fLaneNumber;
   // static bool fInitialisedLaneNums;
};

template <typename Real_v>
TrialStats<Real_v>::TrialStats()
{
   ResetNumActive();
   // InitialiseLaneNumbers();
}

template <typename Real_v>
void TrialStats<Real_v>::Update(vecCore::Mask_v<Real_v> active,
                                vecCore::Mask_v<Real_v> goodStep
      )
{
   vecCore::MaskedAssign( fNumTotalSteps, active,              fNumTotalSteps + 1 );
   vecCore::MaskedAssign( fNumGoodSteps,  goodStep,            fNumGoodSteps  + 1 );
   vecCore::MaskedAssign( fNumBadSteps,   active && !goodStep, fNumBadSteps + 1 );
   vecCore::MaskedAssign( fNumKeepSteps, !active,              fNumKeepSteps + 1 );

   int numActive = countMaskTrue<Real_v>( active );
   fNumActive[ numActive ]++;
   fNumUpdateCalls++;

   // if ....  suppress except in development or check mode
   CheckSums();
}

template <typename Real_v>
void TrialStats<Real_v>::CheckSums() const
{
   using Index_v= vecCore::Index_v<Real_v>;
   
   // Ensure that each lane has the same 'Grand total' 
   //     of active ('Total'=good+bad) + inactive ('keep') steps
   //     summed across the current buffer and the 'banked' sums.
   //
   auto  bigSum = fNumTotalSteps + fNumKeepSteps + fSumTotalSteps + fSumKeepSteps;
   auto  checkLaneSum = ( bigSum == Index_v(fNumUpdateCalls) );
   if( ! vecCore::MaskFull (checkLaneSum ) ) {
      std::cerr << "TS::Update>  Error in lane Sums: Differences in " << !checkLaneSum << std::endl;
      std::cout << "TS::Update>  Error in lane Sums: Differences: " << std::endl;      
      FormattedReporter::ReportRowOfBools<Real_v>(  "TS::Update>  Check: ", checkLaneSum );
      FormattedReporter::ReportRowOfInts<Real_v>(   "TS::Upd> Diff S-exp=", bigSum - fNumUpdateCalls );
   }
   
}

template <typename Real_v>
void TrialStats<Real_v>::ResetLane(unsigned int i)
{
   vecCore::Set( fNumTotalSteps, i, 0 );
   vecCore::Set( fNumBadSteps,   i, 0 );
   vecCore::Set( fNumGoodSteps,  i, 0 );
   vecCore::Set( fNumKeepSteps,  i, 0 );
}



template <typename Real_v>
void TrialStats<Real_v>::SumAndResetLaneOld(unsigned int i)
{
   // fSumTotalSteps[i] +=  vecCore::Get( fNumTotalSteps, i );
   vecCore::Set( fSumTotalSteps,
                 i,
                 vecCore::Get( fSumTotalSteps, i ) +
                 vecCore::Get( fNumTotalSteps, i ) );
   // fSumGoodSteps[i]  +=  vecCore::Get( fNumGoodSteps , i);
   vecCore::Set( fSumGoodSteps,
                 i,
                 vecCore::Get( fSumGoodSteps, i) +
                 +  vecCore::Get( fNumGoodSteps , i) );
   // fSumBadSteps[i]   +=  fNumBadSteps[i];
   vecCore::Set( fSumBadSteps,
                 i,
                 vecCore::Get( fSumBadSteps, i ) +
                 vecCore::Get( fNumBadSteps, i )   );   
   // fSumKeepSteps[i]   +=  fNumKeepSteps[i];
   vecCore::Set( fSumKeepSteps,
                 i,
                 vecCore::Get( fSumKeepSteps, i ) +
                 vecCore::Get( fNumKeepSteps, i )  );
   assert( i <= vecCore::VectorSize<Real_v>() );

   vecCore::Set( fNumTotalSteps, i, 0 );
   vecCore::Set( fNumBadSteps,   i, 0 );
   vecCore::Set( fNumGoodSteps,  i, 0 );
   vecCore::Set( fNumKeepSteps,  i, 0 );
}

template <typename Real_v>
void TrialStats<Real_v>::ResetLanes(vecCore::Mask_v<Real_v> changeMask )
{
   const vecCore::Index_v<Real_v> Zeroes(0);
   
   vecCore::MaskedAssign( fNumTotalSteps, changeMask,  Zeroes ); // Was , 0 );
   vecCore::MaskedAssign( fNumBadSteps,   changeMask,  Zeroes );
   vecCore::MaskedAssign( fNumGoodSteps,  changeMask,  Zeroes );
   vecCore::MaskedAssign( fNumKeepSteps,  changeMask,  Zeroes );
}

template <typename Real_v>
void TrialStats<Real_v>::SumAndReset(vecCore::Mask_v<Real_v> changeMask )
{
   vecCore::MaskedAssign( fSumTotalSteps,  changeMask, fSumTotalSteps + fNumTotalSteps );
   vecCore::MaskedAssign( fSumBadSteps,    changeMask, fSumBadSteps   + fNumBadSteps );
   vecCore::MaskedAssign( fSumGoodSteps,   changeMask, fSumGoodSteps  + fNumGoodSteps );
   vecCore::MaskedAssign( fSumKeepSteps,   changeMask, fSumKeepSteps  + fNumKeepSteps );
   
   ResetLanes( changeMask );
}

template <typename Real_v>
void TrialStats<Real_v>::SumAndResetLane(unsigned int i)
{
   assert( i <= vecCore::VectorSize<Real_v>() );

   // std::cout << "**SumAndResetLane called for lane " << i << std::endl;
   // std::cout << "  - state before: " << std::endl;
   // PrintState();
   
   vecCore::Mask_v<Real_v> changeMask= vecCore::Mask_v<Real_v>(false);

   vecCore::Set( changeMask, i, true );
   
   SumAndReset( changeMask );

   // - For debugging 
   // std::cout << "  - state after: " << std::endl;
   // PrintState();

   // Check the invariant(s)
   CheckSums();
}

template <typename Real_v>
void TrialStats<Real_v>::ResetNumActive()
{
   // std::cout << "Resetting active counts." << std::endl;
   for( int i= 0; i <= VecSize; i++ )
   {
      fNumActive[i]= 0;
   }
}

template <typename Real_v>
void TrialStats<Real_v>::Reset()
{
   fNumTotalSteps= vecCore::Index_v<Real_v>(0);
   fNumBadSteps  = vecCore::Index_v<Real_v>(0);
      
   fNumGoodSteps = vecCore::Index_v<Real_v>(0);
   fNumKeepSteps = vecCore::Index_v<Real_v>(0);
   ResetNumActive();
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
void TrialStats<Real_v>::PrintState() const
{
   std::cout << "-Buffers" << std::endl;
   PrintBuffers();
   std::cout << "-Sums" << std::endl;   
   PrintSums();
   std::cout << "-Summary" << std::endl;   
   PrintSummary();
   std::cout << "-End of State ------------------------------" << std::endl;
}
   
template <typename Real_v>
void TrialStats<Real_v>::PrintBuffers() const
{
   constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();      
   if( VecSize > 1 )
   {
     FormattedReporter::ReportRowOfInts<Real_v>(  "No  Total Steps      ",  fNumTotalSteps );
     FormattedReporter::ReportRowOfInts<Real_v>(  "No  Good  Steps      ",  fNumGoodSteps  );
     FormattedReporter::ReportRowOfInts<Real_v>(  "No  Bad   Steps      ",  fNumBadSteps  );
     FormattedReporter::ReportRowOfInts<Real_v>(  "No  Keep  Steps      ",  fNumKeepSteps );
     FormattedReporter::ReportRowOfInts<Real_v>(  "Add2: Work  (t+k) ",  fNumTotalSteps + fNumKeepSteps );
   } else {
      std::cout << "No(buffers):  Total Steps = " << std::setw(5) << fNumTotalSteps
                << " Good  Steps= " << std::setw(5) <<  fNumGoodSteps
                << " Bad   Steps= " << std::setw(5) <<  fNumBadSteps
                << " Keep  Steps= " << std::setw(5) <<  fNumKeepSteps << std::endl;
   }
}

template <typename Real_v>
void TrialStats<Real_v>::PrintSums() const
{
   constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();      
   if( VecSize > 1 )
   {
     FormattedReporter::ReportRowOfInts<Real_v>(  "Sum Total Steps      ",  fSumTotalSteps );
     FormattedReporter::ReportRowOfInts<Real_v>(  "Sum Good  Steps      ",  fSumGoodSteps  );
     FormattedReporter::ReportRowOfInts<Real_v>(  "Sum Bad   Steps      ",  fSumBadSteps  );
     FormattedReporter::ReportRowOfInts<Real_v>(  "Sum Keep  Steps      ",  fSumKeepSteps );
     FormattedReporter::ReportRowOfInts<Real_v>(  "Add2:  Work  (t+k) ",  fSumTotalSteps + fSumKeepSteps );
     FormattedReporter::ReportRowOfInts<Real_v>(  "Add4:  W/N+S (t+k) ",  fSumTotalSteps + fSumKeepSteps
                                                  + fNumTotalSteps + fNumKeepSteps );
   } else {
      std::cout << "Sums:  Total Steps = " << std::setw(5) << fSumTotalSteps
                << " Good  Steps= " << std::setw(5) <<  fSumGoodSteps
                << " Bad   Steps= " << std::setw(5) <<  fSumBadSteps
                << " Keep  Steps= " << std::setw(5) <<  fSumKeepSteps << std::endl;
   }
   PrintNumActive();
}

template <typename Real_v>
void TrialStats<Real_v>::PrintSummary() const
{
   constexpr unsigned int VecSize = vecCore::VectorSize<Real_v>();   
   unsigned long  GrandTotalSteps(0UL), GrandTotalGood(0UL), GrandTotalBad(0UL), GrandTotalKeep(0UL);
   for (unsigned int i = 0; i < VecSize; ++i)
   {
      GrandTotalSteps += vecCore::Get( fSumTotalSteps, i );
      GrandTotalGood  += vecCore::Get( fSumGoodSteps , i );
      GrandTotalBad   += vecCore::Get( fSumBadSteps,   i );
      GrandTotalKeep  += vecCore::Get( fSumKeepSteps,  i );
   }
   std::cout << "GrandSums:  Total Steps = " << std::setw(5) << GrandTotalSteps
             << " Good  Steps= " << std::setw(5) << GrandTotalGood
             << " Bad   Steps= " << std::setw(5) << GrandTotalBad
             << " Keep  Steps= " << std::setw(5) << GrandTotalKeep
             << " # Update Calls= " << std::setw(5) << fNumUpdateCalls
             << std::endl;
   PrintNumActive();
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

template <typename Real_v>
void TrialStats<Real_v>::PrintNumActive() const
{
   std::cout << "Num of steps with active:";
   for( int i= 0; i <= VecSize; i++ )
   {
      std::cout // << " [" << i << "] = "
                << std::setw(3) << fNumActive[i] << " ";
   }
   std::cout << std::endl;
}

/****
template <typename Real_v>
bool TrialStats<Real_v>::fInitialisedLaneNums = false;

template <typename Real_v>
TrialStats<Real_v>::InitialisedLaneNumbers()
{
   if( ! fInitialisedLaneNums )
   {
      for( unsigned long i= 0; i < VecSize; i++ )
      {
         vecCore::Set( fLaneNumber, i, (vecCore::Index_v<double>) i );
      }
      fInitialisedLaneNums= true;
   }
}
*****/

#endif
