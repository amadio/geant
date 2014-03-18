//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4MultiLevelLocator.cc,v 1.6 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class G4MultiLevelLocator implementation
//
// 27.10.08 - Tatiana Nikitina.
// 04.10.11 - John Apostolakis, revised convergence to use Surface Normal
// ---------------------------------------------------------------------------

//#include <iomanip>
//#include "G4ios.hh"

#include "GPMultiLevelLocator.h"
#include "GPUtils.h"

#include "stdio.h"

FQUALIFIER
void GPMultiLevelLocator_Constructor( GPMultiLevelLocator *This,
				      GPNavigator *theNavigator)
{
  // : G4VIntersectionLocator(theNavigator)
  GPVIntersectionLocator_GPVIntersectionLocator(This,theNavigator);  

  // In case of too slow progress in finding Intersection Point
  // intermediates Points on the Track must be stored.
  // Initialise the array of Pointers [max_depth+1] to do this  
  
  GPThreeVector zeroV = GPThreeVector_create(0.0,0.0,0.0);

  for (G4int idepth=0; idepth<max_depth+1; idepth++ )
  {
    GPFieldTrack aFieldTrack;
    GPFieldTrack_Constructor2(&aFieldTrack,
			      zeroV,zeroV,0.,0.,0.,0.,0.,NULL);
    This->ptrInterMedFT[idepth] = &aFieldTrack;
  }
}      

FQUALIFIER
void GPVIntersectionLocator_GPVIntersectionLocator(GPMultiLevelLocator *This,
                                                   GPNavigator *theNavigator)
{
  This->fUseNormalCorrection = false; 
  This->fiEpsilonStep = -1.0 ;        // Out of range - overridden at each step
  This->fiDeltaIntersection = -1.0;   // Out of range - overridden at each step
  This->fiUseSafety = false;          // Default - overridden at each step

  This->fiChordFinder = NULL;         // Not set - overridden at each step
  This->fiNavigator = theNavigator;
  This->fpTouchable = NULL;           

  // kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  // fVerboseLevel = 0;
  // fHelpingNavigator = new G4Navigator();
  GPNavigator aNavigator;
  GPNavigator_Constructor(&aNavigator);
  This->fHelpingNavigator = &aNavigator;
}    

/*
GPMultiLevelLocator_~GPMultiLevelLocator()
{
  for ( G4int idepth=0; idepth<max_depth+1; idepth++)
  {
    delete ptrInterMedFT[idepth];
  }
}
*/

// --------------------------------------------------------------------------
// G4bool G4PropagatorInField::LocateIntersectionPoint( 
//        const G4FieldTrack&       CurveStartPointVelocity,   // A
//        const G4FieldTrack&       CurveEndPointVelocity,     // B
//        const GPThreeVector&      TrialPoint,                // E
//              G4FieldTrack&       IntersectedOrRecalculated  // Output
//              G4bool&             recalculated )             // Out
// --------------------------------------------------------------------------
//
// Function that returns the intersection of the true path with the surface
// of the current volume (either the external one or the inner one with one
// of the daughters:
//
//     A = Initial point
//     B = another point 
//
// Both A and B are assumed to be on the true path:
//
//     E is the first point of intersection of the chord AB with 
//     a volume other than A  (on the surface of A or of a daughter)
//
// Convention of Use :
//     i) If it returns "true", then IntersectionPointVelocity is set
//       to the approximate intersection point.
//    ii) If it returns "false", no intersection was found.
//          The validity of IntersectedOrRecalculated depends on 'recalculated'
//        a) if latter is false, then IntersectedOrRecalculated is invalid. 
//        b) if latter is true,  then IntersectedOrRecalculated is
//             the new endpoint, due to a re-integration.
// --------------------------------------------------------------------------
// NOTE: implementation taken from G4PropagatorInField
//
FQUALIFIER
G4bool GPMultiLevelLocator_EstimateIntersectionPoint(GPMultiLevelLocator *This, 
		GPFieldTrack&       CurveStartPointVelocity,       // A
                GPFieldTrack&       CurveEndPointVelocity,         // B
                GPThreeVector&      TrialPoint,                    // E
                GPFieldTrack&       IntersectedOrRecalculatedFT,   // Output
                G4bool&             recalculatedEndPoint,          // Out
                G4double&           previousSafety,                // In/Out
                GPThreeVector&      previousSftOrigin)             // In/Out
{
  // Find Intersection Point ( A, B, E )  of true path AB - start at E.

  G4bool found_approximate_intersection = false;
  G4bool there_is_no_intersection       = false;
  
  GPFieldTrack  CurrentA_PointVelocity = CurveStartPointVelocity; 
  GPFieldTrack  CurrentB_PointVelocity = CurveEndPointVelocity;

  GPThreeVector CurrentE_Point = TrialPoint;
  G4bool        validNormalAtE = false;
  GPThreeVector NormalAtEntry;

  //GPFieldTrack ApproxIntersecPointV(CurveEndPointVelocity);//FT-Def-Construct
  GPFieldTrack &ApproxIntersecPointV = CurveEndPointVelocity;

  G4double      NewSafety = 0.0;
  G4bool        last_AF_intersection = false;   

  // G4bool final_section= true;  // Shows whether current section is last
                                  // (i.e. B=full end)
  G4bool first_section = true;
  recalculatedEndPoint = false; 

  G4bool restoredFullEndpoint = false;

  G4int substep_no = 0;

  //  G4int oldprc;   // cout/cerr precision settings

  // Limits for substep number
  //
  const G4int max_substeps=   10000;  // Test 120  (old value 100 )
  const G4int warn_substeps=   1000;  //      100  

  // Statistics for substeps
  //
  //  static G4int max_no_seen= -1; 
  G4int max_no_seen= -1; 

  //--------------------------------------------------------------------------  
  //  Algorithm for the case if progress in founding intersection is too slow.
  //  Process is defined too slow if after N=param_substeps advances on the
  //  path, it will be only 'fraction_done' of the total length.
  //  In this case the remaining length is divided in two half and 
  //  the loop is restarted for each half.  
  //  If progress is still too slow, the division in two halfs continue
  //  until 'max_depth'.
  //--------------------------------------------------------------------------

  const G4int param_substeps=5;  // Test value for the maximum number
                                 // of substeps
  const G4double fraction_done=0.3;

  G4bool Second_half = false;    // First half or second half of divided step

  // We need to know this for the 'final_section':
  // real 'final_section' or first half 'final_section'
  // In algorithm it is considered that the 'Second_half' is true
  // and it becomes false only if we are in the first-half of level
  // depthness or if we are in the first section

  G4int depth=0; // Depth counts how many subdivisions of initial step made

  NormalAtEntry = GPVIntersectionLocator_GetSurfaceNormal(This,
							  CurrentE_Point, 
							  validNormalAtE); 

  // Intermediates Points on the Track = Subdivided Points must be stored.
  // Give the initial values to 'InterMedFt'
  // Important is 'ptrInterMedFT[0]', it saves the 'EndCurvePoint'
  //
  *(This->ptrInterMedFT[0]) = CurveEndPointVelocity;
  for (G4int idepth=1; idepth<max_depth+1; idepth++ )
  {
    *(This->ptrInterMedFT[idepth])=CurveStartPointVelocity;
  }

  // Final_section boolean store
  //
  G4bool fin_section_depth[max_depth];
  for (G4int idepth=0; idepth<max_depth; idepth++ )
  {
    fin_section_depth[idepth]=true;
  }
  // 'SubStartPoint' is needed to calculate the length of the divided step
  //
  GPFieldTrack SubStart_PointVelocity = CurveStartPointVelocity;
   
  do
  {
    G4int substep_no_p = 0;
    G4bool sub_final_section = false; // the same as final_section,
                                      // but for 'sub_section'
    SubStart_PointVelocity = CurrentA_PointVelocity; 

    do // REPEAT param
    {
      GPThreeVector Point_A = GPFieldTrack_GetPosition(&CurrentA_PointVelocity);  
      GPThreeVector Point_B = GPFieldTrack_GetPosition(&CurrentB_PointVelocity);
       
      // F = a point on true AB path close to point E 
      // (the closest if possible)
      //
      ApproxIntersecPointV = GPChordFinder_ApproxCurvePointV(
			     GPVIntersectionLocator_GetChordFinderFor(This),
			     CurrentA_PointVelocity, 
			     CurrentB_PointVelocity, 
                             CurrentE_Point,
                             GPVIntersectionLocator_GetEpsilonStepFor(This));
        // The above method is the key & most intuitive part ...
     
      GPThreeVector CurrentF_Point= GPFieldTrack_GetPosition(&ApproxIntersecPointV);

      // First check whether EF is small - then F is a good approx. point 
      // Calculate the length and direction of the chord AF
      //
      GPThreeVector  ChordEF_Vector = GPThreeVector_sub(CurrentF_Point,CurrentE_Point);

      GPThreeVector  NewMomentumDir= GPFieldTrack_GetMomentumDir(&ApproxIntersecPointV); 
      G4double       MomDir_dot_Norm= GPThreeVector_dot(NewMomentumDir, NormalAtEntry ) ;
      
      G4bool adequate_angle =
             ( MomDir_dot_Norm >= 0.0 ) // Can use ( > -epsilon) instead
          || (! validNormalAtE );       // Invalid, cannot use this criterion
      G4double EF_dist2 = GPThreeVector_mag2(ChordEF_Vector);

      if ( ( EF_dist2 <= sqrt(This->fiDeltaIntersection) && ( adequate_angle ) )
        || ( EF_dist2 <= kCarTolerance*kCarTolerance ) )
      { 
        found_approximate_intersection = true;

        // Create the "point" return value
        //
        IntersectedOrRecalculatedFT = ApproxIntersecPointV;
        GPFieldTrack_SetPosition(&IntersectedOrRecalculatedFT,CurrentE_Point);

        if ( GPVIntersectionLocator_GetAdjustementOfFoundIntersection(This) )
        {
          // Try to Get Correction of IntersectionPoint using SurfaceNormal()
          //  
          GPThreeVector IP;
          GPThreeVector MomentumDir=GPFieldTrack_GetMomentumDirection(&ApproxIntersecPointV);
          G4bool goodCorrection = GPVIntersectionLocator_AdjustmentOfFoundIntersection(This,
                                    Point_A,
                                    CurrentE_Point, CurrentF_Point, MomentumDir,
                                    last_AF_intersection, IP, NewSafety,
                                    previousSafety, previousSftOrigin );
          if ( goodCorrection )
          {
            IntersectedOrRecalculatedFT = ApproxIntersecPointV;
            GPFieldTrack_SetPosition(&IntersectedOrRecalculatedFT,IP);
          }
        }
        // Note: in order to return a point on the boundary, 
        //       we must return E. But it is F on the curve.
        //       So we must "cheat": we are using the position at point E
        //       and the velocity at point F !!!
        //
        // This must limit the length we can allow for displacement!
      }
      else  // E is NOT close enough to the curve (ie point F)
      {
        // Check whether any volumes are encountered by the chord AF
        // ---------------------------------------------------------
        // First relocate to restore any Voxel etc information
        // in the Navigator before calling ComputeStep()
        //
        GPNavigator_LocateGlobalPointWithinVolume(
                    GPVIntersectionLocator_GetNavigatorFor(This), Point_A );

        GPThreeVector PointG;   // Candidate intersection point
        G4double stepLengthAF; 
	G4bool Intersects_AF = GPVIntersectionLocator_IntersectChord( This,
						      Point_A,   CurrentF_Point,
						      NewSafety, previousSafety,
						      previousSftOrigin,
						      stepLengthAF,
						      PointG , NULL);
	last_AF_intersection = Intersects_AF;

        if( Intersects_AF )
        {
          // G is our new Candidate for the intersection point.
          // It replaces  "E" and we will repeat the test to see if
          // it is a good enough approximate point for us.
          //       B    <- F
          //       E    <- G
          //
          CurrentB_PointVelocity = ApproxIntersecPointV;
          CurrentE_Point = PointG;  

          G4bool validNormalLast; 
	  NormalAtEntry  = GPVIntersectionLocator_GetSurfaceNormal(This, 
						  PointG, validNormalLast ); 
          validNormalAtE = validNormalLast; 

          // By moving point B, must take care if current
          // AF has no intersection to try current FB!!
          //
          fin_section_depth[depth]=false;
          
        }
        else  // not Intersects_AF
        {  
          // In this case:
          // There is NO intersection of AF with a volume boundary.
          // We must continue the search in the segment FB!
          //
          GPNavigator_LocateGlobalPointWithinVolume( 
                      GPVIntersectionLocator_GetNavigatorFor(This), 
                      CurrentF_Point );

          G4double stepLengthFB;
          GPThreeVector PointH;

          // Check whether any volumes are encountered by the chord FB
          // ---------------------------------------------------------

          G4bool Intersects_FB = GPVIntersectionLocator_IntersectChord(This, 
                                                 CurrentF_Point, Point_B, 
                                                 NewSafety, previousSafety,
                                                 previousSftOrigin,
                                                 stepLengthFB,
						 PointH, NULL);

          if( Intersects_FB )
          {
            // There is an intersection of FB with a volume boundary
            // H <- First Intersection of Chord FB 

            // H is our new Candidate for the intersection point.
            // It replaces  "E" and we will repeat the test to see if
            // it is a good enough approximate point for us.

            // Note that F must be in volume volA  (the same as A)
            // (otherwise AF would meet a volume boundary!)
            //   A    <- F 
            //   E    <- H
            //
	    CurrentA_PointVelocity = ApproxIntersecPointV;
	    CurrentE_Point = PointH;

	    G4bool validNormalH;
	    NormalAtEntry  = GPVIntersectionLocator_GetSurfaceNormal(This,
						    PointH, validNormalH ); 
	    validNormalAtE = validNormalH; 

          }
          else  // not Intersects_FB
          {
            // There is NO intersection of FB with a volume boundary

            if(fin_section_depth[depth])
            {
              // If B is the original endpoint, this means that whatever
              // volume(s) intersected the original chord, none touch the
              // smaller chords we have used.
              // The value of 'IntersectedOrRecalculatedFT' returned is
              // likely not valid 

              // Check on real final_section or SubEndSection
              //
              if( ((Second_half)&&(depth==0)) || (first_section) )
              {
                there_is_no_intersection = true;   // real final_section
              }
              else
              {
                // end of subsection, not real final section 
                // exit from the and go to the depth-1 level 

		substep_no_p = param_substeps+2;  // exit from the loop

                // but 'Second_half' is still true because we need to find
                // the 'CurrentE_point' for the next loop
                //
                Second_half = true; 
                sub_final_section = true;
              }
            }
            else
            {
              if(depth==0)
              {
                // We must restore the original endpoint
                //
                CurrentA_PointVelocity = CurrentB_PointVelocity;  // Got to B
                CurrentB_PointVelocity = CurveEndPointVelocity;
                SubStart_PointVelocity = CurrentA_PointVelocity;
                restoredFullEndpoint = true;
              }
              else
              {
                // We must restore the depth endpoint
                //
                CurrentA_PointVelocity = CurrentB_PointVelocity;  // Got to B
                CurrentB_PointVelocity =  *(This->ptrInterMedFT[depth]);
                SubStart_PointVelocity = CurrentA_PointVelocity;
                restoredFullEndpoint = true;
              }
            }
          } // Endif (Intersects_FB)
        } // Endif (Intersects_AF)

        // Ensure that the new endpoints are not further apart in space
        // than on the curve due to different errors in the integration
        //
        G4double linDistSq, curveDist; 
        linDistSq = GPThreeVector_mag2(GPThreeVector_sub( 
		    GPFieldTrack_GetPosition(&CurrentB_PointVelocity), 
                    GPFieldTrack_GetPosition(&CurrentA_PointVelocity))); 

        curveDist = GPFieldTrack_GetCurveLength(&CurrentB_PointVelocity)
                  - GPFieldTrack_GetCurveLength(&CurrentA_PointVelocity);

        // Change this condition for very strict parameters of propagation 
        //
        if( curveDist*curveDist*(1+
            2*GPVIntersectionLocator_GetEpsilonStepFor(This)) < linDistSq )
        {  
          // Re-integrate to obtain a new B
          //
          GPFieldTrack newEndPointFT=
	    GPVIntersectionLocator_ReEstimateEndpoint( This,
				   CurrentA_PointVelocity,
				   CurrentB_PointVelocity,
				   linDistSq,    // to avoid recalculation
				   curveDist );
          //GPFieldTrack oldPointVelB = CurrentB_PointVelocity; 
          CurrentB_PointVelocity = newEndPointFT;
           
          if ( (fin_section_depth[depth])           // real final section
             &&( first_section || ((Second_half)&&(depth==0)) ) )
          {
            recalculatedEndPoint = true;
            IntersectedOrRecalculatedFT = newEndPointFT;
              // So that we can return it, if it is the endpoint!
          }
        }

        if( curveDist < 0.0 )
        {
	  //This->fVerboseLevel = 5; // Print out a maximum of information

          if( recalculatedEndPoint )
          {
	    // message << "Recalculation of EndPoint was called with fEpsStep= "
	    //         << GetEpsilonStepFor() << G4endl;
          }
	  //G4Exception("GPMultiLevelLocator_EstimateIntersectionPoint()",
	  //            "GeomNav0003", FatalException, message);
        }
        if(restoredFullEndpoint)
        {
          fin_section_depth[depth] = restoredFullEndpoint;
          restoredFullEndpoint = false;
        }

      } // EndIf ( E is close enough to the curve, ie point F. )
        // tests ChordAF_Vector.mag() <= maximum_lateral_displacement 

      substep_no++; 
      substep_no_p++;

    } while (  ( ! found_approximate_intersection )
            && ( ! there_is_no_intersection )     
            && ( substep_no_p <= param_substeps) );  // UNTIL found or

                                                     // failed param substep
    first_section = false;

    if( (!found_approximate_intersection) && (!there_is_no_intersection) )
    {
      G4double did_len = fabs(GPFieldTrack_GetCurveLength(&CurrentA_PointVelocity)
			    - GPFieldTrack_GetCurveLength(&SubStart_PointVelocity));

      G4double all_len = fabs(GPFieldTrack_GetCurveLength(&CurrentB_PointVelocity)
			    - GPFieldTrack_GetCurveLength(&SubStart_PointVelocity));
   
      G4double stepLengthAB;
      GPThreeVector PointGe;
      // Check if progress is too slow and if it possible to go deeper,
      // then halve the step if so
      //

      if( ( ( did_len )<fraction_done*all_len)
         && (depth<max_depth) && (!sub_final_section) )
      {
        Second_half=false;
        depth++;

        G4double Sub_len = (all_len-did_len)/(2.);
        GPFieldTrack start = CurrentA_PointVelocity;
        GPMagInt_Driver* integrDriver
	  = GPChordFinder_GetIntegrationDriver(
                          GPVIntersectionLocator_GetChordFinderFor(This));

	GPMagInt_Driver_AccurateAdvance(integrDriver,
                                        start, Sub_len, 
					GPVIntersectionLocator_GetEpsilonStepFor(This),0);
        *(This->ptrInterMedFT[depth]) = start;
        CurrentB_PointVelocity = *(This->ptrInterMedFT[depth]);
 
        // Adjust 'SubStartPoint' to calculate the 'did_length' in next loop
        //
        SubStart_PointVelocity = CurrentA_PointVelocity;

        // Find new trial intersection point needed at start of the loop
        //
        GPThreeVector Point_A = GPFieldTrack_GetPosition(&CurrentA_PointVelocity);
        GPThreeVector SubE_point = GPFieldTrack_GetPosition(&CurrentB_PointVelocity);   
     
        GPNavigator_LocateGlobalPointWithinVolume(
		    GPVIntersectionLocator_GetNavigatorFor(This),Point_A);

        G4bool Intersects_AB = GPVIntersectionLocator_IntersectChord(This,
                                              Point_A, SubE_point,
                                              NewSafety, previousSafety,
                                              previousSftOrigin, stepLengthAB,
					      PointGe,NULL);
        if( Intersects_AB )
        {
          last_AF_intersection = Intersects_AB;
          CurrentE_Point = PointGe;
          fin_section_depth[depth]=true;

          G4bool validNormalAB; 
	  NormalAtEntry  = GPVIntersectionLocator_GetSurfaceNormal(This,
	                                          PointGe, validNormalAB ); 
          validNormalAtE = validNormalAB;
        }
        else
        {
          // No intersection found for first part of curve
          // (CurrentA,InterMedPoint[depth]). Go to the second part
          //
          Second_half = true;
        }
      } // if did_len

      if( (Second_half)&&(depth!=0) )
      {
        // Second part of curve (InterMed[depth],Intermed[depth-1])                       ) 
        // On the depth-1 level normally we are on the 'second_half'

        Second_half = true;
        //  Find new trial intersection point needed at start of the loop
        //
        SubStart_PointVelocity = *(This->ptrInterMedFT[depth]);
        CurrentA_PointVelocity = *(This->ptrInterMedFT[depth]);
        CurrentB_PointVelocity = *(This->ptrInterMedFT[depth-1]);
         // Ensure that the new endpoints are not further apart in space
        // than on the curve due to different errors in the integration
        //
        G4double linDistSq, curveDist; 
        linDistSq = GPThreeVector_mag2( GPThreeVector_sub(
                      GPFieldTrack_GetPosition(&CurrentB_PointVelocity), 
                      GPFieldTrack_GetPosition(&CurrentA_PointVelocity))); 

        curveDist = GPFieldTrack_GetCurveLength(&CurrentB_PointVelocity),
		  - GPFieldTrack_GetCurveLength(&CurrentA_PointVelocity);

        if( curveDist*curveDist*(1+
            2*GPVIntersectionLocator_GetEpsilonStepFor(This) ) < linDistSq )
        {
          // Re-integrate to obtain a new B
          //
          GPFieldTrack newEndPointFT= 
	    GPVIntersectionLocator_ReEstimateEndpoint( This,
				   CurrentA_PointVelocity,
				   CurrentB_PointVelocity,
				   linDistSq,    // to avoid recalculation
				   curveDist );
          //GPFieldTrack oldPointVelB = CurrentB_PointVelocity; 
          CurrentB_PointVelocity = newEndPointFT;
          if (depth==1)
          {
            recalculatedEndPoint = true;
            IntersectedOrRecalculatedFT = newEndPointFT;
            // So that we can return it, if it is the endpoint!
          }
        }

        GPThreeVector Point_A    = GPFieldTrack_GetPosition(&CurrentA_PointVelocity);
        GPThreeVector SubE_point = GPFieldTrack_GetPosition(&CurrentB_PointVelocity);   

        GPNavigator_LocateGlobalPointWithinVolume(
                      GPVIntersectionLocator_GetNavigatorFor(This),Point_A);
        G4bool Intersects_AB = GPVIntersectionLocator_IntersectChord(This,
					      Point_A, SubE_point, NewSafety,
                                              previousSafety,
                                              previousSftOrigin, stepLengthAB,
					      PointGe, NULL);
        if( Intersects_AB )
        {
          last_AF_intersection = Intersects_AB;
	  CurrentE_Point = PointGe;

          G4bool validNormalAB = false; 
	  //G4FWP@@@ This call is problematic
	  //NormalAtEntry = GPVIntersectionLocator_GetSurfaceNormal(This,
	  //							    PointGe, 
	  //							    validNormalAB ); 
          validNormalAtE = validNormalAB;
        }       
        depth--;
        fin_section_depth[depth]=true;
      }

    }  // if(!found_aproximate_intersection)

  } while ( ( ! found_approximate_intersection )
            && ( ! there_is_no_intersection )     
            && ( substep_no <= max_substeps) ); // UNTIL found or failed

  if( substep_no > max_no_seen )
  {
    max_no_seen = substep_no; 
  }

  if(  ( substep_no >= max_substeps)
      && !there_is_no_intersection
      && !found_approximate_intersection )
  {
    //    G4Exception("GPMultiLevelLocator_EstimateIntersectionPoint()",
    //                "GeomNav0003", FatalException, message);
  }
  else if( substep_no >= warn_substeps )
  {  
    //    G4Exception("GPMultiLevelLocator_EstimateIntersectionPoint()",
    //                "GeomNav1002", JustWarning, message);
  }


  //  return  !there_is_no_intersection; //  Success or failure

  return false; 
}

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4VIntersectionLocator.cc,v 1.8 2010-07-13 15:59:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class G4VIntersectionLocator implementation
//
// 27.10.08 - John Apostolakis, Tatiana Nikitina.
// ---------------------------------------------------------------------------
 
///////////////////////////////////////////////////////////////////////////
//
// ReEstimateEndPoint.
//
FQUALIFIER
GPFieldTrack 
GPVIntersectionLocator_ReEstimateEndpoint(GPMultiLevelLocator *This, 
					  GPFieldTrack &CurrentStateA,  
					  GPFieldTrack &EstimatedEndStateB,
					  G4double      linearDistSq,
					  G4double      curveDist )
{  
  GPFieldTrack& newEndPoint = CurrentStateA;
  GPMagInt_Driver* integrDriver= GPChordFinder_GetIntegrationDriver(
                     GPVIntersectionLocator_GetChordFinderFor(This)); 

  GPFieldTrack& retEndPoint = CurrentStateA;
  G4bool goodAdvance;
  G4int  itrial=0;
  const G4int no_trials= 20;

  G4double endCurveLen= GPFieldTrack_GetCurveLength(&EstimatedEndStateB);
  do
  {
     G4double currentCurveLen= GPFieldTrack_GetCurveLength(&newEndPoint);
     G4double advanceLength= endCurveLen - currentCurveLen ; 
     if (fabs(advanceLength)<kCarTolerance)
     {
       goodAdvance=true;
     }
     else{
     goodAdvance= 
       GPMagInt_Driver_AccurateAdvance(integrDriver,newEndPoint, advanceLength,
				       GPVIntersectionLocator_GetEpsilonStepFor(This),0);
     }
  }
  while( !goodAdvance && (++itrial < no_trials) );

  if( goodAdvance )
  {
    retEndPoint= newEndPoint; 
  }
  else
  {
    retEndPoint= EstimatedEndStateB; // Could not improve without major work !!
  }

  //  All the work is done
  //  below are some diagnostics only -- before the return!
  // 
  //  static const G4String MethodName("GPVIntersectionLocator_ReEstimateEndpoint");

  // Statistics on the RMS value of the corrections

  //  static 
  G4int    noCorrections=0;
  //  static 
  G4double sumCorrectionsSq = 0;
  noCorrections++; 
  if( goodAdvance )
  {
    sumCorrectionsSq += GPThreeVector_mag2(GPThreeVector_sub(
			GPFieldTrack_GetPosition(&EstimatedEndStateB),
			GPFieldTrack_GetPosition(&newEndPoint)));
  }
  linearDistSq -= curveDist; // To use linearDistSq ... !

  return retEndPoint;
}

///////////////////////////////////////////////////////////////////////////
//
// Method for finding SurfaceNormal of Intersecting Solid 
//
FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetLocalSurfaceNormal(
                      GPMultiLevelLocator *This, 
		      GPThreeVector &CurrentE_Point, G4bool &validNormal)
{
  GPThreeVector NormalAtEntry = GPThreeVector_create( -10. , -10., -10. ); 
  GPVPhysicalVolume* located;

  validNormal = false;
  GPNavigator_SetWorldVolume(This->fHelpingNavigator, 
			     GPNavigator_GetWorldVolume(
			     GPVIntersectionLocator_GetNavigatorFor(This)));

  located = GPNavigator_LocateGlobalPointAndSetup(This->fHelpingNavigator, 
						  CurrentE_Point,NULL,0,0 );

  //  delete fpTouchable;
  This->fpTouchable = NULL;
  //equivalent to GPNavigator_CreateTouchableHistory(This->fHelpingNavigator);
  GPTouchableHistory_Constructor(This->fpTouchable);

  // To check if we can use GetGlobalExitNormal() 
  //
  GPAffineTransform aT = GPNavigationHistory_GetTopTransform(GPTouchableHistory_GetHistory(This->fpTouchable));
  GPThreeVector localPosition = GPAffineTransform_TransformPoint(&aT,CurrentE_Point);

  // Issue: in the case of coincident surfaces, this version does not recognise 
  //        which side you are located onto (can return vector with wrong sign.)
  // TO-DO: use direction (of chord) to identify volume we will be "entering"

  GPThreeVector Normal = GPThreeVector_create(0,0,0);
  if( located != 0)
  { 
    GPLogicalVolume* pLogical= GPVPhysicalVolume_GetLogicalVolume(located); 
    GPVSolid*        pSolid; 

    if((pLogical != 0) && ((pSolid= GPLogicalVolume_GetSolid(pLogical)) !=0 ))
    {
      // G4bool     goodPoint,    nearbyPoint;   
      // G4int   numGoodPoints,   numNearbyPoints;  // --> use for stats
      if ( ( GPVSolid_Inside(pSolid,localPosition)==kSurface )
           || ( GPVSolid_DistanceToOut(pSolid,localPosition) < 1000.0 * kCarTolerance )
         )
      {
        Normal = GPVSolid_SurfaceNormal(pSolid,localPosition);
        validNormal = true;

      }
    }
  }

  return Normal;
}

///////////////////////////////////////////////////////////////////////////
//
// Adjustment of Found Intersection
//
FQUALIFIER
G4bool GPVIntersectionLocator_AdjustmentOfFoundIntersection( 
                               GPMultiLevelLocator *This, 
			       GPThreeVector &CurrentA_Point,
                               GPThreeVector &CurrentE_Point, 
                               GPThreeVector &CurrentF_Point,
                               GPThreeVector &MomentumDir,
                               G4bool         IntersectAF,
			       GPThreeVector &IntersectionPoint,  // I/O
			       G4double      &NewSafety,          // I/O 
			       G4double      &fPreviousSafety,    // I/O
			       GPThreeVector &fPreviousSftOrigin )// I/O
{     
  G4double dist,lambda;
  GPThreeVector Normal, NewPoint, Point_G;
  G4bool goodAdjust=false, Intersects_FP=false, validNormal=false;

  // Get SurfaceNormal of Intersecting Solid
  //
  Normal = GPVIntersectionLocator_GetGlobalSurfaceNormal(This,
				  CurrentE_Point,validNormal);
  if(!validNormal) { return false; }

  // Intersection between Line and Plane
  //
  G4double n_d_m = GPThreeVector_dot(Normal,MomentumDir);
  if ( fabs(n_d_m)>kCarTolerance )
  {
    //    if ( fVerboseLevel>1 )
    //    {
    //        return false;
    //    }
    lambda =- GPThreeVector_dot(Normal,
	      GPThreeVector_sub(CurrentF_Point,CurrentE_Point))/n_d_m;

    // New candidate for Intersection
    //
    NewPoint = GPThreeVector_add(CurrentF_Point,
				 GPThreeVector_mult(MomentumDir,lambda));

    // Distance from CurrentF to Calculated Intersection
    //
    dist = fabs(lambda);

    if ( dist<kCarTolerance*0.001 )  { return false; }

    // Calculation of new intersection point on the path.
    //
    if ( IntersectAF )  //  First part intersects
    {
      G4double stepLengthFP; 
      GPThreeVector Point_P = CurrentA_Point;
      GPNavigator_LocateGlobalPointWithinVolume(
		  GPVIntersectionLocator_GetNavigatorFor(This),Point_P);

      Intersects_FP = GPVIntersectionLocator_IntersectChord(This, 
				     Point_P, NewPoint, NewSafety,
				     fPreviousSafety, fPreviousSftOrigin,
				     stepLengthFP, Point_G,NULL );

    }
    else   // Second part intersects
    {      
      G4double stepLengthFP; 
      GPNavigator_LocateGlobalPointWithinVolume(
		  GPVIntersectionLocator_GetNavigatorFor(This),CurrentF_Point );
      Intersects_FP =  GPVIntersectionLocator_IntersectChord(This,
				     CurrentF_Point, NewPoint, NewSafety,
				     fPreviousSafety, fPreviousSftOrigin,
				     stepLengthFP, Point_G, NULL);
    }
    if ( Intersects_FP )
    {
      goodAdjust = true;
      IntersectionPoint = Point_G;              
    }
  }

  return goodAdjust;
}

FQUALIFIER
GPThreeVector 
GPVIntersectionLocator_GetSurfaceNormal(GPMultiLevelLocator *This, 
					GPThreeVector &CurrentInt_Point,
					G4bool &validNormal) // const
{
  GPThreeVector NormalAtEntry = GPThreeVector_create( -10. , -10., -10. ); 

  GPThreeVector NormalAtEntryLast;
  //  NormalAtEntryGlobal, diffNormals;
  G4bool validNormalLast; 

  // Relies on a call to Navigator::ComputeStep in IntersectChord before this call

  NormalAtEntryLast = GPVIntersectionLocator_GetLastSurfaceNormal(This,
					     CurrentInt_Point, validNormalLast ); 

  // May return valid=false in cases, including
  //  - if the candidate volume was not found (eg exiting world), or
  //  - a replica was involved -- determined the step size.
  // (This list is not complete.) 

  if( validNormalLast ) 
  {
    //    NormalAtEntry=NormalAtEntryLast;  
    GPThreeVector_set(&NormalAtEntry, NormalAtEntryLast.x,
		                      NormalAtEntryLast.y,
		                      NormalAtEntryLast.z);  
    validNormal  = validNormalLast; 
  }
  return NormalAtEntry; 
}

FQUALIFIER
GPThreeVector GPVIntersectionLocator_GetGlobalSurfaceNormal(
				     GPMultiLevelLocator *This, 
		                     GPThreeVector &CurrentE_Point,
		                     G4bool &validNormal)
{
  GPThreeVector     localNormal=
    GPVIntersectionLocator_GetLocalSurfaceNormal(This, 
						 CurrentE_Point, validNormal );
  GPAffineTransform localToGlobal=           //  Must use the same Navigator !!
      GPNavigator_GetLocalToGlobalTransform(This->fHelpingNavigator);
  GPThreeVector     globalNormal =
    GPAffineTransform_TransformAxis(&localToGlobal, localNormal );

  return globalNormal;
}

FQUALIFIER
GPThreeVector 
GPVIntersectionLocator_GetLastSurfaceNormal( GPMultiLevelLocator *This, 
					     GPThreeVector intersectPoint,
					     G4bool &normalIsValid)
{
  GPThreeVector normalVec;
  G4bool        validNorm;
  normalVec    = GPNavigator_GetGlobalExitNormal(This->fiNavigator, 
						 intersectPoint, &validNorm ); 
  normalIsValid= validNorm;

  return normalVec;
}

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VIntersectionLocator.icc,v 1.3 2008-11-14 18:26:35 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class G4VIntersectionLocator inline methods
//
// 27.10.07 - John Apostolakis, Tatiana Nikitina
// ---------------------------------------------------------------------------

FQUALIFIER 
G4double GPVIntersectionLocator_GetDeltaIntersectionFor(GPMultiLevelLocator *This)
{
  return This->fiDeltaIntersection;
} 

FQUALIFIER 
G4double GPVIntersectionLocator_GetEpsilonStepFor(GPMultiLevelLocator *This)
{
  return This->fiEpsilonStep;
}

FQUALIFIER 
GPNavigator* GPVIntersectionLocator_GetNavigatorFor(GPMultiLevelLocator *This)
{
  return This->fiNavigator;
}

FQUALIFIER
GPChordFinder* GPVIntersectionLocator_GetChordFinderFor(GPMultiLevelLocator *This)
{
  return This->fiChordFinder;
}

FQUALIFIER 
G4int GPVIntersectionLocator_GetVerboseFor(GPMultiLevelLocator *This)
{
  //  return This->fVerboseLevel;
  return 0;
}

FQUALIFIER 
G4bool GPVIntersectionLocator_GetAdjustementOfFoundIntersection(GPMultiLevelLocator *This )
{
  return This->fUseNormalCorrection;
}

FQUALIFIER 
void GPVIntersectionLocator_AddAdjustementOfFoundIntersection(
				  GPMultiLevelLocator *This, 
				  G4bool UseCorrection )
{
  This->fUseNormalCorrection=UseCorrection;
}

FQUALIFIER 
void GPVIntersectionLocator_SetEpsilonStepFor(GPMultiLevelLocator *This,
					      G4double EpsilonStep )
{
  This->fiEpsilonStep=EpsilonStep;
}

FQUALIFIER 
void GPVIntersectionLocator_SetDeltaIntersectionFor(GPMultiLevelLocator *This, 
						    G4double deltaIntersection )
{
  This->fiDeltaIntersection=deltaIntersection;
}

FQUALIFIER 
void GPVIntersectionLocator_SetNavigatorFor( GPMultiLevelLocator *This, 
					     GPNavigator *fNavigator )
{
  This->fiNavigator=fNavigator;
}

FQUALIFIER 
void GPVIntersectionLocator_SetChordFinderFor(GPMultiLevelLocator *This, 
					      GPChordFinder *fCFinder )
{
  This->fiChordFinder=fCFinder;
}

FQUALIFIER 
void GPVIntersectionLocator_SetSafetyParametersFor(GPMultiLevelLocator *This,
						   G4bool UseSafety )
{
  This->fiUseSafety=UseSafety;
}

FQUALIFIER 
void GPVIntersectionLocator_SetVerboseFor(GPMultiLevelLocator *This, 
					  G4int fVerbose)
{
  ;
  //  This->fVerboseLevel=fVerbose;
}

FQUALIFIER G4bool
GPVIntersectionLocator_IntersectChord( GPMultiLevelLocator *This, 
				       GPThreeVector  StartPointA, 
				       GPThreeVector  EndPointB,
				       G4double      &NewSafety,
				       G4double      &PreviousSafety,
				       GPThreeVector &PreviousSftOrigin,
				       G4double      &LinearStepLength,
				       GPThreeVector &IntersectionPoint,
				       G4bool        *ptrCalledNavigator )
{
  G4bool CalledNavigator = false; 

  // Calculate the direction and length of the chord AB

  GPThreeVector  ChordAB_Vector = GPThreeVector_sub(EndPointB,StartPointA);
  G4double       ChordAB_Length = GPThreeVector_mag(ChordAB_Vector);  // Magnitude (norm)
  GPThreeVector  ChordAB_Dir =    GPThreeVector_unit(ChordAB_Vector);
  G4bool intersects;
  GPThreeVector OriginShift = GPThreeVector_sub(StartPointA,PreviousSftOrigin);
  G4double      MagSqShift  = GPThreeVector_mag2(OriginShift) ;
  G4double      currentSafety;

  if( MagSqShift >= sqrt(PreviousSafety) )
  {
    currentSafety = 0.0 ;
  }
  else
  {
    currentSafety = PreviousSafety - sqrt(MagSqShift) ;
  }

  if( This->fiUseSafety && (ChordAB_Length <= currentSafety) )
  {
    // The Step is guaranteed to be taken

    LinearStepLength = ChordAB_Length;
    intersects = false;
    NewSafety= currentSafety;
    CalledNavigator= false; 
  }
  else
  {
    // Check whether any volumes are encountered by the chord AB

    LinearStepLength = GPNavigator_ComputeStep(
                       GPVIntersectionLocator_GetNavigatorFor(This), 
		       StartPointA,
                       ChordAB_Dir, ChordAB_Length, &NewSafety );
    intersects = (LinearStepLength <= ChordAB_Length); 

       // G4Navigator contracts to return k_infinity if len==asked
       // and it did not find a surface boundary at that length

    LinearStepLength = GPfmin( LinearStepLength, ChordAB_Length);
    CalledNavigator = true; 

    // Save the last calculated safety!

    PreviousSftOrigin = StartPointA;
    PreviousSafety    = NewSafety;

    if( intersects )
    {
       // Intersection Point of chord AB and either volume A's surface 
       //                                or a daughter volume's surface ..
      IntersectionPoint = GPThreeVector_add(StartPointA,
			  GPThreeVector_mult(ChordAB_Dir,LinearStepLength));
    }
  }
  if( ptrCalledNavigator )
  { 
    *ptrCalledNavigator = CalledNavigator; 
  }

  return intersects;
}
