//

template
<class T_Equation, unsigned int Nvar>
   class StepperClass // : public BaseClass
{
  public:

    // template <typename T>
    //   using Vector3D = vecgeom::Vector3D<T>;
    // using Double_v        = Geant::Double_v;
    using Double_v = double;
    constexpr unsigned int sNstate = 6;
   
  public:
    inline
    StepperClass( T_Equation *EqRhs,
                  unsigned int numStateVariables=0);

    /* virtual */ ~StepperClass();

    template <typename Real_v> struct ScratchSpace; // defined below

#ifdef OUTSIDE
    template <typename Real_v>
    void StepAndEstimate(const Real_v  yInput[],  // Consider __restrict__
                               Real_v  yOut[]
                               //, ScratchSpace<Real_v>* sp
       );
#endif

    // Fixed method, needed for inheritance 
    // GEANT_FORCE_INLINE
    void StepAndEstimate(const Double_v  yInput[],
                               Double_v  yOut[] )  //  override final
    {
       StepAndEstimate<Double_v>(yInput,dydx,charge,hStep,yOut,yErr);
    }
    
    template <typename Real_v>    
    // GEANT_FORCE_INLINE
    void RightHandSideInl(Real_v y[], const Real_v& charge, Real_v dydx[]) 
    { fEquation_Rhs->T_Equation::template RightHandSide<Real_v>(y, charge, dydx); }

    private:
      StepperClass( const StepperClass& ) = delete;
        // No copy c-tor      
      StepperClass& operator=(const StepperClass&) = delete;
        // No assignment operator.

    public:
      
      template <typename Real_v>
      struct ScratchSpace
      {
         // Scratch space (and in future useful state)
         // -------
         Real_v ak2[sNstore];
         Real_v yTemp2[sNstore];  // Separate temporaries per step - to aid compiler

       public:
         ScratchSpace(){}
         ~ScratchSpace(){}         
      };

      template <typename Real_v>
         ScratchSpace<Real_v>* ObtainScratchSpace()
      // Obtain object which can hold the scratch space for integration
      //   ( Should be re-used between calls - preferably long time
         { return new ScratchSpace<Real_v>(); } 
      
      // How to use it:      
      //   auto = stepper->CreatedScratchSpace<Double_v>();
      
    private:
        // 'Invariant' during integration - the pointers must not change
        // -----------
        T_Equation* fEquation_Rhs;
        bool  fDebug= false;

#ifdef OUTSIDE
};
#endif


// template <class Real_v>
// template <class T_Equation, unsigned int Nvar>   
#ifdef OUTSIDE
template <class Real_v, class T_Equation, unsigned int Nvar>
inline void
StepperClass<T_Equation,Nvar>::
  template StepAndEstimate<Real_v>(const Real_v  yInput[],       
#else
   public:                     
         template <typename Real_v>
         void StepAndEstimate(     const Real_v  yInput[],           
#endif   
                                         const Real_v  dydx[],
                                         const Real_v& charge,                        
                                         const Real_v& Step,
                                               Real_v  yOut[],
                        //, StepperClass<T_Equation,Nvar>::template ScratchSpaceCashKarp<Real_v>& sp
     )
{
    // const double a2 = 0.2 , a3 = 0.3 , a4 = 0.6 , a5 = 1.0 , a6 = 0.875;
    StepperClass<T_Equation,Nvar>::template ScratchSpace<Real_v> sp;
    
    unsigned int i;

    const double  b21 = 0.2 ,
           c1 = 37.0/378.0 , c2 = 1- c1;

    //  Saving yInput because yInput and yOut can be aliases for same array - or restrict
    for(i=0;i<Nvar;i++) 
    {
        sp.yIn[i]=yInput[i];
    }
    // RightHandSideInl(yIn, charge,  dydx) ;              // 1st Step

    for(i=0;i<Nvar;i++) 
    {
        sp.yTemp2[i] = sp.yIn[i] + b21*Step*dydx[i] ;
    }
    this->RightHandSideInl(sp.yTemp2, charge,  sp.ak2) ;      // 2nd Step

    for(i=0;i<Nvar;i++)
    {
        // Accumulate increments with correct weights
        yOut[i] = sp.yIn[i] + Step*(c1*dydx[i] + c2*sp.ak3[i] );
    }

    return ;
}

#ifndef OUTSIDE
};   // End of class declaration

//  The remaining functions / methods are defined below
#endif

template <class T_Equation, unsigned int Nvar>
inline
StepperClass<T_Equation,Nvar>::
   StepperClass(  T_Equation   *EqRhs,
              unsigned int  numStateVariables )
   : BaseClass()
{
   if( fDebug ) {
      std::cout<<"\n----Entered constructor of StepperClass "<<std::endl;
      std::cout<<"----In StepperClass constructor, Nvar is: "<<Nvar<<std::endl;
   }
   // assert( dynamic_cast<TemplateVScalarEquationOfMotion<Backend>*>(EqRhs) != 0 );
#if ENABLE_CHORD_DIST
   fLastStepLength= Double_v(0.);
#endif
   assert( (numStateVariables == 0) || (numStateVariables >= Nvar) );
  
   std::cout<<"----end of constructor of StepperClass"<<std::endl;
}

template <class T_Equation, unsigned int Nvar>
void StepperClass<T_Equation,Nvar>::
  SetEquationOfMotion(T_Equation* equation)
{
   fEquation_Rhs= equation;
   // this->GUVVectorIntegrationStepper::SetABCEquationOfMotion(nullptr); // fEquation_Rhs);
}

//  Copy - Constructor
// 
template <class T_Equation,unsigned int Nvar>
inline
StepperClass<T_Equation,Nvar>::
   StepperClass( const StepperClass& right )
   : GUVVectorIntegrationStepper( (GUVVectorEquationOfMotion*) nullptr,
                                              sOrderMethod,
                                              Nvar,
                                              right.GetNumberOfStateVariables() ),
     fEquation_Rhs( (T_Equation*) nullptr ),
     fOwnTheEquation(false)
{
   if( fDebug ) {
      std::cout << "----Entered *copy* constructor of StepperClass " << std::endl;
   }
   SetEquationOfMotion( new T_Equation( *(right.fEquation_Rhs)) );
    // fEquation_Rhs= right.GetEquationOfMotion()->Clone());
   
   // assert( dynamic_cast<GUVVectorEquationOfMotion*>(fEquation_Rhs) != 0 );   // No longer Deriving
   assert( this->GetNumberOfStateVariables() >= Nvar);

#if ENABLE_CHORD_DIST
   fLastStepLength= Double_v(0.);   
#endif

   if( fDebug )
      std::cout << " StepperClass - copy constructor: " << std::endl
                << " Nvar = " << Nvar << " Nstore= " << sNstore 
                << " Own-the-Equation = " << fOwnTheEquation << std::endl;
}

template <class T_Equation,unsigned int Nvar>
GEANT_FORCE_INLINE
StepperClass<T_Equation,Nvar>::~StepperClass()
{
   std::cout<<"----- Vector StepperClass destructor"<<std::endl;
   if( fOwnTheEquation )
      delete fEquation_Rhs; // Expect to own the equation, except if auxiliary (then sharing the equation)

  std::cout<<"----- VectorStepperClass destructor (ended)"<<std::endl;
}


#endif
