ClassNames = [ "TGeoBBox", 
               "TGeoCone", 
               "TGeoConeSeg",
               "TGeoPcon", 
               "TGeoPgon", 
               "TGeoTube", 
               "TGeoTubeSeg", 
               "TGeoCtub", 
               "TGeoEltu", 
               "TGeoHype", 
               "TGeoSphere", 
               "TGeoArb8" ,
               "TGeoTrap" ,
               "TGeoGtra" ,
               "TEveGeoPolyShape" ,
               "TGeoCompositeShape" ,
               "TGeoHalfSpace" ,
               "TGeoPara" ,
               "TGeoParaboloid" ,
               "TGeoScaledShape" ,
               "TGeoShapeAssembly" ,
               "TGeoTorus" ,
               "TGeoTrd1" ,
               "TGeoTrd2" ,
               "TGeoXtru" ]

#virtual Bool_t        Contains(Double_t *point) const;
#virtual Double_t      DistFromInside(Double_t *point, Double_t *dir, Int_t iact=1, Double_t step=TGeoShape::Big(), Double_t *safe=0);
#virtual Double_t      DistFromOutside(Double_t *point, Double_t *dir, Int_t iact=1, Double_t step=TGeoShape::Big(), Double_t *safe=0);
#virtual Double_t      Safety(Double_t *point, Bool_t in=kTRUE);

def EmitCommonStuff( classn ):
    None

def getIndentString( indentlevel ):
    """
    returns a string with the right indentation level
    """
    i=0
    a=""
    while i<indentlevel :
        a = a + "   "
        i = i + 1
    return a

def EmitClassDef( classn ):
    """
    expects a classname as input
    """
    print classn+"_v : public "+  classn + " {"


def EmitClassClosing( ):
    """
    closes the current class definition
    """
    print "};"


def EmitLoopN( indentlevel ):
    """
    writes the loop
    """
    print getIndentString( indentlevel ) + "for(unsigned int k=0;k < vecsize; ++k){"


def EmitContainsDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "virtual void Contains_l( Double_t const *point, Bool_t * isin , Int_t vecsize ) {"


def EmitCallToContains( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "isin[k]= " + classname + "::Contains( (Double_t *) &point[3*k] );"


def EmitSafetyDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "virtual void Safety_l( Double_t const *point, Bool_t inside, Double_t * safe , Int_t vecsize ) {"


def EmitCallToSafety( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "safe[k]= " + classname + "::Safety( (Double_t *) &point[3*k], inside );"


def EmitDistanceFromInsideDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "virtual void DistFromInside_l( Double_t const *point, Double_t const *dir, Int_t iact, Double_t const * step, Double_t *safe , Double_t * dist, Int_t vecsize ) {"


def EmitDistanceFromOutsideDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "virtual void DistFromOutside_l( Double_t const *point, Double_t const *dir, Int_t iact, Double_t const * step, Double_t *safe , Double_t * dist, Int_t vecsize ) {"


def EmitCallToDistanceFromOutside( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "dist[k]= " + classname + "::DistFromOutside( (Double_t *) &point[3*k], (Double_t *) &dist[3*k], 3, step[k] , 0 );"


def EmitCallToDistanceFromInside( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "dist[k]= " + classname + "::DistFromInside( (Double_t *) &point[3*k], (Double_t *) &dist[3*k], 3, step[k] , 0 );"



def EmitClosingParen( indentlevel ):
    None
    print getIndentString( indentlevel ) + "}"


def EmitContainsLoop():
    None

def EmitSafety():
    None

def main():
    for shape in ClassNames:
        EmitCommonStuff( shape )
        EmitClassDef( shape )

    # for the Contains Method
        EmitContainsDecl( 1, shape )
        EmitLoopN( 3 )
        EmitCallToContains( 5 , shape )
        EmitClosingParen( 3 )
        EmitClosingParen( 1 )
        

    # for the Safety Method
        EmitSafetyDecl( 1, shape )
        EmitLoopN( 3 )
        EmitCallToSafety( 5 , shape )
        EmitClosingParen( 3 )
        EmitClosingParen( 1 )
        
    # for the DistanceFromInside Method
        EmitDistanceFromInsideDecl( 1, shape )
        EmitLoopN( 3 )
        EmitCallToDistanceFromInside( 5 , shape )
        EmitClosingParen( 3 )
        EmitClosingParen( 1 )

    # for the DistanceFromOutside Method
        EmitDistanceFromOutsideDecl( 1, shape )
        EmitLoopN( 3 )
        EmitCallToDistanceFromOutside( 5 , shape )
        EmitClosingParen( 3 )
        EmitClosingParen( 1 )
        
    # close the class
        EmitClassClosing(  ) 


if __name__ == "__main__":
    main()
