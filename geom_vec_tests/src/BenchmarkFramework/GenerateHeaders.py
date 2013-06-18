ClassName = [ "TGeoBBox", "TGeoCone", "TGeoPgon" ]
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
    print getIndentString( indentlevel ) + "for(unsigned int k=0;k < vecsize, ++k){"


def EmitContainsDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "void Contains_l( Double_t const *point, Double_t const *dir, Int_t iact=1, Double_t const * step, Double_t * safe=0, Bool_t * isin , Int_t vecsize ) {"


def EmitCallToContains( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "isin[k]= " + classname + "::Contains( &point[3*k] );"


def EmitSafetyDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "void Safety_l( Double_t const *point, Bool_t const *inside, Double_t * safe , Int_t vecsize ) {"


def EmitCallToSafety( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "safe[k]= " + classname + "::Safety( &point[3*k], inside[k] );"



def EmitDistanceFromInsideDecl( indentlevel, classname ):
    print getIndentString( indentlevel ) + "void DistanceFromInside_l( Double_t const *point, Double_t const *dir, Int_t iact=1, Double_t const * step, Double_t *safe=0 , Double_t * dist, Int_t vecsize ) {"


def EmitCallToDistanceFromInside( indentlevel, classname ):
    """
    """
    print getIndentString( indentlevel ) + "dist[k]= " + classname + "::DistanceFromInside( &point[3*k], &dist[3*k], iact, step[k] , safe );"



def EmitClosingParen( indentlevel ):
    None
    print getIndentString( indentlevel ) + "}"


def EmitContainsLoop():
    None

def EmitSafety():
    None

def main():
    EmitCommonStuff( "Foo" )
    EmitClassDef( "Foo" )

    # for the Contains Method
    EmitContainsDecl( 1, "Foo" )
    EmitLoopN( 3 )
    EmitCallToContains( 5 , "Foo" )
    EmitClosingParen( 1 )


    # for the Safety Method
    EmitSafetyDecl( 1, "Foo" )
    EmitLoopN( 3 )
    EmitCallToSafety( 5 , "Foo" )
    EmitClosingParen( 1 )


    # for the DistanceFromInside Method
    EmitDistanceFromInsideDecl( 1, "Foo" )
    EmitLoopN( 3 )
    EmitCallToDistanceFromInside( 5 , "Foo" )
    EmitClosingParen( 1 )

    # close the class
    EmitClassClosing(  ) 

if __name__ == "__main__":
    main()
