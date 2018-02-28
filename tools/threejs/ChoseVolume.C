const char *CheckVolume(TGeoShape *geoShape)
{
  // Chooses the object and method that should be used for processing object
  const char *clsname = geoShape->ClassName();
  Boolean VolFound    = clsname;

  // process different shapes
  if (strcmp(clsname, "TGeoBBox") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoParaboloid") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoSphere") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoConeSeg") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoCone") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoTubeSeg") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoTube") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoPcon") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoTorus") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoPgon") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoHype") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoXtru") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoScaledShape") == 0) {
    cout << "\n---\n-> Geometry of volume: ", << clsname << " \n---\n";
  } else if (strcmp(clsname, "TGeoArb8") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoPara") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoTrap") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoGtra") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoTrd1") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoTrd2") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoCtub") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoEltu") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoCompositeShape") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoUnion") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoIntersection") == 0) {
    VolFound = "";
  } else if (strcmp(clsname, "TGeoSubtraction") == 0) {
    VolFound = "";
  } else {
    VolFound = "";
    cout << "ChooseObject", "ERROR! Solid CANNOT be processed, solid is NOT supported: " << clsname << "\n";
  }
  return VolFound;
}

/*
XMLNodePointer_t TGDMLWrite::CreateCommonBoolN(TGeoCompositeShape *geoShape)
{
// Creates common part of union intersection and subtraction nodes
   XMLNodePointer_t mainN, ndR, ndL, childN;
   TString nodeName = GenName(geoShape->GetName(), TString::Format("%p", geoShape));
   TString lboolType;
   TGeoBoolNode::EGeoBoolType boolType = geoShape->GetBoolNode()->GetBooleanOperator();
   switch (boolType) {
      case TGeoBoolNode::kGeoUnion:
         lboolType = "union";
         break;
      case TGeoBoolNode::kGeoSubtraction:
         lboolType = "subtraction";
         break;
      case TGeoBoolNode::kGeoIntersection:
         lboolType = "intersection";
         break;
   }

   TGDMLWrite::Xyz lrot = GetXYZangles(geoShape->GetBoolNode()->GetLeftMatrix()->Inverse().GetRotationMatrix());
   const Double_t  *ltr = geoShape->GetBoolNode()->GetLeftMatrix()->GetTranslation();
   TGDMLWrite::Xyz rrot = GetXYZangles(geoShape->GetBoolNode()->GetRightMatrix()->Inverse().GetRotationMatrix());
   const Double_t  *rtr = geoShape->GetBoolNode()->GetRightMatrix()->GetTranslation();

   //specific case!
   //Ellipsoid tag preparing
   //if left == TGeoScaledShape AND right  == TGeoBBox
   //   AND if TGeoScaledShape->GetShape == TGeoSphere
   TGeoShape *leftS = geoShape->GetBoolNode()->GetLeftShape();
   TGeoShape *rightS = geoShape->GetBoolNode()->GetRightShape();
   if (strcmp(leftS->ClassName(), "TGeoScaledShape") == 0 &&
       strcmp(rightS->ClassName(), "TGeoBBox") == 0) {
      if (strcmp(((TGeoScaledShape *)leftS)->GetShape()->ClassName(), "TGeoSphere") == 0) {
         if (lboolType == "intersection") {
            mainN = CreateEllipsoidN(geoShape, nodeName);
            return mainN;
         }
      }
   }

   Xyz translL, translR;
   //translation
   translL.x = ltr[0];
   translL.y = ltr[1];
   translL.z = ltr[2];
   translR.x = rtr[0];
   translR.y = rtr[1];
   translR.z = rtr[2];

   //left and right nodes are created here also their names are created
   ndL = ChooseObject(geoShape->GetBoolNode()->GetLeftShape());
   ndR = ChooseObject(geoShape->GetBoolNode()->GetRightShape());

   //retrieve left and right node names by their pointer to make reference
   TString lname = fNameList->fLst[TString::Format("%p", geoShape->GetBoolNode()->GetLeftShape())];
   TString rname = fNameList->fLst[TString::Format("%p", geoShape->GetBoolNode()->GetRightShape())];

   //left and right nodes appended to main structure of nodes (if they are not already there)
   if (ndL != NULL) {
      fGdmlE->AddChild(fSolidsNode, ndL);
      fSolCnt++;
   } else {
      if (lname.Contains("missing_") || lname == "") {
         Info("CreateCommonBoolN", "ERROR! Left node is NULL - Boolean Shape will be skipped");
         return NULL;
      }
   }
   if (ndR != NULL) {
      fGdmlE->AddChild(fSolidsNode, ndR);
      fSolCnt++;
   } else {
      if (rname.Contains("missing_") || rname == "") {
         Info("CreateCommonBoolN", "ERROR! Right node is NULL - Boolean Shape will be skipped");
         return NULL;
      }
   }

   //create union node and its child nodes (or intersection or subtraction)
   /* <union name="...">
    *   <first ref="left name" />
    *   <second ref="right name" />
    *   <firstposition .../>
    *   <firstrotation .../>
    *   <position .../>
    *   <rotation .../>
    * </union>
   */
mainN = fGdmlE->NewChild(0, 0, lboolType.Data(), 0);
fGdmlE->NewAttr(mainN, 0, "name", nodeName);

//<first> (left)
childN = fGdmlE->NewChild(0, 0, "first", 0);
fGdmlE->NewAttr(childN, 0, "ref", lname);
fGdmlE->AddChild(mainN, childN);

//<second> (right)
childN = fGdmlE->NewChild(0, 0, "second", 0);
fGdmlE->NewAttr(childN, 0, "ref", rname);
fGdmlE->AddChild(mainN, childN);

//<firstposition> (left)
if ((translL.x != 0.0) || (translL.y != 0.0) || (translL.z != 0.0)) {
  childN = CreatePositionN((nodeName + lname + "pos").Data(), translL, "firstposition");
  fGdmlE->AddChild(mainN, childN);
}
//<firstrotation> (left)
if ((lrot.x != 0.0) || (lrot.y != 0.0) || (lrot.z != 0.0)) {
  childN = CreateRotationN((nodeName + lname + "rot").Data(), lrot, "firstrotation");
  fGdmlE->AddChild(mainN, childN);
}
//<position> (right)
if ((translR.x != 0.0) || (translR.y != 0.0) || (translR.z != 0.0)) {
  childN = CreatePositionN((nodeName + rname + "pos").Data(), translR, "position");
  fGdmlE->AddChild(mainN, childN);
}
//<rotation> (right)
if ((rrot.x != 0.0) || (rrot.y != 0.0) || (rrot.z != 0.0)) {
  childN = CreateRotationN((nodeName + rname + "rot").Data(), rrot, "rotation");
  fGdmlE->AddChild(mainN, childN);
}

return mainN;
}
* /
