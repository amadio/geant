#include "G4GDMLParser.hh"
#include "G4GDMLWriteStructure.hh"
#include "G4VRangeToEnergyConverter.hh"
#include "G4RToEConvForGamma.hh"
#include "G4RToEConvForElectron.hh"
#include "G4RToEConvForPositron.hh"
#include "G4RToEConvForProton.hh"
#include "G4UserLimits.hh"
#include <xercesc/dom/DOM.hpp>

class MyGDMLWriteStructure: public G4GDMLWriteStructure {
public:
   virtual void AddExtension(xercesc::DOMElement* volumeElement,
			     const G4LogicalVolume* const glv) {
      using CLHEP::GeV;
      using CLHEP::cm;
      using CLHEP::s;
      xercesc::DOMElement* auxiliaryElement = 0;
      std::stringstream ss;
      const char* cutnames[4] = {"pcutg","pcutem","pcutep","pcutp"};
      static G4VRangeToEnergyConverter *converter[4];
      static G4bool ifirst = TRUE;
      if(ifirst) {
	 converter[0] = new G4RToEConvForGamma();
	 converter[1] = new G4RToEConvForElectron();
	 converter[2] = new G4RToEConvForPositron();
	 converter[3] = new G4RToEConvForProton();
	 ifirst = FALSE;
      }
 
      auxiliaryElement = NewElement("auxiliary");
      auxiliaryElement->setAttributeNode(NewAttribute("auxtype","G4Region"));
      auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",glv->GetRegion()->GetName()));
      volumeElement->appendChild(auxiliaryElement);

      auxiliaryElement = NewElement("auxiliary");
      auxiliaryElement->setAttributeNode(NewAttribute("auxtype","pcutunit"));
      auxiliaryElement->setAttributeNode(NewAttribute("auxvalue","GeV"));
      volumeElement->appendChild(auxiliaryElement);

      //     G4cout << "I have been called " << glv->GetName() << " in region " << glv->GetRegion()->GetName() << G4endl;
      G4ProductionCuts *cuts = glv->GetRegion()->GetProductionCuts();

      for(G4int ic=0; ic<4; ++ic) {
	 ss.clear(); ss.str("");
	 ss << converter[ic]->Convert(cuts->GetProductionCut(ic),glv->GetMaterial())/GeV;
	 //	 G4cout << cutnames[ic] << " = " << ss.str() << G4endl;
	 auxiliaryElement = NewElement("auxiliary");
	 auxiliaryElement->setAttributeNode(NewAttribute("auxtype",cutnames[ic]));
	 auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",ss.str()));
	 volumeElement->appendChild(auxiliaryElement);
      }
      G4UserLimits *ulim = glv->GetUserLimits();
      if(ulim) {
	 G4Track gt;
	 G4double value = 0;
	 
	 // Max Allowed Step
	 value = ulim->GetMaxAllowedStep(gt);
	 if(value < 0.99 * DBL_MAX) {
	    auxiliaryElement = NewElement("auxiliary");
	    auxiliaryElement->setAttributeNode(NewAttribute("auxtype","G4MAStep"));
	    auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",value/cm));
	    volumeElement->appendChild(auxiliaryElement);
	 }

	 // Max Track Length
	 value = ulim->GetUserMaxTrackLength(gt);
	 if(value < 0.99 * DBL_MAX) {
	    auxiliaryElement = NewElement("auxiliary");
	    auxiliaryElement->setAttributeNode(NewAttribute("auxtype","G4MTLength"));
	    auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",value/cm));
	    volumeElement->appendChild(auxiliaryElement);
	 }

	 // Max Time
	 value = ulim->GetUserMaxTime(gt);
	 if(value < 0.99 * DBL_MAX) {
	    auxiliaryElement = NewElement("auxiliary");
	    auxiliaryElement->setAttributeNode(NewAttribute("auxtype","G4MTime"));
	    auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",value/s));
	    volumeElement->appendChild(auxiliaryElement);
	 }

	 // Min Kinetic Energy
	 value = ulim->GetUserMinEkine(gt);
	 if(value > 0) {
	    auxiliaryElement = NewElement("auxiliary");
	    auxiliaryElement->setAttributeNode(NewAttribute("auxtype","G4MKEnergy"));
	    auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",value/GeV));
	    volumeElement->appendChild(auxiliaryElement);
	 }

	 // Min Range
	 value = ulim->GetUserMinRange(gt);
	 if(value > 0) {
	    auxiliaryElement = NewElement("auxiliary");
	    auxiliaryElement->setAttributeNode(NewAttribute("auxtype","G4MRange"));
	    auxiliaryElement->setAttributeNode(NewAttribute("auxvalue",value/cm));
	    volumeElement->appendChild(auxiliaryElement);
	 }
      }
   }
};
