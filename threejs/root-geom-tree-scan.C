#include <iostream>
#include "TGeoManager.h"
#include "TGeoNode.h"
void scan() {
// Load some root geometry
  TGeoVolume *top = gGeoManager->GetTopVolume();
  TGeoIterator iter(top);
  TGeoNode *current;
  // Inspect shape properties
  std::cout << "Top volume: " << top->GetName() << std::endl;
  top->InspectShape();
  TString path;
  while ((current=iter.Next())) {
    iter.GetPath(path);
    std::cout << "=IG===================================="
    std::cout << path << std::endl;
    std::cout << "=IG===================================="
    current->InspectNode();
  }
}
 