
#ifndef PHYSICS_LIST1
#define PHYSICS_LIST1


#include "PhysicsList.h"


namespace geantphysics {

class PhysicsList1 : public PhysicsList {

public:
  PhysicsList1(const std::string &name);
 ~PhysicsList1();

  virtual void Initialize();
};


} // namespace geantphysics

#endif // PHYSICS_LIST1
