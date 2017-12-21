
#ifndef USERPHYSICSLIST_H
#define USERPHYSICSLIST_H

#include "PhysicsList.h"
#include <string>

namespace userapplication {

// A user physics list

class UserPhysicsList : public geantphysics::PhysicsList {
public:
  UserPhysicsList(const std::string &name);
 ~UserPhysicsList();
  virtual void Initialize();
};

} // userapplication


#endif // USERPHYSICSLIST_H
