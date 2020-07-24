#pragma once
#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class PhysListEmPolarized : public G4VPhysicsConstructor
{
  public:
    PhysListEmPolarized(const G4String& name = "polarized");
   //~PhysListEmPolarized();

  public:
    // This method is dummy for physics
    virtual void ConstructParticle() {};

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();
};


