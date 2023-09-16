#pragma once
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PhysicsListMessenger;
class G4VPhysicsConstructor;

class PhysicsList: public G4VModularPhysicsList
{
  private:
    G4VPhysicsConstructor*  fEmPhysicsList;
    G4String fEmName;

  public:
    PhysicsList();
    virtual ~PhysicsList();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void AddPhysicsList(const G4String& name);

    void AddStepMax();
};

