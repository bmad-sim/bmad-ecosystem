#include "G4EmStandardPhysics.hh"
#include "G4EmParameters.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4EmProcessOptions.hh"
#include "G4ProcessManager.hh"

#include "PhysListEmPolarized.hpp"
#include "StepMax.hpp"
#include "polarized_physics.hpp"


PhysicsList::PhysicsList()
  : G4VModularPhysicsList(),
    fEmPhysicsList(0),
    fEmName("polarized") {

      G4EmParameters::Instance();
      SetVerboseLevel(0);
      fEmPhysicsList = new PhysListEmPolarized();
}


PhysicsList::~PhysicsList() {}


void PhysicsList::ConstructParticle() {
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}


void PhysicsList::ConstructProcess() {
  // Transportation
  AddTransportation();

  // Electromagnetic physics list
  fEmPhysicsList->ConstructProcess();

  // step limitation (as a full process)
  AddStepMax();
}


void PhysicsList::AddPhysicsList(const G4String& name) {
  // Sets the physics to either standard or polarized
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == fEmName) return;

  if (name == "standard") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();

  } else if (name == "polarized") {

    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmPolarized();

  } else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}


void PhysicsList::AddStepMax() {
  // Step limitation seen as a process
  StepMax* stepMaxProcess = new StepMax();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()){
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (stepMaxProcess->IsApplicable(*particle) && pmanager)
          pmanager ->AddDiscreteProcess(stepMaxProcess);
  }
}

