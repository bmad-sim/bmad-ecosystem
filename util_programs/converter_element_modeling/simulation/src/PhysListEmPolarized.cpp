#include "PhysListEmPolarized.hpp"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4eMultipleScattering.hh"

#include "G4PolarizedCompton.hh"
#include "G4PolarizedGammaConversion.hh"
#include "G4ePolarizedIonisation.hh"
#include "G4ePolarizedBremsstrahlung.hh"
#include "G4eplusPolarizedAnnihilation.hh"
#include "G4PolarizedPhotoElectricEffect.hh"

PhysListEmPolarized::PhysListEmPolarized(const G4String& name)
  : G4VPhysicsConstructor(name) {}


//PhysListEmPolarized::~PhysListEmPolarized() {}

void PhysListEmPolarized::ConstructProcess()
{
  // Add standard EM Processes

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      pmanager->AddDiscreteProcess(new G4PolarizedPhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4PolarizedCompton);
      pmanager->AddDiscreteProcess(new G4PolarizedGammaConversion);

    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4eMultipleScattering,   -1,1,1);
      pmanager->AddProcess(new G4ePolarizedIonisation,  -1,2,2);
      pmanager->AddProcess(new G4ePolarizedBremsstrahlung,      -1,3,3);

    } else if (particleName == "e+") {
      pmanager->AddProcess(new G4eMultipleScattering,  -1, 1,1);
      pmanager->AddProcess(new G4ePolarizedIonisation, -1, 2,2);
      pmanager->AddProcess(new G4ePolarizedBremsstrahlung,    -1, 3,3);
      pmanager->AddProcess(new G4eplusPolarizedAnnihilation,   0,-1,4);
    }
  }
}

