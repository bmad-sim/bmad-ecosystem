#include "StepMax.hpp"

#include "G4ParticleDefinition.hh"
#include "G4Step.hh"


StepMax::StepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),
   fMaxChargedStep(DBL_MAX) { }

StepMax::~StepMax() {}


G4bool StepMax::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0.);
}


void StepMax::SetMaxStep(G4double step) {fMaxChargedStep = step;}


G4double StepMax::PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
                                                  G4double,
                                                  G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

  G4double ProposedStep = DBL_MAX;

  if((fMaxChargedStep > 0.) &&
     (aTrack.GetVolume() != 0) &&
     (aTrack.GetVolume()->GetName() != "World"))
     ProposedStep = fMaxChargedStep;

  return ProposedStep;
}


G4VParticleChange* StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
   // do nothing
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}

