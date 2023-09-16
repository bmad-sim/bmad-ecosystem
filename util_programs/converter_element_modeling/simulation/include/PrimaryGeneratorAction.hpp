#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4MTRunManager.hh"
#endif

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the converter
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class
/// (see the macros provided with this example).

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  private:
    G4ParticleGun*  fParticleGun;

  public:
    PrimaryGeneratorAction(double in_particle_energy, G4ThreeVector in_particle_polarization);
    //PrimaryGeneratorAction(const PrimaryGeneratorAction&) = default;
    //PrimaryGeneratorAction& operator=(const PrimaryGeneratorAction&) = default;
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* event);

    // set methods
    void SetRandomFlag(G4bool value);

};

