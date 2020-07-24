#pragma once
#include "G4VUserActionInitialization.hh"
#include "point_cache.hpp"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4MTRunManager.hh"
#endif

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization
{
  private:
    DetectorConstruction* m_DetConstruction;
    double m_in_particle_energy;
    G4ThreeVector m_particle_polarization;
    PointCache* m_point_cache;

  public:
    ActionInitialization(DetectorConstruction*, double in_particle_energy, G4ThreeVector polarization, PointCache* pc);
    //ActionInitialization(const ActionInitialization&) = default;
    //ActionInitialization& operator=(const ActionInitialization&) = default;
    virtual ~ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

