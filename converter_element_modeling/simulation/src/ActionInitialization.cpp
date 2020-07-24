#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "ActionInitialization.hpp"
#include "PrimaryGeneratorAction.hpp"
#include "RunAction.hpp"
#include "EventAction.hpp"
#include "SteppingAction.hpp"
#include "DetectorConstruction.hpp"
#include "point_cache.hpp"

ActionInitialization::ActionInitialization
                            (DetectorConstruction* detConstruction, double in_particle_energy, G4ThreeVector polarization, PointCache* pc)
 : G4VUserActionInitialization(),
   m_DetConstruction(detConstruction),
   m_in_particle_energy(in_particle_energy),
   m_particle_polarization(polarization),
   m_point_cache(pc)
{}

ActionInitialization::~ActionInitialization() {}


void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}


void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction(m_in_particle_energy, m_particle_polarization));
  SetUserAction(new RunAction);
  EventAction* eventAction = new EventAction;
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(m_DetConstruction,eventAction, m_point_cache));
}

