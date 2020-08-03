#include <string>
#include <thread>

#include "DetectorConstruction.hpp"
#include "ActionInitialization.hpp"
#include "silentSession.hpp"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "polarized_physics.hpp"

#include "Randomize.hh"
#include <ctime>

#include <cmath>
#include <utility>

#include "G4VisExecutive.hh"
#include "GeantMain.hpp"


G4RunManager* Initialize_Geant(void) {
  // Runs the one-time setup for geant and returns a pointer to the resultant G4RunManager

#ifdef G4MULTITHREADED
  int nThreads = std::thread::hardware_concurrency();
  std::cerr << "Note: " << nThreads << " CPU threads detected for use\n";
#endif

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(time(NULL));

  // Construct the run manager
#ifdef G4MULTITHREADED
  static G4MTRunManager *runManager = new G4MTRunManager;
  if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
  }
#else
  static G4RunManager *runManager = new G4RunManager;
#endif

  // Disable excessive printouts from geant
  G4UImanager* UI = G4UImanager::GetUIpointer();
  silentSession * quiet = new silentSession;
  UI->SetCoutDestination(quiet);

  auto *physicsList = new PhysicsList;
  physicsList->AddPhysicsList("polarized");
  runManager->SetUserInitialization(physicsList);


  //G4EventManager::GetEventManager()->GetTrackingManager()->SetStoreTrajectory(1);

  return runManager;
}



int run_simulation(G4RunManager *runManager, const std::string& target_material, double in_energy, double target_thickness, G4ThreeVector polarization, PointCache* point_cache, int run_length) {
  // Alters the properties of the detector and incoming particles if necessary,
  // and runs the simulation for run_length number of incoming particles

  // Store previous parameters for comparison
  static std::string old_material = "";
  static double old_thickness = 0;
  static double old_energy;
  static G4ThreeVector old_polarization;

  DetectorConstruction *detConstruction = nullptr;
  bool new_detector = false, new_gun = false;

  if (target_thickness != old_thickness || target_material != old_material) {
    detConstruction = new DetectorConstruction(target_material, target_thickness);
    runManager->SetUserInitialization(detConstruction);
    old_thickness = target_thickness;
    old_material = target_material;
    new_detector = true;
  }

  if (in_energy != old_energy || polarization != old_polarization) {
    ActionInitialization* actionInitialization
       = new ActionInitialization(detConstruction, in_energy, polarization, point_cache);
    runManager->SetUserInitialization(actionInitialization);
    old_energy = in_energy;
    old_polarization = polarization;
    new_gun = true;
  }

  if (new_detector || new_gun)
    runManager->Initialize();

  runManager->SetPrintProgress(run_length+1);
  runManager->BeamOn(run_length);

  return 0;
}

