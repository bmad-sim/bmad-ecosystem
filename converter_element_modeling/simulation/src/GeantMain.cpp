//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: exampleB4a.cc 95508 2016-02-12 13:52:06Z gcosmo $
//
/// \file exampleB4a.cc
/// \brief Main program of the B4a example

#include <string>

#include "B4DetectorConstruction.hpp"
#include "B4aActionInitialization.hpp"
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

#include "Randomize.hh"
#include <ctime>

#include <cmath>
#include <utility>

#include "G4VisExecutive.hh"
#include "binner.hpp"
#include "GeantMain.hpp"
//#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RunManager* Initialize_Geant(void) {
  // Runs the one-time setup for geant and returns a pointer to the resultant G4RunManager
#ifdef G4MULTITHREADED
  int nThreads = 4;
#endif
  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(time(NULL));

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) {
    static runManager->SetNumberOfThreads(nThreads);
  }
#else
  static G4RunManager * runManager = new G4RunManager;
#endif
  G4UImanager* UI = G4UImanager::GetUIpointer();

  silentSession * quiet = new silentSession;
  UI->SetCoutDestination(quiet);

//  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  runManager->SetUserInitialization(physicsList);

  // Visualization
  //G4VisManager* visManager = new G4VisExecutive;
  //visManager->Initialize();

  G4EventManager::GetEventManager()->GetTrackingManager()->SetStoreTrajectory(1);


  return runManager;
}



int run_simulation(G4RunManager *runManager, const std::string& target_material, double in_energy, double target_thickness, BinnerBase* binner, int run_length) {
  //G4String session;

  // Clear sim_data.txt and write header
  //std::ofstream sim_data;
  //sim_data.open("sim_data.txt");
  //sim_data << "E\tx\ty\tz\tr\tpx\tpy\tpz\tx_adj\tpx_adj\tpz_adj\tdrds\n";
  //sim_data.close();

  // Detect interactive mode (if no macro provided) and define UI session
  //
  //G4UIExecutive* ui = 0;
  //if ( ! macro.size() ) {
  //  ui = new G4UIExecutive(argc, argv, session);
  //}


  // Set mandatory initialization classes
  //
  B4DetectorConstruction* detConstruction = new B4DetectorConstruction(target_material, target_thickness);
  runManager->SetUserInitialization(detConstruction);

  B4aActionInitialization* actionInitialization
     = new B4aActionInitialization(detConstruction, in_energy, binner);
  runManager->SetUserInitialization(actionInitialization);

  runManager->Initialize();

  // Draw the converter
  //G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  //if (pVVisManager) {
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/open OGL 800x800-0+0");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/sceneHandler/create OGL");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/create 800x800-0+0");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/drawVolume");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/flush");
  //}
  //G4UImanager* UI = G4UImanager::GetUIpointer();
  //UI->ApplyCommand("/run/verbose 0");
  //UI->ApplyCommand("/event/verbose 0");
  //UI->ApplyCommand("/tracking/verbose 0");
  //UI->ApplyCommand("/printProgress 10001");
  runManager->SetPrintProgress(run_length+1);
  runManager->BeamOn(run_length);

  // Initialize visualization
  //
  //G4VisManager* visManager = new G4VisExecutive;
  //// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  //// G4VisManager* visManager = new G4VisExecutive("Quiet");
  //visManager->Initialize();

  //// Get the pointer to the User Interface manager
  //G4UImanager* UImanager = G4UImanager::GetUIpointer();

  //// Process macro or start UI session
  ////
  //if ( macro.size() ) {
  //  // batch mode
  //  G4String command = "/control/execute ";
  //  UImanager->ApplyCommand(command+macro);
  //}
  //else  {
  //  // interactive mode : define UI session
  //  UImanager->ApplyCommand("/control/execute init_vis.mac");
  //  if (ui->IsGUI()) {
  //    UImanager->ApplyCommand("/control/execute gui.mac");
  //  }
  //  ui->SessionStart();
  //  delete ui;
  //}

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  //delete runManager;
  //delete detConstruction;
  //delete actionInitialization;
  return 0;
}

