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
// $Id: B4aActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B4aActionInitialization.cc
/// \brief Implementation of the B4aActionInitialization class

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "B4aActionInitialization.hpp"
#include "B4PrimaryGeneratorAction.hpp"
#include "B4RunAction.hpp"
#include "B4aEventAction.hpp"
#include "B4aSteppingAction.hpp"
#include "B4DetectorConstruction.hpp"
#include "point_cache.hpp"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aActionInitialization::B4aActionInitialization
                            (B4DetectorConstruction* detConstruction, double in_particle_energy, PointCache* pc)
 : G4VUserActionInitialization(),
   fDetConstruction(detConstruction),
   m_in_particle_energy(in_particle_energy),
   m_point_cache(pc)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aActionInitialization::~B4aActionInitialization()
{
  //.auto sa = G4RunManager::GetRunManager()->GetUserSteppingAction();
  //.SetUserAction((B4aSteppingAction *) nullptr);
  //.delete sa;

  //.auto ea = G4RunManager::GetRunManager()->GetUserEventAction();
  //.SetUserAction((B4aEventAction *) nullptr);
  //.delete ea;

  //.auto ra = G4RunManager::GetRunManager()->GetUserRunAction();
  //.SetUserAction((B4RunAction *) nullptr);
  //.delete ra;

  //.auto pga = G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  //.SetUserAction((B4PrimaryGeneratorAction *) nullptr);
  //.delete pga;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aActionInitialization::BuildForMaster() const
{
  SetUserAction(new B4RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aActionInitialization::Build() const
{
  //std::cout << "B4aActionInitizalization::Build\n";
  SetUserAction(new B4PrimaryGeneratorAction(m_in_particle_energy));
  //std::cout << "B4PrimaryGeneratorAction constructed\n";
  SetUserAction(new B4RunAction);
  //std::cout << "B4PrimaryGeneratorAction set\n";
  B4aEventAction* eventAction = new B4aEventAction;
  //std::cout << "B4aEventAction constructed\n";
  SetUserAction(eventAction);
  //std::cout << "B4aEventAction set\n";
  SetUserAction(new B4aSteppingAction(fDetConstruction,eventAction, m_point_cache));
  //std::cout << "B4aSteppingAction constructed and set\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
