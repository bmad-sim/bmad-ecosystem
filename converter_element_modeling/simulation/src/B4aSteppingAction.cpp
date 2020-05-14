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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include <cmath>
#include <cstdio>

#include "B4aSteppingAction.hpp"
#include "B4aEventAction.hpp"
#include "B4DetectorConstruction.hpp"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction,
                      PointCache* cache)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    prev_track_id(-1),
    out_energy_min(0),
    out_energy_max(1000),
    in_particle_angle(0),
    point_cache(cache) {
      //std::cout << "Constructor for B4aSteppingAction\n";
      point_cache->Lock();
      data_vec = point_cache->GetVec();
      point_cache->Unlock();
      //std::cout << "data_vec acquired\n";
      data_vec->reserve(40000);
      //std::cout << "data_vec reserved\n";
    }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction() {
  point_cache->Lock();
  point_cache->ReturnVec(data_vec);
  point_cache->Unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
  using namespace std;
  static size_t num_tracks = 0;
// Collect energy and track length step by step
  G4VPhysicalVolume* curPV  = step->GetPreStepPoint()->GetPhysicalVolume();
  G4String name = curPV->GetName();
  name.assign(name,0,12);

  if (name=="detect_layer")
  {
    int track_id=step->GetTrack()->GetTrackID();
    if (track_id != prev_track_id)
    {
      G4String pname=step->GetTrack()->GetDefinition()->GetParticleName();
      if (pname == "e+")
      {
        double px_step=step->GetPreStepPoint()->GetMomentum().x()/MeV;
        double py_step=step->GetPreStepPoint()->GetMomentum().y()/MeV;
        double pz_step=step->GetPreStepPoint()->GetMomentum().z()/MeV;
        double engy=sqrt(px_step*px_step+py_step*py_step+pz_step*pz_step);
        if (engy>out_energy_min)
        if (engy<out_energy_max)
        {
          G4double x=step->GetTrack()->GetPosition().x()/cm;
          G4double y=step->GetTrack()->GetPosition().y()/cm;
          G4double r=sqrt(x*x+y*y);
          double dxp_ds = ((x/r) * px_step + (y/r) * py_step)/pz_step;
          double dyp_ds = ((-y/r) * px_step + (x/r) * py_step)/pz_step;
          data_vec->push_back({engy, r, dxp_ds, dyp_ds});
        }
      }
      prev_track_id=track_id;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
