#include <cmath>
#include <cstdio>
#include <algorithm>

#include "SteppingAction.hpp"
#include "EventAction.hpp"
#include "DetectorConstruction.hpp"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4VisManager.hh"

SteppingAction::SteppingAction(
                      const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction,
                      PointCache* cache)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    prev_track_id(-1),
    out_energy_min(0),
    out_energy_max(1000),
    in_particle_angle(0),
    point_cache(cache),
    data_vec(nullptr) {
      GetVec();
      data_vec->reserve(40000);
    }


SteppingAction::~SteppingAction() {
  ReturnVec();
}

void SteppingAction::GetVec() const {
  if (!data_vec) {
    point_cache->Lock();
    data_vec = point_cache->GetVec();
    point_cache->Unlock();
  }
}

void SteppingAction::ReturnVec() const {
  if (data_vec) {
    point_cache->Lock();
    point_cache->ReturnVec(data_vec);
    point_cache->Unlock();
    data_vec = nullptr;
  }
}


void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4VPhysicalVolume* curPV  = step->GetPreStepPoint()->GetPhysicalVolume();
  G4String name = curPV->GetName();
  name.assign(name,0,12);

  if (name=="detect_layer")
  {
    int track_id=step->GetTrack()->GetTrackID();
    if (track_id != prev_track_id) // prevent double counting
    {
      G4String pname=step->GetTrack()->GetDefinition()->GetParticleName();
      if (pname == "e+")
      {
        double px_step=step->GetPreStepPoint()->GetMomentum().x()/MeV;
        double py_step=step->GetPreStepPoint()->GetMomentum().y()/MeV;
        double pz_step=step->GetPreStepPoint()->GetMomentum().z()/MeV;
        double engy=std::sqrt(px_step*px_step+py_step*py_step+pz_step*pz_step);
        if (engy>out_energy_min)
        if (engy<out_energy_max)
        {
          G4double x=step->GetTrack()->GetPosition().x()/cm;
          G4double y=step->GetTrack()->GetPosition().y()/cm;
          G4double r=std::sqrt(x*x+y*y);
          double theta = std::atan2(y, x);
          double dxp_ds = ((x/r) * px_step + (y/r) * py_step)/pz_step;
          double dyp_ds = ((-y/r) * px_step + (x/r) * py_step)/pz_step;
          auto pol = step->GetTrack()->GetPolarization();

          GeantParticle out;
          out.E = engy;
          out.r = r;
          out.theta = theta;
          out.dxds = dxp_ds;
          out.dyds = dyp_ds;
          out.polx = pol.x();
          out.poly = pol.y();
          out.polz = pol.z();

          data_vec->push_back(out);
        }
      }
      prev_track_id=track_id;
    }
  }
}

