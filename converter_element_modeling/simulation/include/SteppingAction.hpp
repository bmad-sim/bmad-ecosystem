#pragma once
#include "G4UserSteppingAction.hh"
#include "bin.hpp"
#include "point_cache.hpp"

class DetectorConstruction;
class EventAction;

// This class records the properties of positrons that make it
// through the converter, and stores them to a vector obtained
// from the point cache.  These results are read in by the
// single-threaded binning code.

class SteppingAction : public G4UserSteppingAction
{
  private:
    const DetectorConstruction* fDetConstruction;
    EventAction*  fEventAction;
    int prev_track_id;
    double out_energy_min;
    double out_energy_max;
    double in_particle_angle;
    PointCache* point_cache;
    mutable std::vector<GeantParticle>* data_vec;

  public:
    SteppingAction(const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction,
                      PointCache* cache);
    //SteppingAction(const SteppingAction&) = default;
    //SteppingAction& operator=(const SteppingAction&) = default;
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step* step);
    virtual void GetVec() const;
    virtual void ReturnVec() const;
};

