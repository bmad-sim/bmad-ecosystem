#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  private:
    double absorberXY;
    double absorber_thickness;
    std::string material_name;

    G4VPhysicalVolume*   fAbsorberPV; // the absorber physical volume
    G4VPhysicalVolume*   fGapPV;      // the gap physical volume

    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

  public:
    DetectorConstruction(const std::string& target_material, double target_thickness);
    //DetectorConstruction(const DetectorConstruction&) = default;
    //DetectorConstruction& operator=(const DetectorConstruction&) = default;
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    //virtual void ConstructSDandField();

    inline const G4VPhysicalVolume* GetAbsorberPV() const { return fAbsorberPV; };
    inline const G4VPhysicalVolume* GetGapPV() const { return fGapPV; };

};

