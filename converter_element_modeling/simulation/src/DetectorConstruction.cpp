#include <fstream>
#include "DetectorConstruction.hpp"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction(const std::string& target_material, double target_thickness)
 : G4VUserDetectorConstruction(),
   absorberXY(9),
   absorber_thickness(target_thickness),
   material_name(target_material),
   fAbsorberPV(0),
   fGapPV(0),
   fCheckOverlaps(true) {}

DetectorConstruction::~DetectorConstruction() {}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}


void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double calorSizeXY  = absorberXY*cm;
  G4double calorThickness = absorber_thickness*cm;

  G4double worldSizeXY = 2 * absorberXY*cm;
  G4double worldSizeZ  = 2 * absorber_thickness*cm;

  G4double detectorXY  = absorberXY*cm;
  G4double detectorThickness = 1.e-3*mm;

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, pressure, temperature, fractionmass;
  G4String name, symbol;
  G4int nel;


  density = 7.87*g/cm3;
  a = 55.85*g/mole;
  G4Material* Iron = new G4Material(name="Iron", z=26., a, density);

  density = 8.96*g/cm3;
  a = 63.55*g/mole;
  G4Material* Copper = new G4Material(name="Copper", z=29., a, density);

  density = 19.30*g/cm3;
  a = 183.85*g/mole;
  G4Material* Tungsten = new G4Material(name="Tungsten", z=74., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Lead = new G4Material(name="Lead", z=82., a, density);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N2", z=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O2", z=8., a);
  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, 0.8);
  Air->AddElement(elO, 0.2);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* G4vacuum = new G4Material(name="G4vacuum", density, nel=1, kStateGas, temperature, pressure);
  G4vacuum->AddMaterial(Air, fractionmass=1.);


  // Get materials
  G4Material* absorberMaterial;
  if (material_name=="tungsten") {  absorberMaterial = Tungsten;}
  else if (material_name=="copper") {  absorberMaterial = Copper;}
  else if (material_name=="lead") {  absorberMaterial = Lead;}
  else if (material_name=="iron") {  absorberMaterial = Iron;}
  else {
    std::cout << "Material name " << material_name
      << " not recognized, defaulting to tungsten\n";
    absorberMaterial = Tungsten;
  }

  // World
  G4VSolid* worldS = new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2);
  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, G4vacuum, "World");

  G4VPhysicalVolume* worldPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  // Converter
  G4VSolid* calorimeterS = new G4Box("absorber", calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
  G4LogicalVolume* calorLV = new G4LogicalVolume(calorimeterS, absorberMaterial, "absorber");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume
                 "absorber",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  G4VSolid* detectorS = new G4Box("detector", detectorXY/2, detectorXY/2, detectorThickness/2); // its size
  G4LogicalVolume* detLV = new G4LogicalVolume(detectorS, G4vacuum, "det_LV");
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0, 0, calorThickness/2.+detectorThickness),  // at (0,0,0)
                 detLV,            // its logical volume
                 "detect_layer",      // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  // Visualization attributes
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  // Always return the physical World
  return worldPV;
}

//void DetectorConstruction::ConstructSDandField()
//{
//  // Create global magnetic field messenger.
//  // Uniform magnetic field is then created automatically if
//  // the field value is not zero.
//  G4ThreeVector fieldValue = G4ThreeVector();
//  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
//  fMagFieldMessenger->SetVerboseLevel(1);
//
//  // Register the field messenger for deleting
//  G4AutoDelete::Register(fMagFieldMessenger);
//}

