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

#include "B4aSteppingAction.hpp"
#include "B4aEventAction.hpp"
#include "B4DetectorConstruction.hpp"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction,
                      BinnerBase* binner)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    prev_track_id(-1),
    out_energy_min(0),
    out_energy_max(1000),
    in_particle_angle(0),
    binner_ptr(binner)
{
  //using namespace std;
  //prev_track_id=-1;
  //out_energy_max=1.e+50;
  //out_energy_min=-1.;
  //in_particle_angle=0.0;
  //n_part_names=0;

  //string a_line;
  //string word[10];
  //ifstream conf_file("config.txt", ios::in);
  //if ( conf_file.is_open() )
  //{
  //  //read txt file1
  //  while ( ! conf_file.eof() )
  //  {
  //    getline(conf_file, a_line);
  //    if(a_line.size()>0)
  //    {
  //      while (a_line.substr(0,1)==" ") a_line.erase(0,1);
  //      if(a_line.size()>15)
  //      {
  //        istringstream vars(a_line);
  //        if (a_line.substr(0,14)=="out_energy_min")  vars >> word[0] >> word[1] >>out_energy_min;
  //        if (a_line.substr(0,14)=="out_energy_max")  vars >> word[0] >> word[1] >>out_energy_max;
  //        if (a_line.substr(0,17)=="in_particle_angle")  vars >> word[0] >> word[1] >> in_particle_angle;
  //        if (a_line.substr(0,18)=="absorber_thickness")  vars >> word[0] >> word[1] >>absorber_thickness;
  //        if (a_line.substr(0,18)=="out_particle_names")
  //        {
  //          a_line.erase(0,18);
  //          while (a_line.substr(0,1)==" " || a_line.substr(0,1)=="=") a_line.erase(0,1);
  //          while (a_line.size()>1)
  //          {
////              cout<<a_line<<endl;
  //            if (a_line.substr(0,1)=="!") break;
  //            int index=a_line.find(" ");
  //            if (index>0 && index<100)
  //            {
  //              part_names[n_part_names]=a_line.substr(0,index);
  //              n_part_names++;
  //              a_line.erase(0,index);
  //              while (a_line.substr(0,1)==" ") a_line.erase(0,1);
  //            }
  //            else if(index==-1)
  //            {
  //              part_names[n_part_names]=a_line.substr(0,a_line.size());
  //              n_part_names++;
  //              a_line.erase(0,index);
  //              while (a_line.substr(0,1)==" ") a_line.erase(0,1);
  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}
  //conf_file.close();
/*
  for (int i=0; i<n_part_names; i++)
  {
    cout<<"---------------------> "<< i<<"    "<<part_names[i]<<endl;
  }
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
  using namespace std;
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
      //for (int i=0; i<n_part_names; i++)
      if (pname == "e+")
      {
        double px_step=step->GetPreStepPoint()->GetMomentum().x()/MeV;
        double py_step=step->GetPreStepPoint()->GetMomentum().y()/MeV;
        double pz_step=step->GetPreStepPoint()->GetMomentum().z()/MeV;
        double engy=sqrt(px_step*px_step+py_step*py_step+pz_step*pz_step);
        if (engy>out_energy_min)
        if (engy<out_energy_max)
        {
          //G4int evtNb;
          //evtNb = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
          //G4double z=step->GetTrack()->GetPosition().z()/cm;
          G4double x=step->GetTrack()->GetPosition().x()/cm;
          //double x_adj = x-tan(in_particle_angle)*absorber_thickness;
          G4double y=step->GetTrack()->GetPosition().y()/cm;
          //G4double x_vert=step->GetTrack()->GetVertexPosition().x()/cm;
          //G4double y_vert=step->GetTrack()->GetVertexPosition().y()/cm;
          //G4double z_vert=step->GetTrack()->GetVertexPosition().z()/cm;
          G4double r=sqrt(x*x+y*y);
          //double px_adj = cos(in_particle_angle) * px_step - sin(in_particle_angle) * pz_step;
          //double pz_adj = sin(in_particle_angle) * px_step + cos(in_particle_angle) * pz_step;
          //G4double pr_step=sqrt(px_step*px_step + py_step*py_step);
          //G4double drds=pr_step/pz_step;
          double dxp_ds = (x/r) * px_step + (y/r) * py_step;
          double dyp_ds = (-y/r) * px_step + (x/r) * py_step;
          //sim_data.open( "sim_data.txt", ios::out | ios::app );
          //sim_data << engy << '\t' << r << '\t' << drds << '\n';
          binner_ptr->add_point({engy, r, dxp_ds, dyp_ds});
          //sim_data << engy << '\t'
          //  << x << '\t'
          //  << y << '\t'
          //  << z << '\t'
          //  << r << '\t'
          //  << px_step << '\t'
          //  << py_step << '\t'
          //  << pz_step << '\t'
          //  //<< x_adj << '\t'
          //  << px_adj << '\t'
          //  << pz_adj << '\t'
          //  << drds << '\n';
          //sim_data<<setw(10)<<evtNb
//        //      <<setw(10)<<track_id<<setw(10)<<prev_track_id
          //    <<setw(15)<<name<<setw(10)<<pname<<setw(14)<<x<<setw(14)<<y<<setw(14)<<z
          //    <<setw(14)<<x_vert<<setw(14)<<y_vert<<setw(14)<<z_vert
          //    <<setw(14)<<px_step<<setw(14)<<py_step<<setw(14)<<pz_step
          //    <<setw(14)<<step->GetPreStepPoint()->GetMaterial()->GetName()<<endl;
          //sim_data.close();
        }
      }
      prev_track_id=track_id;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
