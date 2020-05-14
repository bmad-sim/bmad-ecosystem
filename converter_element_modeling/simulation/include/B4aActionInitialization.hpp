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
// $Id: B4aActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B4aActionInitialization.hh
/// \brief Definition of the B4aActionInitialization class

#ifndef B4aActionInitialization_h
#define B4aActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "point_cache.hpp"

class B4DetectorConstruction;

/// Action initialization class.
///

class B4aActionInitialization : public G4VUserActionInitialization
{
  public:
    B4aActionInitialization(B4DetectorConstruction*, double in_particle_energy, PointCache* pc);
    B4aActionInitialization(const B4aActionInitialization&) = default;
    B4aActionInitialization& operator=(const B4aActionInitialization&) = default;
    virtual ~B4aActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    B4DetectorConstruction* fDetConstruction;
    double m_in_particle_energy;
    PointCache* m_point_cache;
};

#endif


