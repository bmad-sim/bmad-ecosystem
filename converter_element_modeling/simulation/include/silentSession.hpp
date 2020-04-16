#pragma once
#include "G4UImanager.hh"
#include "G4UIsession.hh"

class silentSession : public G4UIsession {
  public:
    G4int ReceiveG4cout(const G4String&) {return 0;}
    G4int ReceiveG4cerr(const G4String&) {return 0;}
};
