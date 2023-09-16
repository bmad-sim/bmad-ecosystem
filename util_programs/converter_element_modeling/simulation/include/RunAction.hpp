#pragma once
#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

// This class ensures that in single threaded mode, the
// stepping action still acquires and returns a vector
// from the point cache
class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);
};


