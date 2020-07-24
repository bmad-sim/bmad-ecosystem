#pragma once
#include "G4UserEventAction.hh"
#include "globals.hh"

// This class is currently not necessary; it could be removed
// However, I have left it here in case it is needed for some sort of
// event processing in the future

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
};

