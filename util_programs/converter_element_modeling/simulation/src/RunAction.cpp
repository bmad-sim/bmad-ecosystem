#include "RunAction.hpp"
#include "SteppingAction.hpp"

#include "G4Run.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


RunAction::RunAction() : G4UserRunAction() {}
RunAction::~RunAction() {}


void RunAction::BeginOfRunAction(const G4Run* /*run*/) {
#ifndef G4MULTITHREADED
  static_cast<const SteppingAction*>(G4RunManager::GetRunManager()->GetUserSteppingAction())->GetVec();
#endif
}


void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
#ifndef G4MULTITHREADED
  static_cast<const SteppingAction*>(G4RunManager::GetRunManager()->GetUserSteppingAction())->ReturnVec();
#endif
}

