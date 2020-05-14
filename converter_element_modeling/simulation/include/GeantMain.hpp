#pragma once
#include <string>
#include "point_cache.hpp"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#ifdef G4MULTITHREADED
typedef G4RunManager RunManager_t;
#else
typedef G4RunManager RunManager_t;
#endif
RunManager_t* Initialize_Geant(void);
int run_simulation(RunManager_t*, const std::string& , double, double, PointCache*, int);
