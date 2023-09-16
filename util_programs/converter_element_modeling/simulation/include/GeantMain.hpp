#pragma once
#include <string>
#include "point_cache.hpp"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
G4RunManager* Initialize_Geant(unsigned num_threads);
int run_simulation(G4RunManager*, const std::string& , double, double, G4ThreeVector, PointCache*, int);
