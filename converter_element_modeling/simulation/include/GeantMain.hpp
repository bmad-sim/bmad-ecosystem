#pragma once
#include <string>
#include "binner.hpp"
#include "cal_binner.hpp"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
G4RunManager* Initialize_Geant(void);
//int run_simulation(G4RunManager *runManager, const std::string& target_material, double in_energy, double target_thickness, Binner* binner, int run_length);
int run_simulation(G4RunManager*, const std::string& , double, double, BinnerBase*, int);
//int run_simulation(G4RunManager*, const std::string& , double, double, CalibrationBinner*, int);
//std::pair<int, int> calibrate_binner(G4RunManager *runManager, const std::string& target_material, double in_energy, double target_thickness, Binner* binner);
