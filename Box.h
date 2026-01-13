#pragma once
#include <utility>
#include <array>
#include <random>
#include <chrono>
#include <cmath>
#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

class Box {
public:
	static constexpr int ATOM_COUNT = 1000;
	static constexpr int DIMENSION = 3;
	static constexpr float TEMPERATURE = 94.4f; // K
	static constexpr double KB = 1.380649 * 1e-23; // J/K or (m^2*kg)/(s^2*K)
	static constexpr double RHO = 1.374 * 1e-3; // kg/cm^3
	static constexpr double SIGMA = 3.4f; // Angstroms
	static constexpr double EPSILSON = 120 * Box::KB; // J or (m^2*kg)/s^2
	static constexpr double EPSILON_ANGSTROMS = Box::EPSILSON * 1e20; // (A^2*kg)/s^2
	static constexpr double CUTOFF_RANGE = 2.25 * Box::SIGMA; // Angstroms
	static constexpr double ATOM_MASS = 39.95 * 1.6747e-24 * 1e-3; // kg
	static const double BOX_LENGTH; // Angstroms
protected:
	array<string, Box::ATOM_COUNT> names;
	Matrix<double> coordinates;
	Matrix<double> velocities;
	Matrix<double> potentials;
	Matrix<double> forces;
public:
	Box();
	~Box();

	void InitializeNames();
	void InitializeCoordinates();
	void InitializeVelocities();

	void ReflectBoundaryCondition();

	void DisplayCoordinates() const;
	void DisplayVelocities() const;
	void DisplayPotentials() const;

	void UpdatePotentials();
	void UpdatePairPotential(int idx1, int idx2);
	void ComputeForces();
	void ComputeForcesCPU();
	void Integrate(double dt);
	void Equilibrate(int steps, double dt, string dcd_file, int save_freq);

	void WritePDB(string file_name) const;
	void WriteDCDHeader(string file_name, int total_frames, double dt);
	void AppendDCDFrame(string file_name);
};