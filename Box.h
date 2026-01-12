#pragma once
#include <utility>
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
	static const int ATOM_COUNT;
	static const int DIMENSION;
	static const float TEMPERATURE;
	static const double KB;
	static const double RHO;
	static const double SIGMA;
	static const double EPSILSON;
	static const double EPSILON_ANGSTROMS;
	static const double CUTOFF_RANGE;
	static const double ATOM_MASS;
	static const double BOX_LENGTH;
protected:
	Matrix<string>* names;
	Matrix<double>* coordinates;
	Matrix<double>* velocities;
	Matrix<double>* potentials;
	Matrix<double>* forces;
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