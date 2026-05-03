#pragma once
#include <array>
#include <fstream>
#include <random>
#include <chrono>
#include <cmath>

using namespace std;

class Box {
public:
	static constexpr int ATOM_COUNT = 1000;
	static constexpr int DIMENSION = 3;
	static constexpr float TEMPERATURE = 94.4f; // K
	static constexpr double KB = 1.380649 * 1e-23; // J/K or (m^2*kg)/(s^2*K)
	static constexpr double RHO = 1.374 * 1e-3; // kg/cm^3
	static constexpr double SIGMA = 3.4f; // Angstroms
	static constexpr double EPSILON = 120 * Box::KB; // J or (m^2*kg)/s^2
	static constexpr double EPSILON_ANGSTROMS = Box::EPSILON * 1e20; // (A^2*kg)/s^2
	static constexpr double CUTOFF_RANGE = 2.25 * Box::SIGMA; // Angstroms
	static constexpr double ATOM_MASS = 39.95 * 1.6747e-24 * 1e-3; // kg
	static const double BOX_LENGTH; // Angstroms
	static const bool DEFAULT_IS_PERIODIC;
protected:
	bool isPeriodic;
	array<string, Box::ATOM_COUNT> names;
	array<double, 3 * Box::ATOM_COUNT> coordinates;
	array<double, 3 * Box::ATOM_COUNT> initial_coordinates;
	array<double, 3 * Box::ATOM_COUNT> velocities;
	array<double, 3 * Box::ATOM_COUNT> forces;

	double* device_coordinates = nullptr;
	double* device_forces = nullptr;
public:
	Box();
	Box(const bool isPeriodic);
	~Box();

	void InitializeNames();
	void InitializeCoordinates();
	void InitializePotentials();
	void InitializeVelocities();

	double ComputePotential(const double r) const;
	double ComputeForce(const double r, const double dim_diff) const;

	void UpdateCoordinates();
	void UpdateForces();
	void UpdateForcesCPU();
	void Integrate(const double dt);
	void Equilibrate(int steps, double dt, const string& dcd_file, int save_freq, const string& log_file, const string& msd_file);

	void DisplayCoordinates() const;
	void DisplayVelocities() const;
	void DisplayForces() const;

	double ComputeTotalKineticEnergy() const;
	double ComputeTotalPotentialEnergy() const;
	double ComputeTotalEnergy() const;
	double ComputeTemperature() const;
	double ComputePressure() const;
	void StaticStructureFactor() const;
	void Relation() const;
	void MeanSquaredDisplacement(ofstream& out, int step) const;

	void WritePDB(const string& file_name) const;
	void WriteDCDHeader(const string& file_name, int total_frames, double dt);
	void AppendDCDFrame(const string &file_name);
};
