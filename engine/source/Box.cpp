#include "Box.h"
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numbers>

const double Box::BOX_LENGTH = cbrt(Box::ATOM_COUNT * Box::ATOM_MASS / Box::RHO) * 1e8;
const bool Box::DEFAULT_IS_PERIODIC = true;

Box::Box() : isPeriodic(true), coordinates{}, velocities{}, forces {} {}

Box::Box(const bool isPeriodic) : isPeriodic(isPeriodic), coordinates{}, velocities{}, forces{} {
	InitializeNames();
	InitializeCoordinates();
	InitializeVelocities();
	UpdateForces();
}

Box::~Box() {}

void Box::InitializeNames() {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		names[i] = "Ar";
	}
}

void Box::InitializeCoordinates() {
	int atoms_per_side = (int)ceil(cbrt(ATOM_COUNT));
	double spacing = Box::BOX_LENGTH / atoms_per_side;

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		int x = i % atoms_per_side; // cycles 0 to (n-1) for i
		int y = (i / atoms_per_side) % atoms_per_side; // cycles 0 to (n-1) each time x zeroes
		int z = i / (atoms_per_side * atoms_per_side); // cycles 0 to (n-1) each time y zeroes
		coordinates[3 * i + 0] = (x + 0.5) * spacing;
		coordinates[3 * i + 1] = (y + 0.5) * spacing;
		coordinates[3 * i + 2] = (z + 0.5) * spacing;
	}
}

void Box::InitializeVelocities() {
	unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	double mean = 0.0;
	double stdv = sqrt(Box::KB * Box::TEMPERATURE / Box::ATOM_MASS);
	normal_distribution<double> dist(mean, stdv);

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			velocities[3 * i + d]= dist(generator) * 1e10;
		}
	}
}

void Box::CheckCoordinates() {
	if (!isPeriodic) {
		for (int i = 0; i < Box::ATOM_COUNT; ++i) {
			for (int d = 0; d < Box::DIMENSION; ++d) {
				if (coordinates[3 * i + d] < 0) {
					coordinates[3 * i + d] = -coordinates[3 * i + d];
					velocities[3 * i + d] *= -1;
				}
				else if (coordinates[3 * i + d] > Box::BOX_LENGTH) {
					coordinates[3 * i + d] = 2 * Box::BOX_LENGTH - coordinates[3 * i + d];
					velocities[3 * i + d] *= -1;
				}
			}
		}
	}
	if (isPeriodic) {
		for (int i = 0; i < Box::ATOM_COUNT; ++i) {
			for (int d = 0; d < Box::DIMENSION; ++d) {
				if (coordinates[3 * i + d] < 0) { coordinates[3 * i + d] += Box::BOX_LENGTH; }
				if (coordinates[3 * i + d] > Box::BOX_LENGTH) { coordinates[3 * i + d] -= Box::BOX_LENGTH; }
			}
		}
	}
}

void Box::DisplayCoordinates() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			cout << coordinates[3 * i + d] << " ";
		}
		cout << endl;
	}
}

void Box::DisplayVelocities() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			cout << velocities[3 * i + d] << " ";
		}
		cout << endl;
	}
}

void Box::DisplayForces() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			cout << forces[3 * i + d] << " ";
		}
		cout << endl;
	}
}

double Box::ComputeForce(const double r, const double dim_diff) const {
	double u = ((6 * pow(Box::SIGMA, 6)) / pow(r, 8)) - ((12 * pow(Box::SIGMA, 12)) / pow(r, 14)); // Angstroms^-2
	double scalar_f = -4 * Box::EPSILON_ANGSTROMS * u;
	return scalar_f * dim_diff;
}

void Box::UpdateForcesCPU() {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		forces[3 * i + 0] = 0;
		forces[3 * i + 1] = 0;
		forces[3 * i + 2] = 0;
	}

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		double xi = coordinates[3 * i + 0];
		double yi = coordinates[3 * i + 1];
		double zi = coordinates[3 * i + 2];

		for (int j = i + 1; j < Box::ATOM_COUNT; ++j) {
			double dx = xi - coordinates[3 * j + 0];
			double dy = yi - coordinates[3 * j + 1];
			double dz = zi - coordinates[3 * j + 2];

			if (isPeriodic) {
				dx -= Box::BOX_LENGTH * round(dx / Box::BOX_LENGTH);
				dy -= Box::BOX_LENGTH * round(dy / Box::BOX_LENGTH);
				dz -= Box::BOX_LENGTH * round(dz / Box::BOX_LENGTH);
			}

			double r = sqrt(dx * dx + dy * dy + dz * dz);

			if (r < Box::CUTOFF_RANGE) {
				double fx = ComputeForce(r, dx);
				double fy = ComputeForce(r, dy);
				double fz = ComputeForce(r, dz);

				forces[3 * i + 0] += fx;
				forces[3 * i + 1] += fy;
				forces[3 * i + 2] += fz;

				forces[3 * j + 0] -= fx;
				forces[3 * j + 1] -= fy;
				forces[3 * j + 2] -= fz;
			}
		}
	}
}

extern "C" void UpdateForcesGPU(double* coordinates, double* forces, const int ATOM_COUNT,
								const double EPSILON_ANGSTROMS, const double SIGMA,
								const double BOX_LENGTH, const double CUTOFF_RANGE,
								const bool isPeriodic);

void Box::UpdateForces() {
	if (Box::ATOM_COUNT > 999) {
		UpdateForcesGPU(coordinates.data(), forces.data(), Box::ATOM_COUNT,
						Box::EPSILON_ANGSTROMS, Box::SIGMA,
						Box::BOX_LENGTH, Box::CUTOFF_RANGE,
						isPeriodic);
		return;
	}
	UpdateForcesCPU();
}

void Box::Integrate(const double dt) {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			double a = forces[3 * i + d] / Box::ATOM_MASS; // Angstroms/s^2
			coordinates[3 * i + d] = coordinates[3 * i + d] + velocities[3 * i + d] * dt + 0.5 * a * dt * dt;
			velocities[3 * i + d] = velocities[3 * i + d] + 0.5 * a * dt;
		}
	}

	CheckCoordinates();
	UpdateForces();

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			double a = forces[3 * i + d] / Box::ATOM_MASS; // Angstroms/s^2
			velocities[3 * i + d] = velocities[3 * i + d] + 0.5 * a * dt;
		}
	}
}

void Box::Equilibrate(int steps, double dt, const string& dcd_file, int save_freq) {
	WriteDCDHeader(dcd_file, steps / save_freq, dt);

	for (int i = 0; i < steps; ++i) {
		Integrate(dt);
		if (i % 1000 == 0) {
			cout << "Step: " << i << " completed.\n";
		}
		if (i % save_freq == 0) {
			AppendDCDFrame(dcd_file);
		}
	}
}

void Box::WritePDB(const string& file_name) const {
	ofstream fout(file_name);

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		fout << "ATOM";
		fout << setw(7) << right << (i + 1);
		fout << "  ";
		fout << setw(4) << left << "Ar";
		fout << "HET ";
		fout << setw(5) << right << (i + 1);
		fout << setw(12) << right << fixed << setprecision(3) << coordinates[3 * i + 0];
		fout << setw(8) << right << fixed << setprecision(3) << coordinates[3 * i + 1];
		fout << setw(8) << right << fixed << setprecision(3) << coordinates[3 * i + 2];
		fout << setw(6) << right << fixed << setprecision(2) << 1.00;
		fout << setw(6) << right << fixed << setprecision(2) << 0.00;
		fout << setw(12) << right << "HETA";
		fout << endl;
	}
	fout << "END" << endl;
	fout.close();
}

void Box::WriteDCDHeader(const string& file_name, const int total_frames, const double dt) {
	ofstream fout(file_name, ios::binary);

	int block_size = 84;
	fout.write((char*)&block_size, 4);

	char cord[4] = { 'C', 'O', 'R', 'D' };
	fout.write(cord, 4);

	int nframes = total_frames;
	fout.write((char*)&nframes, 4);

	int istart = 0;
	fout.write((char*)&istart, 4);

	int nsavc = 1;
	fout.write((char*)&nsavc, 4);

	int nsteps = total_frames;
	fout.write((char*)&nsteps, 4);

	int padding[5] = { 0, 0, 0, 0, 0 };
	fout.write((char*)padding, 20);

	float timestep = (float)(dt * 1e15 / 48.88821);
	fout.write((char*)&timestep, 4);

	int padding2[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	fout.write((char*)padding2, 36);

	int charmm_version = 24;
	fout.write((char*)&charmm_version, 4);

	fout.write((char*)&block_size, 4);

	int title_block_size = 84;
	fout.write((char*)&title_block_size, 4);

	int ntitle = 1;
	fout.write((char*)&ntitle, 4);

	char title[80] = { 0 };
	strncpy(title, "Created by MD simulation", 80);
	fout.write(title, 80);

	fout.write((char*)&title_block_size, 4);

	block_size = 4;
	fout.write((char*)&block_size, 4);

	int natoms = ATOM_COUNT;
	fout.write((char*)&natoms, 4);

	fout.write((char*)&block_size, 4);

	fout.close();
}

void Box::AppendDCDFrame(const string& file_name) {
	ofstream fout(file_name, ios::binary | ios::app);

	int block_size = ATOM_COUNT * 4;

	fout.write((char*)&block_size, 4);
	for (int i = 0; i < ATOM_COUNT; ++i) {
		float x = (float)coordinates[3 * i + 0];
		fout.write((char*)&x, 4);
	}
	fout.write((char*)&block_size, 4);

	fout.write((char*)&block_size, 4);
	for (int i = 0; i < ATOM_COUNT; ++i) {
		float y = (float)coordinates[3 * i + 1];
		fout.write((char*)&y, 4);
	}
	fout.write((char*)&block_size, 4);

	fout.write((char*)&block_size, 4);
	for (int i = 0; i < ATOM_COUNT; ++i) {
		float z = (float)coordinates[3 * i + 2];
		fout.write((char*)&z, 4);
	}
	fout.write((char*)&block_size, 4);

	fout.close();
}
