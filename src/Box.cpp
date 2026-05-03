#include "Box.h"
#include <iostream>
#include <fstream>
#include <iomanip>

const double Box::BOX_LENGTH = cbrt(Box::ATOM_COUNT * Box::ATOM_MASS / Box::RHO) * 1e8;
const bool Box::DEFAULT_IS_PERIODIC = true;

extern "C" void UpdateForcesGPU(double* device_coordinates, double* device_forces,
								double* host_coordinates, double* host_forces,
								const int ATOM_COUNT,
								const double EPSILON_ANGSTROMS, const double SIGMA,
								const double BOX_LENGTH, const double CUTOFF_RANGE,
								const bool isPeriodic);

extern "C" void AllocateGPUBuffers(double** device_coordinates, double** device_forces,
									const int ATOM_COUNT);

extern "C" void FreeGPUBuffers(double* device_coordinates, double* device_forces);

Box::Box() : isPeriodic(true), coordinates{}, velocities{}, forces{} {}

Box::Box(const bool isPeriodic) : isPeriodic(isPeriodic), coordinates{}, velocities{}, forces{} {
	InitializeNames();
	InitializeCoordinates();
	InitializeVelocities();

	initial_coordinates = coordinates;

	if (Box::ATOM_COUNT > 999) {
		AllocateGPUBuffers(&device_coordinates, &device_forces, Box::ATOM_COUNT);
	}

	UpdateForces();
}

Box::~Box() {
	if (device_coordinates != nullptr) { FreeGPUBuffers(device_coordinates, device_forces); }
}

void Box::InitializeNames() {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		names[i] = "Ar";
	}
}

void Box::InitializeCoordinates() {
	int atoms_per_side = (int)ceil(cbrt(ATOM_COUNT));
	double spacing = Box::BOX_LENGTH / atoms_per_side;

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		int x = i % atoms_per_side;
		int y = (i / atoms_per_side) % atoms_per_side;
		int z = i / (atoms_per_side * atoms_per_side);
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
			velocities[3 * i + d] = dist(generator) * 1e10;
		}
	}

	for (int d = 0; d < Box::DIMENSION; ++d) {
		double sum = 0.0;
		for (int i = 0; i < Box::ATOM_COUNT; ++i) { sum += velocities[3 * i + d]; }
		double avg_vel = sum / Box::ATOM_COUNT;
		for (int i = 0; i < Box::ATOM_COUNT; ++i) { velocities[3 * i + d] -= avg_vel; }
	}
}

double Box::ComputePotential(const double r) const {
	return (4 * Box::EPSILON_ANGSTROMS * (pow((Box::SIGMA / r), 12) - pow((Box::SIGMA / r), 6)));
}

double Box::ComputeForce(const double r, const double dim_diff) const {
	double r2   = r * r;
	double r6   = r2 * r2 * r2;
	double r8   = r6 * r2;
	double r14  = r8 * r6;
	double s6   = pow(Box::SIGMA, 6);
	double s12  = s6 * s6;
	double u    = ((6 * s6) / r8) - ((12 * s12) / r14);
	return -4 * Box::EPSILON_ANGSTROMS * u * dim_diff;
}

void Box::UpdateCoordinates() {
	if (!isPeriodic) {
		for (int i = 0; i < Box::ATOM_COUNT; ++i) {
			for (int d = 0; d < Box::DIMENSION; ++d) {
				if (coordinates[3 * i + d] < 0) { coordinates[3 * i + d] = -coordinates[3 * i + d]; }
				else if (coordinates[3 * i + d] > Box::BOX_LENGTH) { coordinates[3 * i + d] = 2 * Box::BOX_LENGTH - coordinates[3 * i + d]; }
				velocities[3 * i + d] *= -1;
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

void Box::UpdateForces() {
	if (Box::ATOM_COUNT > 999) {
		UpdateForcesGPU(device_coordinates, device_forces,
						coordinates.data(), forces.data(),
						Box::ATOM_COUNT,
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
			double a = forces[3 * i + d] / Box::ATOM_MASS;
			coordinates[3 * i + d] += velocities[3 * i + d] * dt + 0.5 * a * dt * dt;
			velocities[3 * i + d]  += 0.5 * a * dt;
		}
	}

	UpdateCoordinates();
	UpdateForces();

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			double a = forces[3 * i + d] / Box::ATOM_MASS;
			velocities[3 * i + d] += 0.5 * a * dt;
		}
	}
}

void Box::Equilibrate(int steps, double dt, const string& dcd_file, int save_freq, const string& log_file, const string& msd_file) {
	WriteDCDHeader(dcd_file, steps / save_freq, dt);

	ofstream log(log_file);
	ofstream msd_out(msd_file);
	log << "step,KE,PE,E_total,Temperature\n";
	msd_out << "step,MSD\n";

	for (int i = 0; i < steps; ++i) {
		Integrate(dt);

		if (i % 1000 == 0) {
			cout << "Step: " << i << " completed.\n";
		}
		if (i % save_freq == 0) {
			AppendDCDFrame(dcd_file);

			double KE = ComputeTotalKineticEnergy();
			double PE = ComputeTotalPotentialEnergy();
			double ETOT = KE + PE;
			double T = ComputeTemperature();
			log << i << "," << KE << "," << PE << "," << ETOT << "," << T << "\n";

			MeanSquaredDisplacement(msd_out, i);
		}
	}

	log.close();
	msd_out.close();
}

void Box::DisplayCoordinates() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) { cout << coordinates[3 * i + d] << " "; }
		cout << endl;
	}
}

void Box::DisplayVelocities() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) { cout << velocities[3 * i + d] << " "; }
		cout << endl;
	}
}

void Box::DisplayForces() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) { cout << forces[3 * i + d] << " "; }
		cout << endl;
	}
}

double Box::ComputeTotalKineticEnergy() const {
	double KE = 0;
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		double vx = velocities[3 * i + 0];
		double vy = velocities[3 * i + 1];
		double vz = velocities[3 * i + 2];
		KE += 0.5 * Box::ATOM_MASS * (vx * vx + vy * vy + vz * vz);
	}
	return KE;
}

double Box::ComputeTotalPotentialEnergy() const {
	double PE = 0.0;
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
			if (r < Box::CUTOFF_RANGE) { PE += ComputePotential(r); }
		}
	}
	return PE;
}

double Box::ComputeTotalEnergy() const {
	return (ComputeTotalKineticEnergy() + ComputeTotalPotentialEnergy());
}

double Box::ComputeTemperature() const {
	return (2.0 * ComputeTotalKineticEnergy() * 1e-20) / (3.0 * Box::ATOM_COUNT * Box::KB);
}

double Box::ComputePressure() const {
	double virial = 0.0;

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
				// r · F_ij = dx*Fx + dy*Fy + dz*Fz
				virial += dx * ComputeForce(r, dx)
						+ dy * ComputeForce(r, dy)
						+ dz * ComputeForce(r, dz);
			}
		}
	}

	double volume_angstroms3 = Box::BOX_LENGTH * Box::BOX_LENGTH * Box::BOX_LENGTH;
	double NkBT = Box::ATOM_COUNT * Box::KB * ComputeTemperature();
	double volume_m3 = volume_angstroms3 * 1e-30;

	return (NkBT / volume_m3) + (virial * 1e-20) / (3.0 * volume_m3); // Convert Angstroms^3 to m^3 for SI pressure (Pa)
}

void Box::MeanSquaredDisplacement(ofstream& out, int step) const {
	double msd = 0.0;
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			double delta = coordinates[3 * i + d] - initial_coordinates[3 * i + d];
			msd += delta * delta;
		}
	}
	msd /= Box::ATOM_COUNT;
	out << step << "," << msd << "\n";
}

void Box::StaticStructureFactor() const {}
void Box::Relation() const {}

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
		fout << setw(8)  << right << fixed << setprecision(3) << coordinates[3 * i + 1];
		fout << setw(8)  << right << fixed << setprecision(3) << coordinates[3 * i + 2];
		fout << setw(6)  << right << fixed << setprecision(2) << 1.00;
		fout << setw(6)  << right << fixed << setprecision(2) << 0.00;
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