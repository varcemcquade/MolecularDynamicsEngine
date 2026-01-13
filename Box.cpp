#include "Box.h"

const double Box::BOX_LENGTH = cbrt(Box::ATOM_COUNT * Box::ATOM_MASS / Box::RHO) * 1e8;

Box::Box() :
	coordinates(Matrix<double>(Box::ATOM_COUNT, Box::DIMENSION)), 
	velocities(Matrix<double>(Box::ATOM_COUNT, Box::DIMENSION)), 
	potentials(Matrix<double>(Box::ATOM_COUNT, Box::ATOM_COUNT)), 
	forces(Matrix<double>(Box::ATOM_COUNT, Box::DIMENSION)) {
	InitializeNames();
	InitializeCoordinates();
	InitializeVelocities();
	UpdatePotentials();
	ComputeForces();
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

	int idx = 0;
	for (int x = 0; x < atoms_per_side && idx < ATOM_COUNT; ++x) {
		for (int y = 0; y < atoms_per_side && idx < ATOM_COUNT; ++y) {
			for (int z = 0; z < atoms_per_side && idx < ATOM_COUNT; ++z) {
				coordinates(idx, 0) = (x + 0.5) * spacing;
				coordinates(idx, 1) = (y + 0.5) * spacing;
				coordinates(idx, 2) = (z + 0.5) * spacing;
				++idx;
			}
		}
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
			velocities(i, d) = dist(generator) * 1e10;
		}
	}
}

void Box::ReflectBoundaryCondition() {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			if (coordinates(i, d) < 0) {
				coordinates(i, d) = -coordinates(i, d);
				velocities(i, d) *= -1;
			}
			else if (coordinates(i, d) > Box::BOX_LENGTH) {
				coordinates(i, d) = 2 * Box::BOX_LENGTH - coordinates(i, d);
				velocities(i, d) *= -1;
			}
		}
	}
}

void Box::DisplayCoordinates() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int j = 0; j < Box::DIMENSION; ++j) {
			cout << coordinates(i, j) << " ";
		}
		cout << endl;
	}
}

void Box::DisplayVelocities() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int j = 0; j < Box::DIMENSION; ++j) {
			cout << velocities(i, j) << " ";
		}
		cout << endl;
	}
}

void Box::DisplayPotentials() const {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int j = 0; j < Box::ATOM_COUNT; ++j) {
			cout << potentials(i, j) << " ";
		}
		cout << endl;
	}
}

void Box::UpdatePotentials() {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int j = i + 1; j < Box::ATOM_COUNT; ++j) {
			UpdatePairPotential(i, j);
		}
	}
}

void Box::UpdatePairPotential(int idx1, int idx2) {
	double x1 = coordinates(idx1, 0);
	double y1 = coordinates(idx1, 1);
	double z1 = coordinates(idx1, 2);

	double x2 = coordinates(idx2, 0);
	double y2 = coordinates(idx2, 1);
	double z2 = coordinates(idx2, 2);

	double r = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2)); // Angstroms
	potentials(idx1, idx2) = 4 * Box::EPSILON_ANGSTROMS * (pow((Box::SIGMA / r), 12) - pow((Box::SIGMA / r), 6)); // (kg*A^2)/s^2
}

void Box::ComputeForcesCPU() {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		forces(i, 0) = 0;
		forces(i, 1) = 0;
		forces(i, 2) = 0;
	}

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		double xi = coordinates(i, 0);
		double yi = coordinates(i, 1);
		double zi = coordinates(i, 2);

		for (int j = i + 1; j < Box::ATOM_COUNT; ++j) {
			double dx = xi - coordinates(j, 0);
			double dy = yi - coordinates(j, 1);
			double dz = zi - coordinates(j, 2);

			double r2 = dx * dx + dy * dy + dz * dz;
			double r = sqrt(r2);

			double u = ((6 * pow(Box::SIGMA, 6)) / pow(r, 8)) - ((12 * pow(Box::SIGMA, 12)) / pow(r, 14)); // Angstroms^-2
			double scalar_f = -4 * Box::EPSILON_ANGSTROMS * u;

			forces(i, 0) += scalar_f * dx;
			forces(i, 1) += scalar_f * dy;
			forces(i, 2) += scalar_f * dz;

			forces(j, 0) += -scalar_f * dx;
			forces(j, 1) += -scalar_f * dy;
			forces(j, 2) += -scalar_f * dz;
		}
	}
}

extern "C" void ComputeForcesGPU(double* coordinates, double* forces, int n, double epsilon, double sigma);

void Box::ComputeForces() {
	if (Box::ATOM_COUNT > 999) {
		ComputeForcesGPU(coordinates.data_ptr(), forces.data_ptr(), Box::ATOM_COUNT, Box::EPSILON_ANGSTROMS, Box::SIGMA);
	}
	else {
		cout << "Initializing CPU.\n";
		ComputeForcesCPU();
	}
}

void Box::Integrate(double dt) {
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			double a = forces(i, d) / Box::ATOM_MASS; // Angstroms/s^2
			coordinates(i, d) = coordinates(i, d) + velocities(i, d) * dt + 0.5 * a * dt * dt;
			velocities(i, d) = velocities(i, d) + 0.5 * a * dt;
		}
	}

	ReflectBoundaryCondition();
	ComputeForces();

	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		for (int d = 0; d < Box::DIMENSION; ++d) {
			double a = forces(i, d) / Box::ATOM_MASS; // Angstroms/s^2
			velocities(i, d) = velocities(i, d) + 0.5 * a * dt;
		}
	}
}

void Box::Equilibrate(int steps, double dt, string dcd_file, int save_freq) {
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

void Box::WritePDB(string file_name) const {
	ofstream fout(file_name);
	for (int i = 0; i < Box::ATOM_COUNT; ++i) {
		fout << "ATOM";
		fout << setw(7) << right << (i + 1);
		fout << "  ";
		fout << setw(4) << left << "Ar";
		fout << "HET ";
		fout << setw(5) << right << (i + 1);
		fout << setw(12) << right << fixed << setprecision(3) << coordinates(i, 0);
		fout << setw(8) << right << fixed << setprecision(3) << coordinates(i, 1);
		fout << setw(8) << right << fixed << setprecision(3) << coordinates(i, 2);
		fout << setw(6) << right << fixed << setprecision(2) << 1.00;
		fout << setw(6) << right << fixed << setprecision(2) << 0.00;
		fout << setw(12) << right << "HETA";
		fout << endl;
	}
	fout << "END" << endl;
	fout.close();
}

void Box::WriteDCDHeader(string file_name, int total_frames, double dt) {
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

void Box::AppendDCDFrame(string file_name) {
	ofstream fout(file_name, ios::binary | ios::app);

	int block_size = ATOM_COUNT * 4;

	fout.write((char*)&block_size, 4);
	for (int i = 0; i < ATOM_COUNT; ++i) {
		float x = (float)coordinates(i, 0);
		fout.write((char*)&x, 4);
	}
	fout.write((char*)&block_size, 4);

	fout.write((char*)&block_size, 4);
	for (int i = 0; i < ATOM_COUNT; ++i) {
		float y = (float)coordinates(i, 1);
		fout.write((char*)&y, 4);
	}
	fout.write((char*)&block_size, 4);

	fout.write((char*)&block_size, 4);
	for (int i = 0; i < ATOM_COUNT; ++i) {
		float z = (float)coordinates(i, 2);
		fout.write((char*)&z, 4);
	}
	fout.write((char*)&block_size, 4);

	fout.close();
}