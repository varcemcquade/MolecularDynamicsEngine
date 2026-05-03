#include "Box.h"
#include <filesystem>
#include <iostream>

using namespace std;

int main() {
	cout << "Working directory: " << filesystem::current_path() << endl;
	cout << "Box length: " << Box::BOX_LENGTH << " Angstroms" << endl;

	bool isPeriodic = true;

	Box* box = new Box(isPeriodic);

	double dt = 1e-15;
	int steps = 100000;
	int nsavc = 100;

	box->WritePDB("initial.pdb");
	box->Equilibrate(steps, dt, "trajectory.dcd", nsavc, "energies.csv", "msd.csv");
	box->WritePDB("final.pdb");

	cout << "DCD file size:     " << filesystem::file_size("trajectory.dcd") << " bytes" << endl;
	cout << "PDB file size:     " << filesystem::file_size("final.pdb")      << " bytes" << endl;
	cout << "Energy log size:   " << filesystem::file_size("energies.csv")   << " bytes" << endl;

	delete box;
	return 0;
}