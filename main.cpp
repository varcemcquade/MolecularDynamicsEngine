#include "Matrix.h"
#include "Box.h"
#include <filesystem>
#include <iostream>

using namespace std;

int main() {
	cout << "Working directory: " << filesystem::current_path() << endl;
	cout << "Box size length is: " << Box::BOX_LENGTH << endl;

	Box* box = new Box();

	double dt = 1e-15;
	int steps = 10000;
	int nsavc = 100;

	box->DisplayCoordinates();
	box->Equilibrate(steps, dt, "trajectory.dcd", nsavc);


	box->DisplayCoordinates();
	box->WritePDB("final.pdb");

	cout << "DCD file size: " << filesystem::file_size("trajectory.dcd") << " bytes" << endl;
	cout << "PDB file size: " << filesystem::file_size("final.pdb") << " bytes" << endl;

	delete box;
}