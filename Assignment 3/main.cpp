#include <iostream>
#include <chrono>

#include "armadillo"
#include "utilities.h"

using namespace std;

void hamiltonianWithBarrier();

int main() {
	
	//Setup and solve the system with barrier
	hamiltonianWithBarrier();
	
	cout << "Program finished, enter any character... ";
	cin.get();
	return 0;
}

void hamiltonianWithBarrier() {
	//Initialize discretization environment
	const int disc = 3000;
	const double deltaX = 1.0 / (double)(disc - 1);
	const arma::vec discretization = arma::linspace(0, 1, disc);
	
	const double v0 = 0; //Potential barrier height. 0 gives no-barrier-case, naturally

	
	cout << "Setting up " << disc << "x" << disc << " matrix...";
	arma::mat m = setUpWithBarrier(disc, discretization, deltaX, v0);
	cout << " matrix setup complete!" << endl;
	//m.print();
	cout << "Solving for eigenvalues and eigenvectors..." << endl;
	auto start = chrono::high_resolution_clock::now();
	//Solve the problem for eigenvalues and -vectors
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, m);
	auto stop = chrono::high_resolution_clock::now();
	cout << "System solved in " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms." << endl;

	//eigenvectors.print("Eigenvectors");
	//potentialProfile(discretization);

	saveEigenData(eigenvalues, eigenvectors, discretization, v0);
}
