#include <iostream>
#include <chrono>

#include "armadillo"
#include "utilities.h"

using namespace std;

void createAndFindEigensnacks();
void hamiltonianWithBarrier();

int main() {
	
	//createAndFindEigensnacks();

	hamiltonianWithBarrier();
	
	//Automagically run the python script for plotting (32-bit only, also does not work with creating GIFs)
	system(std::string{ "python main.py" }.c_str());
	cout << "Program finished, enter any character... ";
	cin.get();
	return 0;
}

void hamiltonianWithBarrier() {
	//Initialize discretization environment
	const int disc = 3000;
	const double deltaX = 1.0 / (double)(disc - 1);
	const arma::vec discretization = arma::linspace(0, 1, disc);
	
	const double v0 = 0; //Potential barrier height

	
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

	//Task 3.7
	if (true) {
		arma::vec Psi_0 = eigenvectors.col(0); //Initial condition
		arma::vec PSI; //The result, to be used later on

	}

}

//Simply for bundling this stuff together
void createAndFindEigensnacks() {
	
	const int disc = 900; //Discretization of the space
	const double deltaX = 1.0 / (double)(disc -1); //Spacing of the space :P
	const arma::vec discretization = arma::linspace(0, 1, disc);
	const double v0 = 0;
	arma::mat m = initCentralDiff2(disc, deltaX);
	//discretization.print();
	
	//Create the central difference matrix, including the dX factor
	

	//Find eigenvalues and store them in vector
	arma::vec eigenvalues;
	//Find eigenvectors and store them as column vectors in a matrix
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, m);

	//eigenvectors.print();
	//m.print();

	//eigenvalues.print("Eigenvalues: ");
	//eigenvectors.print("Eigenvectors: ");
	saveEigenData(eigenvalues, eigenvectors, discretization, v0);
}