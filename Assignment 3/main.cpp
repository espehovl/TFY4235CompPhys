#include <iostream>

#include "armadillo"
#include "utilities.h"

using namespace std;

void createAndFindEigensnacks();
void hamiltonianWithBarrier();
int main() {
	
	//createAndFindEigensnacks();

	hamiltonianWithBarrier();
	
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
	
	const double v0 = 1000;

	arma::mat m = setUpWithBarrier(disc, discretization, deltaX, v0);
	//m.print();

	//Solve the problem for eigenvalues and -vectors
	arma::vec eigenvalues;
	arma::mat eigenvectors;
	arma::eig_sym(eigenvalues, eigenvectors, m);

	//eigenvectors.print("Eigenvectors");
	//potentialProfile(discretization);

	saveEigenData(eigenvalues, eigenvectors, discretization, v0);
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