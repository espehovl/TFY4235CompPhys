#include "armadillo"
#include "utilities.h"
#include <fstream>

using std::pow;

//Initialize finite difference matrix of problem with barrier
arma::mat setUpWithBarrier(int discretization, arma::vec x_space, double dX, double v0){
	//Set up matrix
	arma::mat m(discretization, discretization);
	m.fill(0.0);

	//Fill subdiagonals (no magic going on here)
	m.diag(-1).fill(-1.0 / pow(dX, 2));
	m.diag( 1).fill(-1.0 / pow(dX, 2));

	//Fill main diagonal (slightly more magic going on here now)
	for (int i = 0; i < m.n_cols; i++) {
		m(i, i) = 2/pow(dX,2) + v0 * v(x_space(i));
	}
	return m; //We are done
}

//Returns the prefactor of the potential (i.e. 1 if we are in the barrier, 0 if we are outside)
double v(double x_prime) {
	if (1.0 / 3.0 <= x_prime && x_prime < 2.0 / 3.0) {
		return 1.0;
	}
	else return 0.0;
}

//For plotting the potential over the graphs later on
arma::vec potentialProfile(arma::vec x_space){
	arma::vec profile(x_space.n_rows);
	for (int i = 0; i < x_space.n_rows; i++) {
		profile(i) = v(x_space(i));
	}
	return profile;
}

// Save eigenvectors and values to file, for plotting purposes
void saveEigenData(const arma::vec & vals, const arma::mat & vecs, const arma::vec& discretization, double v0, int howMany, std::string filename){
	if (howMany == 0) howMany = (int) vals.n_rows; //Save all eigenvalues and vectors

	arma::vec potential = potentialProfile(discretization);

	std::cout << "Saving data in combo file..." << std::endl;
	std::ofstream ofs{ filename, std::ios::out | std::ios::trunc };
	ofs << "#Discretization\t" << std::endl;
	for (auto it = discretization.begin(); it != discretization.end(); it++) {
		ofs << *it << "\t";
	}
	ofs << "\n#v0" << std::endl
		<< v0 << std::endl;

	ofs << "#Potential Profile" <<std::endl;
	for (auto it = potential.begin(); it != potential.end(); it++) {
		ofs << *it << "\t";
	}

	ofs << "\n#\t" << "Lambda \t" << "Eigenvector -->\t" << std::endl;
	
	for (int i = 0; i < howMany; i++) {
		ofs << vals(i) << "\t";
		for (int j = 0; j < vecs.n_rows; j++) {
			ofs << vecs(j, i) << "\t";
		}
		ofs << std::endl;
	}
	std::cout << "Eigendata saved in combo file!" << std::endl;

	std::cout << "Saving eigenvalues in separate file..." << std::endl;
	vals.save("eigenvalues.txt", arma::raw_ascii);
	ofs.close();

}
//For studying the error for different discretizations
void compareDiscretizations(){

	using namespace std;
	string dir{"results/vectors/"};
	for (int i = 10; i < 10000; i += 100) {
		const int disc = i;
		const double deltaX = 1.0 / (double)(disc - 1);
		const arma::vec discretization = arma::linspace(0, 1, disc);
		cout << "Setting up " << disc << "x" << disc << " matrix...";
		arma::mat m = setUpWithBarrier(disc, discretization, deltaX, 0);
		cout << " matrix setup complete!" << endl;

		ofstream ofs{ dir + "vector" + to_string(i) + ".tsv",std::ios::out | std::ios::trunc };
		for (auto x : discretization) { //Write the x-space to file
			ofs << x << "\t";
		}
		ofs << endl;

		cout << "Solving for eigenvalues and eigenvectors..." << endl;
		auto start = chrono::high_resolution_clock::now();
		//Solve the problem for eigenvalues and -vectors
		arma::vec eigenvalues;
		arma::mat eigenvectors;
		arma::eig_sym(eigenvalues, eigenvectors, m);
		auto stop = chrono::high_resolution_clock::now();
		cout << "System solved in " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms." << endl;

		for (int col = 0; col < 10; col++) { //For 10 eigenvectors...
			for (int row = 0; row < eigenvectors.n_rows; row++) { //For each element in the vector...
				ofs << eigenvectors(row, col) << "\t";
			}
			ofs << endl;
		}

	}

}


