
#include "armadillo"
#include "utilities.h"
#include <assert.h>

using std::pow;

arma::mat initCentralDiff2(int discretization, double dX){
	//Set up matrix
	arma::mat m(discretization, discretization);
	m.fill(0.0);

	//Fill diagonals (negative of a "normal" 2nd order central diff scheme)
	m.diag(0).fill(2.0 / pow(dX, 2));
	m.diag(-1).fill(-1.0 / pow(dX, 2)); //Values are divided by dX^2 here
	m.diag(1).fill(-1.0 / pow(dX, 2));

	return m; //We are done
}

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

double v(double x_prime) {
	if (1.0 / 3.0 <= x_prime && x_prime < 2.0 / 3.0) {
		return 1.0;
	}
	else return 0.0;
}

arma::vec potentialProfile(arma::vec x_space){
	arma::vec profile(x_space.n_rows);
	for (int i = 0; i < x_space.n_rows; i++) {
		profile(i) = v(x_space(i));
	}
	return profile;
}

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


