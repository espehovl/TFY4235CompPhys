#pragma once

//Initialize a finite difference matrix of second order
arma::mat initCentralDiff2(int discretization, double dX);

//Initialize finite difference matrix of problem with barrier
arma::mat setUpWithBarrier(int discretization, arma::vec x_space, double dX, double v0);

//Returns the prefactor of the potential (i.e. 1 if we are in the barrier, 0 if we are outside)
double v(double x_prime);
//For plotting the potential over the graphs later on
arma::vec potentialProfile(arma::vec x_space);

//Save eigenvectors and values to file, for plotting purposes
//dX corrects for the discretization spacing, such that the eigenvalues are "correct"
void saveEigenData(const arma::vec& vals, const arma::mat& vecs, const arma::vec& discretization, double v0, int howMany = 0, std::string filename = "eigenData.tsv");

