#include "Source.h"
#include <random>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <iomanip>

using namespace std;

double gauss(double stdDev, double mean, default_random_engine& gen) {
	normal_distribution<double> dist(mean, stdDev);
	return dist(gen); //return random number
}
void testGauss() {
	ofstream out{ "normaldistribution.txt", ios::out | ios::trunc };

	map<double, int> stats; // A map is used as the histogram, which is later plotted in Python
	for (unsigned int i = 0; i < 1000000; i++) {
		double number = gauss();
		if (stats.find(number) != stats.end()) stats[number]++;
		else stats[number] = 1;
	}
	for (auto m : stats) { //Write histogram to file
		out << m.first << " " << m.second << endl;
	}
	out.close();
}


double flasher(double t_hat, double tau, double omega) {
	
	t_hat = fmod(t_hat, tau*omega); //Deperiodicize (float modulo function)

	if (t_hat >= 0 && t_hat < 3.0 * tau / 4.0*omega) return 0.0;
	else if (t_hat >= 3.0*tau / 4.0*omega && t_hat < tau*omega) return 1.0;
}

//Reduced units force
double dU_hat(double x_hat, double alpha, double t_hat, double tau, double omega) {
	x_hat = fmod(x_hat, 1.0); //Deperiodicize
	while (x_hat < 0) x_hat += 1.0; //Correct if x is negative (return to the region of interest)

	if (x_hat >= 0 && x_hat < alpha) return 1.0 / alpha * flasher(t_hat, tau, omega);
	else if (x_hat >= alpha && x_hat < 1.0) return - 1.0 / (1.0 - alpha) * flasher(t_hat, tau, omega);
	else {
		cout << "Fuckup in dU_hat!" << endl; //something went wrong
		return 0.0;
	}
}

//Reduced units potential, NOTE: x must be reduced when sent into the function!
double U_hat(double x, double L, double alpha, double t_hat, double tau, double omega, bool active) {
	x = fmod(x, 1.0); //Deperiodicize
	while (x < 0) x += 1.0;

	double f = flasher(t_hat, tau, omega);
	if (!active) f = 1.0;

	if (x >= 0 && x < alpha) return (x / alpha) * f;
	else if (x >= alpha && x < 1.0) return ((1.0 - x) / (1.0 - alpha)) * f;
	else {
		cout << "Fuckup in U_hat!" << endl; //Something went wrong
		return 0.0;
	}
}

//Reduced units Euler scheme
double EulerScheme(double x_hat, double tn, double dt, double D_hat, double alpha, double tau, double omega) {
	//Quantities are separated for debugging purposes only
	double F = dU_hat(x_hat, alpha, tn, tau, omega)*dt;
	double stoch = sqrt(2 * D_hat*dt)*gauss();
	x_hat = x_hat - F + stoch;
	return x_hat;
}

//Returns max of force, given alpha
double dU_MAX(double alpha) {
	if (alpha <= 0.5) return 1 / alpha;
	else return 1 / (1 - alpha);
}

//Finds dt according to criterion
double find_dt(double D_hat, double alpha, double omega) {
	double dt = 1; //Reduced units
	double factor = 10; //How much bigger do you define the mathematical >> as?
	while (dU_MAX(alpha)*dt + 4 * sqrt(2 * D_hat * dt)*factor > alpha) dt /= 10.0;
	return dt; //Reduced units
}

//Obsolete!
void saveMultisysPos(const vector<Simulation*>& sim) {
	cout << "Saving multiple positions..." << endl;
	ofstream out{ "multisys.txt", ios::out | ios::trunc };
	vector<double> T = sim[0]->getTimes();
	vector<vector<double>> X_hist;
	for (auto S : sim) {
		X_hist.push_back(S->getPos());
	}
	for (unsigned int i = 0; i < T.size(); i++) {
		out << setw(15) << T[i];
		for (unsigned int j = 0; j < X_hist.size(); j++) {
			out << setw(15) << X_hist[j][i];
		}
		out << endl;
	}
	out.close();
	cout << "Multiple positions saved" << endl; 
	
}

//Better, use this! vvv
void saveMultisysPos2(const vector<Simulation*>& sim) {
	cout << "Saving multiple positions 2..." << endl;
	
	//Some 2D array magic - Those turned out to be waayyyfaster than std::vectors
	double** T = new double*[sim.size()];
	double** X_hist = new double*[sim.size()];
	vector<double> xtemp, Ttemp;

	cout << "Generating array..." << endl;
	
	for (unsigned int i = 0; i < sim.size();i++) {
		if (i % 10 == 0)cout << "\t " << i << endl;

		xtemp = (sim[i]->getPos());
		Ttemp = (sim[i]->getTimes());
		X_hist[i] = new double[(int)xtemp.size() + 1]{ (double)xtemp.size() };
		T[i] = new double[(int)xtemp.size() + 1]{ (double)xtemp.size() };
		for (int j = 1; j < xtemp.size()+1; j++) {
			X_hist[i][j] = xtemp[j-1];
			T[i][j] = Ttemp[j-1];
		}
	}

	ofstream out{ "multisys2.txt", ios::out | ios::trunc };
	cout << "Array generated!\nWriting to file..." << endl;
	for (unsigned int t = 0; t < sim.size(); t++) {
		cout << t << endl;
		for (unsigned int i = 1; i < T[t][0]+1; i+=5) {
			out << setw(15) << T[t][i] << setw(15) << X_hist[t][i] << endl;
		}
		out << "#" << endl;
	}
	out.close();
	cout << "Multiple positions saved" << endl;

    //Clean up after the fun is over
	delete[] T;
	delete[] X_hist;
}
