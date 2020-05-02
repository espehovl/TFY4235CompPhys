#pragma once
#include <random>
#include <ctime>
#include "Simulate.h"
using namespace std;
static default_random_engine generator(42);// initialize the random number generator with a fixed seed, for reproducibility purposes


//Trial no.2
//Gaussian random number generator, made slightly complicated, really
double gauss(double stdDev = 1.0, double mean = 0.0, default_random_engine& gen = generator);
void testGauss();

//The flashing function of the ratchet potential
double flasher(double t_hat, double tau, double omega);
										  
//Reduced units force
double dU_hat(double x_hat, double alpha, double t_hat, double tau, double omega);
										  
//Reduced units potential
double U_hat(double x, double L, double alpha, double t_hat, double tau, double omega, bool active = true);
										  
//Reduced units Euler scheme
double EulerScheme(double x_hat, double tn, double dt, double D_hat, double alpha, double tau, double omega);

//Returns max of force, given alpha
double dU_MAX(double alpha);

//Finds dt according to criterion
double find_dt(double D_hat, double alpha, double omega);

//Save the trajectories of multiple simulations into a single file
//Requires equal time steps for all particles
void saveMultisysPos(const vector<Simulation*>& sim);
//more flexible in terms of different numbers of iterations
void saveMultisysPos2(const vector<Simulation*>& sim);
