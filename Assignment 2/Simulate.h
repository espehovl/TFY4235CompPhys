#pragma once
#include <vector>
#include <map>
#include <string>


using namespace std;

class Simulation {
private:
	//Units are SI units, in meter
	double L, alpha, etta, kBT, deltaU, gamma, omega, r, D_hat, dt, time;
	double tau; //[s]
	double elapseTime; //[s] How long time should the particles tracks be simulated for?
	int steps;
	bool hasRun = false;
	bool useTimeLimit; //Are we using time limit or step no. limit?
	
	vector<double> Xs; //Reduced units
	vector<double> Ts; //Reduced units
public:
	Simulation( double _L, double _alpha, double _etta, 
		double _kBT, double _deltaU, double _tau, double _r, int _steps = 0,
		double _elapseTime = 0.0);
	
	//Performs steps iterations, and saves the data to file
	void perform(bool writeToFile = false);

	//Saves the potential profile as a .txt file. 
	//Adapts to the actual span of x-values for niceness
	void savePotProfile() const;

	//Save the entire potential history along the track to a .txt file in [eV]
	void savePotentialOfTrack() const;

	//Find the final position of the particle after a run
	double getEndPos() const; //reduced units!

	double getEndPotential() const; //"Real" units [J]

	//Find the maxest position of the trail
	vector<double> getMaxPos() const; //reduced units!

	//Calculate (naively) average drift velocity
	double getDriftVelocity() const ; //[m/s]

	//Save metadata to .txt file
	void saveMeta()const ;

	//Useful get-functions
	vector<double> getTimes() const; //[s]
	vector<double> getPos() const; //[µm]

	void saveMultisysMetadata(int ctr) const;

};