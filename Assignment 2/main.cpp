#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <iomanip>
#include <random>
#include <fstream>
#include <cmath>
#include <chrono>

#include "Simulate.h"
#include "Source.h"


using namespace std;

int main() {
	//The optimal tau value for the smaller particle (calculated earlier on)
	constexpr double tau_op = 0.46;

	//Task 6 - Laying down the parameters - in SI units
	constexpr double r = 12.0e-9;				//m
	constexpr double L = 20e-6;					//m
	constexpr double alpha = 0.2;				//dimensionless
	constexpr double etta = 1.0e-3;				//Pa*s 
	constexpr double kBT = 26.0e-3*1.609e-19;	//J
	constexpr double deltaU = 80.0*1.609e-19;	//J
	constexpr double tau = tau_op;				//s

	//The radius of the larger particles (often referred to as DNA)
	constexpr double r2 = 3 * r; //Task 11

	//Generate a vector of taus to test
	vector<double> TAU;
	for (double d = 0.01; d < 5.0; d += 0.05) TAU.push_back(d); //Fine res. 
	//for (double d = 1.0; d < 5.0; d += 0.1) TAU.push_back(d); //Coarse res.

	//Single particle simulation or multiple realizations of same particle
	if (true) {
		
		//How long should the experiment go on?
		constexpr double elapseTime = 28.0; //[s]
		
		/*
		cout << "Single particle coming up..." << endl;
		Simulation sim{ L, alpha, etta, kBT, deltaU, tau, r ,1000, elapseTime };
		sim.perform(true);
		sim.savePotProfile();
		sim.saveMeta();
		sim.savePotentialOfTrack();
		*/

		//Compare with DNA
		Simulation DNA{ L, alpha, etta, kBT, deltaU, tau, r, 1000, elapseTime };
		DNA.perform(false);


		//Prepare writing to file!
		ofstream metadata{ "MetadataMultisys.txt", ios::out | ios::trunc };
		metadata << "particle #\t"
			<< "dt(s)\t"
			<< "tau(s)\t"
			<< "L(µm)\t"
			<< "steps\t"
			<< "r(nm)\t"
			<< "deltaU(eV)\t"
			<< "kBT(eV)\t"
			<< "End position(µm) \t"
			<< "Avg drift vel. (µm/s) \t"
			<< "Degeneracy: " << 1 << "\t"
			<< endl;
		metadata.close();

		//sim.saveMultisysMetadata(1);
		//DNA.saveMultisysMetadata(2);
		//Make a vector of pointers to the different particle objets
		vector<Simulation*> DNAs{&DNA};

		cout << "Performing simulations..." << endl;
		for (unsigned int i = 0; i < 1000; i++) {
			DNAs.push_back(new Simulation{ L, alpha, etta, kBT, deltaU, tau, r2, 100000, elapseTime });
			DNAs.back()->perform();
			//DNAs.back()->saveMultisysMetadata(i + 1);
		}
		cout << "Simulations performed!" << endl;
		saveMultisysPos2(DNAs);
	}
	//Plot a few particles with different values for tau
	if (false) {
		cout << "A few particles coming up..." << endl;
		//Create the appropriate vectors
		vector<Simulation*> sims;
		vector<double> fewTau{ 0.25, 0.38, 0.46, 0.70, 1.0 };

		//prepare the metadata file
		ofstream metadata{ "MetadataMultisys.txt", ios::out | ios::trunc };
		metadata << "particle #\t"
			<< "dt(s)\t"
			<< "tau(s)\t"
			<< "L(µm)\t"
			<< "steps\t"
			<< "r(nm)\t"
			<< "deltaU(eV)\t"
			<< "kBT(eV)\t"
			<< "End position(µm) \t"
			<< "Avg drift vel. (µm/s) \t"
			<< "Degeneracy: " << 1 << "\t"
			<< endl;
		metadata.close();

		//Run simulations
		cout << "Running simulations..." << endl;
		for (unsigned int i = 0; i < fewTau.size(); i++) {
			sims.push_back(new Simulation{ L, alpha, etta, kBT, deltaU, fewTau[i], r,100000 });
			sims.back()->perform(false);
			sims.back()->saveMultisysMetadata(i+1);
		}
		
		//Save the metadata
		saveMultisysPos(sims);
	}



	//Multiple realizations (to find tau_op) 
	//Takes a lot of time!
	if (false) { //Multiple realizations
		cout << "Creating multiple simulations..." << endl;
		//Run multiple realizations
		
		constexpr int degeneracy = 1000; //How many trials at each tau?
		const unsigned int runs = (unsigned int)TAU.size(); //How many taus?
		constexpr int steps = 100000; //How many iterations for each particle?

		//Prepare the metadata file
		ofstream metadata{ "MetadataMultisys.txt", ios::out | ios::trunc };
		metadata << "particle #\t"
			<< "dt(s)\t"
			<< "tau(s)\t"
			<< "L(µm)\t"
			<< "steps\t"
			<< "r(nm)\t"
			<< "deltaU(eV)\t"
			<< "kBT(eV)\t"
			<< "End position(µm) \t"
			<< "Avg drift vel. (µm/s) \t"
			<< "Degeneracy: " << degeneracy << "\t"
			<< endl;
		metadata.close();

		Simulation* S = nullptr;
		int ctr = 0;
		auto start = chrono::high_resolution_clock::now(); //For recording the time it takes for each simulation
		for (unsigned int i = 0; i < runs; i++) { //For each tau...
			for (unsigned int j = 0; j < degeneracy; j++) { // For each trial at each tau...
				ctr++;
				S = new Simulation{ L, alpha, etta, kBT, deltaU, TAU[i], r, steps }; //Create a new Simulation object
				S->perform( false);
				S->saveMultisysMetadata(ctr);
				if (ctr % 100 == 0) {
					auto stop = chrono::high_resolution_clock::now();
					auto dur = chrono::duration_cast<chrono::milliseconds>(stop - start);
					cout << setw(6) << ctr << "/" << runs * degeneracy << "\t" << setw(10) << dur.count() / 1000.0 << " sec/100 iterations" << endl;
					start = chrono::high_resolution_clock::now();
				}
				delete S; //Delete and reset the pointer
				S = nullptr;
			}
		}
		cout << "Realizations complete!" << endl;

	}

	//Automagically plot things in Python
	cout << "Calling Python plotting..." << endl;
	//system(string{ "python posTracking.py" }.c_str());
	//system(string{ "python boltzmann.py" }.c_str());
	//system(string{ "python MultisysReader.py" }.c_str());
	//system(string{ "python diffusion.py" }.c_str());


	//Keep window open
	cout << "Program finished, enter any character... ";
	cin.get();
	return 0;
}



