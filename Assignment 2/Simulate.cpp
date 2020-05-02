#include "Simulate.h"
#include "Source.h"

#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;

//Construct the object and initialize variables
Simulation::Simulation(double _L, double _alpha, double _etta, double _kBT,
	double _deltaU, double _tau, double _r, int _steps, double _elapseTime)
	: L{ _L }, 
	alpha{ _alpha }, 
	etta{ _etta }, 
	kBT{ _kBT }, 
	gamma{ 6 * static_cast<double>(_Pi)*_etta*_r },
	omega{ _deltaU / (gamma*_L*_L) }, 
	r{ _r }, 
	D_hat{ _kBT / _deltaU }, 
	deltaU{_deltaU}, 
	tau{_tau},
	steps{_steps},
	elapseTime{_elapseTime},
	Xs(steps), //fix lengths of vectors, saves a LOT of time!
	Ts(steps)
{
	if (elapseTime != 0.0) {
		useTimeLimit = true;
	}
	else {
		useTimeLimit = false;
	}
}

//Perform the actual simulations
void Simulation::perform(bool writeToFile) {
	if (hasRun) { //For avoiding unintentional extra runs
		cout << "Simulation has already been run!" << endl;
		return;

	}
	//Reduced units!
	double x_0 = 0.0;
	double t = 0.0;

	//Find the timestep size
	dt = find_dt(D_hat, alpha, omega); //Reduced unit!

	//Initialize position and time vectors
	Xs[0] = x_0;
	Ts[0] = t;

	ofstream ofs;
	if (writeToFile) {
		ofs.open("pos_data.txt", ios::out | ios::trunc);
		ofs << setw(15) << "x" << setw(20) << "t" << endl;
		ofs << setw(15) << Xs[0] * L*1e6 << setw(20) << t / omega << endl;
	}
	if (!useTimeLimit) {
		for (int i = 1; i < steps; i++) {
			t += dt;
			Xs[i] = (EulerScheme(Xs[i - 1], t, dt, D_hat, alpha, tau, omega));
			Ts[i] = t;
			if (writeToFile) ofs << setw(15) << Xs[i] * L*1e6 << setw(20) << t / omega << endl;
		}
	}
	else if (useTimeLimit) {
		int i = 1;
		//Adjust Xs and Ts according to elapsed time and not number of iterations
		Xs.resize(static_cast<int>(elapseTime / (dt / omega)) + 1); //+1 is to adjust for int conversion
		Ts.resize(static_cast<int>(elapseTime / (dt / omega)) + 1);
		while (i<Xs.size()) { //Perform the iterations
			t += dt;
			Xs[i] = (EulerScheme(Xs[i - 1], t, dt, D_hat, alpha, tau, omega));
			Ts[i] = t;
			if (writeToFile) ofs << setw(15) << Xs[i] * L*1e6 << setw(20) << t / omega << endl;
			i++;
		}
	}
	if (writeToFile) {
		ofs << "# " << dt/omega << " " << tau << endl;
		ofs.close();
	}
	time = t / omega; // in seconds
	//cout << "Simulation performed with " << steps << " steps." << endl;
	hasRun = true;
}


//Save potential profile, adapts to the distance travelled by the particle, for niceness
void Simulation::savePotProfile() const{
	cout << "Saving potential profile..." << endl;
	ofstream ofs{ "pot.txt", ios::out | ios::trunc };
	
	double lower = -1.0, upper = 1.0;
	while (getMaxPos()[1] > upper) upper += 1.0; //Right max
	while (getMaxPos()[0] < lower) lower -= 1.0; //Left max
	for (double x = lower; x <= upper; x += 0.01) {
		ofs << U_hat(x , L, alpha, 0.1*omega, tau, omega, false) << "   " << x*L*1e6 << endl;
	}
	ofs.close();
}

//Save the potential history of the particle through the entire track
void Simulation::savePotentialOfTrack() const {
	if (hasRun) {
		cout << "Saving potential history..." << endl;
		ofstream ofs{ "potentialtracks.txt", ios::out | ios::trunc };
		for (unsigned int i = 0; i < Xs.size(); i++) {
			ofs << U_hat(Xs[i], L, alpha, Ts[i], tau, omega)*deltaU/1.609e-19 << endl; // [eV]
		}
		ofs.close();
	}
	else {
		cout << "Simulation has not run yet!" << endl;
		return;
	}
}

double Simulation::getEndPos()const{
	if (hasRun) return Xs.back(); //Returns the reduced unit end position!
	else {
		cout << "Simulation has not run yet!" << endl;
		return 0.0;
	}
}

double Simulation::getEndPotential()const{ //Returns the "real" end potential
	if (hasRun)	return U_hat(getEndPos(), L, alpha, 0.1*tau*omega, tau, omega, false)*deltaU; //[J]
	else {
		cout << "Simulation has not run yet!" << endl;
		return 0.0;
	}
}

vector<double> Simulation::getMaxPos() const{
	vector< double> MAX{ 0.0,0.0 };
	if (hasRun) {
		for (double x : Xs) {
			if (x > MAX[1]) MAX[1] = x; //Right max (pos)
			if (x < MAX[0]) MAX[0] = x; //Left max (neg)
		}
	}
	else cout << "Simulation has not run yet!" << endl;
	return MAX;
}

double Simulation::getDriftVelocity() const {
	if(hasRun) return getEndPos()*L / time; // [m/s]
	else {
		cout << "Simulation has not run yet!" << endl;
		return 0.0;
	}
}
//Save metadata of the particle
void Simulation::saveMeta() const {
	ofstream ofs{ "metadata.txt" };
	ofs << "dt(s)\t" << "tau(s)\t" << "L(µm)\t"
		<< "steps\t" << "r(nm)\t" << "deltaU(eV)\t" << "kBT(eV)\t" << endl;

	ofs << dt/omega << "\t" << tau << "\t" << L / (1e-6) << "\t"
		<< steps << "\t" << r / (1e-9) << "\t" << deltaU / (1.609e-19) << "\t" 
		<< kBT / (1.609e-19) << "\t" << endl;
	cout << "Metadata saved!" << endl;
	ofs.close();
}

vector<double> Simulation::getTimes() const {
	vector<double> T;
	for (double t : Ts) {
		T.push_back(t / omega);
	}
	return T; //In actual units [s]
}

vector<double> Simulation::getPos() const {
	vector<double> X;
	for (double x : Xs) {
		X.push_back(x*L*1e6);
	}
	return X; //In actual units [µm]
}


void Simulation::saveMultisysMetadata(int ctr) const {
	//cout << "Saving multiple metadata..." << endl;
	ofstream out{ "MetadataMultisys.txt", ios::out | ios::app };
		out << ctr << "\t"
			<< dt / omega << "\t"
			<< tau << "\t"
			<< L *1e6 << "\t"
			<< steps << "\t"
			<< r*1e9 << "\t"
			<< deltaU / 1.609e-19 << "\t"
			<< kBT / 1.609e-19 << "\t"
			<< getEndPos()*L *1e6 << "\t"
			<< getEndPos()*L *1e6 / time << "\t"
			<< endl;
	//cout << "Multiple metadata saved!" << endl;
	out.close();
}