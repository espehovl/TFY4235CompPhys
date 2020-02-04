#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <fstream>


using namespace std;

int main() {
	srand(static_cast<unsigned int>(time(nullptr)));
	int sz = 64; //Size of population

	ofstream file;
	file.open("baksnepperdata.txt");
	
	vector<int> species(sz); //List of available species in the simulation, default fitness is set to 0, number of species is sz
	vector<int> counts(species.size()); // A list of all counts (ages) of all species. This will be incremented as we go

	for (unsigned int i = 0; i < species.size(); i++) {
		species[i] = rand() % 100; //Assigns a random fitness factor [0-99] to each species
		//cout << setw(4) << species[i]; 
		file << species[i] << " ";
	}
	//cout << endl;
	file << ", ";
	for (unsigned int i = 0; i < counts.size(); i++) {
		file << counts[i] << " ";
	}
	file << "\n";
	
	
	int ctr = 10000;

	//Select and destroy the lowest fitness
	do {
		//vector<int> prev{ species }; //previous step, for comparison and aging
		int min = 100;
		int index = -1;
		for (unsigned int i = 0; i < species.size(); i++) {
			if (species[i] < min) {
				min = species[i];
				index = i;
			}
		}
		
		//Destroy, and replace the neighbours. Keep track of counts.
		
		species[index] = rand() % 100;
		
		if (index + 1 < sz) species[index + 1] = rand() % 100;
		else if (index + 1 == sz) species[0] = rand() % 100;

		if (index - 1 >= 0) species[index - 1] = rand() % 100;
		else if (index == 0) species[sz - 1] = rand() % 100;

		for (int i = 0; i < species.size(); i++) { 
			//std::cout << setw(4) << species[i]; //Print for debugging
			file << species[i] << " ";
			counts[i]++; //Increase ages(counts)
		}
		file << ", ";

		//Reset dead species counts
		counts[index] = 0;

		if (index + 1 < sz) counts[index + 1] = 0;
		else if (index + 1 == sz) counts[0] = 0;

		if (index - 1 >= 0) counts[index - 1] = 0;
		else if (index == 0) counts[sz - 1] = 0;


		ctr--;

		//cout << "\n";
		for (int i = 0; i < counts.size(); i++) {
			//cout << setw(4) << counts[i];
			file << counts[i] << " ";
		}
		file << "\n";
		//cout << "\n\n";
	} while (ctr > 0);

	file.close();

	std::cout << "\n\nProgram complete, press any key to continue...";
	//std::cin.get(); //Keeps window open
	return 0;
}

