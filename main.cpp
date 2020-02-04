#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <map>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <chrono>

using namespace std;

void addToMap(map<long long, int> &m, int step);

int main()
{
	srand(static_cast<unsigned int>(time(nullptr)));

	map<long long,int> stats; //Percentages in int form pairs: percentage-number of hits

	auto start = chrono::high_resolution_clock::now();
	int overflows = 0;
	unsigned int limit = (int)1e6;
	for (unsigned int i = 0; i < limit; i++){
		long long ctr = 1;
		double dist = 0;
		bool right; //true == first step was to the right, false == first step was to the left
		//dist += (rand() % 100 + 1) / static_cast<double>(100); //First step is to the right
		
		double step;
		//Step length and direction is based on uniform distributions, with 0.01 resolution
		step = (rand() % 100 + 1) / static_cast<double>(100) * pow((-1), rand() % 2);
		right = ((ctr == 1) && (step > 0)) ? true : false;
		dist += step;
		//cout << dist << ", ";
		while (true)
		{
			ctr++;
			step = (rand() % 100 + 1) / static_cast<double>(100) * pow((-1), rand() % 2); //Randomly picks a step between [-1,1] ("not discrete", i.e. 0.01 resolution)
			dist += step;
			//cout << dist << ", ";
			if (ctr > 1000) {
				//cout << "Overflow, not interesting" << endl;
				overflows++;
				break;
			}

			if (dist < 0 && ctr > 1 && right)
			{
				//cout << "\nReturn at step(initial right): " << ctr << endl;
				addToMap(stats, ctr);
				break;
			}
			else if (dist > 0 && ctr > 1 && !right){
				//cout << "\nReturn at step(initial left): " << ctr << endl;
				addToMap(stats, ctr);
				break;
			}
		}
		if((i%(int)(limit/static_cast<double>(10)))==0){
		 	std::cout << i*100/limit << "% - time elapsed: ";
			auto stop = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
			cout << duration.count() << "s." << endl;
		}
	}

	//Open file for writing
	ofstream file;
	file.open("plotdata.txt");


	//Print a "beautiful" histogram
	cout << "\n--------HISTOGRAM--------" << endl;
	cout << setw(6) << "#Steps" << setw(8) << "freq" << endl;
	int thresh = 20;
	for (auto e : stats)
	{
		if(e.first < thresh){
			cout << setw(6) << e.first << "| " << setw(6) << (int)(e.second * (100 / static_cast<double>(limit-overflows))) << "| ";
			for (int i = 0; i < (int)(e.second * (100 / static_cast<double>(limit-overflows))); i++)
				cout << "*";
			cout << endl;
			file << e.first << "," << e.second << "\n"; //Write to file
		}
		else if (e.first > thresh) break;
	}
	
	file.close();
	cout << "No. overflows: " << overflows << endl;
	std::cin.get();
	return 0;
}

void addToMap(map<long long,int>& m,int step){
	map<long long, int>::iterator it = m.find(step);
	 if (it != m.end()) it->second++;
	 else m.insert(make_pair(step, 1));
	 return;
}