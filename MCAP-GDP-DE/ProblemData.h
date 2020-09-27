#pragma once
#include <iostream>;  // for  cout cin
#include <vector>; 
#include<random>;
#include <string>;
#include <fstream>; // to read and write from/on files
#include <cstdio>;  // for puts and printf
#include <cstdlib>; // for string converstion, rand srand, dynamic memory management (C++ header)
#include <iterator>;
#include<math.h>;
#include<chrono>;   // for high resolution timing (elapsed time in microseconds and so on)
#include<array>;
#include<numeric>;  // for iota member in std which is used to initialize an array with 1 to n.
#include<map>;  // To define dictionary data types
using namespace std;



struct Bundle
{
	vector<int> ID;
	int NodesNumber = this->ID.size();
	vector<int> LeafNodes;
	int root;
	double weight;
	double obj;	
	double* qHat;
	Bundle(vector<int> ID0, int* Nodes0, vector<int> leaf0, int root0) {
		this->ID = ID0;
		this->LeafNodes = leaf0;
		this->root = root0;
	}
	Bundle() {}
	~Bundle() {}
};


class ProblemData
{


public:
	static int T;      // Number of Stages
	static int r;      // Number of Realizations
	static int R;      // Number of Resources
	static int Ts;     // Number of Tasks
	static int N;     // Number of Nodes
	static int nS;     // Number of Scenarios
	static unsigned seed;
	static double Kappa;      // Upper limit of expansion

	static std::vector<int> Pred;         // Predecessor array
	static double* Pr;
	static double** f;
	static double** g;
	static double** s;
	static double** d;
	static double*** c;

	static void PopulateParameters();

	ProblemData();
	~ProblemData();
};



