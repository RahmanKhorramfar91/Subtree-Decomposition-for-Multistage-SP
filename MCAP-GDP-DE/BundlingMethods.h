#pragma once
#include "Functions.h";
#include"ProblemData.h";

 class BundlingMethods
{
public:
	static int T;
	static int R;
	static	int r;
	static int Ts;
	static	double Kappa;
	static	int N;
	static	int nS;
	static	double** f;
	static	double** g;
	static	double** s;
	static	double** d;
	static	double*** c;
	static vector<int> Pred;

	static double* Pr;

	static vector<int> Nodes_in_Period_t(int t);
	static vector<int> Bundle_Nodes(int root, vector<int> leaves);

	static vector<Bundle> DiameterExpansion(int D, int L);

	static void InitiateData();
};

