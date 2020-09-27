
#include "Functions.h";

vector<int> GetRho(int m, int n) {

	std::vector<int> Pred = ProblemData::Pred;
	std::vector<int> rho;
	if (n == 0) { rho.push_back(m); return rho; }
	if (n == m) { rho.push_back(m); return rho; }

	rho.push_back(m); rho.push_back(n);
	int prec = Pred[n];
	while (prec != m)
	{
		rho.push_back(prec);
		prec = Pred[prec];
		if (prec == -1)
		{
			rho.clear();
			return rho;
		}

	}
	return rho;
}


vector<int> Nodes_in_Period_t(int t0)
{
	/*
	argument: t0, period t0
	output: nodes of the scenario tree in period t0
	*/
	int r = ProblemData::r;
	int N = ProblemData::N;
	int N0 = (int)(pow(r, t0) - 1) / (r - 1);
	int S0 = (int)pow(r, t0 - 1);

	vector<int> nodesT;
	int ss = 1; int j = 0;
	for (int i = S0; i > 0; i--)
	{
		nodesT.push_back(N0 - i);
		j++;
	}
	return nodesT;
}


map<int, double> Reference_bundle_distribution_factor(double* DFs) 
{
	map<int, double> V1;
	int T = ProblemData::T;
	for (int t = 0; t < T; t++)
	{
		vector<int> nodesT1 = Nodes_in_Period_t(t);
		for (int h : nodesT1)
		{
			V1[h] = DFs[t];
		}
		nodesT1.end();
	}
	return V1;
}