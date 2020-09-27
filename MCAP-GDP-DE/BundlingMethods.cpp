#include "BundlingMethods.h"

#pragma region Fetch Data
int  BundlingMethods::T = ProblemData::T;
int  BundlingMethods::r = ProblemData::r;
int  BundlingMethods::N = ProblemData::N;
double* BundlingMethods::Pr = ProblemData::Pr;
#pragma endregion
void BundlingMethods::InitiateData() {
	//ProblemData::PopulateParameters();

	BundlingMethods::T = ProblemData::T;
	BundlingMethods::r = ProblemData::r;
	BundlingMethods::N = ProblemData::N;
	BundlingMethods::Pr = ProblemData::Pr;
}
vector<int> BundlingMethods::Nodes_in_Period_t(int t0)
{
	/*
	argument: t0, period t0
	output: nodes of the scenario tree in period t0
	*/
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

vector<int>  BundlingMethods::Bundle_Nodes(int root, vector<int> leaves)
{
	/*
	argument: "root" node of a bundle, and "leaves" of a bundle
	output: nodes of the original tree that are part of the bundle
	*/
	vector<int> Ns;
	for (int i : leaves)
	{
		vector<int> S1 = GetRho(root, i);
		for (int n : S1) { Ns.push_back(n); }

	}

	sort(Ns.begin(), Ns.end()); // sort the vecor
	Ns.erase(unique(Ns.begin(), Ns.end()), Ns.end()); // remove duplicate elements
	return Ns;
}

vector<Bundle>  BundlingMethods::DiameterExpansion(int D, int L)
{
	BundlingMethods::InitiateData();

	/*
	this function divided the scenario tree according to the diameter decomposition.
	The first bundle's root is 0, its leaves consist of nodes int the Dth period.
	nodes at D+1 period make up the root nodes for another set of bundles. Leaf nodes lie in the (D+1+L) period.
	This process continues until no node left.

	Arguments: D: diameter size, L= length of bunlde for bundles other than Bundle0.
	*/
	vector<Bundle> B;

	vector<int> SectionsBegin;
	vector<int> SectionsEnd;
	SectionsBegin.push_back(1);
	SectionsEnd.push_back(D);

	double rem = max(0, T - D);
	int s1 = std::ceil(rem / (double)L);

	for (int i = 0; i < s1; i++)
	{
		int l0 = D + 1 + i * L;
		int l1 = min(D + 1 + i * L + L - 1, T);
		SectionsBegin.push_back(l0);
		SectionsEnd.push_back(l1);
	}
	for (int i = 0; i < SectionsEnd.size(); i++)
	{
		vector<int> fN = BundlingMethods::Nodes_in_Period_t(SectionsBegin[i]);
		vector<int> lN = BundlingMethods::Nodes_in_Period_t(SectionsEnd[i]);

		int step = lN.size() / fN.size();

		for (int n : fN)
		{
			Bundle bundle; // create a new bundle and fill the members
			bundle.root = n;
			vector<int> leafB;
			//double BundleWeight = 0;
			for (int j = 0; j < step; j++)
			{
				leafB.push_back(lN[0]);
			//	BundleWeight += Pr[lN[0]];
				lN.erase(lN.begin());
			}
			bundle.LeafNodes = leafB;
			//bundle.weight = BundleWeight;

			vector<int> nodesB = BundlingMethods::Bundle_Nodes(n, leafB);
			bundle.ID = nodesB;
			bundle.qHat = new double[nodesB.size()]();
			B.push_back(bundle);  // add the new bundle to the list of bundles.
			bundle.~Bundle();
			leafB.end();
			nodesB.end();

		}
		fN.end();
		lN.end();
	}
	SectionsBegin.end();
	SectionsEnd.end();
	return B;
}


