
#include "ProblemData.h";
//#include "Functions.h";
#include "BundlingMethods.h";
#include <ilcplex/ilocplex.h>;

using namespace std;

typedef IloArray<IloNumVarArray> NumVar2D; // to define 2-D decision variables
typedef IloArray<NumVar2D> NumVar3D;	   // to define 3-D decision variables
typedef IloArray<NumVar3D> NumVar4D; 	   // to define 4-D decision variables

/*
HOW TO SET UP CPLEX FOR C++  :
https://bzdww.com/article/134619/
https://www.youtube.com/watch?v=Hbn1pGWLeaA
1) go to IBM directory in the program files in driver C
2) Go to concert folder -> include and when you see the "ilconcert" folder, copy the directory somewhere
3) Goto cplex folder -> include and when you see the "ilcplex" folder, copy the directory somewhere
4) In the solution Explorer tab, click on the project name and select properties
5) Go to C/C++ general ->" additional include directories" -> paste the two directories:
C:\Program Files\IBM\ILOG\CPLEX_Studio129\concert\include
C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\include

6) Go to C/C++ general ->"Preprocessors" and add these words:
WIN32
_CONSOLE
IL_STD
_CRT_SECURE_NO_WARNINGS


Or
NDEBUG
_CONSOLE
IL_STD


7)  In the Project1 property page, select: "c/c++" - "code generation" - "runtime library",
set to "multithreaded DLL (/MD)". determine.

8) In the Project1 property page, select: "Linker" - "Input" - "Additional Dependencies",
and then enter the path of the following three files:
C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\lib\x64_windows_vs2017\stat_mda\cplex1280.lib
C:\Program Files\IBM\ILOG\CPLEX_Studio128\cplex\lib\x64_windows_vs2017\stat_mda\ilocplex.lib
C:\Program Files\IBM\ILOG\CPLEX_Studio128\concert\lib\x64_windows_vs2017\stat_mda\concert.lib


9) if you're using visual studio 2017 with cplex 12.8, you may encounter an error which you then may
follow this link: https://www-01.ibm.com/support/docview.wss?uid=ibm10718671

*/

int main() {
	auto start = chrono::high_resolution_clock::now();

	/*int* ptr;
	int var = 7; int foo = 21;
	ptr = &var;
	ptr = &foo;
	int& ref = var;
	cout << ref << endl;
	cout << &foo << endl;
	cout << &ptr << "\t" << ptr << "\t" << *ptr<< endl;
*/

	ProblemData::PopulateParameters();
#pragma region Fetch Data
	int T = ProblemData::T;
	int R = ProblemData::R;
	int r = ProblemData::r;
	int Ts = ProblemData::Ts;
	double Kappa = ProblemData::Kappa;
	int N = ProblemData::N;
	int nS = ProblemData::nS;
	double** f = ProblemData::f;
	double** g = ProblemData::g;
	double** s = ProblemData::s;
	double** d = ProblemData::d;
	double*** c = ProblemData::c;
	vector<int> Pred = ProblemData::Pred;

	double* Pr = ProblemData::Pr;
#pragma endregion


	vector<Bundle> B = BundlingMethods::DiameterExpansion(4, 4);


#pragma region Get shared and succeding bundles of each node (Sn, Rn)
	map<int, vector<int>> Rn;
	map<int, vector<int>> Sn;
	for (int b = 0; b < B.size(); b++)
	{
		vector<int> rho = GetRho(0, B[b].root);
		for (int h : rho)
		{
			if (h == B[b].root) { continue; }
			//if (Rn.find(h) == Rn.end())
			{
				// the key does not exist in the map
				Rn[h].push_back(B[b].root);
			}
		}
		for (int id : B[b].ID) { Sn[id].push_back(b); }
		rho.end();
	}
#pragma endregion


#pragma region	 Get the distribution factor for the reference nodes. Adjust probs inside a bundle
	
	//double* DFs = new double[T] { 1,1,1,0.95,0.95,0.95,0.7,0.7,0.7,1,1,1,1,1 };
	double* DFs = new double[T] { 1, 1, 1, 0.9,0.9,0.9, 1,1,1, 1, 1, 1, 1, 1 };
	//map<int, double> V = Reference_bundle_distribution_factor(DFs);  // root node's distribution factor at reference
	map<int, double> V;
	for (int b = 0; b < B.size(); b++)
	{		
		for (int n : B[b].ID)
		{
			if (Rn[n].size() >= 1) { V[n] =1; }
			else { V[n] = 1; }
		}

		// adjust probabilitites inside a bundle
		double ss = 0;
		for (int l : B[b].LeafNodes) { ss += Pr[l] / Sn[l].size(); }
		B[b].weight = ss;
		for (int l : B[b].LeafNodes)
		{
			auto ind = std::find(B[b].ID.begin(), B[b].ID.end(), l) - B[b].ID.begin();
			B[b].qHat[ind] = Pr[l] / ss;
		}

		for (int i = B[b].ID.size()-1; i>0; i--)
		{
			int nn = Pred[B[b].ID[i]];
			//vector<int>::iterator ind = std::find(B[b].ID.begin(), B[b].ID.end(), nn);
			auto  ind = std::find(B[b].ID.begin(), B[b].ID.end(), nn)- B[b].ID.begin();
			B[b].qHat[ind] += B[b].qHat[i];
		}
	}

	//for (int n:B[0].ID)
	//{
	//	cout << "V[" << n << "] = " << V[n] << endl;
	//	//cout << "R[" << n << "] = " << Rn[n][0] << endl;
	//	cout << "q hat [" << n << "] = " << B[0].qHat[n] << endl;
	//}

#pragma endregion
	//cout << V[0] << endl;

#pragma region Solve each bundle and aggregate the contribution toward the lower bound
	double TotalObj = 0;
	for (int b = 0; b < B.size(); b++)
	{
		double obj = MCAP(B[b], Rn, V,b);
		TotalObj += obj * B[b].weight;
	}

#pragma endregion
	
	auto end = chrono::high_resolution_clock::now();
	auto Elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
	std::cout << "\n\t Elapsed Time: " << Elapsed.count() << endl;
	std::cout << "DE Lower Bound : " << TotalObj << endl;
	std::system("pause");
	return 0;
}