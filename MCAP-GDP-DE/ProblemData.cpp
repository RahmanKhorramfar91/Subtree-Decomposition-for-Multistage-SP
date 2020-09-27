#include "ProblemData.h";

unsigned ProblemData::seed = 50;
std::default_random_engine Gen(ProblemData::seed);
ProblemData::ProblemData() {}
ProblemData::~ProblemData() {}
using namespace std;

#pragma region Random Functions
double randint(int a, int b) {
	std::uniform_int_distribution<int> randint(a, b);
	return randint(Gen);
}

double unifrnd(double a, double b) {
	std::uniform_real_distribution<double>  unifrnd(a, b);
	return unifrnd(Gen);
}

#pragma endregion

#pragma region Basic Parameters
int ProblemData::T = 8;
int ProblemData::r = 2;

int ProblemData::R = 3;
int ProblemData::Ts = 5;
double ProblemData::Kappa = 1;
int ProblemData::N = (pow(ProblemData::r, ProblemData::T) - 1) / (ProblemData::r - 1);
int ProblemData::nS = pow(ProblemData::r, ProblemData::T - 1);

vector<int> ProblemData::Pred;
double* ProblemData::Pr = new double[ProblemData::N]();
double** ProblemData::f = new double* [ProblemData::R]();
double** ProblemData::g = new double* [ProblemData::R]();
double** ProblemData::s = new double* [ProblemData::Ts]();
double** ProblemData::d = new double* [ProblemData::Ts]();
double*** ProblemData::c = new double** [ProblemData::N]();

#pragma endregion


void ProblemData::PopulateParameters() {


	// populate f and g matrices
	for (int i = 0; i < R; i++)
	{
		f[i] = new double[N]();
		g[i] = new double[N]();
		for (int n = 0; n < N; n++)
		{
			f[i][n] = unifrnd(25.0, 50.0);// uniform[25,50]
			g[i][n] = unifrnd(5.0, 10.0); //// uniform[5,10]
		}
	}

	// Populate s and d
	for (int j = 0; j < Ts; j++)
	{
		s[j] = new double[N]();
		d[j] = new double[N]();
		for (int n = 0; n < N; n++)
		{
			s[j][n] = randint(500, 1000); //uniform[500,1000]
			d[j][n] = unifrnd(0.5, 1.5); //uniform[0.5,1.5]
		}
	}

	for (int n = 0; n < N; n++)
	{
		c[n] = new double* [R]();
		for (int i = 0; i < R; i++)
		{
			c[n][i] = new double[Ts]();
			for (int j = 0; j < Ts; j++)
			{
				c[n][i][j] = unifrnd(5, 10);
			}
		}
	}

	// get the pred vector, indicating the predecessor of each node
	int tr = N - nS;
	//vector<int> Pred;
	Pred.push_back(-1);
	int prs = 0;
	for (int i = 0; i < tr; i++)
	{
		for (int j = 0; j < r; j++)
		{
			Pred.push_back(prs);
		}
		prs++;
	}



	// Populate the probabitlity vector
	std::vector<int> nodes(N);
	std::iota(std::begin(nodes) + 1, std::end(nodes), 1); // iota member requires numeric library #include<numeric>
	double* Leaf_Node_Prob = new double[nS]();
	double ss = 0;  // store sum of probabilties at leaf nodes
	for (int i = 0; i < nS; i++) { Leaf_Node_Prob[i] = unifrnd(0, 1); ss += Leaf_Node_Prob[i]; }
	for (int i = 0; i < nS; i++) { Leaf_Node_Prob[i] = Leaf_Node_Prob[i] / ss; }
	for (int i = N - nS; i < N; i++) { Pr[i] = Leaf_Node_Prob[i - N + nS]; }

	double ssum = 0;
	for (int i = 0; i < nS; i++)
	{
		ssum += Leaf_Node_Prob[i];
	}
	while (nodes.size() >= r)
	{
		int last = nodes[nodes.size() - 1]; ss = 0;
		for (int i = 0; i < r; i++)
		{
			ss += Pr[nodes[nodes.size() - 1]];
			nodes.pop_back();
		}
		Pr[Pred[last]] = ss;
	}

}



