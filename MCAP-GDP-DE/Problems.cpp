#include"Functions.h";
#include<ilcplex/ilocplex.h>;


typedef IloArray<IloNumVarArray> NumVar2D; // to define 2-D decision variables
typedef IloArray<NumVar2D> NumVar3D;	   // to define 3-D decision variables
typedef IloArray<NumVar3D> NumVar4D; 	   // to define 4-D decision variables

double MCAP(Bundle Bdl, map<int, vector<int>> Rn, map<int, double> V, int bundleNum)
{
	//ProblemData::PopulateParameters();
#pragma region Fetch Data
	int T = ProblemData::T;
	int R = ProblemData::R;
	int r = ProblemData::r;
	int Ts = ProblemData::Ts;
	double Kappa = ProblemData::Kappa;
	//int N = ProblemData::N;
	int nS = ProblemData::nS;
	double** f = ProblemData::f;
	double** g = ProblemData::g;
	double** s = ProblemData::s;
	double** d = ProblemData::d;
	double*** c = ProblemData::c;
	vector<int> Pred = ProblemData::Pred;

	double* Pr = ProblemData::Pr;
#pragma endregion

	int Nb = Bdl.ID.size();
#pragma region Define Decision Variables
	IloEnv env;
	IloModel Model(env);
	int root = Bdl.root;
	NumVar2D Xp(env, root);
	NumVar2D Up(env, root);

	NumVar2D X(env, Nb);
	NumVar2D U(env, Nb);

	NumVar3D Y(env, Nb);
	NumVar2D Z(env, Nb);
	bool Linear_Relaxation = false;
	for (int n = 0; n < Nb; n++)
	{		
		X[n] = IloNumVarArray(env, R, 0, Kappa, ILOFLOAT);
		U[n] = IloNumVarArray(env, R, 0, 1, ILOBOOL);
		Z[n] = IloNumVarArray(env, Ts, 0, 1, ILOBOOL);

		Y[n] = NumVar2D(env, R);
		for (int i = 0; i < R; i++)
		{
			Y[n][i] = IloNumVarArray(env, Ts, 0, 1, ILOBOOL);
		}

	}

	for (int n = 0; n < Bdl.root; n++)
	{
		Xp[n] = IloNumVarArray(env, R, 0, Kappa, ILOFLOAT);
		Up[n] = IloNumVarArray(env, R, 0, 1, ILOBOOL);
	}
#pragma endregion



#pragma region Objective Function
	IloExpr ex0(env);
	for (int n = 0; n < Nb; n++)
	{
		for (int i = 0; i < R; i++)
		{
			double Vs = V[Bdl.ID[n]];
			ex0 += Vs * Bdl.qHat[n] * f[i][Bdl.ID[n]] * U[n][i] + Vs * Bdl.qHat[n] * g[i][Bdl.ID[n]] * X[n][i];
		}
		for (int j = 0; j < Ts; j++)
		{
			ex0 += Bdl.qHat[n] * s[j][Bdl.ID[n]] * Z[n][j];
			for (int i = 0; i < R; i++)
			{
				ex0 += Bdl.qHat[n] * c[Bdl.ID[n]][i][j] * Y[n][i][j];
			}
		}
	}

	vector<int> rho = GetRho(0, Bdl.root);
	for (int h : rho)
	{
		if (h == Bdl.root) { continue; }
		// probability of the root node is always 1 after adjustemnt.
		// Therefore, all precedeeding nodes would also have probability of 1. 
		// Think about an extended tree whose X variables only are considered before 
		// the bundle's root node. 
		double Vnr = Pr[h] * (1 - V[h]) / (Pr[h] * Rn[h].size());
		for (int i = 0; i < R; i++)
		{
			ex0 += Vnr*f[i][h] * Up[h][i] +Vnr* g[i][h] * Xp[h][i];
		}

	}
	Model.add(IloMinimize(env, ex0));
	ex0.end();
#pragma endregion


#pragma region Cons1: Upper limit restriction on X, Xp variables
	for (int n = 0; n < Nb; n++)
	{
		for (int i = 0; i < R; i++)
		{
			Model.add(X[n][i] <= Kappa * U[n][i]);
		}
	}

	for (int h : rho)
	{
		if (h == Bdl.root) { continue; }
		for (int i = 0; i < R; i++)
		{
			Model.add(Xp[h][i] <= Kappa * Up[h][i]);
		}
	}
#pragma endregion


#pragma region Cons 2: Capacity assignment for each task
	for (int n = 0; n < Nb; n++)
	{
		for (int i = 0; i < R; i++)
		{
			IloExpr ex2(env);
			for (int j = 0; j < Ts; j++)
			{
				ex2 += d[j][Bdl.ID[n]] * Y[n][i][j];
			}

			vector<int> path1 = GetRho(Bdl.root, Bdl.ID[n]);
			for (int h : path1)
			{
				auto kk = std::find(Bdl.ID.begin(), Bdl.ID.end(), h) - Bdl.ID.begin();
				ex2 -= X[kk][i];
			}

			vector<int> path2 = GetRho(0, Bdl.root);
			for (int p : path2)
			{
				if (p == Bdl.root) { continue; }
				ex2 -= Xp[p][i];
			}
			Model.add(ex2 <= 0);
			ex2.end();

		}

	}
#pragma endregion


#pragma region Cons 3: Task completion
	for (int n = 0; n < Nb; n++)
	{
		for (int j = 0; j < Ts; j++)
		{
			IloExpr ex3(env);
			for (int i = 0; i < R; i++)
			{
				ex3 += Y[n][i][j];
			}
			ex3 += Z[n][j];
			Model.add(ex3 == 1);
		}
	}
#pragma endregion

	IloCplex cplex(Model);
	//cplex.exportModel("Bunlde0.lp");
	//cplex.setOut(env.getNullStream()); // turn off cplex's logging to console
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.001);
	if (bundleNum==0)
	{
		cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.001);
	}
	
	if (!cplex.solve()) {
		env.error() << "Failed to optimize  the problem " << endl;
		throw(-1);
	}
	double obj = cplex.getObjValue();
	//env.out() << "\t ENG Model Solution Status: " << cplex.getStatus() << endl;
	cout << "\t Bundle:  " <<bundleNum<< " Objective Value: " << obj << endl;
	//cplex.writeSolutions("ENGSol.txt");
	//cplex.writeSolution("ENGSol2.txt");

	//for (int n = 0; n < Nb; n++)
	//{
	//	for (int i = 0; i < R; i++)
	//	{
	//		double sl = cplex.getValue(X[n][i]);
	//		if (sl > 0.1)
	//		{
	//			cout << "X[" << n << "][" << i << "] = " << sl << endl;
	//		}
	//	}
	//}
	//std::system("pause");

	Model.end();
	cplex.end();
	Bdl.obj = obj;
	return obj;
}

