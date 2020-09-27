#pragma once
#include "ProblemData.h";

vector<int> GetRho(int m, int n);
double MCAP(Bundle B, map<int, vector<int>> Rn, map<int, double> V, int bundleNum);
vector<int> Nodes_in_Period_t(int t0);
map<int, double> Reference_bundle_distribution_factor(double* DFs);