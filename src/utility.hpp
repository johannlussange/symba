#pragma once
#include <vector>
#include <string>
#include <gsl/gsl_matrix.h>

using namespace std;


/// Initialize a brand new RNG state, with unique seed
gsl_rng* make_rng();

/// Produce a vector with N elements [0, N), randomly shuffled.
vector<int> Shuffle(int n);

double BinaryProjection(gsl_matrix* ReflexiveValues, int t, int Tool, int Lag, int Future);

// Trunks any double number with a number "Digits" of Significant digits
double DigitTrunk (double x, int Digits, string FloorOrCeil);

// Generating random variables for specific time delays
vector<double> STLRandom (int Length, string S);

// Generating random variables for specific time delays
vector<int> STLRandomInt (int Length, string S);

vector<vector<double>> GSLMatrixToSTLMatrix (gsl_matrix* M);