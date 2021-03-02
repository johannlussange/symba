#pragma once
#include <vector>

using namespace std;


// A lognormal random walk simulating an asset: dS/S= mu.dt + sig.N(0,1).sqrt(dt) + \mathcal P(k)
vector<double> RandomWalk(double InitialValue, double Drift, double Volatility, int TimeSteps, double Seed);

// Mu is the poisson-length between each jump, and Sig the intensity of the jump
vector<double> PoissonRandomWalk (double S0, int Mu, int Sig, int Time, double Seed);

/// A lognormal random walk that is cointegrated to another
vector<double> CointegratedWalk(vector<double> Master, double Leash, double LeashVolatility, double Accuracy);