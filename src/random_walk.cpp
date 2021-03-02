#include <ctime>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "random_walk.hpp"
#include "utility.hpp"


// This function returns a lognormal random walk simulating an asset: dS/S= mu.dt + sig.N(0,1).sqrt(dt) + \mathcal P(k)
vector<double> RandomWalk(double InitialValue, double Drift, double Volatility, int TimeSteps, double Seed) {
    double Skewness=0.00;
    double FundamentalValue=InitialValue/5; // Asset cannot go below 20% of its fundamental value
    Drift=0.00;
    vector<double> S,S2, StandardNormal, Jumps;
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(Seed)); // gsl_rng_set(const gsl_rng* r, unsigned long int s)
    for (int t=0; t<TimeSteps; t++) {Jumps.push_back(gsl_ran_poisson (r, 10));};
    for (int t=0; t<TimeSteps; t++) {StandardNormal.push_back(Skewness+gsl_ran_gaussian (r, 1));};
    S.push_back(InitialValue); // The +20 is there to ensure we have an entry price that is not negligeable
    S2.push_back(S[0]+FundamentalValue);
    double dt=1;
    int Condition=0;
    for (int t=1; t<TimeSteps; t++) {
        if ((Jumps[t]) > 15) {Condition=1;};
        if ((StandardNormal[t]) >= 0) {Condition *= 1;};
        if ((StandardNormal[t]) < 0) {Condition *= (-1);};
        //S.push_back(floor(S[t-1] + Volatility*(S[t-1])*(StandardNormal[t]) + (S[t-1])*Condition*(Jumps[t])/100));
        S.push_back(ceil((S[t-1])*(1 + Drift*dt + Volatility*(StandardNormal[t])*sqrt(dt) + Condition*(Jumps[t])/100)));
        Condition=0;
        S2.push_back(S[t]+FundamentalValue);
    };
    return S2;
}

// Mu is the poisson-length between each jump, and Sig the intensity of the jump
vector<double> PoissonRandomWalk (double S0, int Mu, int Sig, int Time, double Seed) {
    vector<double> S, VSig; vector<int> VMu;
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(Seed)); // gsl_rng_set(const gsl_rng* r, unsigned long int s)
    for (int t=0; t<Time; t++) {
        S.push_back(S0);
        VMu.push_back(1+gsl_ran_poisson (r, Mu));
        VSig.push_back((1+gsl_ran_poisson (r, Sig)) * (gsl_ran_gaussian (r, 1)));
    };
    //PlotSTLInt(VMu, "VMu.txt"); PlotSTL(VSig, "VSig.txt", "NoXL");
    int P=VMu[0];
    for (int t=1; t<Time; t++) {
        if (P<=0) {P=VMu[t]; S[t] = S[t-1] + VSig[t];}
        else {P-=1; S[t] = S[t-1];};
        if (S[t]<=0.2*S[0]) {S[t] = S[t-1];}; // Stock cannot go below treshold of 20% of initial value
    };
    return S;
}


// This function returns a lognormal random walk that is cointegrated to another
vector<double> CointegratedWalk(vector<double> Master, double Leash, double LeashVolatility, double Accuracy) {
    double Val;
    vector<double> S, L;
    vector<double> StandardNormal = STLRandom((int(Master.size())), "StandardNormal2"); // Time series of size Master.size() made of values gsl_ran_gaussian (r, 1)
    S.push_back(Master[0]);
    L.push_back(Leash);
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng* r, unsigned long int s)
    for (int t=1; t<(int(Master.size())); t++) {
        L.push_back((L[t-1])*(1 + LeashVolatility*(gsl_ran_gaussian (r, 1)))); // The leash is elastic with its own vol
        Val = S[t-1] + LeashVolatility*(S[t-1])*(StandardNormal[t]);
        // Asset comes back if too far
        if ((Master[t] - Val) > ((L[t])*(Master[t]))) {Val=S[t-1] + LeashVolatility*(S[t-1])*(abs(StandardNormal[t]));};
        if ((Val - Master[t]) > ((L[t])*(Master[t]))) {Val=S[t-1] + LeashVolatility*(S[t-1])*(-abs(StandardNormal[t]));};
        if ((gsl_ran_poisson (r, 1))*(gsl_ran_poisson (r, 1)) >= Accuracy) {Val=Master[t]*(1 + 0.10*StandardNormal[t]);}; //This poisson event is computed as happening 0.2, 0.75, 1, 2 per year for Accuracy=13,11,10,9 on average and brings back the cointegrated walk to the Master time series (after a news for example that reveals TruePresentValue)
        S.push_back(Val);
    };
    return S;
}