#include <iostream>
#include <numeric>
#include <random>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utility.hpp"

using namespace std;


gsl_rng* make_rng() {
    gsl_rng_default_seed = static_cast<unsigned long>(time(NULL)) + rand();
    const gsl_rng_type* T = gsl_rng_ranlux389;
    gsl_rng* r = gsl_rng_alloc (T);
    return r;
}


vector<int> Shuffle(int n) {
    static random_device rd;
    static mt19937 generator(rd());

    vector<int> result(n);
    iota(result.begin(), result.end(), 0);
    shuffle(result.begin(), result.end(), generator);

    return result;
}


double BinaryProjection(gsl_matrix* ReflexiveValues, int t, int Tool, int Lag, int Future) {
    // Fail for Tool=2 and Lag=0 because Start=
    double Result=0; int Start=t-(Lag+1)*Future; if (Start<0) {Start=0;}; int Past=t-Start;
    double Mean0=0; for (int i=Start; i<t-Past/2; i++) {Mean0+=gsl_matrix_get(ReflexiveValues, 0, i)/(Past/2);};
    double Mean1=0; for (int i=t-Past/2; i<t; i++) {Mean1+=gsl_matrix_get(ReflexiveValues, 0, i)/(Past/2);};
    if (Tool==0) {Result=gsl_matrix_get(ReflexiveValues, 0, t)-(Mean1-Mean0);} // counter-trend following
    else if (Tool==1) {Result=0.5*(Mean1+Mean0);} // mean reversion
    else if (Tool==2) {Result=gsl_matrix_get(ReflexiveValues, 0, t)+(Mean1-Mean0);} // trend following
    if (Result<0.001) {Result=0.001;}; // Failsafe
    //if (Result<1.0) {cout << "Pt=" << gsl_matrix_get(ReflexiveValues, 0, t) << ", Mean0=" << Mean0 << ", Mean1=" << Mean1 << ", Tool=" << Tool << ", Future=" << Future << ", Lag=" << Lag << ", Result=" << Result << endl;}; // Failsafe
    return Result;
}


/// Trunks any double number with a number "Digits" of Significant digits
double DigitTrunk(double x, int Digits, string FloorOrCeil) {
    double res=0;
    //ofstream outputDigitTrunk(cli::args::output_dir / "DigitTrunk.txt", ofstream::app);
    /*
     double x=0; double y=0;
     gsl_rng* r = make_rng();
     gsl_rng_set(r, static_cast<unsigned long int>(time(0)/10)); // gsl_rng_set(const gsl_rng* r, unsigned long int s)
     x=100*gsl_rng_uniform(r);
     y=100*gsl_rng_uniform(r);
     outputV << "x=" << x << ", y=" << y << endl;
     outputV.close();
     */
    // Trunk is the interger version of x, and Residue is its decimal rest
    int Trunk=0; double Residue=0;
    if (FloorOrCeil=="Floor") {
        Trunk=int(floor(x));
        Residue=x-Trunk;
        Residue*=pow(10, Digits);
        Residue=floor(Residue);
        Residue/=pow(10, Digits);
        res=Trunk+Residue;
    }
    else if (FloorOrCeil=="Ceil") {
        Trunk=int(ceil(x));
        Residue=x-(Trunk-1);
        Residue*=pow(10, Digits);
        Residue=ceil(Residue);
        Residue/=pow(10, Digits);
        res=Trunk+Residue;
    }
    else if (FloorOrCeil=="Int") {
        Trunk=int(x);
        Residue=x-(Trunk-1);
        Residue*=pow(10, Digits);
        Residue=double(int(Residue));
        Residue/=pow(10, Digits);
        res=Trunk+Residue;
    };
    //outputDigitTrunk.precision(16);
    //outputDigitTrunk << "x=" << x << ", Digits=" << Digits << ", FloorOrCeil=" << FloorOrCeil <<" : result=" << res << endl;
    //outputDigitTrunk.close();
    if (abs(res-ceil(res))>0) {cout << "DIGITTRUNK() ISSUE : res=" << res << endl;};
    return res;
}


vector<double> STLRandom (int Length, string S) {
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng* r, unsigned long int s)
    //gsl_rng* r2 = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r2 = make_rng();
    gsl_rng_set(r2, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng* r, unsigned long int s)
    vector<double> V;
    if (S=="Poisson") {for (int t=0; t<Length; t++) {V.push_back((gsl_ran_poisson (r, 1))*(gsl_ran_poisson (r, 1)));};};
    if (S=="PoissonOriginal") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_poisson (r, 1));};};
    if (S=="PoissonNew") {for (int t=0; t<Length; t++) {V.push_back(abs(10*(1+gsl_ran_gaussian (r, 1)+gsl_ran_poisson (r, 1))));};};
    if (S=="PoissonJump") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_poisson (r, 10));};};
    if (S=="StandardNormal") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_gaussian (r, 1));};};
    if (S=="StandardNormal2") {for (int t=0; t<Length; t++) {V.push_back(gsl_ran_gaussian (r2, 1));};};
    if (S=="NormalNew") {for (int t=0; t<Length; t++) {double x = 0.5+0.2*gsl_ran_gaussian (r, 1); if (x<0.05) {x=gsl_rng_uniform (r);}; if (x>0.95) {x=gsl_rng_uniform (r);}; V.push_back(x);};}; // This is our NEBLearningRate
    if (S=="Uniform") {for (int t=0; t<Length; t++) {V.push_back(gsl_rng_uniform (r));};};
    if (S=="UniformPercent") {for (int t=0; t<Length; t++) {V.push_back(100*(gsl_rng_uniform (r)));};};
    if (S=="UniformPercentFloor") {for (int t=0; t<Length; t++) {V.push_back(floor(100*(gsl_rng_uniform (r))));};};
    if (S=="C") {for (int t=0; t<Length; t++) {V.push_back(rand());};};
    
    gsl_rng_free (r);
    return V;
}


vector<int> STLRandomInt (int Length, string S) {
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); //gsl_rng_set(const gsl_rng* r, unsigned long int s)
    vector<int> V;
    srand (unsigned (time(0)));
    
    if (S=="Uniform") {for (int t=0; t<Length; t++) {V.push_back(int(1000*gsl_rng_uniform (r)));};};
    if (S=="Rand") {for (int t=0; t<Length; t++) {V.push_back(int(rand()));};};
    
    gsl_rng_free (r);
    return V;
};


vector<vector<double>> GSLMatrixToSTLMatrix (gsl_matrix* M) {
    vector<vector<double>> STLM;
    vector<double> TempLine;
    for (int i=0; i<int(M->size1); i++) {
        for (int j=0; j<int(M->size2); j++) {
            //vector<double> Temp; V.push_back(Temp); double x=gsl_matrix_get(M,i,j); if (x!=0) {(V[i]).push_back(x);};
            double x=gsl_matrix_get(M,i,j);
            TempLine.push_back(x);
        };
        STLM.push_back(TempLine);
        TempLine.clear();
    };
    return STLM;
};