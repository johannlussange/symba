//#include <stdafx.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <cstdlib>
#include <numeric>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rstat.h>

#include "global.hpp"
#include "utility.hpp"
#include "agent.hpp"
#include "order_book.hpp"
#include "share.hpp"
#include "random_walk.hpp"
#include "cli.hpp"

using namespace std;


// Function plotting an STL vector
void PlotSTL(vector<double> V, string Name, string XL) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform"))[2])); B2=to_string(B);};
    A2+=B2;
    const char* Title2 = A2.c_str();
    ofstream outputV(Title2, ofstream::app);
    for (int t=0; t<int(V.size()); t++) {
        double x=V[t];
        if (XL=="XL") {x=int(x);};
        if (t==(int(V.size()))-1) {outputV << x; break;};
        outputV << x << ",";
    };
    outputV << endl << endl;
    outputV.close();
};




// Function plotting an STL vector
void PlotSTLInt(vector<int> V, string Name) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform"))[2])); B2=to_string(B);};
    A2+=B2;
    const char* Title2 = A2.c_str();
    ofstream outputV(Title2, ofstream::app);
    for (int t=0; t<(int(V.size())); t++) {
        if (t==(int(V.size()))-1) {outputV << V[t]; break;};
        outputV << V[t] << endl;
    };
    outputV << endl << endl;
    outputV.close();
};




// Function plotting an STL vector
void PlotSTLMarket(vector<double> V, string Name) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform"))[2])); B2=to_string(B);};
    A2+=B2;
    const char* Title2 = A2.c_str();
    ofstream outputV(Title2, ofstream::app);
    for (int t=0; t<(int(V.size())); t++) {
        int x = int(V[t]);
        if (t==(int(V.size()))-1) {outputV << x; break;};
        outputV << x << endl;
    };
    outputV << endl << endl;
    outputV.close();
};




// Function plotting a GSL matrix
void PlotGSLMatrix(gsl_matrix* M, string Name, int F) {
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform"))[2])); B2=to_string(B);};
    A2+=B2;
    const char* Title2 = A2.c_str();
    ofstream outputM(Title2, ofstream::app);
    for (int t=0; t<int(M->size2); t++) {
        for (int j=0; j<int(M->size1); j++) {
            //int x = int(floor(F*gsl_matrix_get(M,j,t)));
            double x = F*gsl_matrix_get(M,j,t);
            if (j==(int(M->size1))-1) {outputM << fixed << x << endl; continue;};
            outputM << fixed << x << "\t";
        };
    };
    outputM << endl << endl;
    outputM.close();
};




// This compute the mean and variance from an STL vector
vector<double> StaticMoments (vector<double> V, int Start, int End, int Delta) {
    vector<double> Result;
    double Mean=0;
    double Variance=0;
    double Skewness=0;
    double Kurtosis=0;
    double RealizedVariance=0;
    double RealizedVolatility=0;
    double SSXlag=0; double SSXX=0; double SSXX2=0;
    //int SeriesSize=int(V.size()); // Original
    //if (End>=SeriesSize) {End=SeriesSize-1;}; // New
    if (Start>=End) {Start=End;}; // Original
    //if (Start<0) {Start=0;}; // New
    
    for (int i=Start; i<=End; i++) {
        Mean+=V[i];
        if (i-Delta>=0) {
            RealizedVariance+=(log(V[i]/V[i-Delta]))*(log(V[i]/V[i-Delta]));
            RealizedVolatility+=(log(V[i]/V[i-Delta]))*(log(V[i]/V[i-Delta]));
        };
    };
    Mean=Mean/(End-Start+1); // Mean computation
    RealizedVolatility=sqrt(252*RealizedVolatility/(End-Start+1));
    for (int i=Start; i<=End; i++) {
        double x=V[i] - Mean;
        Variance+=x*x;
        Skewness+=x*x*x;
        Kurtosis+=x*x*x*x;
    };
    Variance/=End-Start+1; // Variance computation
    Skewness=(Skewness/(End-Start+1))*(1/(Variance*sqrt(Variance)));
    Kurtosis=(Kurtosis/(End-Start+1))*(1/(Variance*Variance));
    
    for (int i=Start; i<=End-Delta; i++) {
        double x=V[i] - Mean;
        double y=V[i+Delta] - Mean;
        SSXlag+=x*y;
        SSXX+=x*x;
        SSXX2+=y*y;
    };
    double AutoCorr=SSXlag/((sqrt(SSXX))*(sqrt(SSXX2))); // Function returning the p-lag autocorrelation of a gsl_matrix, between time steps begin and end
    
    Result.push_back (Mean);
    Result.push_back (Variance);
    Result.push_back (Skewness);
    Result.push_back (Kurtosis); // Distributions with kurtosis<3 (>3) are platykurtic (leptokurtic)
    Result.push_back (sqrt(Variance));
    Result.push_back (RealizedVolatility);
    Result.push_back (AutoCorr);
    return Result;
};





// This computes the mean and variance for a given lag at each time step of a given time series
vector<vector<double>> DynamicMoments (vector<double> V, int Start, int End, int Lag, int Delta) {
    vector<vector<double>> Result;
    vector<double> ResultMean, ResultVariance, ResultSkewness, ResultKurtosis, ResultStdDev, TimeSeriesPlus3Sigma, TimeSeriesMinus3Sigma, ResultRealizedVolatility, ResultAutoCorr, ResultSquarredLogReturns1, ResultSquarredLogReturns2, ResultSquarredLogReturns3;
    if (End >= int(V.size())) {End=int(V.size()) -1;}; // Original
    Lag-=1;
    //if (Start>=End) {Start=End;}; // New
    //if (Start<0) {Start=0;}; // New
    
    double SquarredLogReturns1=0;
    double SquarredLogReturns2=0;
    double SquarredLogReturns3=0;
    for (int i=End; i>=Start; i--) {
        if (i<Start+Lag) {Lag-=1;};
        vector<double> SM = StaticMoments(V, i-Lag, i, Delta);
        double m = SM[0];
        double v = SM[1];
        double s = SM[2];
        double k = SM[3];
        double sd = SM[4];
        double rv = SM[5];
        double ac = SM[6];
        if (i-1>=0) {SquarredLogReturns1=(log(V[i]/V[i-1]))*(log(V[i]/V[i-1]));};
        if (i-Week>=0) {SquarredLogReturns2=(log(V[i]/V[i-Week]))*(log(V[i]/V[i-Week]));};
        if (i-2*Week>=0) {SquarredLogReturns3=(log(V[i]/V[i-2*Week]))*(log(V[i]/V[i-2*Week]));};
        ResultMean.push_back(m);
        ResultVariance.push_back(v);
        ResultSkewness.push_back(s);
        ResultKurtosis.push_back(k);
        ResultStdDev.push_back(sd);
        ResultRealizedVolatility.push_back(rv);
        ResultAutoCorr.push_back(ac);
        ResultSquarredLogReturns1.push_back(SquarredLogReturns1);
        ResultSquarredLogReturns2.push_back(SquarredLogReturns2);
        ResultSquarredLogReturns3.push_back(SquarredLogReturns3);
    };
    reverse(ResultMean.begin(), ResultMean.end());
    reverse(ResultVariance.begin(), ResultVariance.end());
    reverse(ResultSkewness.begin(), ResultSkewness.end());
    reverse(ResultKurtosis.begin(), ResultKurtosis.end());
    reverse(ResultStdDev.begin(), ResultStdDev.end());
    reverse(ResultRealizedVolatility.begin(), ResultRealizedVolatility.end());
    reverse(ResultAutoCorr.begin(), ResultAutoCorr.end());
    reverse(ResultSquarredLogReturns1.begin(), ResultSquarredLogReturns1.end());
    reverse(ResultSquarredLogReturns2.begin(), ResultSquarredLogReturns2.end());
    reverse(ResultSquarredLogReturns3.begin(), ResultSquarredLogReturns3.end());
    
    for (int i=Start; i<=End; i++) {
        TimeSeriesPlus3Sigma.push_back(V[i]+3*ResultStdDev[i-Start]);
        TimeSeriesMinus3Sigma.push_back(V[i]-3*ResultStdDev[i-Start]);
    };
    
    Result.push_back(ResultMean);
    Result.push_back(ResultVariance);
    Result.push_back(ResultSkewness);
    Result.push_back(ResultKurtosis);
    Result.push_back(ResultStdDev);
    Result.push_back(TimeSeriesPlus3Sigma);
    Result.push_back(TimeSeriesMinus3Sigma);
    Result.push_back(ResultRealizedVolatility);
    Result.push_back(ResultAutoCorr);
    Result.push_back(ResultSquarredLogReturns1);
    Result.push_back(ResultSquarredLogReturns2);
    Result.push_back(ResultSquarredLogReturns3);
    
    return Result;
};




//  AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={3, w, 2w, 3w, m}
double MALAutoCorrelation (gsl_matrix* X, int t, int p, int Stock) {
    double Mean1=0; double Mean2=0; double SSXY=0; double SSXX=0; double SSYY=0;
    int Start1=max(t-p, 0); int Denominator1=max(t-Start1+1, 1);
    int Start2=max(t-2*p, 0); int Denominator2=max(Start1-Start2+1, 1);
    for (int i=Start1; i<=t; i++) {Mean1+=gsl_matrix_get (X, Stock, i)/Denominator1;};
    for (int i=Start2; i<=Start1; i++) {Mean2+=gsl_matrix_get (X, Stock, i)/Denominator2;};
    for (int i=Start1; i<=t; i++) {
        double x=gsl_matrix_get (X, Stock, i) - Mean1;
        double y=gsl_matrix_get (X, Stock, max(i-p, 0)) - Mean2;
        SSXX+=x*x;
        SSYY+=y*y;
        SSXY+=x*y;
    };
    double Result=SSXY/((sqrt(SSXX))*(sqrt(SSYY))); // Autocorrelation between intervals [t-p, t] and [t-2p, t-p]
    if (sqrt(SSXX)*sqrt(SSYY)==0) {Result=0;};
    if ((t-p<0) || (!Result)) {Result=0;}; // Failsafe
    return Result;
};


//  AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 3, w, 2w, 3w, m}
double MALAutoCorrelationBlend (gsl_matrix* X, int t, int p, int Shift, int Stock) {
    double Mean1=0; double Mean2=0; double SSXY=0; double SSXX=0; double SSYY=0;
    int Start1=max(t-p, 0); int Start2=max(t-p-Shift, 0);
    int Denominator=max(t-Start1+1, 1);
    for (int i=Start1; i<=t; i++) {Mean1+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start2; i<=max(t-Shift, 0); i++) {Mean2+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start1; i<=t; i++) {
        double x=gsl_matrix_get (X, Stock, i) - Mean1;
        double y=gsl_matrix_get (X, Stock, max(i-Shift, 0)) - Mean2;
        SSXX+=x*x;
        SSYY+=y*y;
        SSXY+=x*y;
    };
    double Result=SSXY/((sqrt(SSXX))*(sqrt(SSYY))); // Autocorrelation between intervals [t-p, t] and [t-2p, t-p]
    if (sqrt(SSXX)*sqrt(SSYY)==0) {Result=0;};
    if ((t-p<0) || (!Result)) {Result=0;}; // Failsafe
    return Result;
};


double MALVolatility (gsl_matrix* X, int t, int p, int Stock) {
    double Mean=0; double Variance=0;
    int Start=max(t-p, 0); int Denominator=max(t-Start+1, 1);
    for (int i=Start; i<=t; i++) {Mean+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start; i<=t; i++) {double x=gsl_matrix_get (X, Stock, i) - Mean; Variance+=x*x;};
    double Result=sqrt(Variance/Denominator)/Mean; // Standard deviation on the interval [t-p, t]
    if ((t-p<0) || (!Result)) {Result=0;}; // Failsafe
    //if (p==2*Week) {cout << "Denominator=" << Denominator << endl;};
    return Result;
};


double MALVolatility2 (gsl_matrix* X, int t1, int t2, int Stock) {
    double Mean=0; double Variance=0;
    int Start=max(t1, 0); int Denominator=max(t2-Start+1, 1);
    for (int i=Start; i<=t2; i++) {Mean+=gsl_matrix_get (X, Stock, i)/Denominator;};
    for (int i=Start; i<=t2; i++) {double x=gsl_matrix_get (X, Stock, i) - Mean; Variance+=x*x;};
    double Result=sqrt(Variance/Denominator)/Mean; // Standard deviation on the interval [t-p, t]
    if ((t1<0) || (!Result)) {Result=0;}; // Failsafe
    //if (p==2*Week) {cout << "Denominator=" << Denominator << endl;};
    return Result;
};




double StaticAutoCorrelation (gsl_matrix* X, int Start, int End, int p, int Stock, int j) { // BBB
    double Mean=0; double SSXlag=0; double SSXX=0; double SSXX2=0;
    for (int i=Start; i<=End; i++) {Mean+=gsl_matrix_get (X, Stock+8*j, i);}; Mean=Mean/(End-Start+1); // Mean of MarketBidAskTrue!!!
    for (int i=Start; i<=End-p; i++) {
        double x=gsl_matrix_get (X, Stock+8*j, i) - Mean;
        double y=gsl_matrix_get (X, Stock+8*j, i+p) - Mean;
        SSXlag+=x*y;
        SSXX+=x*x;
        SSXX2+=y*y;
    };
    double Result=SSXlag/((sqrt(SSXX))*(sqrt(SSXX2))); // p-lag autocorrelation of a gsl_matrix, between Start and End
    return Result;
};






gsl_matrix* DynamicAutoCorrelation (gsl_matrix* X, int Lag, int Stock, int j) {
    gsl_matrix* Result = gsl_matrix_alloc (5, int(X->size2));
    for (int i=0; i<int(X->size2)-Lag; i++) {
        gsl_matrix_set (Result, 0+6*j, i, StaticAutoCorrelation (X, i, i+Lag, DayTicks, Stock, j));
        gsl_matrix_set (Result, 1+6*j, i, StaticAutoCorrelation (X, i, i+Lag, 2*DayTicks, Stock, j));
        gsl_matrix_set (Result, 2+6*j, i, StaticAutoCorrelation (X, i, i+Lag, Week, Stock, j));
        gsl_matrix_set (Result, 3+6*j, i, StaticAutoCorrelation (X, i, i+Lag, 2*Week, Stock, j));
        gsl_matrix_set (Result, 4+6*j, i, StaticAutoCorrelation (X, i, i+Lag, Month, Stock, j));
    }
    return Result;
};






gsl_matrix* ProgressiveAutoCorrelation (gsl_matrix* X, int Stock, int j) {
    gsl_matrix* Result = gsl_matrix_calloc (5, int(X->size2));
    for (int i=1; i<int(X->size2); i++) {
        if (i>DayTicks) {gsl_matrix_set (Result, 0+6*j, i, StaticAutoCorrelation (X, 0, i, DayTicks, Stock, j));};
        if (i>2*DayTicks) {gsl_matrix_set (Result, 1+6*j, i, StaticAutoCorrelation (X, 0, i, 2*DayTicks, Stock, j));};
        if (i>Week) {gsl_matrix_set (Result, 2+6*j, i, StaticAutoCorrelation (X, 0, i, Week, Stock, j));};
        if (i>2*Week) {gsl_matrix_set (Result, 3+6*j, i, StaticAutoCorrelation (X, 0, i, 2*Week, Stock, j));};
        if (i>Month) {gsl_matrix_set (Result, 4+6*j, i, StaticAutoCorrelation (X, 0, i, Month, Stock, j));};
    }
    return Result;
};






gsl_matrix* ExtensiveAutoCorrelation (gsl_matrix* X, int Stock, int j) {
    gsl_matrix* Result = gsl_matrix_calloc (5, int(X->size2));
    for (int i=1; i<int(X->size2); i++) {
        gsl_matrix_set (Result, 0+6*j, i, StaticAutoCorrelation (X, 0, int(X->size2)-1, i, Stock, j));
    }
    return Result;
};






double StaticCorrelation (gsl_matrix* X, gsl_matrix* Y, int XStart, int XEnd, int p, int Stock, int j) {
    double XMean=0; double YMean=0; double SSXY=0; double SSXX=0; double SSYY=0;
    for (int i=XStart; i<=XEnd; i++) {XMean+=gsl_matrix_get (X, Stock+6*j, i)/(XEnd-XStart+1);}; // Mean of MarketBidAskTrue!!!
    for (int i=XStart+p; i<=XEnd+p; i++) {YMean+=gsl_matrix_get (Y, Stock+6*j, i)/(XEnd-XStart+1);}; // Mean of MarketBidAskTrue!!!
    for (int i=XStart; i<=XEnd; i++) {
        double x=gsl_matrix_get (X, Stock+6*j, i) - XMean;
        double y=gsl_matrix_get (Y, Stock+6*j, i+p) - YMean;
        SSXY+=x*y;
        SSXX+=x*x;
        SSYY+=y*y;
    };
    return SSXY/((sqrt(SSXX))*(sqrt(SSYY)));
};





gsl_matrix* CFMCorrelation (gsl_matrix* X, gsl_matrix* Y, int XStart, int XEnd, int Lag, int Stock, int j) {
    gsl_matrix* Result = gsl_matrix_calloc (1, Lag);
    for (int i=0; i<Lag; i++) {gsl_matrix_set (Result, 0, i, StaticCorrelation (X, Y, XStart, XEnd, i, Stock, j));} // Be careful that the Lag+XEnd <= int(X->size2) !!!
    return Result;
};






// Function returning the p-lag autocorrelation of a gsl_matrix, between time steps begin and end
gsl_matrix* AutocorrelationLagP (gsl_matrix* X) {
    gsl_matrix* Res=gsl_matrix_alloc (int(X->size1), int(X->size2));
    for (int j=0; j<int(X->size1); j++) {
        for (int p=0; p<int(X->size2); p++) {
            double SSXlag=0; double SSXX1=0; double SSXX2=0; double Mean=0;
            
            for (int t=0; t<int(X->size2); t++) {Mean+=gsl_matrix_get (X, j, t);};
            Mean/=int(X->size2); // Mean computation
            
            for (int t=0; t<int(X->size2)-p; t++) {
                double x=gsl_matrix_get (X, j, t) - Mean;
                double y=gsl_matrix_get (X, j, t+p) - Mean;
                SSXlag+=x*y;
                SSXX1+=x*x;
                SSXX2+=y*y;
            }; // closes t loop
            gsl_matrix_set (Res, j, p, SSXlag/((sqrt(SSXX1))*(sqrt(SSXX2)))); // Function returning the p-lag autocorrelation
        }; // closes p loop
    }; // closes j loop
    
    return Res;
}








// Function returning the Pearson correlation between two vectors
double StaticCorrelation (vector<double> X, vector<double> Y, int Start, int End) {
    double Res; double SSXY=0; double SSXX=0; double SSYY=0; double MeanX=0; double MeanY=0;
    int Size=min(int(X.size()), int(Y.size()));
    if (Start<0) {Start=0;};
    if (End>Size) {End=Size;};
    if (End==0) {End=Size;};
    for (int i=Start; i<=End; i++) {MeanX+=X[i]; MeanY+=Y[i];}; MeanX/=End-Start+1; MeanY/=End-Start+1; // Mean
    
    for (int i=Start; i<=End; i++) {SSXY+=(X[i] - MeanX)*(Y[i] - MeanY);}
    for (int i=Start; i<=End; i++) {SSXX+=(X[i] - MeanX)*(X[i] - MeanX);}
    for (int i=Start; i<=End; i++) {SSYY+=(Y[i] - MeanY)*(Y[i] - MeanY);}
    Res=SSXY/(sqrt(SSXX*SSYY));
    
    return Res;
};








vector<double> DynamicCorrelation (vector<double> X, vector<double> Y, int Lag) {
    vector<double> Res;
    int Size=min(int(X.size()), int(Y.size()));
    for (int k=0; k<Lag; k++) {Res.push_back(0);};
    for (int k=0; k<Size-Lag; k++) {
        Res.push_back(StaticCorrelation(X,Y,k,k+Lag));
    };
    
    return Res;
};

// KEVIN
// Fonction A(i) dans la log-likelihood
void function_A(vector<double>* t ,double beta, vector<double>& A ){
    size_t n = (*t).size();
    A[0] = 0;
    for (size_t i = 1; i < n; i++) {A[i] = A[i-1]*exp(-beta*( (*t)[i] - (*t)[i-1] )) + exp(-beta*( (*t)[i] - (*t)[i-1] ));}
};

// Fonction log-likelihood
double function_f(const gsl_vector* v, void* params){
    double alpha = gsl_vector_get(v,0);
    double beta = gsl_vector_get(v,1);
    double mu = gsl_vector_get(v,2);
    vector<double>* t = (vector<double>*) params;
    size_t n = (*t).size();
    vector<double> A (n);
    function_A(t,beta,A);
    const double b = alpha/beta;
    double r = mu*(*t)[n-1]-b*(A[n-1] + 1);
    for (size_t i = 0; i < n; i++) {r = r -log(mu + alpha*A[i]) + b;}
    return r;
};

// Fonction B(i) dans la log-likelihood
void function_B( vector<double>* t,double beta, vector<double>& B){
    size_t n = (*t).size();
    B[0] = 0;
    for (size_t i = 1; i < n; i++) {
        double tmp = 0;
        for (size_t j = 0; j < i; j++) {tmp = tmp + ( (*t)[i] - (*t)[j])*exp(-beta*((*t)[i]-(*t)[j]));}
        B[i] = tmp;
    };
};

// Gradient de la log-likelihood
void function_df(const gsl_vector* v, void* params, gsl_vector* df){
    double alpha = gsl_vector_get(v,0);
    double beta = gsl_vector_get(v,1);
    double mu = gsl_vector_get(v,2);
    vector<double>* t = (vector<double>*) params;
    size_t n = (*t).size();
    vector<double> A (n);
    vector<double> B (n);
    function_A(t, beta, A);
    function_B(t, beta, B);
    double dmu = (*t)[n-1];
    double dalpha = 0;
    double dbeta = 0;
    for (size_t i = 0; i < n; i++) {dmu = dmu - 1/( mu + alpha*A[i]);}
    for (size_t i = 0; i < n; i++) {dalpha = dalpha - (A[i]/(mu + alpha*A[i])) + (1/beta);}
    dalpha = dalpha - (1/beta)*(A[n-1] + 1 );
    const double C1 = alpha/beta;
    const double C2 = alpha/pow(beta,2);
    for (size_t i = 0; i < n; i++) {dbeta = dbeta + (alpha*B[i]/(mu + alpha*A[i])) - C2 ;}
    dbeta = dbeta + C1*B[n-1] + C2*A[n-1]+C2;
    gsl_vector_set(df,0,dalpha);
    gsl_vector_set(df,1,dbeta);
    gsl_vector_set(df, 2,dmu);
}

// Fonction fdf, calcule la fonction f et le gradient de f
void function_fdf(const gsl_vector* v, void* params, double* f, gsl_vector* df) {
    *f = function_f(v,params);
    function_df(v, params, df );
}

// Minimizer de la log-likelihood
int EstimationOzaki(std::vector<double>* t, double* Init, const double stepsize, const size_t NbIteration, double* R){
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = 0; // NULL
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    size_t iter = 0;
    int status;
    double size;
    // Starting point
    x = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, Init[0]);
    gsl_vector_set (x, 1, Init[1]);
    gsl_vector_set (x, 2, Init[2]);
    // Set initial step sizes to 1
    ss = gsl_vector_alloc (3);
    gsl_vector_set_all (ss, stepsize);
    // Initialize method and iterate
    minex_func.n = 3;
    minex_func.f = &function_f;
    minex_func.params =  t;
    s = gsl_multimin_fminimizer_alloc (T, 3);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 0.0003);
    }
    while ((status == GSL_CONTINUE) && (iter < NbIteration));
    R[0] = gsl_vector_get(s->x,0);
    R[1] = gsl_vector_get(s->x,1);
    R[2] = gsl_vector_get(s->x,2);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    return status;
}

// Calcul de lambda(t)
double ExpIntensityState(const double mu, const double alpha, const double beta,std::vector<double>* t, const double rescale ,const size_t n){
    double res = mu;
    const double T = (double) n/rescale;
    for (vector<double>::iterator it = (*t).begin(); (it != (*t).end()) && (*it <= T); it++) {res = res + alpha*exp(-beta*(T-*it));};
    return res;
}

// Fonction d'integration
void IntegrateExpIntensity(const double mu, const double alpha, const double beta, std::vector<double>* t, vector<double>& I){
    size_t n = (*t).size();
    vector<double> A(n);
    function_A(t,beta,A);
    vector<double> tmp;
    I.push_back(0);
    for (size_t i = 1; i < n; i++) {
        double res = mu*(t[0])[i];
        res = res - (alpha/beta)*A[i] + i;
        tmp.push_back(res);
    };
    for (size_t i = 1; i < n; i++) {
        I.push_back(tmp[i]-tmp[i-1]);
    };
};

// Fonction retournant l'intensité du processus de hawkes; t = les closing prices dans l'ordre chronologique ; pct = seuil (0.1); NbIteration = le nombre d'iteration du minimiseur; R = les parametres du processus de Hawkes; I = goodness of fit
vector<double>* Intensity(std::vector<double>* t, const double rescale, const double pct , const size_t NbIteration, double* R, vector<double>& I){
    size_t n = (*t).size();
    std::vector<double>* res = new vector<double>[2]; // JJJ10
    std::vector<double> jump;
    
    vector<double> logr;
    logr.reserve(n-1);

    for (size_t i = 0; i < n-1; i++) {
        logr.push_back(((*t)[i+1]/(*t)[i]));
    }

    gsl_rstat_quantile_workspace* quantile1 =  gsl_rstat_quantile_alloc(pct);
    gsl_rstat_quantile_workspace* quantile2 = gsl_rstat_quantile_alloc(1-pct);
    for (size_t i = 0; i < n-1; i++) {
        gsl_rstat_quantile_add(logr[i], quantile1);
        gsl_rstat_quantile_add(logr[i], quantile2);
    };
    double bsup = gsl_rstat_quantile_get(quantile2);
    double binf = gsl_rstat_quantile_get(quantile1);
    for (size_t i = 1; i < n-1; i++) {
        double tmp = logr[i];
        if ( tmp > bsup ) {jump.push_back((i-1)/rescale);};
        if ( tmp < binf) {jump.push_back((i-1)/rescale);};
    };
    double Init[3];
    Init[0] = 0.6;
    Init[1] = 0.4;
    Init[2] = 0.5;
    double stepsize = 0.001;
    EstimationOzaki(&jump, Init, stepsize, NbIteration, R);
    for (size_t i = 0; i < n-1; i++) {
        double tmp = ExpIntensityState( R[2], R[0], R[1], &jump, rescale, i );
        (res[0]).push_back(tmp);
        (res[1]).push_back(i/rescale);
    };
    IntegrateExpIntensity(R[2], R[0], R[1], &jump, I);
    return res;
}

// Hawkes process historical intensity
gsl_matrix* HistoricalIntensity (gsl_matrix* TimeSeries, int k) {
    int T = int(TimeSeries->size2);
    gsl_matrix* Res = gsl_matrix_calloc (1, T);
    vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (TimeSeries, k, t)); GoF.push_back(1);};
    vector<double>* pV= &V;
    const double Rescale=10; // Rescale factor
    const double Threshold=0.1; // Threshold (0.1)
    const size_t NumIteration=100000; // Number of iterations of the minimizer (stops before once a minimum is found)
    double R; double* pR = &R; // Hawkes process parameters (does not need to be initialized)
    vector<double>& pGoF = GoF; // Goodness of fit
    vector<double>* pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
    for (int t=0; t<T; t++) {gsl_matrix_set (Res, 0, t, (*pRes)[t]);}; // Hawkes process intensity
    gsl_matrix_set (Res, 0, T-1, gsl_matrix_get (Res, 0, T-2));
    delete[] pRes; // JJJ10 void gsl_matrix_free(gsl_matrix* m)
    return Res;
}; // closes HistoricalIntensity()

// Hawkes process spot intensity
gsl_matrix* SpotIntensity (gsl_matrix* TimeSeries, int k) {
    gsl_matrix* Res = gsl_matrix_calloc (1, int(TimeSeries->size2));
    for (int T=5; T<int(TimeSeries->size2); T++) {
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (TimeSeries, k, t)); GoF.push_back(1);};
        vector<double>* pV= &V;
        const double Rescale=10; // Rescale factor
        const double Threshold=0.1; // Threshold (0.1)
        const size_t NumIteration=100000; // Number of iterations of the minimizer (stops before once a minimum is found)
        double R; double* pR = &R; // Hawkes process parameters (does not need to be initialized)
        vector<double>& pGoF = GoF; // Goodness of fit
        vector<double>* pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
        gsl_matrix_set (Res, 0, T, (*pRes)[T-2]); // Hawkes process intensity
    }; // closes T loop
    return Res;
}; // closes SpotIntensity()

// Hawkes process spot intensity with first 100 time steps at mean
gsl_matrix* SpotIntensity2 (gsl_matrix* TimeSeries, int k) {
    gsl_matrix* Res = gsl_matrix_calloc (1, int(TimeSeries->size2));
    for (int T=5; T<int(TimeSeries->size2); T++) {
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (TimeSeries, k, t)); GoF.push_back(1);};
        vector<double>* pV= &V;
        const double Rescale=10; // Rescale factor
        const double Threshold=0.1; // Threshold (0.1)
        const size_t NumIteration=1000000; // Number of iterations of the minimizer (stops before once a minimum is found)
        double R; double* pR = &R; // Hawkes process parameters (does not need to be initialized)
        vector<double>& pGoF = GoF; // Goodness of fit
        vector<double>* pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
        gsl_matrix_set (Res, 0, T, (*pRes)[T-2]); // Hawkes process intensity
    }; // closes T loop
    double Mean=0;
    for (int T=100; T<int(TimeSeries->size2); T++) {Mean+=gsl_matrix_get (Res, 0, T)/(int(TimeSeries->size2)-100);};
    for (int T=0; T<100; T++) {gsl_matrix_set (Res, 0, T, Mean);};
    return Res;
}; // closes SpotIntensity2()



// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***
// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***
// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***
// ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE *** ULM ALGORITHM CODE ***


// Takes as input the file constructed by ImportMacro (). Then the file must be updated by hand as such:
// LIBOR rates taken from the FED in StLouis https://fred.stlouisfed.org/categories/33003/downloaddata
// Exchange rates taken from the Bank of England at http://www.bankofengland.co.uk/boeapps/iadb/Rates.asp?TD=12&TM=Dec&TY=2017&
gsl_matrix* Macroeconomics () {
    double FXFee=1-0.3*0.01; // FX fees from IG
    string Path=Machine + "Symba/CSV/MACROECONOMICS.csv";
    ifstream FileInput(Path); string Line; int LineNb=0;
    vector<double> V1, V2, V3, V4, V5;
    while (getline (FileInput, Line)) {LineNb++;
        istringstream LineStream(Line); string Item; int ItemNb=0;
        while (getline (LineStream, Item, ',')) {ItemNb++;
            if ((ItemNb==1) && (LineNb>1)) {V1.push_back(atof(Item.c_str()));};
            if ((ItemNb==2) && (LineNb>1)) {V2.push_back(atof(Item.c_str()));};
            if ((ItemNb==3) && (LineNb>1)) {V3.push_back(atof(Item.c_str()));};
            if ((ItemNb==4) && (LineNb>1)) {V4.push_back(atof(Item.c_str()));};
            if ((ItemNb==5) && (LineNb>1)) {V5.push_back(atof(Item.c_str()));};
        };
    };
    int Size=int(V1.size()); FileInput.close();
    gsl_matrix* Result = gsl_matrix_calloc (6, Size);
    for (int i=0; i<Size; i++) {
        gsl_matrix_set (Result, 0, i, V1[i]); // Dates
        gsl_matrix_set (Result, 1, i, V2[i]); // LIBOR in percent and based on USD
        gsl_matrix_set (Result, 2, i, V3[i]); // LIBOR in percent and based on GBP
        gsl_matrix_set (Result, 3, i, V4[i]*FXFee); // 1 EUR in USD
        gsl_matrix_set (Result, 4, i, V5[i]*FXFee); // 1 EUR in GBP
        gsl_matrix_set (Result, 5, i, 1); // No currency conversion
    };
    /*
     for (int i=1; i<Size; i++) {
     if ((abs(gsl_matrix_get(Result, 1, i))>0.1) || (abs(gsl_matrix_get(Result, 1, i))==0)) {gsl_matrix_set(Result, 1, i, gsl_matrix_get(Result, 1, i-1));};
     if ((abs(gsl_matrix_get(Result, 2, i))>0.1) || (abs(gsl_matrix_get(Result, 2, i))==0)) {gsl_matrix_set(Result, 2, i, gsl_matrix_get(Result, 2, i-1));};
     if ((abs(gsl_matrix_get(Result, 3, i))>2) || (abs(gsl_matrix_get(Result, 3, i))==0)) {gsl_matrix_set(Result, 3, i, gsl_matrix_get(Result, 3, i-1));};
     if ((abs(gsl_matrix_get(Result, 4, i))>2) || (abs(gsl_matrix_get(Result, 4, i))==0)) {gsl_matrix_set(Result, 4, i, gsl_matrix_get(Result, 4, i-1));};
     };
     */
    //PlotGSLMatrix (Result, "Macroeconomics.csv", 1);
    return Result;
};





// Barclay Equity L/S Index between 1.1.2007 and 1.1.2018 (132 months, one month being 21.75=22 days)
void StockIndices () {
    int Size=2872;
    string Path=Machine + "Symba/CSV/BarclayHFI.csv"; vector<double> V;
    ifstream FileInput(Path); string Line; int LineNb=0;
    while (getline (FileInput, Line)) {LineNb++;
        istringstream LineStream(Line); string Item; int ItemNb=0;
        while (getline (LineStream, Item, '\r')) {ItemNb++; V.push_back(atof(Item.c_str()));};
    };
    // Output on 132 months
    gsl_matrix* Result = gsl_matrix_calloc (1, 132); gsl_matrix_set (Result, 0, 0, 1);
    for (int i=1; i<132; i++) {gsl_matrix_set (Result, 0, i, gsl_matrix_get (Result, 0, i-1)*(V[i-1]+100)/100.0);}
    // Output on 2520 days
    gsl_matrix* Result2 = gsl_matrix_calloc (2, Size);
    for (int i=0; i<Size; i++) {gsl_matrix_set (Result2, 0, i, gsl_matrix_get (Result, 0, i/22));}
    
    Path=Machine + "Symba/CSV/NasdaqComposite.csv"; vector<double> W;
    ifstream FileInput2(Path); string Line2; LineNb=0;
    while (getline (FileInput2, Line2)) {LineNb++;
        istringstream LineStream(Line2); string Item; int ItemNb=0;
        while (getline (LineStream, Item, '\r')) {ItemNb++; W.push_back(atof(Item.c_str()));};
    };
    vector<int> J = Shuffle(int(W.size())-10);
    int Diff=Size-int(W.size());
    vector<double> N;
    for (int i=0; i<int(W.size()); i++) {
        N.push_back(W[i]);
        for (int k=0; k<Diff; k++) {
            if (J[k]==i) {N.push_back(W[i]);};
        };
    }
    for (int i=0; i<Size; i++) {
        gsl_matrix_set (Result2, 1, i, N[i]/N[0]);
    }
    PlotGSLMatrix (Result2, "StockIndices.csv", 1); // Column 0 is Barclay, column 1 is Nasdaq
};


vector<Share> PortfolioGenerator (string Name, int FirstDate) {
    string ProcessedData="Yes";
    vector<Share> PF, PFFinal; int TickStart=int(time(0));
    string DPath=Machine + "Symba/CSV/" + Name + ".txt";
    ifstream DFileInput(DPath);string DLine; vector<string> VExchanges, VCountries, VCurrencies;
    while (getline (DFileInput, DLine, '\r')) {
        istringstream LineStream(DLine); string DItem; int DItemNb=0;
        while (getline (LineStream, DItem, ',')) {DItemNb++;
            if (DItemNb==1) {VExchanges.push_back(DItem);}
            else if (DItemNb==2) {VCountries.push_back(DItem);}
            else if (DItemNb==3) {VCurrencies.push_back(DItem);};
        };
    };
    int ExchangeSize=int(VExchanges.size()); DFileInput.close();
    gsl_matrix* Macro = Macroeconomics();
    for (int m=0; m<ExchangeSize; m++) {
        string FPath=Machine + "Symba/CSV/" + VExchanges[m] + "new.txt";
        vector<string> VFiles, VSymbol, VTitles;
        if (ProcessedData=="No") {FPath=Machine + "Symba/CSV/" + VExchanges[m] + ".txt";};
        ifstream FFileInput(FPath); string FLine;
        if (ProcessedData=="No") {
            while (getline (FFileInput, FLine, '\r')) {VFiles.push_back(FLine);};
        }
        else if (ProcessedData=="Yes") {
            while (getline (FFileInput, FLine, '\r')) {
                istringstream LineStream(FLine); string DItem; int DItemNb=0;
                while (getline (LineStream, DItem, ',')) {DItemNb++;
                    if (DItemNb==1) {VFiles.push_back(DItem);}
                    else if (DItemNb==2) {VSymbol.push_back(DItem);}
                    else if (DItemNb==3) {VTitles.push_back(DItem);};
                };
            };
        };
        FFileInput.close();
        
        int StockSize=int(VFiles.size());
        for (int s=0; s<StockSize; s++) {
            Share X; X.Exchange=VExchanges[m]; X.File=VFiles[s]; X.Country=VCountries[m]; X.Currency=VCurrencies[m];
            if (ProcessedData=="Yes") {X.Symbol=VSymbol[s]; X.Title=VTitles[s];}
            else if (ProcessedData=="No") {X.ProcessBrokerData();};
            //cout << X.InPF << endl;
            if ((X.InPF=="In") && (ProcessedData=="No")) {
                X.Gen(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            if (ProcessedData=="Yes") {
                X.Gen(Machine + "Symba/CSV/" + VExchanges[m] + "/" + VFiles[s], Macro, FirstDate);
                PF.push_back(X);
                cout << VExchanges[m] << " (" << floor(s*100.0/StockSize) << "%): " << X.Symbol << ", " << X.Title << " (" << X.Country << ", " << X.Currency << ")" << endl;
            };
            
        }; // closes s loop
    }; // closes m loop
    
    // Considering older shares since Jan 2007 only
    for (int s=0; s<int(PF.size()); s++) {
        int Count=0;
        for (int t=0; t<int(PF[s].Data->size2); t++) {
            if (gsl_matrix_get(PF[s].Data, 2, t)<0.0000001) {Count+=1;};
        };
        if (Count<Month) {PFFinal.push_back(PF[s]);};
    };
    int PFSize=int(PFFinal.size());
    //int PFSize=int(PF.size());
    cout << "****************" << endl << "PF.size()=" << PFSize << endl << "Loading time: " << int(time(0)) - TickStart << "s" << endl << "****************" << endl << endl;
    
    return PFFinal;
}; // closes PortfolioGenerator()




gsl_matrix* PortfolioOutput (vector<Share> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    string A = Machine+"PF.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix* Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2-381);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Symbol << " (" << PF[J[i-1]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2-381; t++) {gsl_matrix_set(Result, 0, t, gsl_matrix_get(PF[0].Data, 0, t+381));}; // Dates
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2-381; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t+381)); // Close prices
        };
    };
    PlotGSLMatrix(Result, "PF.csv", 1);
    return Result;
};



vector<vector<Share>> PFTrainingTesting (vector<Share> PF, int TrainSize) {
    int Size=int(PF.size());
    vector<int> J=Shuffle(Size);
    vector<vector<Share>> PFVec;
    vector<Share> PFTraining, PFTesting;
    for (int k=0; k<TrainSize; k++) {PFTraining.push_back(PF[J[k]]);}; PFVec.push_back(PFTraining);
    for (int k=TrainSize; k<Size; k++) {PFTesting.push_back(PF[J[k]]);}; PFVec.push_back(PFTesting);
    return PFVec;
};



vector<vector<gsl_matrix*>> PortfolioMultiOutput (vector<Share> PF, int Size) {
    if (Size==0) {Size=int(PF.size());};
    string A = Machine+"PF.csv"; ofstream outputSpecs(A.c_str(), ofstream::app);
    vector<int> J=Shuffle(int(PF.size()));
    gsl_matrix* Result = gsl_matrix_calloc (Size+1, PF[0].Data->size2-381);
    gsl_matrix* Result2 = gsl_matrix_calloc (Size+1, PF[0].Data->size2-381);
    outputSpecs << fixed << "Date (yyyymmdd),"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Symbol << " (" << PF[J[i-1]].Exchange << "),";}; outputSpecs << endl;
    //outputSpecs << fixed << "Date,"; for (int i=1; i<Size+1; i++) {outputSpecs << PF[J[i-1]].Exchange <<",";}; outputSpecs << endl;
    outputSpecs.close();
    for (int t=0; t<PF[0].Data->size2-381; t++) {// Dates
        gsl_matrix_set(Result, 0, t, gsl_matrix_get(PF[0].Data, 0, t+381));
        gsl_matrix_set(Result2, 0, t, gsl_matrix_get(PF[0].Data, 0, t+381));
    };
    for (int i=1; i<Size+1; i++) {
        for (int t=0; t<PF[0].Data->size2-381; t++) {
            gsl_matrix_set(Result, i, t, gsl_matrix_get(PF[J[i-1]].Data, 2, t+381)); // Close prices
            gsl_matrix_set(Result2, i, t, gsl_matrix_get(PF[J[i-1]].Data, 3, t+381)); // Volumes
        };
        //cout << PF[J[i-1]].Exchange << ": " << PF[J[i-1]].Symbol << ", " << PF[J[i-1]].Title << " selected" << endl;
    };
    PlotGSLMatrix(Result, "PF.csv", 1);
    //outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << endl;
    vector<vector<gsl_matrix*>> FinalResult;
    int T=int(PF[0].Data->size2-381);
    for (int m=1; m<Size+1; m++) {
        vector<gsl_matrix*> Temp;
        gsl_matrix* Moments = gsl_matrix_calloc (64, T);
        gsl_matrix* Prices = gsl_matrix_calloc (1, T); for (int t=0; t<T; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (Result, m, t));};
        //PlotGSLMatrix(Prices, "Prices.csv", 1);
        for (int t=1; t<T; t++) {
            if (gsl_matrix_get(Result, m, t-1)!=0) {gsl_matrix_set(Moments, 0, t, log(gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)));} // Log-return
            else {gsl_matrix_set(Moments, 0, t, 0);}; // Log-return failsafe
            gsl_matrix_set(Moments, 7, t, abs(gsl_matrix_get(Moments, 0, t))); // Abs log-return
            //if (gsl_matrix_get(Result2, m, t-1)!=0) {gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)/gsl_matrix_get(Result2, m, t-1));} // Volumes-return
            //else {gsl_matrix_set(Moments, 26, t, 1);}; // Volumes-return failsafe
            gsl_matrix_set(Moments, 26, t, gsl_matrix_get(Result2, m, t)); // Volumes
        };
        double Crash=0; double TrendCounter=0;
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 1, t, MALAutoCorrelation (Moments, t, Week, 0)); // AC of log-returns at lag of Week
            gsl_matrix_set(Moments, 2, t, MALAutoCorrelation (Moments, t, 2*Week, 0)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 3, t, MALAutoCorrelation (Moments, t, Month, 0)); // AC of log-returns at lag of Month
            gsl_matrix_set(Moments, 4, t, MALAutoCorrelation (Moments, t, 3*Month, 0)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 5, t, MALAutoCorrelation (Moments, t, 6*Month, 0)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 6, t, MALAutoCorrelation (Moments, t, Year, 0)); // AC of log-returns at lag of Year
            gsl_matrix_set(Moments, 8, t, abs(MALAutoCorrelation (Moments, t, Week, 7))); // AC of abs log-returns at lag of Week
            gsl_matrix_set(Moments, 9, t, abs(MALAutoCorrelation (Moments, t, 2*Week, 7))); // AC of abs log-returns at lag of 2*Week
            gsl_matrix_set(Moments, 10, t, abs(MALAutoCorrelation (Moments, t, Month, 7))); // AC of abs log-returns at lag of Month
            gsl_matrix_set(Moments, 11, t, abs(MALAutoCorrelation (Moments, t, 3*Month, 7))); // AC of abs log-returns at lag of 3*Month
            gsl_matrix_set(Moments, 12, t, abs(MALAutoCorrelation (Moments, t, 6*Month, 7))); // AC of abs log-returns at lag of 6*Month
            gsl_matrix_set(Moments, 13, t, abs(MALAutoCorrelation (Moments, t, Year, 7))); // AC of abs log-returns at lag of Year
            gsl_matrix_set(Moments, 14, t, MALVolatility (Result, t, Week, m)); // Volatility at lag of Week
            gsl_matrix_set(Moments, 16, t, MALVolatility (Result, t, 2*Week, m)); // Volatility at lag of 2*Week
            gsl_matrix_set(Moments, 18, t, MALVolatility (Result, t, Month, m)); // Volatility at lag of Month
            gsl_matrix_set(Moments, 20, t, MALVolatility (Result, t, 3*Month, m)); // Volatility at lag of 3*Month
            gsl_matrix_set(Moments, 22, t, MALVolatility (Result, t, 6*Month, m)); // Volatility at lag of 6*Month
            gsl_matrix_set(Moments, 24, t, MALVolatility (Result, t, Year, m)); // Volatility at lag of Year
            gsl_matrix_set(Moments, 27, t, MALAutoCorrelation (Moments, t, Week, 26)); // AC of volumes at lag of Week
            gsl_matrix_set(Moments, 28, t, MALAutoCorrelation (Moments, t, 2*Week, 26)); // AC of volumes at lag of 2*Week
            gsl_matrix_set(Moments, 29, t, MALAutoCorrelation (Moments, t, Month, 26)); // AC of volumes at lag of Month
            gsl_matrix_set(Moments, 30, t, MALAutoCorrelation (Moments, t, 3*Month, 26)); // AC of volumes at lag of 3*Month
            gsl_matrix_set(Moments, 31, t, MALAutoCorrelation (Moments, t, 6*Month, 26)); // AC of volumes at lag of 6*Month
            gsl_matrix_set(Moments, 32, t, MALAutoCorrelation (Moments, t, Year, 26)); // AC of volumes at lag of Year
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33, t, MALAutoCorrelation (Moments, t, 3*Week, 0)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0)); // AC of log-returns at lag of Month shifted by Month
            
            //Systemics
            gsl_matrix_set(Moments, 54, t, gsl_matrix_get(Result, m, t)); // Market prices
            gsl_matrix_set(Moments, 55, t, -100+100*gsl_matrix_get(Result, m, t)/gsl_matrix_get(Result, m, t-1)); // Percentage of returns -100+100*P(t)/P(t-1)
            gsl_matrix_set(Moments, 56, t, gsl_matrix_get(Result2, m, t)); // Volumes
            gsl_matrix_set(Moments, 57, t, gsl_matrix_get(Moments, 14, t)); // w-volatility
            gsl_matrix_set(Moments, 58, t, gsl_matrix_get(Moments, 18, t)); // m-volatility
            gsl_matrix_set(Moments, 59, t, gsl_matrix_get(Moments, 22, t)); // 6m-volatility
            if (t>2*Week) {
                double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (Result, m, k)/Week;};
                double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (Result, m, k)/Week;};
                gsl_matrix_set (Moments, 60, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            };
            if (gsl_matrix_get (Moments, 60, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (Moments, 61, t, Crash); // Crash
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter=0;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)>0) && (gsl_matrix_get (Moments, 60, t)>0)) {TrendCounter+=1;}; // Trend reversion
            if ((gsl_matrix_get (Moments, 60, t-1)<0) && (gsl_matrix_get (Moments, 60, t)<0)) {TrendCounter-=1;}; // Trend reversion
            gsl_matrix_set (Moments, 62, t, TrendCounter);
            
        };
        for (int t=1; t<T; t++) {
            gsl_matrix_set(Moments, 15, t, MALAutoCorrelation (Moments, t, Week, 14)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set(Moments, 17, t, MALAutoCorrelation (Moments, t, 2*Week, 16)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set(Moments, 19, t, MALAutoCorrelation (Moments, t, Month, 18)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set(Moments, 21, t, MALAutoCorrelation (Moments, t, 3*Month, 20)); // AC of volatility at lag of 3*Month for lag of 3*Month
            gsl_matrix_set(Moments, 23, t, MALAutoCorrelation (Moments, t, 6*Month, 22)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set(Moments, 25, t, MALAutoCorrelation (Moments, t, Year, 24)); // AC of volatility at lag of Year for lag of Year
        };
        
        vector<double> V, GoF; for (int t=0; t<T; t++) {V.push_back(gsl_matrix_get (Prices, 0, t)); GoF.push_back(1);}; // FFF2
        vector<double>* pV= &V;
        const double Rescale=10; // Rescale factor
        const double Threshold=0.1; // Threshold (0.1)
        const size_t NumIteration=10000000; // Number of iterations of the minimizer (stops before once a minimum is found)
        double R; double* pR = &R; // Hawkes process parameters (does not need to be initialized)
        vector<double>& pGoF = GoF; // Goodness of fit
        vector<double>* pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
        for (int t=0; t<T; t++) {gsl_matrix_set (Moments, 63, t, (*pRes)[t]);}; // Hawkes process intensity in column 9
        for (int t=0; t<T; t++) { // Hawkes process intensity in column 9
            if ((*pRes)[t]>10) { // Maximum value of Hawkes at 10
                gsl_matrix_set (Moments, 63, t, 10);
            }
            else if ((*pRes)[t]<0) { // Minimum value of Hawkes at 0
                gsl_matrix_set (Moments, 63, t, 0);
            };
        };
        //gsl_matrix* Hawkes = HistoricalIntensity (Prices, 0);
        //for (int t=1; t<T; t++) {gsl_matrix_set(Moments, 42, t, gsl_matrix_get (Hawkes, 0, t));};
        
        //delete pV; delete pR; delete pRes;
        //pV=0; pR=0; pRes=0;
        
        PlotGSLMatrix(Moments, "MomentsReal.xls", 1);
        ofstream outputLog(Machine+"MomentsReal.xls", ofstream::app); outputLog << endl;
        Temp.push_back(Moments);
        FinalResult.push_back(Temp);
        //Temp.erase(Temp.begin(), Temp.end());
    };
    
    
    return FinalResult;
};

// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***
// END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE *** END OF ULM ALGORITHM CODE ***


// This plots the distribution of an STL vector taken as parameter in Mathematica. Setting Xmin=Xmax=0 finds these bounds automatically
vector<vector<double>> Distribution (vector<double> X, int Precision, double Xmin, double Xmax, string Name, string Plot) {
    vector<vector<double>> DistRes;
    // Randomization of output file name
    string A2 = Machine;
    string B2=Name;
    if (Name=="000") {int B = int(1000*((STLRandom(5, "Uniform"))[2])); B2=to_string(B);};
    A2+=B2;
    const char* Title = A2.c_str();
    ofstream output5(Title);
    
    int SeriesSize = int(X.size());
    double res;
    vector<double> X1, X2, Xbin, Xgrid; X1=X; X2=X;
    sort (X1.begin(), X1.begin() + SeriesSize); sort (X2.begin(), X2.begin() + SeriesSize);
    // After sorting the series, we take the first and last element as min and max resp.
    if ((Xmin==0) && (Xmax==0)) {Xmin = X1[0]; Xmax = X1[SeriesSize-1];};
    
    // Now generating the distribution itself in Excel
    for (int i=0; i<Precision; i+=1) {
        res=0;
        double BinStep = (Xmax-Xmin)/Precision;
        Xgrid.push_back(Xmin + i*BinStep/Precision);
        if (Plot=="On") {output5 << Xmin + i*BinStep << "\t";};
        for (int j=0; j<SeriesSize; j+=1) {
            if((X2[j] >= Xmin + i*BinStep) && (X2[j] < Xmin + (i+1)*BinStep)) {res=res+1;};
        }
        Xbin.push_back(res);
        if (Plot=="On") {output5 << Xbin[i] << endl;}
    }
    if (Plot=="On") {output5.close();}
    
    DistRes.push_back(Xgrid);
    DistRes.push_back(Xbin);
    return DistRes;
}




// This plots the distribution of an STL vector taken as parameter in Mathematica. Setting Xmin=Xmax=0 finds these bounds automatically
gsl_matrix* GSLDistribution (gsl_matrix* M, int Precision) {
    int SeriesSize1 = int(M->size1);
    int SeriesSize2 = int(M->size2);
    gsl_matrix* D = gsl_matrix_calloc (SeriesSize1*2, Precision);
    double Xmin, Xmax;
    // Determining the min and max of each time series j from gsl_matrix M
    vector<double> X;
    for (int j=0; j<SeriesSize1; j++) {
        for (int t=0; t<SeriesSize2; t++) {X.push_back(gsl_matrix_get(M,j,t));};
        sort (X.begin(), X.begin() + SeriesSize2);
        Xmin = X[0]; Xmax = X[SeriesSize2-1];
        double BinStep = (Xmax-Xmin+1)/Precision; // The +1 term is what fixed the whole range issue
        // Generating the distribution
        for (int k=0; k<Precision; k+=1) {
            double res=0;
            gsl_matrix_set(D, 2*j, k, Xmin+k*BinStep);
            for (int t=0; t<SeriesSize2; t+=1) {
                if ((X[t] >= Xmin + k*BinStep) && (X[t] < Xmin + (k+1)*BinStep)) {res+=1;};
            }; // closes t loop
            gsl_matrix_set(D, 2*j+1, k, res);
        }; // closes k loop
        //gsl_matrix_set(D, 2*j+1, Precision-1, gsl_matrix_get(D, 2*j+1, Precision-1)+1); // Adding one bin count at Xmax
        X.clear();
    }; // closes j loop
    return D;
}


double Compare (vector<double> X1, vector<double> X2) {
    int Size = int(min(X1.size(),X2.size()));
    double Result=0;
    for (int i=0; i<Size; i++) {
        Result += abs(X1[i]-X2[i]);
    }; // closes i loop
    return 100*Result/Size;
}; // closes Compare()






// This function returns the Future parameters for a NumberOfAgents and a simulation of time Time, according to a quasi-hyperbolic discount function F with parameters Beta and Delta: F(0)=1, F(t)=Beta*pow(Delta, t) with 0 < Beta,Delta < 1
vector<int> Quasihyperbolic(int NumberOfAgents, int Time, double Beta, double Delta, double Seed) {
    vector<int> Res;
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(Seed)); // gsl_rng_set(const gsl_rng* r, unsigned long int s)
    for (int i=0; i<NumberOfAgents+10; i++) {
        while (i>0) {
            double x = gsl_rng_uniform (r);
            double y = gsl_rng_uniform (r);
            if (y<=Beta*pow(Delta, floor(10*x))) {Res.push_back(x*Time); break;};
        };
    };
    return Res;
};









// Outputs the agent parameters for the stitch into one file: for each agent, the first line are integers, the second doubles, the third vector<int> Accuracy, the fourth vector<double> Leash, the fifth vector<double> LeashVol, the sixth vector<vector<int>> FPi[j][i], and the seventh vector<vector<int>> TPi[j][i]
void AgentParametersOutput (vector<Agent>& Market, int NumberOfAgents, int NumberOfStocks, string OutputName) {
    string RootName = Machine; RootName+=OutputName;
    ofstream outputParameters(RootName.c_str(), ofstream::app);
    for (int i=0; i<NumberOfAgents ; i++) {
        // int Future, History, TradingWindow, LiquidationFloor, Exploration, NEB, NEB1, NEBLossAversion, NEB3, NEB4, NEBPositivity, NEB6
        outputParameters << Market[i].Future << ", " << Market[i].History << ", " << Market[i].TradingWindow << ", " << Market[i].LiquidationFloor << ", " << Market[i].Exploration << ", " << Market[i].NEB << ", " << Market[i].NEBLossAversion << ", " << Market[i].NEBPositivity << endl;
        
        // double Reflexivity, Versatility, Gesture, Epsilon, Human, NEB1p, NEB3p, NEB4p, NEBPositivity, NEBLearningRate, NEB8, NEB9
        outputParameters << Market[i].Reflexivity << ", " << Market[i].Versatility << ", " << Market[i].Gesture << ", " << Market[i].Epsilon << ", " << Market[i].Human << ", " << Market[i].NEBPositivity << ", " << Market[i].NEBLearningRate << endl;
        
        // vector<int> Accuracy
        for (int j=0; j<NumberOfStocks ; j++) {
            if (j<NumberOfStocks-1) {outputParameters << Market[i].Accuracy[j] << ", ";}
            else {outputParameters << Market[i].Accuracy[j] << endl;};
        }; // closes j loop
        
        // vector<double> Leash
        for (int j=0; j<NumberOfStocks ; j++) {
            if (j<NumberOfStocks-1) {outputParameters << Market[i].Leash[j] << ", ";}
            else {outputParameters << Market[i].Leash[j] << endl;};
        }; // closes j loop
        
        // vector<double> LeashVol
        for (int j=0; j<NumberOfStocks ; j++) {
            if (j<NumberOfStocks-1) {outputParameters << Market[i].LeashVol[j] << ", ";}
            else {outputParameters << Market[i].LeashVol[j] << endl;};
        }; // closes j loop
        
        // vector<vector<int>> FPi[j][i]
        for (int j=0; j<NumberOfStocks ; j++) {
            for (int k=0; k < Market[0].FS * Market[0].FA; k++) {
                if (j<Market[0].FS * Market[0].FA * NumberOfStocks - 1) {outputParameters << Market[i].FPi[j][k] << ", ";}
                else {outputParameters << Market[i].FPi[j][k] << endl;};
            }; // closes k loop
        }; // closes j loop
        
        // vector<vector<int>> TPi[j][i]
        for (int j=0; j<NumberOfStocks ; j++) {
            for (int k=0; k < Market[0].TS * Market[0].TA; k++) {
                if (j<Market[0].TS * Market[0].TA * NumberOfStocks - 1) {outputParameters << Market[i].TPi[j][k] << ", ";}
                else {outputParameters << Market[i].TPi[j][k] << endl;};
            }; // closes k loop
        }; // closes j loop
        
    }; // closes i loop
    
    
    outputParameters.close();
}; // closes OutputParameters()










// Outputs the agent parameters of agents in a line
void AgentParameters (vector<Agent>& Market, int NumberOfAgents, int NumberOfAgentsFull, int NumberOfStocks, int Time) {
    ofstream outputParameters(Machine+"AgentParameters.txt", ofstream::app);
    vector<int> J = Shuffle(NumberOfAgentsFull);
    outputParameters << "I=" << NumberOfAgentsFull << ", J=" << NumberOfStocks << ", T=" << Time << endl;
    outputParameters << "***********************" << endl;
    outputParameters << "Selected agents: i={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << J[i] << ", ";} else {outputParameters << J[i] << "}" << endl;};};
    outputParameters << "Future={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Future << ", ";} else {outputParameters << Market[J[i]].Future << "}" << endl;};};
    outputParameters << "History={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].History << ", ";} else {outputParameters << Market[J[i]].History << "}" << endl;};};
    outputParameters << "TradingWindow={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].TradingWindow << ", ";} else {outputParameters << Market[J[i]].TradingWindow << "}" << endl;};};
    outputParameters << "Exploration={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Exploration << ", ";} else {outputParameters << Market[J[i]].Exploration << "}" << endl;};};
    outputParameters << "Epsilon={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Epsilon << ", ";} else {outputParameters << Market[J[i]].Epsilon << "}" << endl;};};
    outputParameters << "Gesture={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Gesture << ", ";} else {outputParameters << Market[J[i]].Gesture << "}" << endl;};};
    outputParameters << "Versatility={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Versatility << ", ";} else {outputParameters << Market[J[i]].Versatility << "}" << endl;};};
    outputParameters << "Reflexivity={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Reflexivity << ", ";} else {outputParameters << Market[J[i]].Reflexivity << "}" << endl;};};
    //outputParameters << "DiscountFactor={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].DiscountFactor << ", ";} else {outputParameters << Market[J[i]].DiscountFactor << "}" << endl;};};
    outputParameters << "RFA={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].RFAFirst << ", ";} else {outputParameters << Market[J[i]].RFAFirst << "}" << endl;};};
    outputParameters << "RBA={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].RBAFirst << ", ";} else {outputParameters << Market[J[i]].RBAFirst << "}" << endl;};};
    outputParameters << "Bankruptcy={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Bankruptcy << ", ";} else {outputParameters << Market[J[i]].Bankruptcy << "}" << endl;};};
    outputParameters << "LiquidationFloor={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].LiquidationFloor << ", ";} else {outputParameters << Market[J[i]].LiquidationFloor << "}" << endl;};};
    outputParameters << "Human={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Human << ", ";} else {outputParameters << Market[J[i]].Human << "}" << endl;};};
    outputParameters << "NEBLossAversion={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBLossAversion << ", ";} else {outputParameters << Market[J[i]].NEBLossAversion << "}" << endl;};};
    outputParameters << "NEBPositivity={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBPositivity << ", ";} else {outputParameters << Market[J[i]].NEBPositivity << "}" << endl;};};
    outputParameters << "NEBLearningRate={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBLearningRate << ", ";} else {outputParameters << Market[J[i]].NEBLearningRate << "}" << endl;};};
    outputParameters << "NEBPositivity={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].NEBPositivity << ", ";} else {outputParameters << Market[J[i]].NEBPositivity << "}" << endl;};};
    outputParameters << "***********************" << endl;
    outputParameters.close();
}; // closes OutputParameters()






// Outputs the RL variables (mostly vectors) of agents in a line, for a given stock j
void AgentVariables (vector<Agent>& Market, int NumberOfAgents, int NumberOfAgentsFull, int NumberOfStocks, int Time, gsl_matrix* ReflexiveValues, int j) {
    vector<int> J = Shuffle(NumberOfAgentsFull);
    string A = Machine+"AgentVariables_j";
    string B=to_string(j);
    string C= ".txt";
    A+=B+C;
    const char* Title=A.c_str();
    ofstream outputParameters(Title, ofstream::app);
    outputParameters << "***********************" << "j=" << j << "***********************" << endl;
    outputParameters << "I=" << NumberOfAgentsFull << ", J=" << NumberOfStocks << ", T=" << Time << endl;
    outputParameters << "Selected agents: i={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << J[i] << ", ";} else {outputParameters << J[i] << "}" << endl;};};
    outputParameters << "Q={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Stocks[j].StockQuantity << ", ";} else {outputParameters << Market[J[i]].Stocks[j].StockQuantity << "}" << endl;};};
    outputParameters << "Bid={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Stocks[j].StockBidValue << ", ";} else {outputParameters << Market[J[i]].Stocks[j].StockBidValue << "}" << endl;};};
    outputParameters << "Ask={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Stocks[j].StockAskValue << ", ";} else {outputParameters << Market[J[i]].Stocks[j].StockAskValue << "}" << endl;};};
    outputParameters << "Leash={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Leash[j] << ", ";} else {outputParameters << Market[J[i]].Leash[j] << "}" << endl;};};
    outputParameters << "LeashVol={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].LeashVol[j] << ", ";} else {outputParameters << Market[J[i]].LeashVol[j] << "}" << endl;};};
    outputParameters << "Accuracy={"; for (int i=0; i<NumberOfAgents ; i++) {if (i!=NumberOfAgents-1)  {outputParameters << Market[J[i]].Accuracy[j] << ", ";} else {outputParameters << Market[J[i]].Accuracy[j] << "}" << endl;};};
    outputParameters << "***********************" << endl;
    outputParameters << "ReflexiveValues={"; for (int t=0; t<Time ; t++) {if (t!=Time-1)  {outputParameters << gsl_matrix_get (ReflexiveValues, j, t) << ", ";} else {outputParameters << gsl_matrix_get (ReflexiveValues, j, t) << "}" << endl;};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FSIndex[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FSIndex[j].size()) ; u++) {if (u!=int(Market[J[i]].FSIndex[j].size())-1)  {outputParameters << Market[J[i]].FSIndex[j][u] << ", ";} else {outputParameters << Market[J[i]].FSIndex[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dRoughVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dRoughVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dRoughVec[j].size())-1)  {outputParameters << Market[J[i]].dRoughVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dRoughVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dSmoothVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dSmoothVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dSmoothVec[j].size())-1)  {outputParameters << Market[J[i]].dSmoothVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dSmoothVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Lvol[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LvolLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LvolLog[j].size())-1)  {outputParameters << Market[J[i]].LvolLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LvolLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Svol[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].SvolLog[j].size()) ; u++) {if (u!=int(Market[J[i]].SvolLog[j].size())-1)  {outputParameters << Market[J[i]].SvolLog[j][u] << ", ";} else {outputParameters << Market[J[i]].SvolLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Rough[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].RoughLog[j].size()) ; u++) {if (u!=int(Market[J[i]].RoughLog[j].size())-1)  {outputParameters << Market[J[i]].RoughLog[j][u] << ", ";} else {outputParameters << Market[J[i]].RoughLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Smooth[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].SmoothLog[j].size()) ; u++) {if (u!=int(Market[J[i]].SmoothLog[j].size())-1)  {outputParameters << Market[J[i]].SmoothLog[j][u] << ", ";} else {outputParameters << Market[J[i]].SmoothLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Reflexive[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ReflexiveLog[j].size()) ; u++) {if (u!=int(Market[J[i]].ReflexiveLog[j].size())-1)  {outputParameters << Market[J[i]].ReflexiveLog[j][u] << ", ";} else {outputParameters << Market[J[i]].ReflexiveLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FAIndexReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FAIndexReal[j].size()) ; u++) {if (u!=int(Market[J[i]].FAIndexReal[j].size())-1)  {outputParameters << Market[J[i]].FAIndexReal[j][u] << ", ";} else {outputParameters << Market[J[i]].FAIndexReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ToolReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ToolRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].ToolRealLog[j].size())-1)  {outputParameters << Market[J[i]].ToolRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].ToolRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "LagReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LagRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LagRealLog[j].size())-1)  {outputParameters << Market[J[i]].LagRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LagRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "WeightReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].WeightRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].WeightRealLog[j].size())-1)  {outputParameters << Market[J[i]].WeightRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].WeightRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FAIndexVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FAIndexVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].FAIndexVirtual[j].size())-1)  {outputParameters << Market[J[i]].FAIndexVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].FAIndexVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ToolVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ToolVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].ToolVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].ToolVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].ToolVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "LagVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LagVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LagVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].LagVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LagVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "WeightVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].WeightVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].WeightVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].WeightVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].WeightVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ForecastReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ForecastReal[j].size()) ; u++) {if (u!=int(Market[J[i]].ForecastReal[j].size())-1)  {outputParameters << Market[J[i]].ForecastReal[j][u] << ", ";} else {outputParameters << Market[J[i]].ForecastReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dFResultVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dFResultVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dFResultVec[j].size())-1)  {outputParameters << Market[J[i]].dFResultVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dFResultVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FResultDisReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FResultDisReal[j].size()) ; u++) {if (u!=int(Market[J[i]].FResultDisReal[j].size())-1)  {outputParameters << Market[J[i]].FResultDisReal[j][u] << ", ";} else {outputParameters << Market[J[i]].FResultDisReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "ForecastVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].ForecastVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].ForecastVirtual[j].size())-1)  {outputParameters << Market[J[i]].ForecastVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].ForecastVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "FResultDisVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].FResultDisVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].FResultDisVirtual[j].size())-1)  {outputParameters << Market[J[i]].FResultDisVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].FResultDisVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TSIndex[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TSIndex[j].size()) ; u++) {if (u!=int(Market[J[i]].TSIndex[j].size())-1)  {outputParameters << Market[J[i]].TSIndex[j][u] << ", ";} else {outputParameters << Market[J[i]].TSIndex[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dMuPosVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dMuPosVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dMuPosVec[j].size())-1)  {outputParameters << Market[J[i]].dMuPosVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dMuPosVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dMuNegVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dMuNegVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dMuNegVec[j].size())-1)  {outputParameters << Market[J[i]].dMuNegVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dMuNegVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dSigVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dSigVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dSigVec[j].size())-1)  {outputParameters << Market[J[i]].dSigVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dSigVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "LiquidPercentVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LiquidPercentVec[j].size()) ; u++) {if (u!=int(Market[J[i]].LiquidPercentVec[j].size())-1)  {outputParameters << Market[J[i]].LiquidPercentVec[j][u] << ", ";} else {outputParameters << Market[J[i]].LiquidPercentVec[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Mu[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].MuLog[j].size()) ; u++) {if (u!=int(Market[J[i]].MuLog[j].size())-1)  {outputParameters << Market[J[i]].MuLog[j][u] << ", ";} else {outputParameters << Market[J[i]].MuLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Sig[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].SigLog[j].size()) ; u++) {if (u!=int(Market[J[i]].SigLog[j].size())-1)  {outputParameters << Market[J[i]].SigLog[j][u] << ", ";} else {outputParameters << Market[J[i]].SigLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "RFALevel[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].RFALog[j].size()) ; u++) {if (u!=int(Market[J[i]].RFALog[j].size())-1)  {outputParameters << Market[J[i]].RFALog[j][u] << ", ";} else {outputParameters << Market[J[i]].RFALog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "RBALevel[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].RBALog[j].size()) ; u++) {if (u!=int(Market[J[i]].RBALog[j].size())-1)  {outputParameters << Market[J[i]].RBALog[j][u] << ", ";} else {outputParameters << Market[J[i]].RBALog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "Liquidity[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].LiquidityLog[j].size()) ; u++) {if (u!=int(Market[J[i]].LiquidityLog[j].size())-1)  {outputParameters << Market[J[i]].LiquidityLog[j][u] << ", ";} else {outputParameters << Market[J[i]].LiquidityLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TAIndexReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TAIndexReal[j].size()) ; u++) {if (u!=int(Market[J[i]].TAIndexReal[j].size())-1)  {outputParameters << Market[J[i]].TAIndexReal[j][u] << ", ";} else {outputParameters << Market[J[i]].TAIndexReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantRealLog[j].size())-1)  {outputParameters << Market[J[i]].QuantRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "PinchReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].PinchRealLog[j].size()) ; u++) {if (u!=int(Market[J[i]].PinchRealLog[j].size())-1)  {outputParameters << Market[J[i]].PinchRealLog[j][u] << ", ";} else {outputParameters << Market[J[i]].PinchRealLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TAIndexVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TAIndexVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].TAIndexVirtual[j].size())-1)  {outputParameters << Market[J[i]].TAIndexVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].TAIndexVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].QuantVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "PinchVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].PinchVirtualLog[j].size()) ; u++) {if (u!=int(Market[J[i]].PinchVirtualLog[j].size())-1)  {outputParameters << Market[J[i]].PinchVirtualLog[j][u] << ", ";} else {outputParameters << Market[J[i]].PinchVirtualLog[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantitiesReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantitiesReal[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantitiesReal[j].size())-1)  {outputParameters << Market[J[i]].QuantitiesReal[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantitiesReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "QuantitiesVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].QuantitiesVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].QuantitiesVirtual[j].size())-1)  {outputParameters << Market[J[i]].QuantitiesVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].QuantitiesVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TransactionPriceReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TransactionPriceReal[j].size()) ; u++) {if (u!=int(Market[J[i]].TransactionPriceReal[j].size())-1)  {outputParameters << Market[J[i]].TransactionPriceReal[j][u] << ", ";} else {outputParameters << Market[J[i]].TransactionPriceReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TransactionPriceVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TransactionPriceVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].TransactionPriceVirtual[j].size())-1)  {outputParameters << Market[J[i]].TransactionPriceVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].TransactionPriceVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    double dTRmean=0;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "dTResultVec[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].dTResultVec[j].size()) ; u++) {if (u!=int(Market[J[i]].dTResultVec[j].size())-1)  {outputParameters << Market[J[i]].dTResultVec[j][u] << ", ";} else {outputParameters << Market[J[i]].dTResultVec[j][u] << "}";}; dTRmean+=Market[J[i]].dTResultVec[j][u]/int(Market[J[i]].dTResultVec[j].size());}; outputParameters << " => <" << dTRmean << ">" << endl;}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TResultDisReal[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TResultDisReal[j].size()) ; u++) {if (u!=int(Market[J[i]].TResultDisReal[j].size())-1)  {outputParameters << Market[J[i]].TResultDisReal[j][u] << ", ";} else {outputParameters << Market[J[i]].TResultDisReal[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;
    for (int i=0; i<NumberOfAgents ; i++) {outputParameters << "TResultDisVirtual[j](i=" << J[i] << ")={"; for (int u=0; u<int(Market[J[i]].TResultDisVirtual[j].size()) ; u++) {if (u!=int(Market[J[i]].TResultDisVirtual[j].size())-1)  {outputParameters << Market[J[i]].TResultDisVirtual[j][u] << ", ";} else {outputParameters << Market[J[i]].TResultDisVirtual[j][u] << "}" << endl;};};}; outputParameters << "***********************" << endl;     outputParameters.close();
    
}; // closes OutputParameters()



gsl_matrix*  StockMoments (gsl_matrix* ReflexiveValues, gsl_matrix* Spread, gsl_matrix* TotalStockQuantityTraded, int Time, int Lag, int Delta, int CorrLag, int j) {
    gsl_matrix* StockIndexMoments = gsl_matrix_calloc (12, Time);
    gsl_matrix* StocksDistributions = gsl_matrix_calloc (24, 1000);
    gsl_matrix* AutoCorrLagP = gsl_matrix_alloc (24, Time);
    vector<double> RealizedVolatility, SpreadV, Volume, Stock;
    for (int t=0; t<Time; t++) {Stock.push_back(gsl_matrix_get(ReflexiveValues, j, t));};
    vector<vector<double>> MomentsStockIndex = DynamicMoments (Stock, 0, int(Stock.size()), Lag, Delta);
    for (int t=0; t<Time; t++) {
        gsl_matrix_set(StockIndexMoments, 0, t, Stock[t]);                              // Raw Stock
        gsl_matrix_set(StockIndexMoments, 1, t, MomentsStockIndex[0][t]); // Mean of Stock
        gsl_matrix_set(StockIndexMoments, 2, t, MomentsStockIndex[4][t]); // StdDev of Stock
        gsl_matrix_set(StockIndexMoments, 3, t, MomentsStockIndex[2][t]); // Skewness of Stock
        gsl_matrix_set(StockIndexMoments, 4, t, MomentsStockIndex[3][t]); // Kurtosis of Stock
        gsl_matrix_set(StockIndexMoments, 5, t, MomentsStockIndex[5][t]); // Stock + 3 Sigma
        gsl_matrix_set(StockIndexMoments, 6, t, MomentsStockIndex[6][t]); // Stock - 3 Sigma
        gsl_matrix_set(StockIndexMoments, 7, t, MomentsStockIndex[7][t]); // Realized volatility of Stock
        gsl_matrix_set(StockIndexMoments, 8, t, MomentsStockIndex[8][t]); // Autocorrelations of lag Delta in bsp for Stock
        gsl_matrix_set(StockIndexMoments, 9, t, MomentsStockIndex[9][t]); // Squarred of daily log returns
        gsl_matrix_set(StockIndexMoments, 10, t, MomentsStockIndex[10][t]); // Squarred of weekly log returns
        gsl_matrix_set(StockIndexMoments, 11, t, MomentsStockIndex[11][t]); // Squarred of bi-weekly log returns
        RealizedVolatility.push_back(gsl_matrix_get(StockIndexMoments, 7, t));
        SpreadV.push_back(gsl_matrix_get(Spread, j, t));
        Volume.push_back(gsl_matrix_get(TotalStockQuantityTraded, j, t));
    };
    StocksDistributions = GSLDistribution(StockIndexMoments, 1000);
    AutoCorrLagP = AutocorrelationLagP (StockIndexMoments);
    string B=to_string(j); string C= ".xls";
    string A = "Moments_j"; A+=B+C;
    PlotGSLMatrix(StockIndexMoments, A.c_str(), 1);
    A = "Distributions_j"; A+=B+C;
    PlotGSLMatrix(StocksDistributions, A.c_str(), 1);
    A = "AutoCorrLagP_j"; A+=B+C;
    PlotGSLMatrix(AutoCorrLagP, A.c_str(), 1);
    
    // Outputing correlations
    vector<double> CorrRealizedVolatilityVolume = DynamicCorrelation (RealizedVolatility, Volume, CorrLag);
    vector<double> CorrRealizedVolatilitySpread = DynamicCorrelation (RealizedVolatility, SpreadV, CorrLag);
    vector<double> CorrSpreadVolume = DynamicCorrelation (SpreadV, Volume, CorrLag);
    int Minsize = min(int(CorrRealizedVolatilityVolume.size()), min(int(CorrRealizedVolatilitySpread.size()), int(CorrSpreadVolume.size())));
    gsl_matrix* Correlations = gsl_matrix_alloc (3, Minsize);
    for (int t=0; t<Minsize; t++) {
        gsl_matrix_set (Correlations, 0, t, CorrRealizedVolatilityVolume[t]);
        gsl_matrix_set (Correlations, 1, t, CorrRealizedVolatilitySpread[t]);
        gsl_matrix_set (Correlations, 2, t, CorrSpreadVolume[t]);
    };
    A = "Correlations_j"; A+=B+C;
    PlotGSLMatrix(Correlations, A.c_str(), 1);
    
    return StockIndexMoments;
}; // closes function StockMoments()










// Policy distances (NON-MULTIVARIATE!)
// The simulation can output matrices of dimension IxI of the average of the absolute differences between each value P(s,a) of two agents. This is what we call Policy Distance (PD): FPolicyDistances, TPolicyDistances, MostSuccessfulFPolicyDistances, MostSuccessfulTPolicyDistances. But the outputs are 3d plots that are hard to read. So we should create a measure of policy distance not between an agent i1 and an agent i2, but between an agent i1 and many other agents (either all others or just a subgroup, such as the most successful percentile) as the average all its policy distances with all other agents, and then we should rank these agents by their increasing Average Policy Distances (APD). Then we know that an agent with a large APD compared to the market (or subset thereof) Average of other APD's (AAPD) has a largely different policy than the rest of the market or a subset thereof. Then we can extract features about the parameters of these most different or most similar agents to market or submarket.
gsl_matrix* PolicyDistances (vector<Agent> Market, int NumberOfAgents, int NumberOfStocks, int Time, int Percentile, int t) {
    // Most successful and less successful agents
    vector<double> CapitalEvolution;
    for (int i=0; i<NumberOfAgents ; i++) {
        double x = Market[i].Capital();
        double y = Market[i].RFAFirst + Market[i].RBAFirst;
        CapitalEvolution.push_back(100*x/y);
    }; // closes i loop
    vector<int> DescendingRank;
    for (int k=0; k<NumberOfAgents; k++) {
        double Bar=0; int BarX=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if ((CapitalEvolution[i]>=Bar) && (CapitalEvolution[i]>=0)) {Bar=CapitalEvolution[i]; BarX=i;};
        }; // closes i loop
        CapitalEvolution[BarX]=-1;
        DescendingRank.push_back(BarX);
    }; // closes k loop
    int AgentNumberPercentile=Percentile*NumberOfAgents/100;
    if (AgentNumberPercentile<=5) {AgentNumberPercentile=5;};
    
    gsl_matrix* Res = gsl_matrix_calloc (10, AgentNumberPercentile);
    gsl_matrix* FPolicyDistancesMostSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* TPolicyDistancesMostSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* FPolicyDistancesLessSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* TPolicyDistancesLessSuccessfulToMarket = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* FPolicyDistancesMostSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* TPolicyDistancesMostSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* FPolicyDistancesLessSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* TPolicyDistancesLessSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* FPolicyDistancesMostSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* TPolicyDistancesMostSuccessfulToLessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* FPolicyDistancesLessSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* TPolicyDistancesLessSuccessfulToMostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, NumberOfStocks);
    gsl_matrix* FPolicyDistancesM = gsl_matrix_calloc (NumberOfAgents, NumberOfAgents);
    gsl_matrix* TPolicyDistancesM = gsl_matrix_calloc (NumberOfAgents, NumberOfAgents);
   
    for (int j=0; j<NumberOfStocks; j++) {
        // Computing FPolicyDistancesM and TPolicyDistancesM as the IxI matrices of all PD of each agent with each agent
        int SizeFPi = int(Market[0].FPi[j].size());
        int SizeTPi = int(Market[0].TPi[j].size());
        for (int i1=0; i1<NumberOfAgents; i1++) {
            for (int i2=i1+1; i2<NumberOfAgents; i2++) {
                cout << "****************************" << endl;
                cout << "** STOCK MARKET SIMULATOR **" << endl;
                cout << "****************************" << endl;
                cout << "I=" << NumberOfAgents << ", J=" << NumberOfStocks << ", T=" << Time << endl  << endl;
                cout.setf(ios::fixed);
                cout << "PD step             : " << i1*NumberOfAgents + i2 << "/" << NumberOfAgents*NumberOfAgents << " (" << 100.0*(i1*NumberOfAgents + i2)/(NumberOfAgents*NumberOfAgents) << "%)" << endl;
                
                double ResultF=0; double ResultT=0;
                for (int k=0; k<SizeFPi; k++) {ResultF += abs(Market[i1].FPi[j][k] - Market[i2].FPi[j][k]);};
                gsl_matrix_set (FPolicyDistancesM, i1, i2, 100*ResultF/SizeFPi);
                gsl_matrix_set (FPolicyDistancesM, i2, i1, 100*ResultF/SizeFPi);
                for (int k=0; k<SizeTPi; k++) {ResultT += abs(Market[i1].TPi[j][k] - Market[i2].TPi[j][k]);};
                gsl_matrix_set (TPolicyDistancesM, i1, i2, 100*ResultT/SizeTPi);
                gsl_matrix_set (TPolicyDistancesM, i2, i1, 100*ResultT/SizeTPi);
            }; // closes i2 loop
        }; // closes i1 loop
        
        // Computing FPolicyDistancesMostSuccessfulToMarket and TPolicyDistancesMostSuccessfulToMarket as matrices of j APD of most successful agents with all agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<NumberOfAgents; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1]].AgentName, i2);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1]].AgentName, i2);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesMostSuccessfulToMarket, i1, j, ResultF/(NumberOfAgents-1));
            gsl_matrix_set (TPolicyDistancesMostSuccessfulToMarket, i1, j, ResultT/(NumberOfAgents-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of less successful agents with all agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<NumberOfAgents; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, i2);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, i2);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesLessSuccessfulToMarket, i1, j, ResultF/(NumberOfAgents-1));
            gsl_matrix_set (TPolicyDistancesLessSuccessfulToMarket, i1, j, ResultT/(NumberOfAgents-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of most successful agents with most successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesMostSuccessfulToMostSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesMostSuccessfulToMostSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of less successful agents with less successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesLessSuccessfulToLessSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesLessSuccessfulToLessSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of most successful agents with less successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1]].AgentName, Market[DescendingRank[i2+NumberOfAgents-AgentNumberPercentile]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesMostSuccessfulToLessSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesMostSuccessfulToLessSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
        // Computing FPolicyDistancesLessSuccessfulToMarket and TPolicyDistancesLessSuccessfulToMarket as matrices of j APD of less successful agents with most successful agents
        for (int i1=0; i1<AgentNumberPercentile; i1++) {
            double ResultF=0; double ResultT=0;
            for (int i2=0; i2<AgentNumberPercentile; i2++) {
                ResultF += gsl_matrix_get (FPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2]].AgentName);
                ResultT += gsl_matrix_get (TPolicyDistancesM, Market[DescendingRank[i1+NumberOfAgents-AgentNumberPercentile]].AgentName, Market[DescendingRank[i2]].AgentName);
            }; // closes i2 loop
            gsl_matrix_set (FPolicyDistancesLessSuccessfulToMostSuccessful, i1, j, ResultF/(AgentNumberPercentile-1));
            gsl_matrix_set (TPolicyDistancesLessSuccessfulToMostSuccessful, i1, j, ResultT/(AgentNumberPercentile-1));
        }; // closes i1 loop
        
    }; // closes j loop
    
    for (int k=0; k<AgentNumberPercentile; k++) {
        gsl_matrix_set (Res, 0, k, gsl_matrix_get (FPolicyDistancesMostSuccessfulToMostSuccessful, k, 0));
        gsl_matrix_set (Res, 1, k, gsl_matrix_get (FPolicyDistancesMostSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 2, k, gsl_matrix_get (FPolicyDistancesMostSuccessfulToLessSuccessful, k, 0));
        gsl_matrix_set (Res, 3, k, gsl_matrix_get (FPolicyDistancesLessSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 4, k, gsl_matrix_get (FPolicyDistancesLessSuccessfulToLessSuccessful, k, 0));
        gsl_matrix_set (Res, 5, k, gsl_matrix_get (TPolicyDistancesMostSuccessfulToMostSuccessful, k, 0));
        gsl_matrix_set (Res, 6, k, gsl_matrix_get (TPolicyDistancesMostSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 7, k, gsl_matrix_get (TPolicyDistancesMostSuccessfulToLessSuccessful, k, 0));
        gsl_matrix_set (Res, 8, k, gsl_matrix_get (TPolicyDistancesLessSuccessfulToMarket, k, 0));
        gsl_matrix_set (Res, 9, k, gsl_matrix_get (TPolicyDistancesLessSuccessfulToLessSuccessful, k, 0));
    }; // closes k loop
    PlotGSLMatrix(Res, "PDistances.xls", 1);
    
    // Matrix memory freeing
    gsl_matrix_free(FPolicyDistancesM);
    gsl_matrix_free(FPolicyDistancesMostSuccessfulToMostSuccessful);
    gsl_matrix_free(FPolicyDistancesMostSuccessfulToMarket);
    gsl_matrix_free(FPolicyDistancesMostSuccessfulToLessSuccessful);
    gsl_matrix_free(FPolicyDistancesLessSuccessfulToMarket);
    gsl_matrix_free(FPolicyDistancesLessSuccessfulToLessSuccessful);
    gsl_matrix_free(TPolicyDistancesM);
    gsl_matrix_free(TPolicyDistancesMostSuccessfulToMostSuccessful);
    gsl_matrix_free(TPolicyDistancesMostSuccessfulToMarket);
    gsl_matrix_free(TPolicyDistancesMostSuccessfulToLessSuccessful);
    gsl_matrix_free(TPolicyDistancesLessSuccessfulToMarket);
    gsl_matrix_free(TPolicyDistancesLessSuccessfulToLessSuccessful);
    
    return Res;

};




void PoliciesAverages (vector<Agent> Market, int NumberOfAgents, int NumberOfStocks, int Percentile) {
    vector<gsl_matrix*> Res;
    // Most successful and less successful agents
    vector<double> CapitalEvolution;
    for (int i=0; i<NumberOfAgents ; i++) {
        double x = Market[i].Capital();
        double y = Market[i].RFAFirst + Market[i].RBAFirst;
        CapitalEvolution.push_back(100*x/y);
    }; // closes i loop
    vector<int> DescendingRank;
    for (int k=0; k<NumberOfAgents; k++) {
        double Bar=0; int BarX=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if ((CapitalEvolution[i]>=Bar) && (CapitalEvolution[i]>=0)) {Bar=CapitalEvolution[i]; BarX=i;};
        }; // closes i loop
        CapitalEvolution[BarX]=-1;
        DescendingRank.push_back(BarX);
    }; // closes k loop
    int AgentNumberPercentile=Percentile*NumberOfAgents/100;
    if (AgentNumberPercentile<=5) {AgentNumberPercentile=5;};
    
    gsl_matrix* FPoliciesAverageBest = gsl_matrix_calloc (Market[0].FA, Market[0].FS); // Heat maps
    gsl_matrix* FPoliciesAverageWorst = gsl_matrix_calloc (Market[0].FA, Market[0].FS); // Heat maps
    gsl_matrix* TPoliciesAverageBest = gsl_matrix_calloc (Market[0].TA, Market[0].TS); // Heat maps
    gsl_matrix* TPoliciesAverageWorst = gsl_matrix_calloc (Market[0].TA, Market[0].TS); // Heat maps

    for (int j=0; j<NumberOfStocks; j++) {
        // Computing the F() policies averages as heat maps for Best and Worst populations
        for (int n=0; n<AgentNumberPercentile; n++) {
            for (int i=0; i<(Market[0].FS)*(Market[0].FA); i++) {
                gsl_matrix_set (FPoliciesAverageBest, i%Market[0].FA, i/Market[0].FA, Market[DescendingRank[n]].FPi[j][i]/AgentNumberPercentile + gsl_matrix_get (FPoliciesAverageBest, i%Market[0].FA, i/Market[0].FA));
            }; // closes i loop
            for (int i=0; i<(Market[0].FS)*(Market[0].FA); i++) {
                gsl_matrix_set (FPoliciesAverageWorst, i%Market[0].FA, i/Market[0].FA, Market[DescendingRank[n+NumberOfAgents-AgentNumberPercentile]].FPi[j][i]/AgentNumberPercentile + gsl_matrix_get (FPoliciesAverageWorst, i%Market[0].FA, i/Market[0].FA));
            }; // closes i loop
        }; // closes n loop

        // Computing the T() policies averages as heat maps for Best and Worst populations
        for (int n=0; n<AgentNumberPercentile; n++) {
                for (int i=0; i<(Market[0].TS)*(Market[0].TA); i++) {
                        gsl_matrix_set (TPoliciesAverageBest, i%Market[0].TA, i/Market[0].TA, Market[DescendingRank[n]].TPi[j][i]/AgentNumberPercentile + gsl_matrix_get (TPoliciesAverageBest, i%Market[0].TA, i/Market[0].TA));
                }; // closes i loop
                for (int i=0; i<(Market[0].TS)*(Market[0].TA); i++) {
                        gsl_matrix_set (TPoliciesAverageWorst, i%Market[0].TA, i/Market[0].TA, Market[DescendingRank[n+NumberOfAgents-AgentNumberPercentile]].TPi[j][i]/AgentNumberPercentile + gsl_matrix_get (TPoliciesAverageWorst, i%Market[0].TA, i/Market[0].TA));
                }; // closes i loop
        }; // closes n loop

    }; // closes j loop
    
    PlotGSLMatrix(FPoliciesAverageBest, "FPoliciesAverageBest.xls", 1);
    PlotGSLMatrix(FPoliciesAverageWorst, "FPoliciesAverageWorst", 1);
    PlotGSLMatrix(TPoliciesAverageBest, "TPoliciesAverageBest.xls", 1);
    PlotGSLMatrix(TPoliciesAverageWorst, "TPoliciesAverageWorst", 1);
    
    // Matrix memory freeing
    gsl_matrix_free(FPoliciesAverageBest);
    gsl_matrix_free(FPoliciesAverageWorst);
    gsl_matrix_free(TPoliciesAverageBest);
    gsl_matrix_free(TPoliciesAverageWorst);
    
};



// Function to Simulate a given Market over a given number of time steps
//pNEBLossAversion24 => pNEBLossAversion, pNEBDelayDiscounting56=>pNEBPositivity, pNEB89=>pNEBDelayDiscounting
vector<gsl_matrix*> MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int pNEB, double LearningRateScale, int s, int Trunk, int MetaorderImpact) {
    
    int pBots=100; int pNEBLossAversion=0; int pNEBPositivity=0; int pNEBNegativity=0; int pNEBDelayDiscounting=0; int pNEBFear=0; int pNEBGreed=0; int pNEBLearningRate=0;
    if (TypeNEB=="Classic") {pBots=80; pNEBLossAversion=3; pNEBPositivity=3; pNEBNegativity=3; pNEBDelayDiscounting=3; pNEBFear=3; pNEBGreed=3; pNEBLearningRate=2;}
    else if (TypeNEB=="Algorithmic") {pBots=100; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Human") {pBots=100-pNEB; pNEBLossAversion=pNEB-(6*pNEB/7); pNEBPositivity=pNEB/7; pNEBNegativity=pNEB/7; pNEBDelayDiscounting=pNEB/7; pNEBFear=pNEB/7; pNEBGreed=pNEB/7; pNEBLearningRate=pNEB/7;}
    else if (TypeNEB=="LossAversion") {pBots=100-pNEB; pNEBLossAversion=pNEB; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Positivity") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=pNEB; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Negativity") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=pNEB; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="DelayDiscounting") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=pNEB; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Fear") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=pNEB; pNEBGreed=0; pNEBLearningRate=0;}
    else if (TypeNEB=="Greed") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=pNEB; pNEBLearningRate=0;}
    else if (TypeNEB=="LearningRate") {pBots=100-pNEB; pNEBLossAversion=0; pNEBPositivity=0; pNEBNegativity=0; pNEBDelayDiscounting=0; pNEBFear=0; pNEBGreed=0; pNEBLearningRate=pNEB;};
    
    
    
    string VersatilityCondition="Off"; string TradingFrequencyCond="On"; int Percentile=10; string OutputName="OutputCondOff";
    // ofstream outputDebug("/Users/admin/Documents/GNT/SYMBA/Debug.txt", ofstream::app);
    ofstream outputLog(Machine+"SimLog.txt", ofstream::app);
    ofstream outputMetalog(Machine+"MetaorderLog.txt", ofstream::app);
    ofstream outputClustered(Machine+"ClusteredBestWorstAgent.csv", ofstream::app);
    if (Plot=="On") {
        outputLog << "**********SUM UP OF SIMULATION**********" << endl;
        outputLog << "SS1: Initialization at t=0" << endl;
        outputLog << "SS2: Start of t-loop t>=1" << endl;
        outputLog << "SS3: Computing ReflexiveValues" << endl;
        //outputLog << "SS4: Computing ReferenceValues" << endl;
        outputLog << "SS5: Computing interest rates" << endl;
        outputLog << "SS6: RL()" << endl;
        outputLog << "SS7: MarketEvolution.push_back(Market)" << endl;
        outputLog << "SS8: Filling of OB" << endl;
        outputLog << "SS9: Sort()" << endl;
        outputLog << "SS10: Output of OB" << endl;
        outputLog << "SS11: Clear()" << endl;
        outputLog << "SS12: Market indicators" << endl;
        outputLog << "SS13: End of t-loop" << endl;
        outputLog << "*********************************************" << endl << endl;
    };
    
    /*
     int t=0;
     string N1="/Users/admin/Documents/GNT/SYMBA/SimLog";
     string N2=to_string(t);
     string N3=".txt";
     N1+=N2+N3;
     const char* Title = N1.c_str();
     ofstream outputLog(Title, ofstream::app);
     */
    // INITIALIZING CLOCK (SS1)
    int TickStart=int(time(0)); // Starting the clock to measure time of computation
    int TickEnd=TickStart;
    int CountBrankrupcies=0;
    double DividendYield=2.276*0.01; // According to http://indexarb.com/dividendYieldSorteddj.html
    double ClusteredBestAgent=0; // Number of  best agents that are within the ClusterLimit population
    double ClusteredWorstAgent=0; // Number of  best agents that are within the ClusterLimit population
    int SharesOutstanding=0; // For stock j=0 only (study of metaorder impact)
    int LastMetaorder=0; // Last time step a metaorder was injected in the OB
    vector<int> MetaorderInjectionTimes; // Times at which metaorders are injected
    //vector<int> MetaorderInjectionAmplitudes; // Effective volume transacted at metaorder injection
    
    //system("CLS");
    cout << "INITIALIZED CLOCK..." << endl;
    
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0))); // gsl_rng_set(const gsl_rng* r, unsigned long int s)
    vector<double> StockVlm, StockVlmInit, StockIndex, StockLOD, TotalStockQuantities;
    //string TimeScope="NoTimeScope";
    
    // INITIALIZING EACH AGENT
    gsl_matrix* ReflexiveValues = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all bids
    gsl_matrix* AverageBids = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all bids
    gsl_matrix* AverageAsks = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all asks
    gsl_matrix* AverageTopBids = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all bids in buisness
    gsl_matrix* AverageTopAsks = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the average of all asks in buisness
    gsl_matrix* MedianTopBids = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the median of all bids in buisness
    gsl_matrix* MedianTopAsks = gsl_matrix_calloc (NumberOfStocks, Time); // Computed as the median of all asks in buisness
    gsl_matrix* HighestBid = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix* LowestAsk = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix* GSLTrueValues = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix* Spread = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix* TotalStockQuantityTraded = gsl_matrix_calloc (NumberOfStocks, Time);
    gsl_matrix* VolumesDistributions = gsl_matrix_alloc (2*NumberOfStocks, Time);
    gsl_matrix* ReturnsDistributions = gsl_matrix_alloc (2*NumberOfStocks, Time);
    gsl_matrix* AutoCorrLagPVolumes = gsl_matrix_alloc (2*NumberOfStocks, Time);
    gsl_matrix* CapitalEvolution = gsl_matrix_calloc (NumberOfAgents, Time);
    gsl_matrix* Bankruptcies = gsl_matrix_calloc (2, Time); // Bankruptcies
    gsl_matrix_set(Bankruptcies, 1, 0, 100);
    gsl_matrix* BiasedDelta = gsl_matrix_alloc (NumberOfStocks, Time);
    gsl_matrix* TrueDelta = gsl_matrix_alloc (NumberOfStocks, Time);
    vector<gsl_matrix*> PDs;
    vector<gsl_matrix*> PAverages;
    
    // Generating the Future parameter of each agent according to quasi-hyperbolic discounting
    //vector<double> Gen = STLRandom (NumberOfStocks, "Uniform", "NoPlot");
    int LearningPhase=1000;
    /*
     int Tau=Time;
     if (Time>2*Year) {Tau=8*Month;}
     else if ((Time<=2*Year) && (Time>1*Year)) {Tau=4*Month;}
     else if ((Time<=1*Year) && (Time>6*Month)) {Tau=3*Month;}
     else if ((Time<=6*Month) && (Time>1*Month)) {Tau=2*Month;}
     else if (Time<=1*Month) {Tau=Week;}
     //Tau/=2; // AAA
     vector<int> QH = Quasihyperbolic(NumberOfAgents, Tau, 1, 0.8, 1000*(gsl_rng_uniform (r)+0.0001)); // Beta=1, Delta=0.8 for the Hyperbolic discounting model
     vector<int> Reflex = Quasihyperbolic(NumberOfAgents, 100, 1, 0.75, 1000*(gsl_rng_uniform (r)+0.0001)); // e.g Beta=1, Delta=0.75 yields ~75% of agents with reflexivity below 0.5 // GGG
     //gsl_matrix* QHPlot = gsl_matrix_calloc (1, int(QH.size())); for (int i=0; i<int(QH.size()); i++) {gsl_matrix_set(QHPlot, 0, i, QH[i]);}; PlotGSLMatrix(QHPlot, "QHPlot.xls", 1); exit(0);
     //vector<int> Quasihyperbolic(int NumberOfAgents, int Time, double Beta, double Delta, double Seed) {
     */
    
    vector<Agent> Market;
    int Counter=1;
    for (int j=0; j<NumberOfStocks; j++) {
        StockVlmInit.push_back(0);
    }; // closes j loop
    for (int i=0; i<NumberOfAgents; i++) {
        //outputLog << "Agent " << i << ", ";
        Agent TempAgent;
        vector<Stock> TempStocks;
        TempAgent.AgentName=i;
        TempAgent.RFA=1000*abs((gsl_ran_gaussian (r, 10)));
        TempAgent.RFAFirst=TempAgent.RFA; // RFA at time t=0
        //TempAgent.DiscountFactor = gsl_rng_uniform (r);
        //TempAgent.Future=QH[i]+Week;
        TempAgent.Future=ceil(6*Month*gsl_rng_uniform (r))+Week;
        //if (Time<=1*Month) {TempAgent.Future=int(ceil(Week*gsl_rng_uniform (r)));}; // For very small Time, we go for a basic uniform distribution
        //if ((TempAgent.Future<2) || (TempAgent.Future>Tau)) { // Failsafe
        //outputLog << "Agent " << i << ": Future=" << TempAgent.Future << " (QH[" << i << "]=" << QH[i] << ") replaced to Week" << endl;
        //TempAgent.Future=Week;
        //};
        TempAgent.History=abs(Week+int(gsl_rng_uniform (r)*(Time-TempAgent.Future-1-(2*Week)))); //Week to Time-Future-1
        TempAgent.Reflexivity = gsl_rng_uniform (r);
        
        //TempAgent.Reflexivity*=0.75; // GGG
        //TempAgent.Reflexivity=0.01*Reflex[i]; // GGG
        //TempAgent.TradingWindow=int(0.25*TempAgent.Future + (gsl_rng_uniform(r) * ((1.5*TempAgent.Future)-(0.25*TempAgent.Future))));
        //TempAgent.TradingWindow=int(TempAgent.Future*(0.33 + 0.66*gsl_rng_uniform(r))); // EEE
        //if (TempAgent.TradingWindow<=2) {TempAgent.TradingWindow=3;}; // Lower bound
        //if (TempAgent.TradingWindow>TempAgent.Future) {TempAgent.TradingWindow=TempAgent.Future;}; // Higher bound
        TempAgent.TradingWindow=5+gsl_rng_uniform(r)*TempAgent.Future;
        TempAgent.Bankruptcy = 1; // Bankruptcy = 0 is effective bankruptcy
        //TempAgent.LiquidationFloor = 80 + 10*gsl_rng_uniform (r);
        TempAgent.LiquidationFloor = LiquidationFloor + 10*gsl_rng_uniform (r);
        TempAgent.BiasedValues = gsl_matrix_calloc (NumberOfStocks, Time);
        TempAgent.Versatility=1+floor(10*abs(gsl_rng_uniform (r)));
        //TempAgent.Gesture=floor(100*(4*gsl_rng_uniform (r)/5 + 0.2))/100; // Commercial gesture of the agent in percent of the spread
        TempAgent.Gesture=0.2+0.8*gsl_rng_uniform (r);
        //TempAgent.Gesture*=2;
        if (HPGesture==0) {TempAgent.Gesture*=1;}
        else if (HPGesture==1) {TempAgent.Gesture*=1.5;}
        else if (HPGesture==2) {TempAgent.Gesture*=2;}
        else if (HPGesture==3) {TempAgent.Gesture*=3;};
        //TempAgent.Exploration=int((1000*(gsl_rng_uniform (r))))%4 + 1;
        TempAgent.Exploration=1;
        TempAgent.OffPolicy=1;
        //TempAgent.Epsilon=0.1+0.15*gsl_rng_uniform (r); // 0.3r: The epsilon factor must be small
        TempAgent.Epsilon=0;
        TempAgent.Quant=0;
        TempAgent.Pinch=0;
        TempAgent.FS=27; // Number of states for Forecast()
        TempAgent.FA=27; // Number of actions for Forecast()
        TempAgent.TS=108; // Number of states for Trade()
        TempAgent.TA=9; // Number of actions for Trade()
        TempAgent.FSDim[0]=3; TempAgent.FSDim[1]=3; TempAgent.FSDim[2]=3;
        TempAgent.FADim[0]=3; TempAgent.FADim[1]=3; TempAgent.FADim[2]=3;
        TempAgent.FQDim[0]=3; TempAgent.FQDim[1]=3; TempAgent.FQDim[2]=3; TempAgent.FQDim[3]=3; TempAgent.FQDim[4]=3; TempAgent.FQDim[5]=3;
        TempAgent.FPiDim[0]=3; TempAgent.FPiDim[1]=3; TempAgent.FPiDim[2]=3; TempAgent.FPiDim[3]=3; TempAgent.FPiDim[4]=3; TempAgent.FPiDim[5]=3;
        //TempAgent.FPiDim[0]=4; TempAgent.FPiDim[1]=4; TempAgent.FPiDim[2]=3; TempAgent.FPiDim[3]=3; // DDD2
        TempAgent.TSDim[0]=3; TempAgent.TSDim[1]=3; TempAgent.TSDim[2]=2; TempAgent.TSDim[3]=2; TempAgent.TSDim[4]=3;
        TempAgent.TADim[0]=3; TempAgent.TADim[1]=3;
        TempAgent.TQDim[0]=3; TempAgent.TQDim[1]=3; TempAgent.TQDim[2]=2; TempAgent.TQDim[3]=2; TempAgent.TQDim[4]=3; TempAgent.TQDim[5]=3; TempAgent.TQDim[6]=3;
        TempAgent.TPiDim[0]=3; TempAgent.TPiDim[1]=3; TempAgent.TPiDim[2]=2; TempAgent.TPiDim[3]=2; TempAgent.TPiDim[4]=3; TempAgent.TPiDim[5]=3; TempAgent.TPiDim[6]=3;
        // RL() and LOG of vectors of vectors
        for (int j=0; j<NumberOfStocks; j++) {
            TempAgent.TradingWindowClock.push_back(0);
            vector<double> V; TempAgent.dRoughVec.push_back(V); TempAgent.dSmoothVec.push_back(V); TempAgent.dMuPosVec.push_back(V); TempAgent.dMuNegVec.push_back(V); TempAgent.dSigVec.push_back(V); TempAgent.LiquidPercentVec.push_back(V); TempAgent.dFResultVec.push_back(V); TempAgent.dTResultVec.push_back(V); TempAgent.FPi.push_back(V); TempAgent.TPi.push_back(V); TempAgent.FQ.push_back(V); TempAgent.TQ.push_back(V); TempAgent.ForecastReal.push_back(V); TempAgent.ForecastVirtual.push_back(V); TempAgent.Forecast5.push_back(V); TempAgent.TransactionPriceReal.push_back(V); TempAgent.TransactionPriceVirtual.push_back(V); TempAgent.Qs.push_back(V); // RL
            for (int u=0; u<TempAgent.TS; u++) {TempAgent.Qs[j].push_back(0);}; // Initializatioon of model-based RL framework
            TempAgent.LvolLog.push_back(V); TempAgent.SvolLog.push_back(V); // LOG
            vector<int> U; TempAgent.FSIndex.push_back(U); TempAgent.TSIndex.push_back(U); TempAgent.FAIndexReal.push_back(U); TempAgent.FAIndexVirtual.push_back(U); TempAgent.TAIndexReal.push_back(U); TempAgent.TAIndexVirtual.push_back(U); TempAgent.QuantitiesReal.push_back(U); TempAgent.QuantitiesVirtual.push_back(U); TempAgent.FNumberA.push_back(U); TempAgent.TNumberA.push_back(U); // RL
            TempAgent.RoughLog.push_back(U); TempAgent.SmoothLog.push_back(U); TempAgent.ReflexiveLog.push_back(U); TempAgent.ToolRealLog.push_back(U); TempAgent.LagRealLog.push_back(U); TempAgent.WeightRealLog.push_back(U); TempAgent.ToolVirtualLog.push_back(U); TempAgent.LagVirtualLog.push_back(U); TempAgent.WeightVirtualLog.push_back(U); TempAgent.MuLog.push_back(U); TempAgent.SigLog.push_back(U); TempAgent.RFALog.push_back(U); TempAgent.RBALog.push_back(U); TempAgent.LiquidityLog.push_back(U); TempAgent.PinchRealLog.push_back(U); TempAgent.QuantRealLog.push_back(U); TempAgent.PinchVirtualLog.push_back(U); TempAgent.QuantVirtualLog.push_back(U); TempAgent.FResultDisReal.push_back(U); TempAgent.FResultDisVirtual.push_back(U); TempAgent.TResultDisReal.push_back(U); TempAgent.TResultDisVirtual.push_back(U); // LOG
        }; // closes j loop
        
        
        // Attribution of biases (for all: if NEBp>r then NEB is activated)
        TempAgent.Human=100*gsl_rng_uniform (r); // If it is smaller than pBots then it is a bot
        TempAgent.NEBLossAversion=0.5*gsl_rng_uniform (r);
        TempAgent.NEBPositivity=0.5*gsl_rng_uniform (r);
        TempAgent.NEBNegativity=0.5*gsl_rng_uniform (r);
        // DelayDiscounting is set below
        TempAgent.NEBLearningRate=0.05+0.2*gsl_rng_uniform (r);
        TempAgent.NEBFear=0.2*gsl_rng_uniform (r);
        TempAgent.NEBGreed=0.2*gsl_rng_uniform (r);

        
        // 50% have nothing, 10% PROFILE1 (NEB124), 10% PROFILE2 (NEB356), 5% PROFILE3 (NEB89), 5% PROFILE12 (NEB123456), 5% PROFILE13 (NEB12489), 5% PROFILE23 (NEB35689), 10% everything
        //pNEBLossAversion24 => pNEBLossAversion, pNEBDelayDiscounting56=>pNEBPositivity, pNEB89=>pNEBDelayDiscounting
        if (TempAgent.Human<pBots) {
            TempAgent.NEB=0; // Bot
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots) && (TempAgent.Human<pBots + pNEBLossAversion)) {
            TempAgent.NEB=1; // Loss aversion
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity)) {
            TempAgent.NEB=2; // Positivity
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity)) {
            TempAgent.NEB=3; // Negativity
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting)) {
            TempAgent.NEB=4; // Delay discounting
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBNegativity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
            TempAgent.Future=ceil(Week*gsl_rng_uniform (r))+Week;
            TempAgent.History=abs(Week+int(gsl_rng_uniform (r)*(Time-TempAgent.Future-1-(2*Week)))); //2*Week to Time-Future-1
            TempAgent.TradingWindow=5+gsl_rng_uniform(r)*TempAgent.Future;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear)) {
            TempAgent.NEB=5; // Fear
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBGreed=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear + pNEBGreed)) {
            TempAgent.NEB=6; // Greed
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
        }
        else if ((TempAgent.Human>=pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear + pNEBGreed) && (TempAgent.Human<pBots + pNEBLossAversion + pNEBPositivity + pNEBNegativity + pNEBDelayDiscounting + pNEBFear + pNEBGreed + pNEBLearningRate)) {
            TempAgent.NEB=7; // Learning
            TempAgent.NEBLossAversion=-1;
            TempAgent.NEBPositivity=-1;
            TempAgent.NEBFear=-1;
            TempAgent.NEBGreed=-1;
            TempAgent.NEBLearningRate*=LearningRateScale;
        };
       
        for (int j=0; j<NumberOfStocks; j++) {
            TempAgent.Leash.push_back(0.001*gsl_rng_uniform (r));
            TempAgent.LeashVol.push_back(5*0.01*gsl_rng_uniform (r)); // TTT4
            //TempAgent.Accuracy.push_back(7+int((1000*(gsl_rng_uniform (r))))%5); // Prob. of bubble burst, on average 0.2,0.75,1,2 per year for Accuracy=13,11,10,9 resp.
            //TempAgent.Accuracy.push_back(10+gsl_rng_uniform (r)); // TTT3
            TempAgent.Accuracy.push_back(HPAccuracy+gsl_rng_uniform (r)); // TTT3
            Stock TempStock;
            TempStock.StockName = j;
            TempStock.StockQuantity = abs(int(gsl_ran_gaussian (r, 100)));
            TempStock.StockQuantityFirst = TempStock.StockQuantity;
            TempStocks.push_back(TempStock);
            Counter += 1;
        }; // closes j loop
        TempAgent.NAVOnJanFirst=TempAgent.Capital(); // NAV on Jan 1st
        TempAgent.Stocks = TempStocks;
        //TempStocks.erase(TempStocks.begin(),TempStocks.end());
        Market.push_back(TempAgent);
        SharesOutstanding+=TempAgent.Stocks[0].StockQuantity;
    }; // closes i loop
    //outputLog << "INITIALIZED AGENTS..." << endl;
    cout << "INITIALIZED AGENTS..." << endl;
    
    
    // GENERATING THE STOCKS TRUE VALUES
    vector<double> Gen2 = STLRandom (NumberOfStocks, "Uniform");
    //double HPTrueMu=0.1;
    for (int j=0; j<NumberOfStocks; j++) {
        //vector<double> S = PoissonRandomWalk (100, Week+3*Month*Gen2[j], int(0.1*Time/250+0.9*Gen2[j]*Time/250), Time, 1000*Gen2[j]);
        vector<double> S = PoissonRandomWalk (100, Week+HPTrueMu*3*Month*Gen2[j], int(0.1*Time/250+0.9*Gen2[j]*Time/250), Time, 1000*Gen2[j]);
        for (int t=0; t<Time; t++) {gsl_matrix_set (GSLTrueValues, j, t, S[t]);};
    };
    vector<vector<double>> TrueValues = GSLMatrixToSTLMatrix(GSLTrueValues);
    cout << "INITIALIZED STOCKS TRUE VALUES..." << endl;
    
    
    // INITIALIZING THE AGENTS BIASED VALUES
    for (int i=0; i<NumberOfAgents; i++) {
        for (int j=0; j<NumberOfStocks; j++) {
            vector<double> C = CointegratedWalk(TrueValues[j], Market[i].Leash[j], Market[i].LeashVol[j], Market[i].Accuracy[j]); // Master, Leash, LeashVolatility, Accuracy
            // vector<double> C = CointegratedWalk(TrueValues[j], 0.0005*2*gsl_rng_uniform (r), 0.005*2*gsl_rng_uniform (r), Accuracy); // Master, Leash, LeashVolatility, Accuracy
            //if ((j==0) || (j== NumberOfStocks-1)) {PlotSTL(C, "000");};
            for (int t=0; t<Time; t++) {
                gsl_matrix_set(Market[i].BiasedValues,j,t, C[t]);
            }; // closes t loop
            Market[i].Stocks[j].StockBidValue = gsl_matrix_get(Market[i].BiasedValues,j,0) + abs(gsl_ran_gaussian (r, 1));
            //Market[i].Stocks[j].StockAskValue = (Market[i].Stocks[j].StockBidValue)*(1 + Market[i].ProfitMargin + Market[i].RiskPremium);
            Market[i].Stocks[j].StockAskValue = (Market[i].Stocks[j].StockBidValue)*(1 + 0.02*gsl_rng_uniform (r));
            C.clear();
        }; // closes j loop
        Market[i].RBAFirst=Market[i].StockHoldings(); // RBA at time t=0
    }; // clloses i loop
    cout << "INITIALIZED AGENTS BIASED VALUES..." << endl;
    
    
    
    // Checking on distribution of agents Future parameter
    //gsl_matrix* M=gsl_matrix_alloc (1, NumberOfAgents);
    //for (int i=0; i<NumberOfAgents; i++) {gsl_matrix_set (M, 0, i, Market[i].Future);};
    //gsl_matrix* D = gsl_matrix_alloc (2, NumberOfAgents);
    //D = GSLDistribution(M, 10);
    //PlotGSLMatrix(M, "MM.xls", 1);
    //PlotGSLMatrix(D, "DD.xls", 1);
    
    
    
    // INITILAZING THE MARKET SIMULATION
    // Initializing reflexive values as averages at t=0, and reference values = f(biasedtrue,reflexive)
    vector<double> VSpread;
    for (int j=0; j<NumberOfStocks; j++) {
        double InitValBid=0;
        for (int i=0; i<NumberOfAgents; i++) {
            InitValBid+=Market[i].Stocks[j].StockBidValue;
        }; // closes i loop
        gsl_matrix_set(ReflexiveValues, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageBids, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageAsks, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageTopBids, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(AverageTopAsks, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(MedianTopBids, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(MedianTopAsks, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(HighestBid, j, 0, InitValBid/NumberOfAgents);
        gsl_matrix_set(LowestAsk, j, 0, InitValBid/NumberOfAgents);
        VSpread.push_back(0);
    }; // closes j loop
    
    // Initializing TotalStockQuantities, the vector of total quantity of each stock
    for (int j=0; j<NumberOfStocks; j++) {
        double TempStockQuantity=0;
        for (int i=0; i<NumberOfAgents; i++) {
            TempStockQuantity+=Market[i].Stocks[j].StockQuantity;
        }; // closes i loop
        TotalStockQuantities.push_back(TempStockQuantity);
    }; // closes j loop
    // Initializing Market Indicators
    for (int i=0; i<NumberOfAgents; i++) {gsl_matrix_set(CapitalEvolution, i, 0, 100);};
    cout << "INITIALIZED MARKET INDICATORS..." << endl;
    
    // INITIALIZATION OF MARKETEVOLUTION
    vector<vector<Agent>> MarketEvolution;
    MarketEvolution.push_back(Market);
    cout << "INITIALIZED MARKET EVOLUTION..." << endl;
    
    
    // INITIALIZATION OF THE ORDER BOOK
    vector<vector<StockOrderBook>> FullOrderBook; // This is our final order book with many pages t
    vector<StockOrderBook> TempOrderBook; // This is the order book at time t, composed of many order books (one for each stock)
    for (int j=0; j<NumberOfStocks; j++) {
        StockOrderBook TempStockOrderBook;
        vector<int> J = Shuffle(NumberOfAgents);
        for (int i=0; i<NumberOfAgents ; i++) {
            Order TempOrderBid;
            TempOrderBid.Nature=0; // 0 for Bid, 1 for Ask
            TempOrderBid.Agent=J[i];
            TempOrderBid.Stock=j;
            TempOrderBid.Q=MarketEvolution[0][J[i]].Stocks[j].StockQuantity; // CHECK!!! For now the buyer buys all the seller has !!!
            TempOrderBid.Price=MarketEvolution[0][J[i]].Stocks[j].StockBidValue;
            Order TempOrderAsk;
            TempOrderAsk.Nature=1; // 0 for Bid, 1 for Ask
            TempOrderAsk.Agent=J[i];
            TempOrderAsk.Stock=j;
            TempOrderAsk.Q=MarketEvolution[0][J[i]].Stocks[j].StockQuantity; // CHECK!!! For now the seller sells all he has !!!
            TempOrderAsk.Price=MarketEvolution[0][J[i]].Stocks[j].StockAskValue;
            // Now filling the final order book of stock j
            TempStockOrderBook.Bids.push_back(TempOrderBid);
            TempStockOrderBook.Asks.push_back(TempOrderAsk);
        }; // closes i loop
        TempOrderBook.push_back(TempStockOrderBook);
        TempStockOrderBook.Bids.clear();
        TempStockOrderBook.Asks.clear();
    }; // closes j loop
    FullOrderBook.push_back(TempOrderBook); // FullOrderBook is filled on first page at t=0
    cout << "INITIALIZED ORDER BOOK..." << endl;
    
    
    // INITIALIZATION t=0 FINISHES HERE
    int TickMid=(int(time(0))-TickStart); // Checking the clock for time of computation of initialization
    //ofstream outputT2("/Users/admin/Documents/GNT/SYMBA/Details.txt", ofstream::app);
    system("CLS");
    cout << "****************************" << endl;
    cout << "** STOCK MARKET SIMULATOR **" << endl;
    cout << "****************************" << endl;
    cout << "I=" << NumberOfAgents << ", J=" << NumberOfStocks << ", T=" << Time << endl << endl;
    cout.setf(ios::fixed);
    cout << "Initialization : " << TickMid << "s" << endl;
    if (Plot=="On") {outputLog << "INITIALIZATION t=0 FINISHED..." << endl << endl;};
    
    
    
    
    // ***** MARKET SIMULATION STARTS HERE WITH t-LOOP ***** (SS2)
    if (Plot=="On") {outputLog << "***** MARKET SIMULATION WITH t-LOOP STARTS HERE *****" << endl;};
    for (int t=1; t<Time; t++) {
        for (int j=0; j<NumberOfStocks; j++) {
            // REFLEXIVE VALUES = CLEARING PRICE OF ORDER BOOK (SS3)
            gsl_matrix_set(ReflexiveValues, j, t, gsl_matrix_get(ReflexiveValues, j, t-1));
            if (Plot=="On") {outputLog << "t=" << t << ", j=" << j << ": ReflexiveValues computed..." << endl;};
        }; // closes j loop
        // outputDebug << "t=" << t << ": ReflexiveValues computed..." << endl;
        
        
        // INTEREST RATES (SS5)
        for (int i=0; i<NumberOfAgents; i++) {
            Market[i].RFA*=1+(Rate/Year); // Saving account banking rate, modeled w/o inflation
            Market[i].RFA+=Market[i].StockHoldings()*DividendYield/Year; // Modeled as compounded dividend yield
        }; // closes i loop
        if (Plot=="On") {outputLog << "t=" << t << ", Interest rates computed..." << endl << endl;};
        // outputDebug << "t=" << t << ": Interest rates computed..." << endl;
        
        // REINITIALIZATION OF THE AGENTS' PORTFOLIO AFTER LEARNING PHASE
        for (int i=0; i<NumberOfAgents ; i++) {
            if (t==LearningPhase) {
                Market[i].Bankruptcy=1;
                Market[i].RFA=Market[i].RFAFirst;
                for (int j=0; j<NumberOfStocks; j++) {
                    Market[i].Stocks[j].StockQuantity=Market[i].Stocks[j].StockQuantityFirst;
                };
            };
        };
        
        // RECORD OF LEADER ACTIONS FOR CLUSTER TRADING
        int LeaderAgent=0; // "Random" agent is the leader for cluster trading
        double MinCap=9999999; double MaxCap=-9999999; int WorstAgent=0; int BestAgent=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if (Market[i].Capital()<MinCap) {MinCap=Market[i].Capital(); WorstAgent=i;};
            if (Market[i].Capital()>MaxCap) {MaxCap=Market[i].Capital(); BestAgent=i;};
        };
        if (LeaderType=="Worst") {LeaderAgent=WorstAgent;};
        if (LeaderType=="Best") {LeaderAgent=BestAgent;};
        if (LeaderType=="Static") {LeaderAgent=0;};
        if ((BestAgent<ClusterLimit) && (t>=LearningPhase)) {ClusteredBestAgent+=1;}; // Nb of LeaderAgent that are part of the ClusterLimit
        if ((WorstAgent<ClusterLimit) && (t>=LearningPhase)) {ClusteredWorstAgent+=1;}; // Nb of LeaderAgent that are part of the ClusterLimit
        
        if ((Market[LeaderAgent].Bankruptcy==0) && (ClusterLimit>1)) { // Agent leader cannot be bankrupt
            Market[LeaderAgent].RFA=Market[LeaderAgent].RFAFirst;
            Market[LeaderAgent].Bankruptcy=1;
            for (int j=0; j<NumberOfStocks; j++) {Market[LeaderAgent].Stocks[j].StockQuantity=Market[LeaderAgent].Stocks[j].StockQuantityFirst;};
        };
        int LeaderQuant=Market[LeaderAgent].Quant;
        int LeaderPinch=Market[LeaderAgent].Pinch;
        if (t<=2) {LeaderQuant=0; LeaderPinch=0;};
        
        
        
        
        // REINFORCEMENT LEARNING RL() (SS6)
        gsl_matrix_set(Bankruptcies, 1, t, gsl_matrix_get(Bankruptcies, 1, t-1));
        for (int j=0; j<NumberOfStocks; j++) {
            for (int i=0; i<NumberOfAgents; i++) {
                string NewPlot=Plot;
                if (LeaderType=="Noise") {
                    LeaderQuant=int((1000*(gsl_rng_uniform (r))))%3; //0, 1, 2
                    LeaderPinch=int((1000*(gsl_rng_uniform (r))))%3; //0, 1, 2
                };
                if (Market[i].Bankruptcy!=0) {Market[i].RL(j, t, Rate, ReflexiveValues, VSpread[j], gsl_matrix_get(TotalStockQuantityTraded, j, t-1), Time, NumberOfStocks, TradingFrequencyCond, NewPlot, VersatilityCondition, gsl_matrix_get(Bankruptcies, 1, t), t-(t-LearningPhase)%Year, LearningPhase, LeaderAgent, LeaderQuant, LeaderPinch, ClusterLimit, Trunk);}; // XYZ1
                
            }; // closes i loop
        }; // closes j loop
        if (Plot=="On") {outputLog << "t=" << t << ", RL() method ran in all agents..." << endl << endl;};
        // outputDebug << "t=" << t << ": RL() computed..." << endl;
        
        
        // FILLING THE ORDER BOOK (SS8)
        TempOrderBook.clear();
        for (int j=0; j<NumberOfStocks; j++) {
            StockOrderBook TempStockOrderBook;
            vector<int> J = Shuffle(NumberOfAgents);
            for (int i=0; i<NumberOfAgents ; i++) {
                bool B1 = Market[(J[i])].Bankruptcy != 1; // Agent is bankrupt
                bool B2 = Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size())-1]==0; // Agent is not trading at this time step
                if (Plot=="On") {
                    outputLog << "Filling TempOrderBook of Agent " << J[i] << ", Stock " << j << ", at t=" << t << ": ";
                    outputLog << "(Bankruptcy=" << B1 << ", Waiting=" << B2 << "): ";
                };
                if (B1 || B2) {if (Plot=="On") {outputLog << "order cancelled" << endl;}; continue;}
                else if ((Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1]>0) && !B1 && !B2) {
                    Order TempOrderBid;
                    TempOrderBid.Nature=0; // 0 for Bid, 1 for Ask
                    TempOrderBid.Agent=Market[(J[i])].AgentName;
                    TempOrderBid.Stock=j;
                    TempOrderBid.Price=Market[(J[i])].Stocks[j].StockBidValue;
                    //if (Trunk>-1) {TempOrderBid.Price=DigitTrunk (TempOrderBid.Price, Trunk, "Ceil");}; // Trunk according to the number of significant digits
                    TempOrderBid.Q=Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1];
                    TempStockOrderBook.Bids.push_back(TempOrderBid);
                    if (Plot=="On") {outputLog << "Bid Q=" << TempOrderBid.Q << " at $" << TempOrderBid.Price << " (RFA=$" << Market[(J[i])].RFA << ")" << endl;};
                } // closes else if
                else if ((Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1]<0) && !B1 && !B2) {
                    Order TempOrderAsk;
                    TempOrderAsk.Nature=1; // 0 for Bid, 1 for Ask
                    TempOrderAsk.Agent=Market[(J[i])].AgentName;
                    TempOrderAsk.Stock=j;
                    TempOrderAsk.Price=Market[(J[i])].Stocks[j].StockAskValue;
                    //if (Trunk>-1) {TempOrderAsk.Price=DigitTrunk (TempOrderAsk.Price, Trunk, "Ceil");}; // Trunk according to the number of significant digits
                    TempOrderAsk.Q=abs(Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1]);
                    TempStockOrderBook.Asks.push_back(TempOrderAsk);
                    if (Plot=="On") {outputLog << "Ask Q=" << TempOrderAsk.Q << " at $" << TempOrderAsk.Price << " (QuantitiesReal[j]=" << Market[(J[i])].QuantitiesReal[j][int(Market[(J[i])].QuantitiesReal[j].size()) - 1] << ")" << endl;};
                }; // closes else if
            }; // closes i loop
            TempStockOrderBook.Sort(); // (SS9)
            if (Plot=="On") {outputLog << "t=" << t << ", Sorted OB[j=" << j << "]" << endl << endl;};
            
            TempOrderBook.push_back(TempStockOrderBook);
            TempStockOrderBook.Bids.clear();
            TempStockOrderBook.Asks.clear();
        }; // closes j loop
        FullOrderBook.push_back(TempOrderBook); // FullOrderBook is filled at page t
        if (Plot=="On") {outputLog << "t=" << t << ", Filled the OB..." << endl;};
        // outputDebug << "t=" << t << ": OB filled..." << endl;
        
        
        // RECORDING AVERAGE BIDS AND ASKS
        // The market price is taken as the lowest of clearing prices (not the median), and the bid and ask are taken as the average (not median) of all bid and ask prices of orders that will be transacted: hence the bid and ask being sometimes not below and above of market price
        for (int j=0; j<NumberOfStocks; j++) {
            double MeanBid=0; double MeanAsk=0; double MeanTopBid=0; double MeanTopAsk=0; int Level=-1;
            int OBSize=int(min(FullOrderBook[t][j].Asks.size(), FullOrderBook[t][j].Bids.size()));
            if (OBSize==0) {
                gsl_matrix_set(AverageTopBids, j, t, gsl_matrix_get (AverageTopBids, j, t-1));
                gsl_matrix_set(AverageTopAsks, j, t, gsl_matrix_get (AverageTopAsks, j, t-1));
                gsl_matrix_set(MedianTopBids, j, t, gsl_matrix_get (MedianTopBids, j, t-1));
                gsl_matrix_set(MedianTopAsks, j, t, gsl_matrix_get (MedianTopAsks, j, t-1));
                gsl_matrix_set(AverageBids, j, t, gsl_matrix_get (AverageBids, j, t-1));
                gsl_matrix_set(AverageAsks, j, t, gsl_matrix_get (AverageAsks, j, t-1));
                gsl_matrix_set(HighestBid, j, t, gsl_matrix_get (HighestBid, j, t-1));
                gsl_matrix_set(LowestAsk, j, t, gsl_matrix_get (LowestAsk, j, t-1));
            }
            else {
                for (int k=0; k<OBSize; k++) {if (FullOrderBook[t][j].Bids[k].Price>=FullOrderBook[t][j].Asks[k].Price) {Level+=1;}}; // Level
                gsl_matrix_set(AverageTopBids, j, t, gsl_matrix_get (AverageTopBids, j, t-1));
                gsl_matrix_set(AverageTopAsks, j, t, gsl_matrix_get (AverageTopAsks, j, t-1));
                gsl_matrix_set(HighestBid, j, t, gsl_matrix_get (HighestBid, j, t-1));
                gsl_matrix_set(LowestAsk, j, t, gsl_matrix_get (LowestAsk, j, t-1));
                if (Level>-1) {
                    for (int k=0; k<Level+1; k++) {
                        MeanTopBid+=FullOrderBook[t][j].Bids[k].Price/(Level+1);
                        MeanTopAsk+=FullOrderBook[t][j].Asks[k].Price/(Level+1);
                    }; // closes k loop
                    gsl_matrix_set(AverageTopBids, j, t, MeanTopBid);
                    gsl_matrix_set(AverageTopAsks, j, t, MeanTopAsk);
                    gsl_matrix_set(MedianTopBids, j, t, FullOrderBook[t][j].Bids[(Level+1)/2].Price);
                    gsl_matrix_set(MedianTopAsks, j, t, FullOrderBook[t][j].Asks[(Level+1)/2].Price);
                    gsl_matrix_set(HighestBid, j, t, FullOrderBook[t][j].Bids[0].Price);
                    gsl_matrix_set(LowestAsk, j, t, FullOrderBook[t][j].Asks[0].Price);
                }; // closes if
                for (int k=0; k<OBSize; k++) {
                    MeanBid+=FullOrderBook[t][j].Bids[k].Price/OBSize;
                    MeanAsk+=FullOrderBook[t][j].Asks[k].Price/OBSize;
                }; // closes k loop
                gsl_matrix_set(AverageBids, j, t, MeanBid);
                gsl_matrix_set(AverageAsks, j, t, MeanAsk);
            }; // closes else
            VSpread.at(j)=abs(gsl_matrix_get(AverageBids, j, t) - gsl_matrix_get(AverageAsks, j, t));
            /*
             if (HPSpread==0) {VSpread.at(j)=abs(gsl_matrix_get(AverageBids, j, t) - gsl_matrix_get(AverageAsks, j, t));} // whole average definition of spread XXX3 // Small volatility
             else if (HPSpread==1) {VSpread.at(j)=abs(gsl_matrix_get(AverageTopBids, j, t) - gsl_matrix_get(AverageTopAsks, j, t));} // average definition of spread XXX1 // Too smooth
             else if (HPSpread==2) {VSpread.at(j)=abs(gsl_matrix_get(MedianTopBids, j, t) - gsl_matrix_get(MedianTopAsks, j, t));} // median definition of spread XXX2 // Too smooth
             else {VSpread.at(j)=gsl_matrix_get(HighestBid, j, t) - gsl_matrix_get(LowestAsk, j, t); if (VSpread[j]<0) {VSpread.at(j)=0;}}; // Formal definition (difference between the highest bid and lowest ask) XXX4 // Too volatile
             */
            double FormalDefinitionSpread=gsl_matrix_get(HighestBid, j, t) - gsl_matrix_get(LowestAsk, j, t); if (FormalDefinitionSpread<0) {FormalDefinitionSpread=0;};
            if (Plot=="On") {
                outputLog << "gsl_matrix_get(AverageTopBids, j, t)=" << gsl_matrix_get(AverageTopBids, j, t) << ", gsl_matrix_get(AverageTopAsks, j, t)=" << gsl_matrix_get(AverageTopAsks, j, t) << ", Type-I VSpread[j]=" << abs(gsl_matrix_get(AverageTopBids, j, t) - gsl_matrix_get(AverageTopAsks, j, t)) << endl;
                outputLog << "gsl_matrix_get(MedianTopBids, j, t)=" << gsl_matrix_get(MedianTopBids, j, t) << ", gsl_matrix_get(MedianTopAsks, j, t)=" << gsl_matrix_get(MedianTopAsks, j, t) << ", Type-II VSpread[j]=" << abs(gsl_matrix_get(MedianTopBids, j, t) - gsl_matrix_get(MedianTopAsks, j, t)) << endl;
                outputLog << "gsl_matrix_get(AverageBids, j, t)=" << gsl_matrix_get(AverageBids, j, t) << ", gsl_matrix_get(AverageAsks, j, t)=" << gsl_matrix_get(AverageAsks, j, t) << ", Type-III VSpread[j]=" << abs(gsl_matrix_get(AverageBids, j, t) - gsl_matrix_get(AverageAsks, j, t)) << endl;
                outputLog << "gsl_matrix_get(HighestBid, j, t)=" << gsl_matrix_get(HighestBid, j, t) << ", gsl_matrix_get(LowestAsk, j, t)=" << gsl_matrix_get(LowestAsk, j, t) << ", Type-IV VSpread[j]=" << FormalDefinitionSpread << endl;
                outputLog << "VSpread[j]=" << VSpread[j] << endl;
            };
        }; // closes j loop
        // outputDebug << "t=" << t << ": Average bids and asks computed..." << endl;
        
        
        // LOG OF CREDIT EVENTS
        // Bankruptcy = 0 is effective bankruptcy
        int LearningPhase2=10;
        double NAVOfNonBankruptAgents=0; int NbOfNonBankruptAgents=0; double MarketPerformanceJan1st=0;
        if ((t>=LearningPhase2) && ((t-LearningPhase2)%Year2==0)) {
            for (int i=0; i<NumberOfAgents ; i++) {
                Market[i].NAVOnJanFirst=Market[i].Capital(); // Recording all Jan 1st NAVs
            };
            MarketPerformanceJan1st=gsl_matrix_get(ReflexiveValues, 0, t);
        };
        if ((t>LearningPhase2) && ((t-LearningPhase2)%Year2>Week)) {
            for (int i=0; i<NumberOfAgents ; i++) {
                if (Market[i].Bankruptcy!=0) {NAVOfNonBankruptAgents+=100*Market[i].Capital()/Market[i].NAVOnJanFirst; NbOfNonBankruptAgents+=1;};
            };
            NAVOfNonBankruptAgents/=NbOfNonBankruptAgents; // Mean NAV ratio (%) since Jan 1st of all non-bankrupt agents
            if (NbOfNonBankruptAgents==0) {NAVOfNonBankruptAgents=100;}; // Failsafe
            gsl_matrix_set(Bankruptcies, 1, t, NAVOfNonBankruptAgents);
            outputLog << "At t=" << t << ", NAVOfNonBankruptAgents=" << NAVOfNonBankruptAgents << ", NbOfNonBankruptAgents=" << NbOfNonBankruptAgents << ", Pt(Jan1st)=" << MarketPerformanceJan1st << ", Pt=" << gsl_matrix_get(ReflexiveValues, 0, t) << endl;
        };
        
        int Discount=max(0, 100-int(NAVOfNonBankruptAgents));
        for (int i=0; i<NumberOfAgents ; i++) {
            if ((Market[i].Bankruptcy!=0) && (t>LearningPhase2) && ((t-LearningPhase2)%Year2>Week) && (100*Market[i].Capital()/Market[i].NAVOnJanFirst+Discount<=Market[i].LiquidationFloor)) {
                CountBrankrupcies+=1;
                outputLog << "At t=" << t << ", NAV of agent " << i << " reached " << 100*Market[i].Capital()/Market[i].NAVOnJanFirst << "% (while average agent performance was " << NAVOfNonBankruptAgents << "% since Jan 1st at t=" << t-(t-LearningPhase2)%Year2 << "), and got bankrupt (from NAV=$" << Market[i].Capital() << ") due to agent LiquidationFloor=" << Market[i].LiquidationFloor << "%" << endl;
                Market[i].Liquidation(i, Market);
            };
        };
        gsl_matrix_set(Bankruptcies, 0, t, 100.0*CountBrankrupcies/NumberOfAgents);
        
        
        
        if (Plot=="On") {outputLog << "t=" << t << ", Log of credit events completed..." << endl;};
        // outputDebug << "t=" << t << ": Credit events computed..." << endl;
        
        
        
        
        // // METAORDERS
        // // Counting the OB business levels
        // //if (t>1000) {Plot="On";};
        // int OBLevelSize=0;
        // int CountOB=0; int OBLevel=0;
        // for (int k=0; k<int(min(FullOrderBook[t][0].Bids.size(), FullOrderBook[t][0].Asks.size())); k++) {
        //     if ((FullOrderBook[t][0].Bids[k].Price < FullOrderBook[t][0].Asks[k].Price) && (OBLevel==0)) {OBLevel=1; OBLevelSize=k-1;};
        //     CountOB=1;
        // }; // closes k loop
        // //ofstream outputOBLevels(Machine+"OBlevels.log", ofstream::app); outputOBLevels << " At t=" << t << " OB size =" << OBLevelSize << endl;;
        // // Meta order injection
        // if ((t>500) && (t-LastMetaorder>=6*Month) && (MetaorderImpact>0)) {
        //     int Res=FullOrderBook[t][0].MetaorderInjection(Market, SharesOutstanding, MetaorderImpact, OBLevelSize);
        //     if (Res>-1) {
        //         LastMetaorder=t;
        //         outputMetalog << "t=" << t << ", metaorder successfully injected at OB level " << Res << ", last metaorder was at t=" << LastMetaorder << endl;
        //         if ((t-LearningPhase>1000) && (Time-t>7*Month)) {
        //             MetaorderInjectionTimes.push_back(LastMetaorder-LearningPhase);
        //             outputMetalog << "And metaorder recorded in MetaorderInjectionTimes[]..." << endl;
        //         };
        //     };
        // };
        
        
        // OUTPUTING THE ORDER BOOK AT TIME t (SS10)
        if (Plot=="On") {
            for (int j=0; j<NumberOfStocks; j++) {
                int CountOB=0;
                int OBLevel=0;
                outputLog << endl << "ORDER BOOK OF STOCK " << j << " AT TIME t=" << t << endl;
                if (int(min(FullOrderBook[t][j].Bids.size(), FullOrderBook[t][j].Asks.size()))==0) {outputLog << "t=" << t << ", the Order Book has zero size..." << endl;};
                for (int k=0; k<int(min(FullOrderBook[t][j].Bids.size(), FullOrderBook[t][j].Asks.size())); k++) { // ABCD
                    if ((FullOrderBook[t][j].Bids[k].Price < FullOrderBook[t][j].Asks[k].Price) && (OBLevel==0)) {
                        outputLog << "============================================================================================" << endl;
                        OBLevel=1;
                    };
                    CountOB=1;
                    outputLog << "Level " << k << ": Agent " << FullOrderBook[t][j].Bids[k].Agent << " bids " << FullOrderBook[t][j].Bids[k].Q << " Stocks " << FullOrderBook[t][j].Bids[k].Stock << " at $" << FullOrderBook[t][j].Bids[k].Price;
                    outputLog << " || ";
                    outputLog << "Agent " << FullOrderBook[t][j].Asks[k].Agent << " asks " << FullOrderBook[t][j].Asks[k].Q << " Stocks " << FullOrderBook[t][j].Asks[k].Stock << " at $" << FullOrderBook[t][j].Asks[k].Price << endl;
                    //}; // closes if
                }; // closes k loop
                if (CountOB==0) {outputLog << "Order book is empty!" << endl;};
                outputLog << endl;
            }; // closes j loop
            outputLog << "t=" << t << ", Outputed the Order Book at time t" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ": OB outputted..." << endl;
        
        
        
        

        
        // ACTUAL TRADING VIA CLEAR() METHOD OF THE ORDER BOOK (SS11)
        for (int j=0; j<NumberOfStocks; j++) {
            // REMOVING METAORDERS EFFECTS (CREDITING THE METAORDER AGENT WITH ITS PREVIOUS Q OR RFA) WHICH MAY BRING THIS AGENT BANKRUPT
            double CreditBeforeTrading=0; double CreditAfterTrading=0;
            // Before trading
            if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==0)) { // 0 for Bid and 1 for Ask
                CreditBeforeTrading=Market[FullOrderBook[t][0].MetaorderLastAgent].RFA;
                outputMetalog << "Before metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had a (meta) RFA of £" << CreditBeforeTrading << endl;
            }
            else if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==1)) { // 0 for Bid and 1 for Ask
                CreditBeforeTrading=double(Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity);
                outputMetalog << "Before metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had a (meta) number of " << int(CreditBeforeTrading) << " stocks" << endl;
            };
            // ACTUAL TRADING
            vector<double> OBClearing = FullOrderBook[t][j].Clear(Market, t, j, Plot); // Actual transations
            
            // After trading
            if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==0)) { // 0 for Bid and 1 for Ask
                CreditAfterTrading=Market[FullOrderBook[t][0].MetaorderLastAgent].RFA;
                Market[FullOrderBook[t][0].MetaorderLastAgent].RFA-=FullOrderBook[t][0].Credit-(CreditBeforeTrading-CreditAfterTrading);
                outputMetalog << "After metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had its unused cash of £" << FullOrderBook[t][0].Credit-(CreditBeforeTrading-CreditAfterTrading) << " removed and thus RFA lowered to £" << Market[FullOrderBook[t][0].MetaorderLastAgent].RFA << endl;
            }
            else if ((t==LastMetaorder) && (MetaorderImpact>0) && (FullOrderBook[t][0].MetaorderNature==1)) { // 0 for Bid and 1 for Ask
                CreditAfterTrading=double(Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity);
                Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity-=int(FullOrderBook[t][0].Credit)-(int(CreditBeforeTrading)-int(CreditAfterTrading));
                outputMetalog << "After metaorder transaction, agent " << FullOrderBook[t][0].MetaorderLastAgent << " had its unused number of " << int(FullOrderBook[t][0].Credit)-(int(CreditBeforeTrading)-int(CreditAfterTrading)) << " stocks removed and thus lowered to " << Market[FullOrderBook[t][0].MetaorderLastAgent].Stocks[0].StockQuantity << endl;
            };
            
            if (OBClearing[2]==0) { // If OB is empty
                gsl_matrix_set(ReflexiveValues, j, t, gsl_matrix_get (ReflexiveValues, j, t-1));
                gsl_matrix_set(AverageBids, j, t, gsl_matrix_get (AverageBids, j, t-1));
                gsl_matrix_set(AverageAsks, j, t, gsl_matrix_get (AverageAsks, j, t-1));
                gsl_matrix_set(TotalStockQuantityTraded, j, t, 0);
            } // closes if
            else {  // If OB is non-empty
                gsl_matrix_set(ReflexiveValues, j, t, OBClearing[0]); // Lowest of clearing prices is the market price, not the median
                // Calculating the average bids and asks at the top levels of the OB
                double MeanBid=0; int CountBid=int(FullOrderBook[t][j].Bids.size());
                for (int k=0; k<CountBid; k++) {MeanBid+=FullOrderBook[t][j].Bids[k].Price/CountBid;};
                double MeanAsk=0; int CountAsk=int(FullOrderBook[t][j].Asks.size());
                for (int k=0; k<CountAsk; k++) {MeanAsk+=FullOrderBook[t][j].Asks[k].Price/CountAsk;};
                if (CountBid==0) {MeanBid=gsl_matrix_get(AverageBids, j, t-1);};
                if (CountAsk==0) {MeanAsk=gsl_matrix_get(AverageAsks, j, t-1);};
                gsl_matrix_set(AverageBids, j, t, MeanBid);
                gsl_matrix_set(AverageAsks, j, t, MeanAsk);
                //gsl_matrix_set(TotalStockQuantityTraded, j, t, 10000*OBClearing[1]/TotalStockQuantities[j]); // Quantity of stock traded at time t in basis points of total quantity of stock RPR
                gsl_matrix_set(TotalStockQuantityTraded, j, t, OBClearing[1]); // Quantity of stock traded at time t in nominal
                if ((Plot=="On") && (Plot=="Mini")) {outputLog << "t=" << t << ", Actual trading via StockOrderBook.Clear()" << endl;};
            }; // closes else
            // FAILSAFE
            if (gsl_matrix_get(ReflexiveValues, j, t)<0.000001) {
                gsl_matrix_set(ReflexiveValues, j, t, 0.5*gsl_matrix_get(ReflexiveValues, j, t-1));
                outputLog << "At t=" << t << " P(t) was zero => now P(t)=P(t-1)/2" << endl;
            };
        }; // closes j loop
        if (Plot=="On") {outputLog << "t=" << t << ", Actual trading via StockOrderBook.Clear()" << endl;};
        // outputDebug << "t=" << t << ": Clear() computed..." << endl;
        
        
        
        
        
        //MetaorderLastAgent=Bids[MetaorderLevel].Agent; // Agent placing the order
        //MetaorderNature=0; // 0 for Bid and 1 for Ask
        //Credit=ceil(Bids[MetaorderLevel].Q*(Bids[MetaorderLevel].Price+Asks[MetaorderLevel].Price)/2); // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
        //Credit=double(Asks[MetaorderLevel].Q+1);
        
        
        
        // MARKET INDICATORS (SS12)
        // Capital evolution
        for (int i=0; i<NumberOfAgents ; i++) {
            double x = Market[i].Capital();
            double y = MarketEvolution[0][i].Capital();
            gsl_matrix_set(CapitalEvolution, i, t, (Market[i].Bankruptcy)*100*x/y);
        }; // closes i loop
        if (Plot=="On") {outputLog << "t=" << t << ", Indicator: CapitalEvolution" << endl;};
        
        // RELATIVE MAX DRAWDOWN
        // The (relative) max drawdown is the largest (adjusted to market) cumulated loss over a given period. We consider each year (after the second year) the lowest value of gsl_matrix* CapitalEvolution, and if this measure is below 100-Drawdown the agent reaches bankruptcy.
        
        
        // Computing the average of the biased values of all agents
        for (int j=0; j<NumberOfStocks; j++) {
            double Mean=0;
            for (int i=0; i<NumberOfAgents; i++) {
                Mean+=abs(-100+100*gsl_matrix_get (Market[i].BiasedValues, j, t)/gsl_matrix_get (ReflexiveValues, j, t)); // Absolute percentage of difference between agent biased and market price
            }; // closes i loop
            gsl_matrix_set (BiasedDelta, j, t, Mean/NumberOfAgents);
            gsl_matrix_set (TrueDelta, j, t, abs(-100+100*gsl_matrix_get (GSLTrueValues, j, t)/gsl_matrix_get (ReflexiveValues, j, t))); // Absolute percentage of difference between true and market price
        }; // closes j loop
        
        
        
        // CONSOLE INFO SCREEN
        TickEnd=(int(time(0))-TickStart); // Closing the clock for time of computation
        //system("CLS");
        cout << "****************************" << endl;
        cout << "** STOCK MARKET SIMULATOR **" << endl;
        cout << "****************************" << endl;
        cout << "I=" << NumberOfAgents << ", J=" << NumberOfStocks << ", T=" << Time << ", s=" << s << endl  << endl;
        cout.setf(ios::fixed);
        cout << "Initialization : " << TickMid << "s" << endl;
        cout << "Timestep       : " << t+1 << "/" << Time << endl;
        cout << "Elapsed time   : " << TickEnd/3600 << " h " << TickEnd/60 << " min " << TickEnd%60 << " s" << endl;
        cout << "Remaining time : " << (Time*TickEnd/t - TickEnd)/3600 << " h " << ((Time*TickEnd/t - TickEnd)%3600)/60 << " min " << (Time*TickEnd/t - TickEnd)%60 << " s" << endl << endl;
        cout << "Progress : " << 100*double(t)/(double(Time)) << " %" << endl;
        if (TickEnd%4==0) {cout << "                                    |" << endl;}
        else if (TickEnd%4==1) {cout << "                                    /" << endl;}
        else if (TickEnd%4==2) {cout << "                                    --" << endl;}
        else if (TickEnd%4==3) {cout << "                                    /" << endl;};
        if ((PDCondition=="PDOn") && (t>=1000) && ((t-1000)%Year2==0)) {outputLog << "Computing PDistances..." << endl;};
            
        // END OF t-LOOP (SS13)
        if (Plot=="On") {
            outputLog << "t=" << t << ", End of t-loop" << endl;
            outputLog << "***************************************************" << endl << endl << endl;
            // outputDebug << "t=" << t << ": End of t-loop" << endl;
        }; // closes Plot condition
    
    
        // Policy distances
        if ((PDCondition=="PDOn") && (t>=1000) && ((t-1000)%Year2==0)) {
            PDs.push_back(PolicyDistances (Market, NumberOfAgents, NumberOfStocks, Time, Percentile, t));
        };
        
        
    }; // closes t loop
    cout << "Closed t-loop..." << endl;
    outputMetalog.close();
    
    
    // Policy heat maps for best and worst populations
    if (PDCondition=="PDOn") {PoliciesAverages (Market, NumberOfAgents, NumberOfStocks, Percentile);};
    
    // Outputing the agents parameters and especially policies for eventual re-use
    if (OutputName!="OutputCondOff") {AgentParametersOutput (Market, NumberOfAgents, NumberOfStocks, OutputName);};
    
    // Outputting the percentage during the simulation of best and worst agents that belonged to the agent population within ClusterLimit
    outputClustered << 100*ClusteredBestAgent/(Time-LearningPhase) << '\t' << 100*ClusteredWorstAgent/(Time-LearningPhase) << endl;
    
    
    
    
    // MarketBidAskTopBidTopAskTrue of j stocks, and Spread of j stocks
    gsl_matrix* MarketBidAskTrue = gsl_matrix_alloc (12*NumberOfStocks, Time);
    gsl_matrix* MarketBidAskTrue2 = gsl_matrix_alloc (12*NumberOfStocks, Time-LearningPhase); // Without LearningPhase
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=12*j;
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (MarketBidAskTrue, 0+Modulo, t, gsl_matrix_get (ReflexiveValues, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 1+Modulo, t, gsl_matrix_get (AverageBids, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 2+Modulo, t, gsl_matrix_get (AverageAsks, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 1+Modulo, t, gsl_matrix_get (AverageTopBids, j, t));
            //gsl_matrix_set (MarketBidAskTrue, 2+Modulo, t, gsl_matrix_get (AverageTopAsks, j, t));
            gsl_matrix_set (MarketBidAskTrue, 1+Modulo, t, gsl_matrix_get (GSLTrueValues, j, t));
            gsl_matrix_set (MarketBidAskTrue, 2+Modulo, t, gsl_matrix_get (HighestBid, j, t));
            gsl_matrix_set (MarketBidAskTrue, 3+Modulo, t, gsl_matrix_get (LowestAsk, j, t));
            if (t==0) {gsl_matrix_set (MarketBidAskTrue, 4+Modulo, t, 1);}
            else {gsl_matrix_set (MarketBidAskTrue, 4+Modulo, t, gsl_matrix_get (ReflexiveValues, j, t)/gsl_matrix_get (ReflexiveValues, j, t-1));}; // Returns P(t)/P(t-1)
            gsl_matrix_set (Spread, j, t, 100*VSpread[j]/gsl_matrix_get (ReflexiveValues, j, t));
            gsl_matrix_set (MarketBidAskTrue, 5+Modulo, t, 100*VSpread[j]/gsl_matrix_get (ReflexiveValues, j, t)); // Spread pct
            //gsl_matrix_set (MarketBidAskTrue, 6+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)); // Quantity of stock traded at time t in basis points of total quantity of stock RPR
            gsl_matrix_set (MarketBidAskTrue, 6+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)); // Quantity of stock traded at time t in nominal
            gsl_matrix_set (MarketBidAskTrue, 7+Modulo, t, gsl_matrix_get(Bankruptcies, 0, t)); // Bankruptcy pct
            gsl_matrix_set (MarketBidAskTrue, 8+Modulo, t, log(gsl_matrix_get(MarketBidAskTrue, 4+Modulo, t))); // Log-returns log[P(t)/P(t-1)]
            gsl_matrix_set (MarketBidAskTrue, 9+Modulo, t, abs(log(gsl_matrix_get(MarketBidAskTrue, 4+Modulo, t)))); // Absolute log-returns abs(log[P(t)/P(t-1)])
            //if ((t==0) || (gsl_matrix_get (TotalStockQuantityTraded, j, t-1)==0)) {gsl_matrix_set (MarketBidAskTrue, 10+Modulo, t, 1);}
            //else {gsl_matrix_set (MarketBidAskTrue, 10+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)/gsl_matrix_get (TotalStockQuantityTraded, j, t-1));}; // Volume returns
            gsl_matrix_set (MarketBidAskTrue, 10+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t)); // Volume
            gsl_matrix_set (MarketBidAskTrue, 11+Modulo, t, gsl_matrix_get(Bankruptcies, 1, t)); // Market performance in % since Jan 1st
        }; // closes t loop
    }; // closes j loop
    
    // Populating MarketBidAskTrue2, which is like MarketBidAskTrue but without the first LearningPhase time steps
    for (int j=0; j<12*NumberOfStocks; j++) {
        for (int t=0; t<Time-LearningPhase; t++) {
            gsl_matrix_set(MarketBidAskTrue2, j, t, gsl_matrix_get(MarketBidAskTrue, j, t+LearningPhase));
        };
    };
    
    
    
    // Metrics: log-returns (& AC), absolute log-returns (& AC), volatilities (& AC), volumes (& AC)
    // Calibration: we can make distributions of stacked real/simulated log-returns, absolute log-returns, volatilities, volumes, but not AC's thereof.
    gsl_matrix* Moments = gsl_matrix_calloc (54*NumberOfStocks, Time);
    gsl_matrix* Moments2 = gsl_matrix_calloc (54*NumberOfStocks, Time-LearningPhase);
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=12*j; // Size1 of MarketBidAskTrue above
        int ModuloMoments=54*j; // Size1 of Moments above
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (Moments, 0+ModuloMoments, t, gsl_matrix_get (MarketBidAskTrue, 8+Modulo, t)); // Log-returns log[P(t)/P(t-1)]
            gsl_matrix_set (Moments, 1+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Week, 8+Modulo)); // AC of log-returns at lag of Week
            gsl_matrix_set (Moments, 2+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 2*Week, 8+Modulo)); // AC of log-returns at lag of 2*Week
            gsl_matrix_set (Moments, 3+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Month, 8+Modulo)); // AC of log-returns at lag of Month
            gsl_matrix_set (Moments, 4+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 3*Month, 8+Modulo)); // AC of log-returns at lag of 3*Month
            gsl_matrix_set (Moments, 5+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 6*Month, 8+Modulo)); // AC of log-returns at lag of 6*Month
            gsl_matrix_set (Moments, 6+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Year, 8+Modulo)); // AC of log-returns at lag of Year
            gsl_matrix_set (Moments, 7+ModuloMoments, t, gsl_matrix_get (MarketBidAskTrue, 9+Modulo, t)); // Absolute log-returns abs(log[P(t)/P(t-1)])
            gsl_matrix_set (Moments, 8+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, Week, 9+Modulo))); // AC of absolute log-returns at lag of Week
            gsl_matrix_set (Moments, 9+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, 2*Week, 9+Modulo))); // AC of absolute log-returns at lag of 2*Week
            gsl_matrix_set (Moments, 10+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, Month, 9+Modulo))); // AC of absolute log-returns at lag of Month
            gsl_matrix_set (Moments, 11+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, 3*Month, 9+Modulo))); // AC of absolute log-returns at lag of 3*Month
            gsl_matrix_set (Moments, 12+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, 6*Month, 9+Modulo))); // AC of absolute log-returns at lag of 6*Month
            gsl_matrix_set (Moments, 13+ModuloMoments, t, abs(MALAutoCorrelation (MarketBidAskTrue, t, Year, 9+Modulo))); // AC of absolute log-returns at lag of Year
            gsl_matrix_set (Moments, 14+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, Week, 0+Modulo)); // Volatility at lag of Week
            gsl_matrix_set (Moments, 16+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, 2*Week, 0+Modulo)); // Volatility at lag of 2*Week
            gsl_matrix_set (Moments, 18+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, Month, 0+Modulo)); // Volatility at lag of Month
            gsl_matrix_set (Moments, 20+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, 3*Month, 0+Modulo)); // Volatility at lag of 3*Month
            gsl_matrix_set (Moments, 22+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, 6*Month, 0+Modulo)); // Volatility at lag of 6*Month
            gsl_matrix_set (Moments, 24+ModuloMoments, t, MALVolatility (MarketBidAskTrue, t, Year, 0+Modulo)); // Volatility at lag of Year
            gsl_matrix_set (Moments, 26+ModuloMoments, t, gsl_matrix_get (MarketBidAskTrue, 10+Modulo, t)); // Volumes returns
            gsl_matrix_set (Moments, 27+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Week, 10+Modulo)); // AC of volumes returns at lag of Week
            gsl_matrix_set (Moments, 28+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 2*Week, 10+Modulo)); // AC of volumes returns at lag of 2*Week
            gsl_matrix_set (Moments, 29+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Month, 10+Modulo)); // AC of volumes returns at lag of Month
            gsl_matrix_set (Moments, 30+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 3*Month, 10+Modulo)); // AC of volumes returns at lag of 3*Month
            gsl_matrix_set (Moments, 31+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, 6*Month, 10+Modulo)); // AC of volumes returns at lag of 6*Month
            gsl_matrix_set (Moments, 32+ModuloMoments, t, MALAutoCorrelation (MarketBidAskTrue, t, Year, 10+Modulo)); // AC of volumes returns at lag of Year
            
            
            
            // AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 33+ModuloMoments, t, MALAutoCorrelation (Moments, t, 3*Week, 0+ModuloMoments)); // AC of log-returns at lag of Week
            // AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 2, 3, 4, w} and ∆={w, 2w, 3w, m}
            gsl_matrix_set(Moments, 34+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 1, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 1
            gsl_matrix_set(Moments, 35+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 2, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 2
            gsl_matrix_set(Moments, 36+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 3, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 3
            gsl_matrix_set(Moments, 37+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, 4, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by 4
            gsl_matrix_set(Moments, 38+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Week, Week, 0+ModuloMoments)); // AC of log-returns at lag of Week shifted by Week
            
            gsl_matrix_set(Moments, 39+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*1, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 2
            gsl_matrix_set(Moments, 40+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*2, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 4
            gsl_matrix_set(Moments, 41+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*3, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 6
            gsl_matrix_set(Moments, 42+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*4, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 8
            gsl_matrix_set(Moments, 43+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 2*Week, 2*Week, 0+ModuloMoments)); // AC of log-returns at lag of 2Week shifted by 2Week
            
            gsl_matrix_set(Moments, 44+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*1, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 3
            gsl_matrix_set(Moments, 45+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*2, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 6
            gsl_matrix_set(Moments, 46+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*3, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 9
            gsl_matrix_set(Moments, 47+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*4, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 12
            gsl_matrix_set(Moments, 48+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, 3*Week, 3*Week, 0+ModuloMoments)); // AC of log-returns at lag of 3Week shifted by 3Week
            
            gsl_matrix_set(Moments, 49+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*1, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 4
            gsl_matrix_set(Moments, 50+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*2, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 8
            gsl_matrix_set(Moments, 51+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*3, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 12
            gsl_matrix_set(Moments, 52+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, 4*4, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by 16
            gsl_matrix_set(Moments, 53+ModuloMoments, t, MALAutoCorrelationBlend (Moments, t, Month, Month, 0+ModuloMoments)); // AC of log-returns at lag of Month shifted by Month
            
            
        }; // closes t loop
    }; // closes j loop
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=12*j; // Size1 of MarketBidAskTrue above
        int ModuloMoments=54*j; // Size1 of Moments above
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (Moments, 15+ModuloMoments, t, MALAutoCorrelation (Moments, t, Week, 14+Modulo)); // AC of volatility at lag of Week for lag of Week
            gsl_matrix_set (Moments, 17+ModuloMoments, t, MALAutoCorrelation (Moments, t, 2*Week, 16+Modulo)); // AC of volatility at lag of 2*Week for lag of 2*Week
            gsl_matrix_set (Moments, 19+ModuloMoments, t, MALAutoCorrelation (Moments, t, Month, 18+Modulo)); // AC of volatility at lag of Month for lag of Month
            gsl_matrix_set (Moments, 21+ModuloMoments, t, MALAutoCorrelation (Moments, t, 3*Month, 20+Modulo)); // AC of volatility at lag of 3*Month for lag of 2*Month
            gsl_matrix_set (Moments, 23+ModuloMoments, t, MALAutoCorrelation (Moments, t, 6*Month, 22+Modulo)); // AC of volatility at lag of 6*Month for lag of 6*Month
            gsl_matrix_set (Moments, 25+ModuloMoments, t, MALAutoCorrelation (Moments, t, Year, 24+Modulo)); // AC of volatility at lag of Year for lag of Year
        }; // closes t loop
    }; // closes j loop
    
    // Populating Moments2, which is like Moments but without the first LearningPhase time steps
    for (int j=0; j<54*NumberOfStocks; j++) {
        for (int t=0; t<Time-LearningPhase; t++) {
            gsl_matrix_set(Moments2, j, t, gsl_matrix_get(Moments, j, t+LearningPhase));
        };
    };
    
    // Meta order impact
    // We want to see at a given time t where a large order is sent (say once per year for each simulation) the metaorder impact on prices i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] (and likewise i2(t) on volumes) as a function of the metaorder size j(t)=Q_order(t)/Q_total(t), sigma being the standard deviation (of prices or volumes) over an interval of size tau, and Q being the stock quantity. For each of the S=20 simulation runs, we record and output the last 5 such tuples (j, i1(tau=w), i1(tau=2w), i1(tau=m), i1(tau=3m), i1(tau=6m), i2(tau=w), i2(tau=2w), i2(tau=m), i2(tau=3m), i2(tau=6m)), disregarding the first 5 because of being not yet learned by the agents. We do 4 such S=20 simulation runs for values of j=5%, 10%, 25%, 50%. The metaorder is sent during each run once every year by a random (yet non-bankrupt) agent whose stock holding is suddenly increased to these values (5%, 10%, 25%, 50%) of the total share outstanding at time t and then revert back to its value at time t+1.
    gsl_matrix* MetaorderImpactResults = gsl_matrix_calloc (16, int(MetaorderInjectionTimes.size()));
    int Tau=Week;
    for (int k=0; k<int(MetaorderInjectionTimes.size()); k++) {
        // Quantities
        gsl_matrix_set (MetaorderImpactResults, 0, k, 100*gsl_matrix_get (MarketBidAskTrue2, 10, MetaorderInjectionTimes[k])/SharesOutstanding); // j(t)=Q_order(t)/Q_total(t)
        //Prices
        Tau=Week; gsl_matrix_set (MetaorderImpactResults, 1, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=we
        Tau=2*Week; gsl_matrix_set (MetaorderImpactResults, 2, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=2we
        Tau=3*Week; gsl_matrix_set (MetaorderImpactResults, 3, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=mo
        Tau=4*Week; gsl_matrix_set (MetaorderImpactResults, 4, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=3mo
        Tau=5*Week; gsl_matrix_set (MetaorderImpactResults, 5, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=6mo
        Tau=6*Week; gsl_matrix_set (MetaorderImpactResults, 6, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=6mo
        Tau=7*Week; gsl_matrix_set (MetaorderImpactResults, 7, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 0)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 0)); // i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] for tau=6mo
        //Volumes
        Tau=Week; gsl_matrix_set (MetaorderImpactResults, 8, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=we
        Tau=2*Week; gsl_matrix_set (MetaorderImpactResults, 9, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=2we
        Tau=3*Week; gsl_matrix_set (MetaorderImpactResults, 10, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=mo
        Tau=4*Week; gsl_matrix_set (MetaorderImpactResults, 11, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=3mo
        Tau=5*Week; gsl_matrix_set (MetaorderImpactResults, 12, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=6mo
        Tau=6*Week; gsl_matrix_set (MetaorderImpactResults, 13, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=6mo
        Tau=7*Week; gsl_matrix_set (MetaorderImpactResults, 14, k, -100+100*MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k], MetaorderInjectionTimes[k]+Tau, 10)/MALVolatility2 (MarketBidAskTrue2, MetaorderInjectionTimes[k]-Tau, MetaorderInjectionTimes[k], 10)); // i1(t)=sigma_volumes[t->t+tau]/sigma_volumes[t-tau->t] for tau=6mo
        // Injection times
        gsl_matrix_set (MetaorderImpactResults, 15, k, MetaorderInjectionTimes[k]); // Times of metaorder injection
    };
    
    

    
    
    // Bull market: -20% then +20% then -20% =>
    // Crash: sudden -20% => return Pt-20% and then no Pt+10% after a month
    // Recession: sudden -20% with no recovery =>
    gsl_matrix* SystemicRisk = gsl_matrix_calloc (15*NumberOfStocks, Time-LearningPhase); // Without LearningPhase
    gsl_matrix* Prices = gsl_matrix_calloc (1, Time-LearningPhase); for (int t=0; t<Time-LearningPhase; t++) {gsl_matrix_set (Prices, 0, t, gsl_matrix_get (ReflexiveValues, 0, t+LearningPhase));};
    for (int j=0; j<NumberOfStocks; j++) {
        int Modulo=15*j; int Modulo2=54*j; double Crash=0; double TrendCounter=0;
        for (int t=0; t<Time-LearningPhase; t++) {
            gsl_matrix_set (SystemicRisk, 0+Modulo, t, gsl_matrix_get (ReflexiveValues, j, t+LearningPhase)); // Price ($)
            gsl_matrix_set (SystemicRisk, 1+Modulo, t, gsl_matrix_get (BiasedDelta, j, t+LearningPhase)); // Percentage of average of differences between biased and market prices (%)
            gsl_matrix_set (SystemicRisk, 2+Modulo, t, gsl_matrix_get (TrueDelta, j, t+LearningPhase)); // Percentage of average of differences between biased and market prices (%)
            double FormalSpread=gsl_matrix_get(HighestBid, j, t+LearningPhase) - gsl_matrix_get(LowestAsk, j, t+LearningPhase); // Formal spread definition ($)
            gsl_matrix_set (SystemicRisk, 3+Modulo, t, 100*FormalSpread/gsl_matrix_get (ReflexiveValues, j, t+LearningPhase)); // Formal spread as a percentage of market price (%)
            if (t==0) {gsl_matrix_set (SystemicRisk, 4+Modulo, t, 0);}
            else {gsl_matrix_set (SystemicRisk, 4+Modulo, t, -100+100*gsl_matrix_get (ReflexiveValues, j, t+LearningPhase)/gsl_matrix_get (ReflexiveValues, j, t-1+LearningPhase));}; // Percentage of returns -100+100*P(t)/P(t-1) in
            gsl_matrix_set (SystemicRisk, 5+Modulo, t, gsl_matrix_get (TotalStockQuantityTraded, j, t+LearningPhase)); // Quantity of stock traded at time t in basis points of total quantity of stock RPR
            gsl_matrix_set (SystemicRisk, 6+Modulo, t, gsl_matrix_get (Bankruptcies, 0, t+LearningPhase)); // Bankruptcy pct
            gsl_matrix_set (SystemicRisk, 7+Modulo, t, gsl_matrix_get (Bankruptcies, 1, t+LearningPhase)); // Market performance in % since Jan 1st
            gsl_matrix_set (SystemicRisk, 8+Modulo, t, gsl_matrix_get (Moments, 14+Modulo2, t+LearningPhase)); // Week normalized volatility
            gsl_matrix_set (SystemicRisk, 9+Modulo, t, gsl_matrix_get (Moments, 18+Modulo2, t+LearningPhase)); // Month normalized volatility
            gsl_matrix_set (SystemicRisk, 10+Modulo, t, gsl_matrix_get (Moments, 22+Modulo2, t+LearningPhase)); // 6 month normalized volatility
            double Mean1=0; for (int k=t-2*Week; k<t-Week; k++) {Mean1+=gsl_matrix_get (ReflexiveValues, j, k+LearningPhase)/Week;};
            double Mean2=0; for (int k=t-Week; k<t; k++) {Mean2+=gsl_matrix_get (ReflexiveValues, j, k+LearningPhase)/Week;};
            gsl_matrix_set (SystemicRisk, 11+Modulo, t, -100+100*Mean2/Mean1); // -100+100*<P(t-w,t)/P(t-2w,t-w)> Allows to see crashes (<-20% in the distribution)
            if (gsl_matrix_get (SystemicRisk, 11+Modulo, t)<=-20) {Crash=1;} else {Crash=0;}; gsl_matrix_set (SystemicRisk, 12+Modulo, t, Crash); // Crash
            if (t>0) {
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)>0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)<0)) {TrendCounter=0;}; // Trend reversion
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)<0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)>0)) {TrendCounter=0;}; // Trend reversion
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)>0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)>0)) {TrendCounter+=1;}; // Trend reversion
                if ((gsl_matrix_get (SystemicRisk, 11+Modulo, t-1)<0) && (gsl_matrix_get (SystemicRisk, 11+Modulo, t)<0)) {TrendCounter-=1;}; // Trend reversion
            };
            gsl_matrix_set (SystemicRisk, 13+Modulo, t, TrendCounter);
        }; // closes t loop
        
        
        PlotGSLMatrix(Prices, "Prices.xls", 1);
        //gsl_matrix* Hawkes = HistoricalIntensity(Prices, 0);
        //for (int t=0; t<Time-LearningPhase; t++) {gsl_matrix_set(SystemicRisk, 14+Modulo, t, gsl_matrix_get (Hawkes, 0, t));};
        
        // vector<double> V, GoF; for (int t=0; t<Time-LearningPhase; t++) {V.push_back(gsl_matrix_get (Prices, 0, t)); GoF.push_back(1);};
        // vector<double>* pV= &V;
        // const double Rescale=10; // Rescale factor
        // const double Threshold=0.1; // Threshold (0.1)
        // const size_t NumIteration=10000000; // Number of iterations of the minimizer (stops before once a minimum is found)
        // double R; double* pR = &R; // Hawkes process parameters (does not need to be initialized)
        // vector<double>& pGoF = GoF; // Goodness of fit
        // vector<double>* pRes = Intensity (pV, Rescale, Threshold, NumIteration, pR, pGoF);
        // for (int t=0; t<Time-LearningPhase; t++) {gsl_matrix_set (SystemicRisk, 14, t, (*pRes)[t]);}; // Hawkes process intensity in column 9
        
        // for (int t=0; t<Time-LearningPhase; t++) { // Hawkes process intensity in column 9
        //     if ((*pRes)[t]>10) { // Maximum value of Hawkes at 10
        //         gsl_matrix_set (SystemicRisk, 14, t, 10);
        //     }
        //     else if ((*pRes)[t]<0) { // Minimum value of Hawkes at 0
        //         gsl_matrix_set (SystemicRisk, 14, t, 0);
        //     };
        // };
        
        // pV=0; pR=0; pRes=0;
        // //delete pV; delete pR; delete pRes;
        
        
        
        
        
    }; // closes j loop
    
    
    
    
    // Volumes of j stocks, VolumesDistributions of j stocks, AutoCorrLagPVolumes of j stocks
    ReturnsDistributions = GSLDistribution(MarketBidAskTrue, 20);
    VolumesDistributions = GSLDistribution(TotalStockQuantityTraded, 1000);
    AutoCorrLagPVolumes = AutocorrelationLagP (TotalStockQuantityTraded);
    
    // Most successful and less successful agents
    vector<double> SummedCapitalEvolution;
    for (int i=0; i<NumberOfAgents; i++) {
        SummedCapitalEvolution.push_back(0);
        //for (int t=Time/2; t<Time; t++) {SummedCapitalEvolution[i]+=gsl_matrix_get(CapitalEvolution, i, t);};
        for (int t=0; t<Time; t++) {SummedCapitalEvolution[i]+=gsl_matrix_get(CapitalEvolution, i, t);};
    }; // closes i loop
    vector<int> DescendingRank;
    for (int k=0; k<NumberOfAgents; k++) {
        double Bar=0; int BarX=0;
        for (int i=0; i<NumberOfAgents; i++) {
            if ((SummedCapitalEvolution[i]>=Bar) && (SummedCapitalEvolution[i]>=0)) {Bar=SummedCapitalEvolution[i]; BarX=i;};
        }; // closes i loop
        SummedCapitalEvolution[BarX]=-1;
        DescendingRank.push_back(BarX);
    }; // closes k loop
    int AgentNumberPercentile=Percentile*NumberOfAgents/100;
    if (AgentNumberPercentile<=5) {AgentNumberPercentile=5;};
    gsl_matrix* MostSuccessful = gsl_matrix_alloc (AgentNumberPercentile, Time);
    gsl_matrix* LessSuccessful = gsl_matrix_alloc (AgentNumberPercentile, Time);
    for (int t=0; t<Time; t++) {
        for (int i=0; i<AgentNumberPercentile; i++) {
            gsl_matrix_set (MostSuccessful, i, t, gsl_matrix_get(CapitalEvolution, DescendingRank[i], t));};
        for (int i=0; i<AgentNumberPercentile; i++) {
            gsl_matrix_set (LessSuccessful, i, t, gsl_matrix_get(CapitalEvolution, DescendingRank[i+NumberOfAgents-AgentNumberPercentile], t));};
    }; // closes t loop
    
    
    //Plot the distribution of each parameter of the AgentNumberPercentile% most/less successful agents (cf. Versatility, Gesture, and all NEBs: NEB1, NEBLossAversion, NEB3, NEB4, NEBPositivity, NEB6, NEB1p, NEB3p, NEB4p, NEBPositivity, NEBLearningRate, NEB8, NEB9).
    gsl_matrix* MostSuccessfulParameters = gsl_matrix_alloc (7, AgentNumberPercentile);
    gsl_matrix* LessSuccessfulParameters = gsl_matrix_alloc (7, AgentNumberPercentile);
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 0, i, Market[DescendingRank[i]].Future);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 0, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].Future);};
    gsl_matrix_set (MostSuccessfulParameters, 1, 0, 0); // To make sure the distributions will be between 0 and 100!! (X1<=X<X2)
    gsl_matrix_set (MostSuccessfulParameters, 1, 1, 100); // To make sure the distributions will be between 0 and 100!!
    gsl_matrix_set (LessSuccessfulParameters, 1, 0, 0); // To make sure the distributions will be between 0 and 100!!
    gsl_matrix_set (LessSuccessfulParameters, 1, 1, 100); // To make sure the distributions will be between 0 and 100!!
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 1, i, Market[DescendingRank[i]].Reflexivity*100);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 1, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].Reflexivity*100);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 2, i, Market[DescendingRank[i]].TradingWindow);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 2, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].TradingWindow);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 3, i, Market[DescendingRank[i]].Gesture);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 3, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].Gesture);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 4, i, Market[DescendingRank[i]].History);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 4, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].History);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 5, i, Market[DescendingRank[i]].NEBLearningRate);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 5, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].NEBLearningRate);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (MostSuccessfulParameters, 6, i, Market[DescendingRank[i]].NEB);};
    for (int i=0; i<AgentNumberPercentile; i++) {gsl_matrix_set (LessSuccessfulParameters, 6, i, Market[DescendingRank[i+NumberOfAgents-AgentNumberPercentile]].NEB);};
    
    
    
    // Ploting the distributions of these parameters
    gsl_matrix* MostSuccessfulParametersDistributions = gsl_matrix_alloc (AgentNumberPercentile, 2*7);
    gsl_matrix* LessSuccessfulParametersDistributions = gsl_matrix_alloc (AgentNumberPercentile, 2*7);
    MostSuccessfulParametersDistributions = GSLDistribution(MostSuccessfulParameters, 10);
    LessSuccessfulParametersDistributions = GSLDistribution(LessSuccessfulParameters, 10);
    //int pBots, int pNEBLossAversion24, int pNEBDelayDiscounting56, int pNEB89, int pNEBLossAversion23456, int pNEBLossAversion2489, int pNEBDelayDiscounting5689
    // Zipf's Law and Pareto's Principle
    gsl_matrix* SortedNAV = gsl_matrix_alloc (1, NumberOfAgents);
    gsl_matrix* ZipfDistribution = gsl_matrix_alloc (2, 50);
    gsl_matrix* Pareto = gsl_matrix_alloc (2, NumberOfAgents); // Percentage of agents who own a percentage of total market capital (80% should own around 20% according to Pareto's Principle)
    double TotalCapital=0;
    for (int i=0; i<NumberOfAgents; i++) {
        TotalCapital+=Market[DescendingRank[i]].Capital();
    };
    double TempTotalCapital=0;
    double ParetoPercentage=0;
    for (int i=0; i<NumberOfAgents; i++) {
        gsl_matrix_set (SortedNAV, 0, i, Market[DescendingRank[i]].Capital());
        TempTotalCapital+=100*(Market[DescendingRank[i]].Capital())/TotalCapital;
        ParetoPercentage=100*i/NumberOfAgents;
        gsl_matrix_set (Pareto, 0, i, ParetoPercentage);
        gsl_matrix_set (Pareto, 1, i, TempTotalCapital);
    };
    ZipfDistribution = GSLDistribution(SortedNAV, 50); // Distribution
    gsl_matrix* MAINDistribution = GSLDistribution(MarketBidAskTrue2, 50); // Distribution of MAIN
    gsl_matrix* MostSuccessfulDistribution = GSLDistribution(MostSuccessfulParameters, 50); // Distribution of MostSuccessfulParameters
    gsl_matrix* LessSuccessfulDistribution = GSLDistribution(LessSuccessfulParameters, 50); // Distribution of LessSuccessfulParameters
    gsl_matrix* MomentsDistribution = GSLDistribution(Moments2, 50); // Distribution of Moments
    
    
    // Output // XYZ
    ofstream outputMain(Machine+"MAIN.xls", ofstream::app);
    outputMain << "Market($)" << '\t' << "True($)" << '\t' << "Bid($)" << '\t' << "Ask($)" << '\t' << "Return" << '\t' << "Spread(%)" << '\t' << "Volume(bsp)" << '\t' << "Banruptcy(%)" << '\t' << "LogReturns" << '\t' << "AbsLogReturns" << '\t' << "Volumes(bsp)" << '\t' << "MarketPerfYTD(%)" << endl;
    PlotGSLMatrix(MarketBidAskTrue2, "MAIN.xls", 1); // Contains info for all j stocks: market, average bids, average asks, top bids, top asks, True, daily returns pct, spread pct, volume bsp
    PlotGSLMatrix(MAINDistribution, "MAINDistribution.xls", 1);
    outputMain.close();
    //PlotGSLMatrix(MarketBidAskTrue, "MarketBidAskTrueReturnspctSpreadpctVolumebspBankrupcypct.xls", 1); // Contains info for all j stocks: market, average bids, average asks, top bids, top asks, True, daily returns pct, spread pct, volume bsp
    //PlotGSLMatrix(Spread, "MSpread.xls", 1); // Contains info for all j stocks
    //PlotGSLMatrix(ReturnsDistributions, "ReturnsDistributions.xls", 1);
    //PlotGSLMatrix(TotalStockQuantityTraded, "MVolumes.xls", 1);
    //PlotGSLMatrix(VolumesDistributions, "VolumesDistributions.xls", 1);
    //PlotGSLMatrix(AutoCorrLagPVolumes, "AutoCorrLagPVolumes.xls", 1);
    //PlotGSLMatrix(MostSuccessful, "MostSuccessful.xls", 1);
    //PlotGSLMatrix(LessSuccessful, "LessSuccessful.xls", 1);
    
    ofstream outputSystemic(Machine+"SystemicRisk.xls", ofstream::app);
    outputSystemic << "Market($)" << '\t' << "Biased/Pt(%)" << '\t' << "True/Pt(%)" << '\t' << "FormalSpread(%)" << '\t' << "Return(%)" << '\t' << "Volume(bsp)" << '\t' << "Banruptcy(%)" << '\t' << "MarketPerfYTD(%)" << '\t' << "w-Volatility" << '\t' << "m-Volatility" << '\t' << "6m-Volatility" << '\t' << "P[t-w,t]/P[t-2w,t-w](%)" << '\t' << "Crash" << '\t' << "TrendCounter" << '\t' << "Hawkes" << endl;
    PlotGSLMatrix(SystemicRisk, "SystemicRisk.xls", 1); // Contains info for all j stocks: market, average bids, average asks, top bids, top asks, True, daily returns pct, spread pct, volume bsp
    outputSystemic.close();
    
    ofstream outputSystemicDistribution(Machine+"SystemicRiskDistribution.xls", ofstream::app);
    outputSystemicDistribution<< "Market($)" << '\t' << " " << '\t' << "Biased/Pt(%)" << '\t' << " " << '\t' << "True/Pt(%)" << '\t' << " " << '\t' << "FormalSpread(%)" << '\t' << " " << '\t' << "Return(%)" << '\t' << " " << '\t' << "Volume(bsp)" << '\t' << " " << '\t' << "Banruptcy(%)" << '\t' << " " << '\t' << "MarketPerfYTD(%)" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "P[t-w,t]/P[t-2w,t-w](%)" << '\t' << " " << '\t' << "Crash" << '\t' << " " << '\t' << "TrendCounter" << '\t' << " " << '\t' << "Hawkes" << endl;
    PlotGSLMatrix(GSLDistribution(SystemicRisk, 50), "SystemicRiskDistribution.xls", 1);
    outputSystemicDistribution.close();
    
    
    
    //MetaorderImpactResults
    //We want to see at a given time t where a large order is sent (say once per year for each simulation) the metaorder impact on prices i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] (and likewise i2(t) on volumes) as a function of the metaorder size j(t)=Q_order(t)/Q_total(t), sigma being the standard deviation (of prices or volumes) over an interval of size tau, and Q being the stock quantity. For each of the S=20 simulation runs, we record and output the last 5 such tuples (j, i1(tau=w), i1(tau=2w), i1(tau=m), i1(tau=3m), i1(tau=6m), i2(tau=w), i2(tau=2w), i2(tau=m), i2(tau=3m), i2(tau=6m)), disregarding the first 5 because of being not yet learned by the agents. We do 4 such S=20 simulation runs for values of j=5%, 10%, 25%, 50%. The metaorder is sent during each run once every year by a random (yet non-bankrupt) agent whose stock holding is suddenly increased to these values (5%, 10%, 25%, 50%) of the total share outstanding at time t and then revert back to its value at time t+1.
    ofstream outputMetaorderImpactResults(Machine+"MetaorderImpact_size" + to_string(MetaorderImpact) + "pct.xls", ofstream::app);
    outputMetaorderImpactResults<< "Q/Qtot(%)" << '\t' << "Psig(t,t+w)/Psig(t-w,t)(%)" << '\t' << "Psig(t,t+2w)/Psig(t-2w,t)(%)" << '\t' << "Psig(t,t+3w)/Psig(t-3w,t)(%)" << '\t' << "Psig(t,t+4w)/Psig(t-4w,t)(%)" << '\t' << "Psig(t,t+5w)/Psig(t-5w,t)(%)" << '\t' << "Psig(t,t+6w)/Psig(t-6w,t)(%)" << '\t' << "Psig(t,t+7w)/Psig(t-7w,t)(%)" << '\t' << "Vsig(t,t+w)/Vsig(t-w,t)(%)" << '\t' << "Vsig(t,t+2w)/Vsig(t-2w,t)(%)" << '\t' << "Vsig(t,t+3w)/Vsig(t-3w,t)(%)" << '\t' << "Vsig(t,t+4w)/Vsig(t-4w,t)(%)" << '\t' << "Vsig(t,t+5w)/Vsig(t-5w,t)(%)" << '\t' << "Vsig(t,t+6w)/Vsig(t-6w,t)(%)" << '\t' << "Vsig(t,t+7w)/Vsig(t-7w,t)(%)" << '\t' << "t" << endl;
    PlotGSLMatrix(MetaorderImpactResults, "MetaorderImpact_size" + to_string(MetaorderImpact) + "pct.xls", 1);
    outputMetaorderImpactResults.close();
    
    
    ofstream outputSuccessful(Machine+"MostLessSuccessfulParameters.xls", ofstream::app);
    outputSuccessful<< "Future" << '\t' << "Reflexivity(%)" << '\t' << "TradingWindow" << '\t' << "Gesture" << '\t' << "History" << '\t' << "NEBLearningRate" << '\t' << "NEB" << endl;
    PlotGSLMatrix(MostSuccessfulParameters, "MostLessSuccessfulParameters.xls", 1);
    PlotGSLMatrix(LessSuccessfulParameters, "MostLessSuccessfulParameters.xls", 1);
    PlotGSLMatrix(MostSuccessfulDistribution, "MostLessSuccessfulParametersDistribution.xls", 1);
    PlotGSLMatrix(LessSuccessfulDistribution, "MostLessSuccessfulParametersDistribution.xls", 1);
    
    PlotGSLMatrix(CapitalEvolution, "CapitalEvolution.xls", 1);
    
    outputSuccessful.close();
    //PlotGSLMatrix(MostSuccessfulParametersDistributions, "MostSuccessfulParametersDistributions.xls", 1);
    //PlotGSLMatrix(LessSuccessfulParametersDistributions, "LessSuccessfulParametersDistributions.xls", 1);
    PlotGSLMatrix(SortedNAV, "Zipf&Distribution.xls", 1);
    PlotGSLMatrix(ZipfDistribution, "Zipf&Distribution.xls", 1);
    ofstream outputMoments(Machine+"Moments.xls", ofstream::app);
    outputMoments << "Log-return" << '\t' << "AC-1w Log-return" << '\t' << "AC-2w Log-return" << '\t' << "AC-m Log-return" << '\t' << "AC-3m Log-return" << '\t' << "AC-6m Log-return" << '\t' << "AC-y Log-return" << '\t' << "Abs-log-return" << '\t' << "AC-1w Abs-log-return" << '\t' << "AC-2w Abs-log-return" << '\t' << "AC-m Abs-log-return" << '\t' << "AC-3m Abs-log-return" << '\t' << "AC-6m Abs-log-return" << '\t' << "AC-y Abs-log-return" << '\t' << "w-Volatility" << '\t' << "AC-w w-Volatility" << '\t' << "2w-Volatility" << '\t' << "AC-2w 2w-Volatility" << '\t' << "m-Volatility" << '\t' << "AC-m m-Volatility" << '\t' << "3m-Volatility" << '\t' << "AC-3m 3m-Volatility" << '\t' << "6m-Volatility" << '\t' << "AC-6m 6m-Volatility" << '\t' << "y-Volatility" << '\t' << "AC-y y-Volatility" << '\t' << "Volumes(bsp)" << '\t' << "AC-1w Volume" << '\t' << "AC-2w Volume" << '\t' << "AC-m Volume" << '\t' << "AC-3m Volume" << '\t' << "AC-6m Volume" << '\t' << "AC-y Volume" << '\t' << "AC-3w Log-return" << '\t' << "b1AC-1w Log-return" << '\t' << "b2AC-1w Log-return" << '\t' << "b3AC-1w Log-return" << '\t' << "b4AC-1w Log-return" << '\t' << "bwAC-1w Log-return" << '\t' << "b2AC-2w Log-return" << '\t' << "b4AC-2w Log-return" << '\t' << "b6AC-2w Log-return" << '\t' << "b8AC-2w Log-return" << '\t' << "b2wAC-2w Log-return" << '\t' << "b3AC-3w Log-return" << '\t' << "b6AC-3w Log-return" << '\t' << "b9AC-3w Log-return" << '\t' << "b12AC-3w Log-return" << '\t' << "b3wAC-3w Log-return" << '\t' << "b4AC-m Log-return" << '\t' << "b8AC-m Log-return" << '\t' << "b12AC-m Log-return" << '\t' << "b16AC-m Log-return" << '\t' << "bmAC-m Log-return" << endl;
    
    PlotGSLMatrix(Moments2, "Moments.xls", 1);
    ofstream outputMomentsDistribution(Machine+"MomentsDistribution.xls", ofstream::app);
    outputMomentsDistribution << "Log-return" << '\t' << " " << '\t' << "AC-1w Log-return" << '\t' << " " << '\t' << "AC-2w Log-return" << '\t' << " " << '\t' << "AC-m Log-return" << '\t' << " " << '\t' << "AC-3m Log-return" << '\t' << " " << '\t' << "AC-6m Log-return" << '\t' << " " << '\t' << "AC-y Log-return" << '\t' << " " << '\t' << "Abs-log-return" << '\t' << " " << '\t' << "AC-1w Abs-log-return" << '\t' << " " << '\t' << "AC-2w Abs-log-return" << '\t' << " " << '\t' << "AC-m Abs-log-return" << '\t' << " " << '\t' << "AC-3m Abs-log-return" << '\t' << " " << '\t' << "AC-6m Abs-log-return" << '\t' << " " << '\t' << "AC-y Abs-log-return" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "AC-w w-Volatility" << '\t' << " " << '\t' << "2w-Volatility" << '\t' << " " << '\t' << "AC-2w 2w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "AC-m m-Volatility" << '\t' << " " << '\t' << "3m-Volatility" << '\t' << " " << '\t' << "AC-3m 3m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "AC-6m 6m-Volatility" << '\t' << " " << '\t' << "y-Volatility" << '\t' << " " << '\t' << "AC-y y-Volatility" << '\t' << " " << '\t' << "Volumes(bsp)" << '\t' << " " << '\t' << "AC-1w Volume" << '\t' << " " << '\t' << "AC-2w Volume" << '\t' << " " << '\t' << "AC-m Volume" << '\t' << " " << '\t' << "AC-3m Volume" << '\t' << " " << '\t' << "AC-6m Volume" << '\t' << " " << '\t' << "AC-y Volume" << '\t' << " " << '\t' << "AC-3w Log-return" << '\t' << " " << '\t' << "b1AC-1w Log-return" << '\t' << " " << '\t' << "b2AC-1w Log-return" << '\t' << " " << '\t' << "b3AC-1w Log-return" << '\t' << " " << '\t' << "b4AC-1w Log-return" << '\t' << " " << '\t' << "bwAC-1w Log-return" << '\t' << " " << '\t' << "b2AC-2w Log-return" << '\t' << " " << '\t' << "b4AC-2w Log-return" << '\t' << " " << '\t' << "b6AC-2w Log-return" << '\t' << " " << '\t' << "b8AC-2w Log-return" << '\t' << " " << '\t' << "b2wAC-2w Log-return" << '\t' << " " << '\t' << "b3AC-3w Log-return" << '\t' << " " << '\t' << "b6AC-3w Log-return" << '\t' << " " << '\t' << "b9AC-3w Log-return" << '\t' << " " << '\t' << "b12AC-3w Log-return" << '\t' << " " << '\t' << "b3wAC-3w Log-return" << '\t' << " " << '\t' << "b4AC-m Log-return" << '\t' << " " << '\t' << "b8AC-m Log-return" << '\t' << " " << '\t' << "b12AC-m Log-return" << '\t' << " " << '\t' << "b16AC-m Log-return" << '\t' << " " << '\t' << "bmAC-m Log-return" << endl;
    PlotGSLMatrix(MomentsDistribution, "MomentsDistribution.xls", 1);
    //PlotGSLMatrix(Pareto, "Pareto.xls", 1);
    //cout << "Produced outputs..." << endl;
    
    
    
    // Memory freeing // JJJ10
    gsl_matrix_free(Moments);
    gsl_matrix_free(MarketBidAskTrue);
    gsl_matrix_free(MarketBidAskTrue2);
    gsl_matrix_free(Prices);
    gsl_matrix_free(ReflexiveValues);
    gsl_matrix_free(AverageBids);
    gsl_matrix_free(AverageAsks);
    gsl_matrix_free(AverageTopBids);
    gsl_matrix_free(AverageTopAsks);
    gsl_matrix_free(MedianTopBids);
    gsl_matrix_free(MedianTopAsks);
    gsl_matrix_free(HighestBid);
    gsl_matrix_free(LowestAsk);
    gsl_matrix_free(GSLTrueValues);
    gsl_matrix_free(Spread);
    gsl_matrix_free(TotalStockQuantityTraded);
    gsl_matrix_free(VolumesDistributions);
    gsl_matrix_free(ReturnsDistributions);
    gsl_matrix_free(AutoCorrLagPVolumes);
    gsl_matrix_free(CapitalEvolution);
    gsl_matrix_free(Bankruptcies);
    gsl_matrix_free(BiasedDelta);
    gsl_matrix_free(TrueDelta);
    gsl_matrix_free(MomentsDistribution);
    gsl_matrix_free(ZipfDistribution);
    gsl_matrix_free(LessSuccessfulDistribution);
    gsl_matrix_free(MostSuccessfulDistribution);
    gsl_matrix_free(MAINDistribution);
    gsl_matrix_free(MetaorderImpactResults);
    
    for (int i=0; i<NumberOfAgents; i++) {gsl_matrix_free(Market[i].BiasedValues);};
    
    
    
    outputLog.close();
    cout << "Outputted .xls files..." << endl;
    
    // BACKUP OF MATRIX RESULTS
    vector<gsl_matrix*> MatrixResults;
    MatrixResults.push_back(Moments2);
    MatrixResults.push_back(MostSuccessfulParameters);
    MatrixResults.push_back(LessSuccessfulParameters);
    MatrixResults.push_back(SortedNAV);
    MatrixResults.push_back(SystemicRisk);
    for (int k=0; k<int(PDs.size()); k++) {MatrixResults.push_back(PDs[k]);}; // PDistances
    
    
    Market.erase (Market.begin(), Market.end()); // JJJ10
    
    cout << "DONE" << endl;
    
    return MatrixResults;
};


/// Each run of MarketSimulator() outputs vector<gsl_matrix*> MatrixResults. We then compute S of these runs and save them in a vector of length S called vector<vector<gsl_matrix*>> MultiSim. Now function Bollinger() gets that MultiSim and accesses the m matrix of all its MatrixResults simulations (m=0 corresponds to MarketBidAskTrue, m=1 to MostSuccessfulParameters, m=2 to LessSuccessfulParameters, m=3 to SortedNAV). It then accesses all the cxr elements of all S matrix m and outputs their mean in a matrix of same dimension than matrix m.
vector<gsl_matrix*> Bollinger (vector<vector<gsl_matrix*>> MultiSim, int m) {
    vector<gsl_matrix* > Res;
    gsl_matrix* ResMean = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix* ResMeanPlusSigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix* ResMeanMinusSigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix* ResMeanPlus2Sigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    gsl_matrix* ResMeanMinus2Sigma = gsl_matrix_calloc (MultiSim[0][m]->size1, MultiSim[0][m]->size2);
    int MSsize=int(MultiSim.size());
    for (int c=0; c<int(MultiSim[0][m]->size1); c++) {
        for (int r=0; r<int(MultiSim[0][m]->size2); r++) {
            double Mean=0; for (int s=0; s<MSsize; s++) {Mean+=gsl_matrix_get (MultiSim[s][m], c, r)/MSsize;};
            gsl_matrix_set (ResMean, c, r, Mean);
            double Variance=0; for (int s=0; s<MSsize; s++) {Variance+=(Mean - gsl_matrix_get (MultiSim[s][m], c, r)) * (Mean - gsl_matrix_get (MultiSim[s][m], c, r))/MSsize;};
            gsl_matrix_set (ResMeanPlusSigma, c, r, Mean+sqrt(Variance));
            gsl_matrix_set (ResMeanMinusSigma, c, r, Mean-sqrt(Variance));
            gsl_matrix_set (ResMeanPlusSigma, c, r, Mean+2*sqrt(Variance));
            gsl_matrix_set (ResMeanMinusSigma, c, r, Mean-2*sqrt(Variance));
        }; // closes r loop
    }; // closes c loop
    string A = "BollingerMean_m"; string B = to_string(m); string C = ".xls"; A+=B+C; PlotGSLMatrix(ResMean, A.c_str(), 1);
    A = "BollingerMeanPlusSigma_m"; B = to_string(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanPlusSigma, A.c_str(), 1);
    A = "BollingerMeanMinusSigma_m"; B = to_string(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanMinusSigma, A.c_str(), 1);
    A = "BollingerMeanPlus2Sigma_m"; B = to_string(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanPlus2Sigma, A.c_str(), 1);
    A = "BollingerMeanMinus2Sigma_m"; B = to_string(m); C = ".xls"; A+=B+C; PlotGSLMatrix(ResMeanMinus2Sigma, A.c_str(), 1);
    Res.push_back(ResMean); Res.push_back(ResMeanPlusSigma); Res.push_back(ResMeanMinusSigma); Res.push_back(ResMeanPlus2Sigma); Res.push_back(ResMeanMinus2Sigma);
    return Res;
}; // closes Bollinger()



/// Computing mean, variance, skewness, kurtosis, median, q5, q25, q75, q95 of every column of a matrix M
gsl_matrix* Momenta (gsl_matrix* M) {
    int Size1=int(M->size1);
    int Size2=int(M->size2);
    gsl_matrix* Res = gsl_matrix_alloc (2*Size1, 9);
    for (int i=0; i<Size1; i++) {
        double Mean=0; double Variance=0; double Stdev=0; double Skewness=0; double Kurtosis=0;
        for (int j=0; j<Size2; j++) {Mean+=gsl_matrix_get (M, i, j)/Size2;};
        for (int j=0; j<Size2; j++) {Variance+=(gsl_matrix_get (M, i, j)-Mean)*(gsl_matrix_get (M, i, j)-Mean)/Size2;}; Stdev=sqrt(Variance);
        for (int j=0; j<Size2; j++) {
            Skewness+=pow(((gsl_matrix_get (M, i, j)-Mean)/Stdev), 3);
            Kurtosis+=pow(((gsl_matrix_get (M, i, j)-Mean)/Stdev), 4);
        }; // closes j loop
        vector<double> V; for (int j=0; j<Size2; j++) {V.push_back(gsl_matrix_get (M, i, j));};
        sort (V.begin(), V.begin()+Size2); // Sorting V ascendingly
        double Median=V[50*Size2/100];
        double q5=V[5*Size2/100];
        double q25=V[25*Size2/100];
        double q75=V[75*Size2/100];
        double q95=V[95*Size2/100];
        // Filling the matrix result Res
        gsl_matrix_set (Res, 1+2*i, 0, Mean);
        gsl_matrix_set (Res, 1+2*i, 1, Variance);
        gsl_matrix_set (Res, 1+2*i, 2, Skewness);
        gsl_matrix_set (Res, 1+2*i, 3, Kurtosis);
        gsl_matrix_set (Res, 1+2*i, 4, Median);
        gsl_matrix_set (Res, 1+2*i, 5, q5);
        gsl_matrix_set (Res, 1+2*i, 6, q25);
        gsl_matrix_set (Res, 1+2*i, 7, q75);
        gsl_matrix_set (Res, 1+2*i, 8, q95);
        V.erase(V.begin(), V.end());
    }; // closes i loop
    return Res;
}; // closes Momenta()



/// Each run of MarketSimulator() outputs vector<gsl_matrix*> MatrixResults. We then compute S of these runs and save them in a vector of length S called vector<vector<gsl_matrix*>> MultiSim. Now function JointDistributions() gets that MultiSim and accesses the m matrix of all its MatrixResults simulations. It then accesses the [c,r]-element of all S matrices m and outputs their joint distributions over all S matrices and output a matrix similar to m but with twice its number of columns (each column has a range for 10 bins between Xmin and Xmax, and another for the number of counts per bin)
void JointDistributions (vector<vector<gsl_matrix*>> MultiSim, int m, int Precision) {
    int MSsize = int(MultiSim.size());
    int Size1=int(MultiSim[0][m]->size1);
    int Size2=int(MultiSim[0][m]->size2);
    gsl_matrix* Res = gsl_matrix_alloc (Size1, (MSsize*Size2));
    for (int s=0; s<MSsize; s++) {
        for (int c=0; c<Size1; c++) {
            for (int r=0; r<Size2; r++) {
                gsl_matrix_set (Res, c, s*Size2+r, gsl_matrix_get (MultiSim[s][m], c, r));
            }; // closes r loop
        }; // closes c loop
    }; // closes s loop
    //PlotGSLMatrix(Res, "StackedRes.xls", 1);
    gsl_matrix* JointDistributions = GSLDistribution(Res, Precision);
    //string A = "JointDistributions_m"; string B = to_string(m); string C = ".xls"; A+=B+C;
    string A = "JointDistributions.xls";
    ofstream outputJointDistribution(Machine+"JointDistributions.xls", ofstream::app);
    if (m==0) {
        outputJointDistribution<< "Moments" << endl;
        outputJointDistribution << "Log-return" << '\t' << " " << '\t' << "AC-1w Log-return" << '\t' << " " << '\t' << "AC-2w Log-return" << '\t' << " " << '\t' << "AC-m Log-return" << '\t' << " " << '\t' << "AC-3m Log-return" << '\t' << " " << '\t' << "AC-6m Log-return" << '\t' << " " << '\t' << "AC-y Log-return" << '\t' << " " << '\t' << "Abs-log-return" << '\t' << " " << '\t' << "AC-1w Abs-log-return" << '\t' << " " << '\t' << "AC-2w Abs-log-return" << '\t' << " " << '\t' << "AC-m Abs-log-return" << '\t' << " " << '\t' << "AC-3m Abs-log-return" << '\t' << " " << '\t' << "AC-6m Abs-log-return" << '\t' << " " << '\t' << "AC-y Abs-log-return" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "AC-w w-Volatility" << '\t' << " " << '\t' << "2w-Volatility" << '\t' << " " << '\t' << "AC-2w 2w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "AC-m m-Volatility" << '\t' << " " << '\t' << "3m-Volatility" << '\t' << " " << '\t' << "AC-3m 3m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "AC-6m 6m-Volatility" << '\t' << " " << '\t' << "y-Volatility" << '\t' << " " << '\t' << "AC-y y-Volatility" << '\t' << " " << '\t' << "Volumes(bsp)" << '\t' << " " << '\t' << "AC-1w Volume" << '\t' << " " << '\t' << "AC-2w Volume" << '\t' << " " << '\t' << "AC-m Volume" << '\t' << " " << '\t' << "AC-3m Volume" << '\t' << " " << '\t' << "AC-6m Volume" << '\t' << " " << '\t' << "AC-y Volume" << '\t' << " " << '\t' << "AC-3w Log-return" << '\t' << " " << '\t' << "b1AC-1w Log-return" << '\t' << " " << '\t' << "b2AC-1w Log-return" << '\t' << " " << '\t' << "b3AC-1w Log-return" << '\t' << " " << '\t' << "b4AC-1w Log-return" << '\t' << " " << '\t' << "bwAC-1w Log-return" << '\t' << " " << '\t' << "b2AC-2w Log-return" << '\t' << " " << '\t' << "b4AC-2w Log-return" << '\t' << " " << '\t' << "b6AC-2w Log-return" << '\t' << " " << '\t' << "b8AC-2w Log-return" << '\t' << " " << '\t' << "b2wAC-2w Log-return" << '\t' << " " << '\t' << "b3AC-3w Log-return" << '\t' << " " << '\t' << "b6AC-3w Log-return" << '\t' << " " << '\t' << "b9AC-3w Log-return" << '\t' << " " << '\t' << "b12AC-3w Log-return" << '\t' << " " << '\t' << "b3wAC-3w Log-return" << '\t' << " " << '\t' << "b4AC-m Log-return" << '\t' << " " << '\t' << "b8AC-m Log-return" << '\t' << " " << '\t' << "b12AC-m Log-return" << '\t' << " " << '\t' << "b16AC-m Log-return" << '\t' << " " << '\t' << "bmAC-m Log-return" << endl;
    }
    else if (m==1) {
        outputJointDistribution<< "MostSuccessfulParameters" << endl;
        outputJointDistribution<< "Future" << '\t' << " " << '\t' << "Reflexivity(%)" << '\t' << " " << '\t' << "TradingWindow" << '\t' << " " << '\t' << "Gesture" << '\t' << " " << '\t' << "History" << '\t' << " " << '\t' << "Epsilon" << '\t' << " " << '\t' << "NEB" << endl;
    }
    else if (m==2) {
        outputJointDistribution<< "LessSuccessfulParameters" << endl;
        outputJointDistribution<< "Future" << '\t' << " " << '\t' << "Reflexivity(%)" << '\t' << " " << '\t' << "TradingWindow" << '\t' << " " << '\t' << "Gesture" << '\t' << " " << '\t' << "History" << '\t' << " " << '\t' << "Epsilon" << '\t' << " " << '\t' << "NEB" << endl;
    }
    else if (m==3) {
        outputJointDistribution<< "Pareto" << endl;
        outputJointDistribution<< "NAV" << '\t' << " " << '\t' << endl;
    }
    else if (m==4) {
        outputJointDistribution << "Systemics" << endl;
        outputJointDistribution<< "Market($)" << '\t' << " " << '\t' << "Biased/Pt(%)" << '\t' << " " << '\t' << "True/Pt(%)" << '\t' << " " << '\t' << "FormalSpread(%)" << '\t' << " " << '\t' << "Return(%)" << '\t' << " " << '\t' << "Volume(bsp)" << '\t' << " " << '\t' << "Banruptcy(%)" << '\t' << " " << '\t' << "MarketPerfYTD(%)" << '\t' << " " << '\t' << "w-Volatility" << '\t' << " " << '\t' << "m-Volatility" << '\t' << " " << '\t' << "6m-Volatility" << '\t' << " " << '\t' << "P[t-w,t]/P[t-2w,t-w](%)" << '\t' << " " << '\t' << "Crash" << '\t' << " " << '\t' << "TrendCounter" << '\t' << " " << '\t' << "Hawkes" << endl;
    }
    else {
        outputJointDistribution << "PDistances_y" << m-5 << endl;
        outputJointDistribution<< "BestToBest" << '\t' << " " << '\t' << "BestToAll" << '\t' << " " << '\t' << "BestToWorst" << '\t' << " " << '\t' << "WorstToAll" << '\t' << " " << '\t' << "WorstToWorst" << endl;
    };
    
    PlotGSLMatrix(JointDistributions, A.c_str(), 1);
    PlotGSLMatrix(Momenta(Res), A.c_str(), 1);

}; // closes JointDistributions()


/// Takes a string name of the csv file as input (N2 can be something like "AAPL.csv") and returns a vector<double> for a matrix of dimension MxN at batch m (in a file of stacked matrices). For Input() to work we must 1- pre-specify the dimensions of the input, 2- the input must be converted from xls to csv (by changing the name file, not via LibreOffice)
gsl_matrix* InputOLD (int M, int N, string S1, int m) {
    gsl_matrix* Res = gsl_matrix_alloc (M, N);
    string S2=Machine; S2+=S1; ifstream FileInput (S2);
    string Line;
    int j = -1 + m*(N+3);
    while (getline (FileInput, Line)) {
        j++;
        istringstream LineStream(Line);
        string Item;
        int i = -1;
        while (getline (LineStream, Item, '\t')) {
            //while (getline (LineStream, Item, ',')) {
            i++;
            if (j<N) {gsl_matrix_set (Res, i, j-m*(N+3), atof(Item.c_str()));};
        }
    }
    string S3="M"; string S3a=to_string(m); string S3b=".txt"; S3+=S3a+S3b;
    PlotGSLMatrix(Res, S3.c_str(), 1);
    FileInput.close();
    return Res;
};




/// Takes a string name of the csv file as input (N2 can be something like "AAPL.csv") and returns a vector<double> for a matrix of dimension MxN at batch m (in a file of stacked matrices). For Input() to work we must 1- pre-specify the dimensions of the input, 2- the input must be converted from xls to csv (by changing the name file, not via LibreOffice)
gsl_matrix* Input (int M, int N, string S1, int m, int Letter, int Opt, int Period) {
    gsl_matrix* Res = gsl_matrix_alloc (M, N);
    string S2=Machine; S2+=S1; ifstream FileInput (S2);
    string Line;
    //for (int j=m*(N+2); j<m*(N+2)+Period; j++) {
    //for (int j=m*(N+2)+N-Period; j<m*(N+2)+N; j++) {
    for (int j=0; j<m*(N+2)+N; j++) {
        getline (FileInput, Line);
        istringstream LineStream(Line);
        string Item;
        int i = -1;
        if (Opt==0) {
            while (getline (LineStream, Item, '\t')) {
                i++;
                //gsl_matrix_set (Res, i, j-(m*(N+2)+N-Period), atof(Item.c_str()));
                if (j>=m*(N+2)) {
                    gsl_matrix_set (Res, i, j-(m*(N+2)), atof(Item.c_str()));
                };
            }
        }
        else {
            while (getline (LineStream, Item, ',')) {
                i++;
                gsl_matrix_set (Res, i, j-(m*(N+2)+N-Period), atof(Item.c_str()));
            }
        }
    }
    gsl_matrix* Res2 = gsl_matrix_alloc (M, Period);
    for (int i=0; i<M; i++) {
        for (int j=N-Period; j<N; j++) {
            gsl_matrix_set (Res2, i, j-N+Period, gsl_matrix_get (Res, i, j));
        };
    };
    string S3="/OUTPUT/M"; string S3a=to_string(m+Letter); string S3b=".txt"; S3+=S3a+S3b;
    PlotGSLMatrix(Res2, S3.c_str(), 1);
    FileInput.close();
    return Res2;
};




/// Takes a string name of the csv file as input (N2 can be something like "AAPL.csv") and returns a vector<double> for a matrix of dimension MxN at batch m (in a file of stacked matrices). For Input() to work we must 1- pre-specify the dimensions of the input, 2- the input must be converted from xls to csv (by changing the name file, not via LibreOffice)
gsl_matrix* CSVInput (int M, int N, string S1, int Opt) {
    gsl_matrix* Res = gsl_matrix_alloc (M, N);
    string S2=Machine; S2+=S1; ifstream FileInput (S2);
    string Line;
    for (int j=0; j<=N; j++) {
        getline (FileInput, Line);
        istringstream LineStream(Line);
        string Item;
        int i = -1;
        if (Opt==0) {
            while (getline (LineStream, Item, '\t')) {
                i++;
                gsl_matrix_set (Res, i, j, atof(Item.c_str()));
            }
        }
        else {
            while (getline (LineStream, Item, ',')) {
                i++;
                if ((j==0) || (i==0)) {continue;};
                gsl_matrix_set (Res, i-1, j-1, atof(Item.c_str()));
            }
        }
    }
    // Inverting the matrix column elements
    gsl_matrix* Res2 = gsl_matrix_alloc (M, N);
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            gsl_matrix_set (Res2, i, j, gsl_matrix_get (Res, i, N-j-1));
        };
    };
    string S3="/OUTPUT/M_"; string S3b=".txt"; S3+=S1+S3b; PlotGSLMatrix(Res, S3.c_str(), 1);
    FileInput.close();
    return Res2;
};


/// Compute the factorial of an integer.
double Factorial (int n) {
    double Res = 1;

    for (int i = 1; i <= n; i++) {
        Res *= i;
    }

    return Res;
}


/**
 * @brief Compute the number of k-combinations of n elements.
 * In a set of n elements, the number of possible subsets of k elements (when non-ordered and without repetition) is
 * given by the number of k-combination C(n,k)=n!/(k!(n-k)!).
 */
double Combination (int n, int k) {
    double Res=Factorial(n)/(Factorial(k) * Factorial(n-k));
    return Res;
}


/// Number of all poissble k-combinations of a set of n elements.
double TotalCombination (int n) {
    double Res=n;
    for (int k=2; k<=n; k++) {Res+=Combination(n,k);}
    return Res;
}


/// Number of all poissble k-combinations of a set of n elements.
void TotalCombinationMatrix (int N) {
    gsl_matrix* Res = gsl_matrix_calloc (1, N+1);
    for (int k=1; k<=N; k++) {
        gsl_matrix_set (Res, 0, k, TotalCombination (k));
        cout << "Combination at N=" << k << "/" << N << endl;
    }
    PlotGSLMatrix(Res, "TotalCombinationMatrix.csv", 1);
};
//TotalCombinationMatrix(100); exit(0);


void HPMarketSimulatorAll (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S) {
    vector<vector<gsl_matrix*>> MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix*> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", "PDOff", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
    JointDistributions (MultiSim, 0, 100); // Moments2
    JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
    JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
    JointDistributions (MultiSim, 3, 100); // SortedNAV
    JointDistributions (MultiSim, 4, 100); // SystemicRisk
    for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
    for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
    for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
    for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
    for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    */
};
void HPMarketSimulatorAll2 (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, string PD, double Rate, int S) {
    vector<vector<gsl_matrix*>> MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix*> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", PD, TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};
void MarketSimulatorTrunk (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S, int Trunk) {
    vector<vector<gsl_matrix*>> MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix*> MSim = MarketSimulator (HPI, J, HPTime, r, "On", "PDOff", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, Trunk, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};

void HPMarketSimulatorAllPD (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S) {
    vector<vector<gsl_matrix*>> MultiSim;
    // int NEB0=50; int NEB1=13; int NEBLossAversion=13; int NEB3=12; int NEB4=12; // resp. bots, conservatism bias, loss aversion, positivity bias
    // HPAccuracy=9, 10, 11, 13 => probability of 0.2, 0.75, 1, 2 per year of bubble burst
    int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
    for (int s=0; s<S; s++) {
        //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
        vector<gsl_matrix*> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", "PDOn", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, 0);
        MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
    };
    for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
    // Memory freeing
    for (int k=0; k<int(MultiSim[0].size()); k++) {
        for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
    };
    MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
    /*
     JointDistributions (MultiSim, 0, 100); // Moments2
     JointDistributions (MultiSim, 1, 100); // MostSuccessfullParameters
     JointDistributions (MultiSim, 2, 100); // LessSuccessfullParameters
     JointDistributions (MultiSim, 3, 100); // SortedNAV
     JointDistributions (MultiSim, 4, 100); // SystemicRisk
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][0]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][1]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][2]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][3]);}; // JJJ10
     for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][4]);}; // JJJ10
     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
     */
};

void HPMarketSimulator (int T, int S) {HPMarketSimulatorAll (500, 0, 0.10, 10, T, 50, "Classic", "None", 1, 0, 1, S);};

// void HPMarketSimulator2 (string TypeNEB, string Leader, int ClusterLimit, int S) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Classic", Leader, ClusterLimit, 0, 1, S);};

// void HPMarketSimulatorNEB (string TypeNEB, int S) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, TypeNEB, "None", 1, 0, 1, S);};

// void HPMarketSimulatorDummy () {HPMarketSimulatorAll (500, 0, 0.10, 10, 1000+1000, 50, "Classic", "None", 1, 0, 1, 3);};





// void HPMarketSimulatorAllOB (int HPI, int HPGesture, double HPTrueMu, int HPAccuracy, int HPTime, int HPLiquidationFloor, string TypeNEB, string Leader, int ClusterLimit, int pNEB, double Rate, int S, int MetaOrderImpact) {
//     vector<vector<gsl_matrix*>> MultiSim;
//     int J=1; double r=0.01; // int T=LearningPhase+2875; int LearningPhase=1000;
//     for (int s=0; s<S; s++) {
//         //MarketSimulator (int NumberOfAgents, int NumberOfStocks, int Time, double Rate, string Plot, string PDCondition, string TypeNEB, int HPGesture, double HPTrueMu, int HPAccuracy, int LiquidationFloor, string LeaderType, int ClusterLimit, int s)
//         vector<gsl_matrix*> MSim = MarketSimulator (HPI, J, HPTime, r, "Off", "PDOff", TypeNEB, HPGesture, HPTrueMu, HPAccuracy, HPLiquidationFloor, Leader, ClusterLimit, pNEB, Rate, s, -1, MetaOrderImpact);
//         MultiSim.push_back(MSim); //delete MSim[0]; MSim[0]=NULL; MSim.erase(MSim.begin(), MSim.end());
//     };
//     for (int k=0; k<int(MultiSim[0].size()); k++) {JointDistributions (MultiSim, k, 100);};
//     // Memory freeing
//     for (int k=0; k<int(MultiSim[0].size()); k++) {
//         for (int i=0; i<S; i++) {gsl_matrix_free(MultiSim[i][k]);};
//     };
//     MultiSim.erase(MultiSim.begin(), MultiSim.end()); // Simulated moments distribution
// };



/*
// Memory freeing check
gsl_matrix* Coucou1 (gsl_matrix* M) {
    gsl_matrix* Res = gsl_matrix_calloc (100, 100);
    for (int i=0; i<int(M->size1); i++) {
        for (int j=0; j<int(M->size2); j++) {gsl_matrix_set (Res, i, j, gsl_matrix_get (M, i, j));};
    };
    return Res;
};
gsl_matrix* Coucou2 (gsl_matrix* M) {
    gsl_matrix* Res = gsl_matrix_calloc (100, 100);
    for (int i=0; i<int(M->size1); i++) {
        for (int j=0; j<int(M->size2); j++) {gsl_matrix_set (Res, i, j, gsl_matrix_get (M, i, j));};
    };
    return Res;
};
gsl_matrix* MemoryFreeingTest () {
    gsl_matrix* A = gsl_matrix_calloc (100, 100);
    for (int i=0; i<int(A->size1); i++) {
        for (int j=0; j<int(A->size2); j++) {gsl_matrix_set (A, i, j, i*j);};
    };
    gsl_matrix* Temp = Coucou1(A); gsl_matrix* Res = Coucou2(Temp);
    gsl_matrix_free(A); gsl_matrix_free(Temp);
    return Res;
};
*/




/*
vector<double> CheckTrueGeneration (double HPTrueMu, double HPLeash) {
    vector<double> Res;
    // INITIALIZING AGENTS
    int NumberOfAgents=10; int NumberOfStocks=1; int Time=1*Year; int HPAccuracy=10;
    gsl_rng* r = make_rng(); gsl_rng_set(r, static_cast<unsigned long int>(time(0)));
    gsl_matrix* GSLTrueValues = gsl_matrix_calloc (NumberOfStocks, Time);
    //gsl_matrix* Spikes = gsl_matrix_calloc (1, NumberOfStocks);
    gsl_matrix* TrueBiased = gsl_matrix_calloc (NumberOfAgents+1, Time); // For j=0 only
    vector<Agent> Market;
    for (int i=0; i<NumberOfAgents; i++) {
        Agent TempAgent; TempAgent.BiasedValues = gsl_matrix_calloc (NumberOfStocks, Time);
        for (int j=0; j<NumberOfStocks; j++) {
            TempAgent.Leash.push_back(HPLeash*0.001*gsl_rng_uniform (r));
            TempAgent.LeashVol.push_back(0.05*gsl_rng_uniform (r));
            TempAgent.Accuracy.push_back(HPAccuracy+gsl_rng_uniform (r));
        };
        Market.push_back(TempAgent);
    };
    // GENERATING TRUE VALUES
    vector<double> Gen = STLRandom (NumberOfStocks, "Uniform");
    double AnnualCounter=0; double AnnualAmplitude=0;
    vector<double> VAnnualCounter, VAnnualAmplitude;
    for (int j=0; j<NumberOfStocks; j++) {
        vector<double> S = PoissonRandomWalk (100, Week+HPTrueMu*3*Month*Gen[j], int(0.1*Time/250+0.9*Gen[j]*Time/250), Time, 1000*Gen[j]);
        for (int t=0; t<Time; t++) {
            gsl_matrix_set (GSLTrueValues, j, t, S[t]);
            gsl_matrix_set (TrueBiased, 0, t, S[t]);
        };
        int Counter=0; double Amplitude=0;
        for (int t=1; t<Time; t++) {
            if (S[t]-S[t-1]>0.001) {Counter+=1; Amplitude+=100*(S[t]-S[t-1])/S[t];};
        };
        AnnualCounter+=Counter/(Time/Year);
        AnnualAmplitude+=Amplitude/Counter;
        VAnnualCounter.push_back(Counter/(Time/Year));
        VAnnualAmplitude.push_back(Amplitude/Counter);
        //cout << "Counter=" << Counter/(Time/Year) << ", Amplitude=" << Amplitude/Counter << endl;
    }; // closes j loop
    AnnualCounter/=NumberOfStocks;
    AnnualAmplitude/=NumberOfStocks;
    cout << "HPTrueMu=" << HPTrueMu << ", HPLeash=" << HPLeash << ": ";
    cout << "AnnualCounter=" << AnnualCounter << ", AnnualAmplitude=" << AnnualAmplitude << "%";
    //PlotGSLMatrix(Spikes, "Spikes.csv", 1);
    vector<vector<double>> TrueValues = GSLMatrixToSTLMatrix(GSLTrueValues);
    // GENERATING BIASED VALUES
    double BiasedDistance=0; vector<double> VBiasedDistance;
    for (int i=0; i<NumberOfAgents; i++) {
        for (int j=0; j<NumberOfStocks; j++) {
            vector<double> C = CointegratedWalk(TrueValues[j], Market[i].Leash[j], Market[i].LeashVol[j], Market[i].Accuracy[j]); // Master, Leash, LeashVolatility, Accuracy
            for (int t=0; t<Time; t++) {
                gsl_matrix_set(Market[i].BiasedValues, j, t, C[t]);
                gsl_matrix_set(TrueBiased, i+1, t, C[t]);
                BiasedDistance+=100*((gsl_matrix_get (GSLTrueValues, j, t)-C[t])/gsl_matrix_get (GSLTrueValues, j, t))/(Time*NumberOfStocks*NumberOfAgents);
                VBiasedDistance.push_back(BiasedDistance);
            };
            C.clear();
        }; // closes j loop
    }; // closes i loop
    cout << ", BiasedDistance=" << BiasedDistance << "%" << endl;
    PlotGSLMatrix(TrueBiased, "TrueBiased.csv", 1);
    sort(VAnnualCounter.begin(), VAnnualCounter.end()); // sort() is by ascending order
    sort(VAnnualAmplitude.begin(), VAnnualAmplitude.end()); // sort() is by ascending order
    sort(VBiasedDistance.begin(), VBiasedDistance.end()); // sort() is by ascending order
    double VarAnnualCounter=0; for (int k=0; k<int(VAnnualCounter.size()); k++) {VarAnnualCounter+=(VAnnualCounter[k]-AnnualCounter)*(VAnnualCounter[k]-AnnualCounter)/int(VAnnualCounter.size());};
    double VarAnnualAmplitude=0; for (int k=0; k<int(VAnnualAmplitude.size()); k++) {VarAnnualAmplitude+=(VAnnualAmplitude[k]-AnnualAmplitude)*(VAnnualAmplitude[k]-AnnualAmplitude)/int(VAnnualAmplitude.size());};
    double VarBiasedDistance=0; for (int k=0; k<int(VBiasedDistance.size()); k++) {VarBiasedDistance+=(VBiasedDistance[k]-BiasedDistance)*(VBiasedDistance[k]-BiasedDistance)/int(VBiasedDistance.size());};
    Res.push_back(AnnualCounter); // Mean
    Res.push_back(AnnualAmplitude); // Mean
    Res.push_back(BiasedDistance); // Mean
    Res.push_back(VAnnualCounter[int(VAnnualCounter.size()/2)]); // Median
    Res.push_back(VAnnualAmplitude[int(VAnnualAmplitude.size()/2)]); // Median
    Res.push_back(VBiasedDistance[int(VBiasedDistance.size()/2)]); // Median
    Res.push_back(sqrt(VarAnnualCounter)); // Stdev
    Res.push_back(sqrt(VarAnnualAmplitude)); // Stdev
    Res.push_back(sqrt(VarBiasedDistance)); // Stdev
    return Res;
}; // closes CheckTrueGeneration()
void CheckTrueGenerationAll () {
    gsl_matrix* Res = gsl_matrix_calloc (9, 1); // For j=0 only
    //double Mu[3] = {0.10, 0.10, 0.10}; double Leash[3] = {1, 1, 1};
    double Mu[1] = {0.10}; double Leash[1] = {1};
    for (int k1=0; k1<1; k1++) {
        for (int k2=0; k2<1; k2++) {
            vector<double> V = CheckTrueGeneration (Mu[k1], Leash[k2]);
            gsl_matrix_set(Res, 0, k1*3+k2, V[0]);
            gsl_matrix_set(Res, 1, k1*3+k2, V[1]);
            gsl_matrix_set(Res, 2, k1*3+k2, V[2]);
            gsl_matrix_set(Res, 3, k1*3+k2, V[3]);
            gsl_matrix_set(Res, 4, k1*3+k2, V[4]);
            gsl_matrix_set(Res, 5, k1*3+k2, V[5]);
            gsl_matrix_set(Res, 6, k1*3+k2, V[6]);
            gsl_matrix_set(Res, 7, k1*3+k2, V[7]);
            gsl_matrix_set(Res, 8, k1*3+k2, V[8]);
        }; // closes k2 loop
    }; // closes k1 loop
    PlotGSLMatrix(Res, "Res.csv", 1);
};
*/


// MAIN PROGRAM - MAIN PROGRAM - MAIN PROGRAM - MAIN PROGRAM - MAIN PROGRAM
// SPECIFY "string Machine" JUSTE BELOW "using namespace std" ACCORDING TO LOCAL PATHWAY
// DOWNLOAD DATA ON THE GOOGLE DRIVE: https://drive.google.com/open?id=1J0ZHXnf_2wg80QyUyQwMN2YaUU9U45Xl
// DOWNLOAD GNU SCIENTIFIC LIBRARY GSL-v2.5 at https://www.gnu.org/software/gsl/
// SOME MEMORY ISSUES MAY REQUIRE REBOOT AND RESUME DURING COMPUTATIONS (work in progress)
int main (int argc, char** argv) {
    CliArgs args(argc, argv);
    Machine = args.output_dir.string() + "/";

    int S=20; // NUMBER OF SIMULATIONS

    // REAL DATA STATISTICS
    //vector<Share> PF=PortfolioGenerator ("ExchangesLSE", 20070102); // LOADS FULL LONDON STOCK EXCHANGE RAW DATA
    //PortfolioOutput (PF, 20); exit(0);
    //JointDistributions (PortfolioMultiOutput (PF, 0), 0, 100); // OUTPUTS FULL LONDON STOCK EXCHANGE CALIBRATION STATISTICS
    //vector<vector<Share>> PFTwin = PFTrainingTesting (PF, int(PF.size()/2)); // LOADS FULL LONDON STOCK EXCHANGE RAW DATA IN TRAINING AND TESTING SETS
    //JointDistributions (PortfolioMultiOutput (PFTwin[0], 0), 0, 100); // OUTPUTS TRAINING SET OF LONDON STOCK EXCHANGE CALIBRATION STATISTICS
    //JointDistributions (PortfolioMultiOutput (PFTwin[1], 0), 0, 100); // OUTPUTS TESTING SET OF LONDON STOCK EXCHANGE CALIBRATION STATISTICS
    //HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Classic", "NoCluster", 1, S); // Best calib
    //HPMarketSimulatorDummy (); exit(0);
    //for (int k=1; k<=8; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "LearningRate", "NoCluster", 1, 100, 0.5*k, S);}; // LearningRateScale={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}
    
    
    // HPMarketSimulatorAllOB (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, 50); // Metaorder impact at 50%
    // HPMarketSimulatorAllOB (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, 50); // Metaorder impact at 50%
    // HPMarketSimulatorAllOB (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, 50); // Metaorder impact at 50%
    // HPMarketSimulatorAllOB (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, 50); // Metaorder impact at 50%
    // HPMarketSimulatorAllOB (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, 50); // Metaorder impact at 50%
    // HPMarketSimulatorAllOB (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, 50); // Metaorder impact at 50%
    
    //HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S); // Classic
    
    //HPMarketSimulatorAllPD (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S); // Ordinary PD
    
    //for (int Trunk=0; Trunk<=6; Trunk++) {MarketSimulatorTrunk (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S, Trunk);};
    HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S);
    
    //for (int k=1; k<5; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "Best", k*500/5, 0, 1, S);};
    //for (int k=1; k<5; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "Worst", k*500/5, 0, 1, S);};
    //for (int k=1; k<5; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "Noise", k*500/5, 0, 1, S);};
    
    /*
    // New profile testing
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "LossAversion", "NoCluster", 1, k, 1, S);}; // LossAversion for p={20, 40, 60, 80, 100}
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Positivity", "NoCluster", 1, k, 1, S);}; // Positivity for p={20, 40, 60, 80, 100}
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Negativity", "NoCluster", 1, k, 1, S);}; // Negativity for p={20, 40, 60, 80, 100}
    double RateScale=2; for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "LearningRate", "NoCluster", 1, k, RateScale, S);}; // RateScale=2 for p={20, 40, 60, 80, 100}
    for (int k=1; k<=8; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "LearningRate", "NoCluster", 1, 100, 0.5*k, S);}; // LearningRateScale={0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0}, learning rate=U(5%,20%)
    
    // Profiles testing
    HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S);
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "LossAversion", "NoCluster", 1, k, 1, S);};
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Positivity", "NoCluster", 1, k, 1, S);};
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Negativity", "NoCluster", 1, k, 1, S);};
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "DelayDiscounting", "NoCluster", 1, k, 1, S);};
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Fear", "NoCluster", 1, k, 1, S);};
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Greed", "NoCluster", 1, k, 1, S);};
    for (int k=20; k<=100; k+=20) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Human", "NoCluster", 1, k, 1, S);};
    
    
    // Agent Learning
    for (int k=1; k<5; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "Best", k*500/5, 0, 1, S);};
    for (int k=1; k<5; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "Worst", k*500/5, 0, 1, S);};
    for (int k=1; k<5; k++) {HPMarketSimulatorAll (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "Noise", k*500/5, 0, 1, S);};
    HPMarketSimulatorAll2 (500, 0, 0.10, 10, 2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, "PDOn", 1, S); // First switch to "PDOn" in HPMarketSimulatorAll()
    
    // Pareto
    HPMarketSimulatorAll (500, 0, 0.10, 10, 1*2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S);
    HPMarketSimulatorAll (500, 0, 0.10, 10, 2*2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S);
    HPMarketSimulatorAll (500, 0, 0.10, 10, 3*2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S);
    HPMarketSimulatorAll (500, 0, 0.10, 10, 4*2875+1000, 50, "Algorithmic", "NoCluster", 1, 0, 1, S);
    */
    
    
    return 0;
}


/*
 NEW
 1- DONE! CFM paper: Find metric to study agent impact of large orders and meta orders (cf. https://arxiv.org/abs/1901.05332) + expand references for tick study => we want to see at a given time t where a large order is sent (say once per year for each simulation) the metaorder impact on prices i1(t)=sigma_prices[t->t+tau]/sigma_prices[t-tau->t] (and likewise i2(t) on volumes) as a function of the metaorder size j(t)=Q_order(t)/Q_total(t), sigma being the standard deviation (of prices or volumes) over an interval of size tau, and Q being the stock quantity. For each of the S=20 simulation runs, we record and output the last 5 such tuples (j, i1(tau=w), i1(tau=2w), i1(tau=m), i1(tau=3m), i1(tau=6m), i2(tau=w), i2(tau=2w), i2(tau=m), i2(tau=3m), i2(tau=6m)), disregarding the first 5 because of being not yet learned by the agents. We do 4 such S=20 simulation runs for values of j=5%, 10%, 25%, 50%. The metaorder is sent during each run once every year by a random (yet non-bankrupt) agent whose stock holding is suddenly increased to these values (5%, 10%, 25%, 50%) of the total share outstanding at time t and then revert back to its value at time t+1. => we won't likely get for example 50% of shares outstanding transacted at a given time step t, but we can simply set MetaorderImpact at 50% and let the sumilations run for a large S, in hope we reach large values of j
 2- DONE! Opinion paper : output basic plots () and submit to https://saiconference.com/IntelliSys2020/CallforPapers#Guidelines
 3- DONE! Calibration paper : submit
 4- DONE! Mesoscale 1 : Focus on Sections 3.1 and 3.2 and especially on PD studies (stress the results of Section IIa) + submit
 5- DONE! Mesoscale 2 : Write intro from psycho paper and put results of Section 4 + submit
 6- DONE! Psycho paper : discuss thursday 6 feb (Set dynamic learning scale? Put learning rate results of Mesoscale in Psycho paper?) => issue of psycho biases not rigorously modelled

 
 
OLDER
 Stefano
 - DONE! Paper #1 at XXX: cite previous ABM papers whose calibration results were below the performance of ours
 - Paper #2 at XXX: stress the results of Section IIa
 - (SBG) Split paper #2 into two: one with the results of Section IIa, and another with the other results
 - (SP) Set dynamic learning scale? => for later
 - (SP) Put learning rate results in Psycho paper? => ask team
 - DONE! (SP) Remove fig 8 of Mesoscale paper? => no not needed unless team asks
 - DONE! Among best and worst agents, check those who belong to reflexive population p.
 - DONE! Increase learning rate to 1? => already done but inconclusive results
 - DONE! Put in Mesoscale paper main figures of the calibration? => in supplementary material

 CFM
 - DONE! Check tick digits => recoded
 - Check impact of meta orders
 - Set agents with very large PF and study their market and OB impact
 - Frédéric Bucci et al. “Slow decay of impact in equity markets: insights from the Ancerno database” 2019 : https://arxiv.org/abs/1901.05332
 - Check AC with [P(t+T)-P(t)]^2
 
 
 OLD
 - Split article (V- calibration, VI.A- DelayDiscounting/Fear/Greed/LossAversion/Negativity/Positivity/Human, VI.B- agent performance, VI.C- agent NAVs, VI.D- Agent PD & noise, VI. Agent reflexivity best & worst) in four papers : i- model calibration (V- calibration), ii- agent performance/trading (VI.B- agent performance, VI.C- agent NAVs, VI.D- Agent PD), iii- market stability (VI.D- Agent noise, VI. Agent reflexivity best & worst) , iv- psychological biases (VI.A- DelayDiscounting/Fear/Greed/LossAversion/Negativity/Positivity/Human)
 - Apply HRL with PF management metrics
 - CFM after arXiv on July 15th
 - Positivitiy and negativity could be merged, under a positivitity chapeau. It could be instantiated as a difference in learning rate rather than a bonus malus in the reward (which is not memory afterall).
 - Work with Stefano on the psychological paper in mid august.
 - DONE! Reset LossAversion, Positivity and Negativity : loss aversion has nothing conclusive (even contradicts previous runs!), negativity shows some spread variation, positivity shows means of AC-2w of volumes
 
 REMARKS
 - Think about studying online survival probability wrt. performance
 - Do multithreading with boost
 - NOTE: Model consistency: i- no long-short strategies, ii- reliance on True, iii- no derivatives, iv- number of agents & percentages of profiles populations, v- realism of biases via agent RL framework, vi- intraday and fixing market activity, vii- lack of asynchronicity between agents learning (asynchronicity via biases), viii- issue of irrealistic PM metrics (daily/weekly/monthly/yearly returns, monthly/yearly drawdowns/bankruptcy, alpha, beta, sharpes, PF turnover, liquidity/volumes, etc.), ix- no leverage, x- no metaorders, xi- no legal nor regulation constraints, xii- execution costs (transaction fees, slippage, taxes, BrokerFee, CostOfCarry, SpreadFee, DividendFee, frais d'emprunts), xiii- dividend yields & inflation, xiv- diversity of orders (limit, market, market-on-close, VWAP), xv- no PF diversification (equity, cross-asset), xvi- parameters of agents (Future, RFAFirst, RBAFirst, History, etc.).
 - NOTE: See experiment prospect theory, description gap, endowement effect + works of Drazen Prelec. Pour installer python, recommandation de Conda avec tous les packages numeriques (NumPy et Panda) mais attention aux conflits avec les installations antérieures.
 - NOTE: On the long term: design limit/market orders & short-selling
 - NOTE: It is true we need to include a short presentation of the different sections in the beginning of the paper. Also we must include in the beginning a brief Introduction of the topic, one or two lines. Finally we should remove the "we", "our", etc.
 F(s)27={Svol3,Lvol3,Reflexive3}, F(a)27={Tool3,Lag3,Weight3}
 T(s)108={Mu3,Lvol3,RBALevel2,RFALevel2,Liquid3}, T(a)9=(Quant3,Pinch3)
 - PROFILE1 LOSS AVERSION NEBLossAversion=0.2r by +Pinch: if ((NEBLossAversion>VRan[13]) && (Pinch!=2)) {Pinch+=1}
 - PROFILE2 POSITIVITY NEBPositivityp=0.2r by +dFResultDisReal,+dTResultDisReal: if ((NEBPositivityp>VRan[17]) && (dFResultDisReal>0)) {dFResultDisReal+=1}, if ((NEBPositivityp>VRan[17]) && (dTResultDisReal>0)) {dTResultDisReal+=1}
 - PROFILE3 DelayDiscounting NEB10=0.2r by Future/2
 - LEARNING RATE NEBLearningRate=0.4r+0.1
 
 OTHERS
 - Voir Scott Meyers "Effective C++" et C++14 aussi avec Scott Meyers "Effective modern C++"
 - http://www.cplusplus.com/doc/tutorial/typecasting/
 - http://www.cplusplus.com/doc/oldtutorial/typecasting/
 - https://www.threadingbuildingblocks.org/
 - thread building blocks parallel_for
 - boost parallel for
 - https://stackoverflow.com/questions/15206446/parallel-tasks-get-better-performances-with-boostthread-than-with-ppl-or-openm
 - Au lieu de static_cast<int>(x) utiliser par exemple int(x) si x est non int
 - Dark Reader pour Add-on Firefox

 OLD
 - DONE! Check and simplify True and Biased => we tried different parameters and checked with CheckTrueGenerationAll() the annual number of jumps in True (AnnualCounter), the mean of those jumps in percentage (AnnualAmplitude), and the mean of difference between biased and true for all agents (BiasedDistance). We found that varying Accuracy or Leash had little effect : however increasing HPTrueMu greatly impacted AnnualCounter : HPTrueMu={0.10, 0.50, 0.75} implied AnnualCounter={13, 6, 4}. Hence we find that for our calibration (HPTrueMu=0.10, Leash=1, Accuracy=10, etc.) and I=500, J=10, T=10y, we have the following means : AnnualCounter=12.70 ± 1.85, AnnualAmplitude=5.90 ± 1.84%, and BiasedDistance=2.37 ± 1.36% (where the ± term correspond to standard deviations).
 - DONE! Redesign new profiles
 - DONE! Check statistics and PD lines of code
 - DONE! S=1200 and send to Ivan Lazarevitch
 - DONE! Update static_cast<int> with int() (the latter is a bit faster)
 - DONE! Check policy distance as a function of time of computation at t=1Y, 2Y, 3Y, etc.., and do it for several simulations and plot the distribution of each PD for best-to-best, best-to-all, etc... for t=1y, 2y, 3y, etc...
 - DONE! Study arbitrage with various distributions at all t: i- AC-logreturns of [t-∆, t] & [t-2∆, t-∆] for ∆={3, w, 2w, 3w, m}, ii- AC-logreturns of [t-∆, t] & [t-∆-∂,t-∂] for ∂={1, 3, w, 2w, 3w, m}, iii- spread (think of right definition),
 - DONE! Output in JointDistributions() medians, variance, Q1, Q2, Q3, Q4, skewness, kurtosis of all mesoscale metrics (returns, bankruptices, NAV, etc.)
 - DONE! Make a S=100 and send to Ivan Lazarevitch
 - DONE! JJJ9 changes to get more frequent biases needs to be reverted to 0.2 when Classic (or else recalibrate), but with JJJ9 it gets interesting for => i- Loss aversion (maybe increasing Pt-Tt), ii- positivity (maybe decreasing volumes), iii- (maybe stability of #crashes)
 - DONE! Recalibrate with Algorithmic only => same parameters as before
 
*/
