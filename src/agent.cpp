#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "agent.hpp"
#include "global.hpp"
#include "utility.hpp"
#include "cli.hpp"

using namespace std;


// Returns the index of a vector<double> Tensor (decomposed from a tensor of dimensions vector<int> Dimensions) corresponding to the vector<int> Indices. This can be used to update the corresponding element in Tensor, or to simply access it
int Agent::get_tensor_index (const int IndSize, int Indices[], int Dimensions[]) {
    int Res=0;
    int ResTemp;
    for (int i=0; i<IndSize; i++) {
        ResTemp=Indices[i];
        for (int j=i+1; j<IndSize; j++) {
            ResTemp*=Dimensions[j];
        };
        Res+=ResTemp;
    };
    return Res;
}; // closes method get_tensor_index()


// Returns the D-dim tensor coordinates corresponding to the Index of the STL vector representation of that tensor
vector<int> Agent::get_tensor_coord (int Index, int DimSize, int Dimensions[]) {
    //ofstream outputLog("/Users/admin/Documents/GNT/SYMBA/SimLog.txt", ofstream::app);
    vector<int> Res;
    int Multi1, Multi2, Sum;
    for (int k=0; k<DimSize; k++) {
        Res.push_back(0);
        Multi1=1; Sum=0;
        for (int n=k+1; n<DimSize; n++) {Multi1*=Dimensions[n];};
        for (int i=0; i<k; i++) {
            if (k==0) {break;};
            Multi2=1;
            for (int j=i+1; j<DimSize; j++) {Multi2*=Dimensions[j];};
            Sum+=Res[i]*Multi2;
        };
        //outputLog << "Index=" << Index << ", Sum=" << Sum << ", Multi1=" << Multi1 << endl;
        Res[k]=(Index - Sum)/Multi1;
    };
    return Res;
}; // closes method get_coord()


double Agent::Capital() { // Computation of total capital amount of portfolio
    // The NAV or Net Asset Value is given by the sum of the RFA and the collection of stocks valued at market price
    double Temp=0;
    for (int i=0; i<int(Stocks.size()); i++) {Temp += ((Stocks[i]).StockAskValue)*((Stocks[i]).StockQuantity);};
    return RFA+Temp;
};


double Agent::StockHoldings() { // Computation of stock holdings
    double Temp=0;
    for (int i=0; i<int(Stocks.size()); i++) {Temp += ((Stocks[i]).StockAskValue)*((Stocks[i]).StockQuantity);};
    return Temp;
};


void Agent::Liquidation(int i, vector<Agent> &Market) {
    Market[i].RFA=0; Market[i].RBA=0; // Clearing all his bond holdings
    for (int j=0; j<int(Market[i].Stocks.size()); j++) {Market[i].Stocks[j].StockQuantity = 0;}; // Clearing all the agent stocks holdings
    Market[i].Bankruptcy=0; // The agent gets entirely liquidated
};


// RL algorithm for forecasting and placing order in the OB (SS6)
void Agent::RL(
    int j, int t, double Rate, gsl_matrix* ReflexiveValues, double VSpread, double LiquidPercent, int Time,
    int NumberOfStocks, string TradingFrequencyCond, string Plot, string VersatilityCondition, double MarketPerformance,
    int TimeSinceJan1st, int LearningPhase, int LeaderAgent, int LeaderQuant, int LeaderPinch, int ClusterLimit,
    int Trunk, const filesystem::path& output_dir) {
    // F(s): Rough={0,1,2,3}, Smooth={0,1,2,3}, Reflexive={0,1,2}
    // F(a): Tool={0,1,2}, Lag={0,1,2}, Weight={0,1,2}
    // T(s): Mu={0,1,2,3,4,5,6}, Sig={0,1,2}, RFA={0,1,2}, RBA={0,1,2}, Liquid={0,1,2,3}
    // T(a): Quant={0,1,2,3,4,5,6}, Pinch={0,1,2}
    ofstream outputLog(output_dir / "SimLog.txt", ofstream::app);
    // ofstream outputDebug("/Users/admin/Documents/GNT/SYMBA/Debug.txt", ofstream::app);
    double Seed=1000*t*(RFA+RBA+AgentName);
    //double Seed=1000*t*gsl_matrix_get(ReflexiveValues, j, t);
    //ERROR RNG
    //#define GSL_DLL
    //#define DGSL_DLL
    //#define DWIN32
    //gsl_rng* r = gsl_rng_alloc (gsl_rng_mt19937); // or gsl_rng_default instead of gsl_rng_mt19937
    gsl_rng* r = make_rng();
    gsl_rng_set(r, static_cast<unsigned long int>(time(0)+Seed)); //gsl_rng_set(const gsl_rng* r, unsigned long int s)
    vector<double> VRan; for (int i=0; i<36; i++) {VRan.push_back(gsl_rng_uniform (r));}; // Generating a vector of random numbers without seeding issues
    //#undef GSL_DLL
    //#undef DGSL_DLL
    //#undef DWIN32
    //delete r; r=0; // Delete and de-allocate pointer
    if (Plot=="On") {
        outputLog << "***************RL() FOR AGENT " << AgentName << " STOCK " << j << " AT TIME t=" << t << "**************" << endl;
        outputLog << "ReflexiveValues[]={"; for (int h=0; h<=t; h++) {outputLog << gsl_matrix_get(ReflexiveValues, j, h); if (h<t) {outputLog << ", ";};};
        outputLog << "}" << endl;
        outputLog << "Seed=" << Seed << endl;
        outputLog << "VRan[]={" ; for (int u=0; u<int(VRan.size()); u++) {if (u<int(VRan.size())-1) {outputLog << VRan[u] << ", ";} else {outputLog << VRan[u];};}; outputLog << "}" << endl;
        outputLog << "At step ABM I.6, the bid and ask values are first defined (as resp. the time-discounted minimum and maximum values between present price and its forecast). Then at step ABM II.3, the gesture on bid-ask spread (''Pinch'') is incorporated in these values. Then at step ABM II.3, the actual quantity to short (QTradeReal<0) or long (QTradeReal>0) is also updated in class agent member (''Quant''). Notice that these quantities are already formated to the real numbers to trade, and as such are ready for transaction (except the sign)." << endl;
        
        
        /*************************************ABM I.1******************************************/
        // Selecting the trading window (Future) according to present time
        outputLog << endl << "    ABM I.1: Selecting the trading window (Future) according to present time" << endl;
        outputLog << "Future=" << Future << ", History=" << History << endl;
        
        
        /*************************************ABM I.2******************************************/
        outputLog << endl << "    ABM I.2: Defining present state s" << endl;
    };// closes Plot condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.1 computed..." << endl;
    // Get the state s
    double Mean=0; double Variance=0; int Start=t-(3*Future); if (Start<0) {Start=0;}; int End=t; // BOLLINGER ISSUE
    for (int p=Start; p<=End; p++) {Mean+=gsl_matrix_get(ReflexiveValues, j, p);}; Mean=Mean/(End-Start+1); // Mean
    for (int p=Start; p<=End; p++) {Variance+=abs(gsl_matrix_get(ReflexiveValues, j, p) - Mean);};
    double Lvol=Variance/(End-Start+1); // Variance Lvol
    if (Lvol<0.001) {Lvol=0.001;}; // Failsafe
    if (Plot=="On") {outputLog << "For t=[" << Start << "," << End << "], Mean=" << Mean << ", Variance=" << Variance << "  " << Lvol << endl;
        LvolLog[j].push_back(Lvol);}; // LLL
    Mean=0; Variance=0; Start=t-(Future)-1; if (Start<0) {Start=0;}; End=t; // BOLLINGER ISSUE
    for (int p=Start; p<=End; p++) {Mean+=gsl_matrix_get(ReflexiveValues, j, p);}; Mean=Mean/(End-Start+1); // Mean
    for (int p=Start; p<=End; p++) {Variance+=abs(gsl_matrix_get(ReflexiveValues, j, p) - Mean);};
    double Svol=Variance/(End-Start+1); // Variance Svol
    if (Svol<0.001) {Svol=0.001;}; // Failsafe
    if (Plot=="On") {outputLog << "For t=[" << Start << "," << End << "], Mean=" << Mean << ", Variance=" << Variance << " => Svol=" << Svol << endl;
        SvolLog[j].push_back(Svol);}; // LLL

    double dRough=Svol/Lvol; double dRoughPercentile=0; int RoughVecSize=int(dRoughVec[j].size());
    if (Lvol==0) {dRough=0;};
    if (t==1) {dRoughVec[j].push_back(dRough); dRoughPercentile=0.5;}
    else { // sorting the vector dRoughVec by ascending values
        for (int i=0; i<RoughVecSize; i++) {
            if (dRough<dRoughVec[j][i]) {dRoughVec[j].insert(dRoughVec[j].begin()+i, dRough*(1 + 0.0001*VRan[29])); dRoughPercentile=i*1.0/RoughVecSize; break;};
            if (i==RoughVecSize-1) {dRoughVec[j].push_back(dRough*(1 + 0.0001*VRan[29])); dRoughPercentile=1; break;};
        }; // closes for loop
    }; // closes else
    
    /*
        if (dRoughPercentile<0.25) {Rough=0;}
        else if ((dRoughPercentile>=0.25) && (dRoughPercentile<0.5)) {Rough=1;}
        else if ((dRoughPercentile>=0.5) && (dRoughPercentile<0.75)) {Rough=2;}
        else if (dRoughPercentile>=0.75) {Rough=3;};
        */
    
    if (dRoughPercentile<0.25) {Rough=0;}
    else if ((dRoughPercentile>=0.25) && (dRoughPercentile<0.75)) {Rough=1;}
    else if (dRoughPercentile>=0.75) {Rough=2;};
    
    if ((Lvol==0) || (t<=3)) {Rough=1;};
    if (int(dRoughVec[j].size()) - Future - History >0) {dRoughVec[j].erase (dRoughVec[j].begin());}; // OOO
    if (Plot=="On") {RoughLog[j].push_back(Rough);}; // LLL
    
    double dSmooth=Lvol/(3*Future); double dSmoothPercentile=0; int SmoothVecSize=int(dSmoothVec[j].size());
    if (3*Future==0) {dSmooth=0;};
    if (t==1) {dSmoothVec[j].push_back(dSmooth); dSmoothPercentile=0.5;}
    else { // sorting the vector SmoothVec by ascending values
        for (int i=0; i<SmoothVecSize; i++) {
            if (dSmooth<dSmoothVec[j][i]) {dSmoothVec[j].insert(dSmoothVec[j].begin()+i, dSmooth*(1 + 0.0001*VRan[29])); dSmoothPercentile=i*1.0/SmoothVecSize; break;};
            if (i==SmoothVecSize-1) {dSmoothVec[j].push_back(dSmooth*(1 + 0.0001*VRan[29])); dSmoothPercentile=1; break;};
        }; // closes for loop
    }; // closes else
    /*
        if (dSmoothPercentile<0.25) {Smooth=0;}
        else if ((dSmoothPercentile>=0.25) && (dSmoothPercentile<0.5)) {Smooth=1;}
        else if ((dSmoothPercentile>=0.5) && (dSmoothPercentile<0.75)) {Smooth=2;}
        else if (dSmoothPercentile>=0.75) {Smooth=3;};
        */
    
    if (dSmoothPercentile<0.25) {Smooth=0;}
    else if ((dSmoothPercentile>=0.25) && (dSmoothPercentile<0.75)) {Smooth=1;}
    else if (dSmoothPercentile>=0.75) {Smooth=2;};
    
    if ((t<=3) || (3*Future==0)) {Smooth=1;};
    if (int(dSmoothVec[j].size()) - Future - History >0) {dSmoothVec[j].erase (dSmoothVec[j].begin());}; // OOO
    if (Plot=="On") {SmoothLog[j].push_back(Smooth);}; // LLL
    
    // DDD
    double dReflexive=0; Start=t-3*Future; if (Start<0) {Start=0;};
    for (int p=Start; p<=t; p++) {dReflexive+=100*abs(gsl_matrix_get(ReflexiveValues, j, p)-gsl_matrix_get(BiasedValues, j, p))/(gsl_matrix_get(ReflexiveValues, j, p)*(t-Start+1));};
    Reflexive=0; // Reflexive={0,1,2}={[0%,10%[, [10%,30%[, [30%,+\infty[} for dReflexive=<100*abs(Pp-Biased)/Pp> over p[t-3*Future,t]
    if ((dReflexive>=0) && (dReflexive<10)) {Reflexive=0;}
    else if ((dReflexive>=10) && (dReflexive<30)) {Reflexive=1;}
    else if (dReflexive>=30) {Reflexive=2;};
    if (Plot=="On") {ReflexiveLog[j].push_back(Reflexive);}; // LLL
    
    if (Plot=="On") {
        outputLog << "dRough=Svol/Lvol=" << dRough << " and Rough=" << Rough << endl;
        outputLog << "dSmooth=Lvol/(3*Future)=" << dSmooth << " and Smooth=" << Smooth << endl;
        outputLog << "dReflexive=mean <100*(P(t)-B(t))/P(t)> over t=[t-3*Future, t]=" << dReflexive << " and Reflexive=" << Reflexive << endl; // DDD
        outputLog << "P(t)=" << gsl_matrix_get(ReflexiveValues, j, t) << ", B(t)=" << gsl_matrix_get(BiasedValues, j, t) << endl;
    }; // closes Plot condition
    FSInd[0]=Rough; FSInd[1]=Smooth; FSInd[2]=Reflexive; // DDD
    int STensorIndex = get_tensor_index (3, FSInd, FSDim); // Vector index in the tensor of all possible RL states // DDD
    FSIndex[j].push_back(STensorIndex); // Recording state s_t (ABM#1 for Forecast)
    if (int(FSIndex[j].size()) - Future - 1 >0) {FSIndex[j].erase (FSIndex[j].begin());}; // Updates by checking FSIndex[int(XXX.size()) - Future - 1] //OOO
    if (Plot=="On") {
        outputLog << "STensorIndex=" << STensorIndex << endl;
        outputLog << "FSIndex[]={" ; for (int u=0; u<int(FSIndex[j].size()); u++) {if (u<int(FSIndex[j].size())-1) {outputLog << FSIndex[j][u] << ", ";} else {outputLog << FSIndex[j][u];};}; outputLog << "}" << endl;
        outputLog << "FSDim[]={" ; for (int u=0; u<3; u++) {if (u<3-1) {outputLog << FSDim[u] << ", ";} else {outputLog << FSDim[u];};}; outputLog << "}" << endl;
        outputLog << "FSInd[]={" ; for (int u=0; u<3; u++) {if (u<3-1) {outputLog << FSInd[u] << ", ";} else {outputLog << FSInd[u];};}; outputLog << "}" << endl;
        
        
        /*************************************ABM I.3******************************************/
        outputLog << endl << "    ABM I.3: Initializing the start policy pi_0" << endl;
    }; // closes Plot condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.2 computed..." << endl;
    // Initializing the start policy $\pi_0$ out of which an an action will be selected (we initialize such that in all states s we take Tool=0, Lag=0)
    if (t==1) {
        FPiInd[0]=0; FPiInd[1]=0; FPiInd[2]=0; FPiInd[3]=0; FPiInd[4]=0; FPiInd[5]=0; // DDD
        for (int i=0; i<FS*FA; i++) {FPi[j].push_back(1.0/FA);};
        if (Plot=="On") {outputLog << "Policy initialized... all " << FA << " possible actions are given equiprobability 1/" << FA << "=" << 1.0/FA << endl;};
    }; // HHH All actions are equiprobable
    if (t>1) {if (Plot=="On") {outputLog << "Policy already initialized..." << endl;};};
    
    
    /*************************************ABM I.4******************************************/
    if (Plot=="On") {outputLog << endl << "    ABM I.4: Using current policy pi to select real action a=(Tool,Lag,Reflexive) in present state s" << endl;};
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.3 computed..." << endl;
    // We now use the current policy to perform an action $a$ in our state $s$. Note: There is no need to use \verb?PiStack?, all we need is when we have a vector of action each with a probability in $[0,1]$, we generate a random number from a uniform distribution: if that number is larger than the first probability it is subtracted from it, and if that new number is still larger than the second probability, we do likewise until we arrive to the probability corresponding to action $a$ to be picked up. This is because the sum of all these probabilities is $1$.
    FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=Reflexive; FTenPiId[3]=0; FTenPiId[4]=0; FTenPiId[5]=0; // DDD
    //FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=0; FTenPiId[3]=0; // DDD2
    int Fid = 0;
    Fid = get_tensor_index (6, FTenPiId, FPiDim); // Vector index of pi-tensor representation corresponding to present s and a=(0,0) // DDD
    //Fid = get_tensor_index (4, FTenPiId, FPiDim); // Vector index of pi-tensor representation corresponding to present s and a=(0,0) // DDD2
    if (Plot=="On") {
        outputLog << "FTenPiId[]={"; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FTenPiId[u] << ", ";} else {outputLog << FTenPiId[u];};}; outputLog << "}" << endl; // DDD
        //outputLog << "FTenPiId[]={"; for (int u=0; u<4; u++) {if (u<4-1) {outputLog << FTenPiId[u] << ", ";} else {outputLog << FTenPiId[u];};}; outputLog << "}" << endl; // DDD2
        outputLog << "Fid=" << Fid << ", FPi[j].size()=" << FPi[j].size() << endl;
    }; // closes Plot condition
    
    double Fix=VRan[0];
    double Fix2=VRan[1];
    double Act=0; int Pos=Fid;
    if (Plot=="On") {outputLog << "Fix=" << Fix << ", Fix2=" << Fix2 << endl;};
    
    if (Exploration==1) { // Epsilon-greedy: Choosing always the action with highest probability, but once in a while (with a small probability epsilon) choosing a (uniform) random one (which is a problem: the worst potential action can be taken as likely as the best, with in some cases devastating returns).
        if (Plot=="On") {outputLog << "Epsilon-greedy method selected" << endl;}; // DDD1
        if (Fix<Epsilon) { // Exploration
            Pos=Fid + int(floor(FA*Fix2));
            Act=FPi[j][Pos];
            if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }
        else { // Exploitation
            for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", FPi[n]=" << floor(100*FPi[j][n]) << "%" << endl;};
                if (FPi[j][n]>=Act) {Pos=n; Act=FPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
            }; // closes for
            vector<int> PosVec; // Indices of equiprobable actions
            for (int n=Fid; n<Fid+FA; n++) {if (FPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
            Pos=PosVec[int(floor(VRan[23]*PosVec.size()))]; // We choose a random indice among these!
            PosVec.erase(PosVec.begin(),PosVec.end());
            if (Plot=="On") {outputLog << "Exploitation epsilon-greedy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else
    } // closes if
    
    else if (Exploration==2) { // Softmax: Same as above except the random choice in now not (uniform) equiprobable for all actions, but graded according to their estimated value. One method (using a Gibbs/Boltzmann distribution $e^{Q_t(a) / \tau} / \sum_{b=1}^n e^{Q_t(b) / \tau}$, with $\tau$ the temperature) allows to even shift this grading all the way back to epsilon greedy methods (when $\tau \rightarrow 0$).
        if (Plot=="On") {outputLog << "Softmax method selected" << endl;};
        if (Fix<Epsilon) { // Exploration
            Pos=Fid + int(floor(FA*Fix2));
            Act=FPi[j][Pos];
            if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }
        else { // Exploitation
            vector<int> N = Shuffle(FA);
            if (Plot=="On") {outputLog << "From n=" << Fid << " to " << Fid+FA-1 << ": Fix2=" << Fix2 << endl;}; // DDD1
            for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                //if (Plot=="On") {outputLog << "For n=" << n << ", FPi[n]=" << floor(100*FPi[n]) << "%" << endl;};
                int Newn = N[n-Fid]+Fid;
                Fix2-=FPi[j][Newn];
                if (Plot=="On") {outputLog << "For n=" << Newn << ", FPi[j][Newn]=" << floor(100*FPi[j][Newn]) << "%: Fix2-=FPi[j][Newn]=" << Fix2 << endl;};
                if (Fix2<=0) {Pos=Newn; Act=FPi[j][Pos]; break;}; // We try and find a new action according to its graded probability
            }; // closes for
            if (Plot=="On") {outputLog << "Exploitation softmax: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else
    } // closes else if
    
    else if (Exploration==3) { // Pursuit: The probability of selecting the greedy action $a_{t+1}=a_{t+1}^{\ast}$ at next time step is incremented a fraction $\beta$ of the way toward $1$, while the probabilities of selecting the other actions are decremented toward $0$, respectively shown by the following two expressions: $\pi_{t+1}(a_{t+1}^{\ast}) = \pi_{t}(a_{t+1}^{\ast}) + \beta \left[ 1- \pi_{t}(a_{t+1}^{\ast}) \right]$ and $\pi_{t+1}(a) = \pi_{t}(a) + \beta \left[ 0- \pi_{t}(a) \right] , \hspace{5mm} \forall a\neq a_{t+1}^{\ast}$.
        if (Plot=="On") {outputLog << "Pursuit method selected" << endl;}; // DDD1
        if (Fix<Epsilon*(1-(double(t)/Time))) { // Exploration
            Pos=Fid + int(floor(FA*Fix2));
            Act=FPi[j][Pos];
            if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }
        else { // Exploitation
            for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", FPi[j][n]=" << floor(100*FPi[j][n]) << "%" << endl;};
                if (FPi[j][n]>=Act) {Pos=n; Act=FPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
            }; // closes for
            vector<int> PosVec; // Indices of equiprobable actions
            for (int n=Fid; n<Fid+FA; n++) {if (FPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
            Pos=PosVec[int(floor(VRan[23]*PosVec.size()))]; // We choose a random indice among these!
            PosVec.erase(PosVec.begin(),PosVec.end());
            if (Plot=="On") {outputLog << "Exploitation pursuit: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else
    } // closes else if
    
    else if (Exploration==4) { // This is for off-policy Watkin's Q-learning
        if (Plot=="On") {outputLog << "Off-policy Watkin's Q-learning method selected" << endl;}; // DDD1
        for (int n=Fid; n<Fid+FA; n++) { // Running over the actions a associated with that state s
            if (Plot=="On") {outputLog << "For n=" << n << ", FPi[j][n]=" << floor(100*FPi[j][n]) << "%" << endl;};
            if (FPi[j][n]>=Act) {Pos=n; Act=FPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
        }; // closes for
        vector<int> PosVec; // Indices of equiprobable actions
        for (int n=Fid; n<Fid+FA; n++) {if (FPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
        Pos=PosVec[int(floor(VRan[23]*PosVec.size()))]; // We choose a random indice among these!
        PosVec.erase(PosVec.begin(),PosVec.end());
        if (Plot=="On") {outputLog << "Exploitation Off-policy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
    }; // closes else if
    
    int FormerWeight=Weight;
    vector<int> TenPiCoo = get_tensor_coord (Pos, 6, FPiDim); // Getting the tensor coordinates from that index
    //vector<int> TenPiCoo = get_tensor_coord (Pos, 4, FPiDim); // Getting the tensor coordinates from that index // DDD2
    Tool=TenPiCoo[3]; // Selection of first action from policy
    Lag=TenPiCoo[4]; // Selection of second action from policy
    Weight=TenPiCoo[5]; // Selection of third action from policy
    
    //if (Tool==0) {Tool+=1;}; // we change Poly to Linear UUU
    
    // Tuning the frequency of changing Weight according to Versatility
    if ((VRan[35] > 1.0/Versatility) && (t>1)) { // Case where we do not want to switch but maintain the previous choices for Weight
        Weight=FormerWeight;
        FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=Reflexive; FTenPiId[3]=Tool; FTenPiId[4]=Lag; FTenPiId[5]=Weight;
        Pos = get_tensor_index (6, FTenPiId, FPiDim); // Index representation of present s and a with former reflexivity weight
        TenPiCoo = get_tensor_coord (Pos, 6, FPiDim); // Getting the tensor coordinates from that index
        if (Plot=="On") {outputLog << "Weight overridden to FormerWeight=" << FormerWeight << " since VRan[35]=" << VRan[35] << ">" << "1.0/Versatility=" << 1.0/Versatility << endl;};
    }; // closes if
    
    // Turning Versatility condition on or off
    if (VersatilityCondition=="Off") {
        Weight=1; // Always fixed at Reflexivity
        FTenPiId[0]=Rough; FTenPiId[1]=Smooth; FTenPiId[2]=Reflexive; FTenPiId[3]=Tool; FTenPiId[4]=Lag; FTenPiId[5]=Weight;
        Pos = get_tensor_index (6, FTenPiId, FPiDim); // Index representation of present s and a with reflexivity weight (Weight=1)
        TenPiCoo = get_tensor_coord (Pos, 6, FPiDim); // Getting the tensor coordinates from that index
        if (Plot=="On") {outputLog << "Weight overridden to 1 (i.e to the natural Reflexivity of the agent) since VersatilityCondition is turned off" << endl;};
    }; // closes if
    
    if (Plot=="On") {outputLog << "TenPiCoo[]={" ; for (int u=0; u<int(TenPiCoo.size()); u++) {if (u<int(TenPiCoo.size())-1) {outputLog << TenPiCoo[u] << ", ";} else {outputLog << TenPiCoo[u];};}; outputLog << "}, ";
        outputLog << "Tool=TenPiCoo[3]=" << Tool << ", ";
        outputLog << "Lag=TenPiCoo[4]=" << Lag << ", ";
        outputLog << "Weight=TenPiCoo[5]=" << Weight << endl;
    }; // closes Plot condition
    // Defining the real action as tensor index to take according to above policy which selected Tool and Lag
    //Lag+=1; // To make sure the value 0 means 1, 1 means 2, and 2 means 3, resp.
    FAIndReal[0]=Tool; FAIndReal[1]=Lag; FAIndReal[2]=Weight;
    int ARealTensorIndex = get_tensor_index (3, FAIndReal, FADim); // Vector index in the tensor of all possible RL actions
    FAIndexReal[j].push_back(ARealTensorIndex); // Recording real action a_t (ABM#2 for Forecast)
    if (int(FAIndexReal[j].size()) - Future - 1 >0) {FAIndexReal[j].erase (FAIndexReal[j].begin());}; // OOO
    if (Plot=="On") {
        outputLog << "FAIndReal[0]=Tool=" << Tool << ", FAIndReal[1]=Lag=" << Lag << ", FAIndReal[2]=Weight=" << Weight << ", ";
        outputLog << "ARealTensorIndex=" << ARealTensorIndex << ", ";
        outputLog << "FAIndexReal[j][]={" ; for (int u=0; u<int(FAIndexReal[j].size()); u++) {if (u<int(FAIndexReal[j].size())-1) {outputLog << FAIndexReal[j][u] << ", ";} else {outputLog << FAIndexReal[j][u];};}; outputLog << "}" << endl;
        if (Plot=="On") {ToolRealLog[j].push_back(Tool); LagRealLog[j].push_back(Lag); WeightRealLog[j].push_back(Weight);}; // LLL
    }; // closes Plot condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.4 computed..." << endl;
    
    /*************************************ABM I.5******************************************/
    int VirtualTool=0; int VirtualLag=0; int VirtualWeight=0;
    if (Exploration==4) {
        if (Plot=="On") {    outputLog << endl << "    ABM I.5: Randomly selecting a virtual action a=(VirtualTool,VirtualLag)" << endl;};
        // Defining the virtual action as tensor index to take (not based on any exploration policy)
        VirtualTool=int(1000*VRan[2])%3; // 0, 1, or 2
        VirtualLag=int(1000*VRan[3])%3; // 1, 2, or 3
        VirtualWeight=int(1000*VRan[33])%3; // 1, 2, or 3
        for (int i=0; i<100; i++) {int VT=int(1000*VRan[4])%3; if (VT!=Tool) {VirtualTool=VT; break;};}; // Making sure virtual action not like real one
        for (int i=0; i<100; i++) {int VL=int(1000*VRan[5])%3; if (VL!=Lag) {VirtualLag=VL; break;};}; // Making sure virtual action not like real one
        for (int i=0; i<100; i++) {int VL=int(1000*VRan[34])%3; if (VL!=Weight) {VirtualWeight=VL; break;};}; // Making sure virtual action not like real one
        FAIndVirtual[0]=VirtualTool; FAIndVirtual[1]=VirtualLag; FAIndVirtual[2]=VirtualWeight; // Turlututu
        int AVirtualTensorIndex = get_tensor_index (3, FAIndVirtual, FADim); // Vector index in the tensor of all possible RL actions
        FAIndexVirtual[j].push_back(AVirtualTensorIndex); // Recording real action a_t (ABM#2 for Forecast)
        if (int(FAIndexVirtual[j].size()) - Future - 1 >0) {FAIndexVirtual[j].erase (FAIndexVirtual[j].begin());}; // OOO
        if (Plot=="On") {
            outputLog << "FAIndVirtual[0]=VirtualTool=" << VirtualTool << ", FAIndVirtual[1]=VirtualLag=" << VirtualLag << ", FAIndVirtual[2]=VirtualWeight=" << VirtualWeight << ", ";
            outputLog << "AVirtualTensorIndex=" << AVirtualTensorIndex << endl;
            outputLog << "FAIndexVirtual[j][]={" ; for (int u=0; u<int(FAIndexVirtual[j].size()); u++) {if (u<int(FAIndexVirtual[j].size())-1) {outputLog << FAIndexVirtual[j][u] << ", ";} else {outputLog << FAIndexVirtual[j][u];};}; outputLog << "}" << endl;
            if (Plot=="On") {ToolVirtualLog[j].push_back(VirtualTool); LagVirtualLog[j].push_back(VirtualLag); WeightVirtualLog[j].push_back(VirtualWeight);}; // LLL
        }; // closes Plot condition
    }; // closes Exploration==4 condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.5 computed..." << endl;
    
    
    /*************************************ABM I.6******************************************/
    if (Plot=="On") {outputLog << endl << "    ABM I.6: Forecasting the time series Future time steps ahead" << endl;};
    // F() must work on Reflexive not Reference, and the result of this is backed up and used for evaluation of the result of the whole F(). Then once this is done, we weight-average it with the value at time t+Future (not a forecast thereof!) of BiasedValues the cointegrated time series, and produce Mini and Maxi that will be input to T().
    double dForecastReal=gsl_matrix_get(ReflexiveValues, j, t); double Arma=0; double Linear=0; double Polynomial=0;
    double MeanGen=0; double MeanStep=0; double y0=0; double y1=0; double Poly=1; double DiscountSum=0; int Count=0;
    vector<double> Vec;
    if (Plot=="On") {outputLog << "We are in the case Tool=" << Tool;};
    if (t<(Lag+1)*Future) { // If not enough steps backwards everyone (Arma, Linear, Polynomial) is given the Mean of past reference values
        for (int n=0; n<=t; n++) {MeanGen+=gsl_matrix_get(ReflexiveValues, j, n);};
        dForecastReal=MeanGen/(t+1);
        if (Plot=="On") {outputLog << " and t=" << t << "<" << "(Lag+1)*Future-1=" << (Lag+1)*Future-1 << endl;};
        if (Plot=="On") {outputLog << "MeanGen=" << MeanGen << ", dForecastReal=" << dForecastReal << endl;};
    } // closes if
    else { // If enough steps backwards computing the mean over every (Lag+1) past array
        if (Plot=="On") {outputLog << " and t=" << t << ">=" << "(Lag+1)*Future-1=" << (Lag+1)*Future-1 << endl;};
        dForecastReal=BinaryProjection (ReflexiveValues, t, Tool, Lag, Future);
    }; // closes else
    //if (dForecastReal<0) {dForecastReal=0.01;}; // Just in case because of Polynomial easily having exaggerated expectations
    if (dForecastReal<=0) {dForecastReal=gsl_matrix_get(ReflexiveValues, j, t);}; // Failsafe for Polynomial exaggerated expectations
    if (Plot=="On") {outputLog << "dForecastReal=" << dForecastReal << endl;};
    //dForecastReal*=1/(1-(Future)*Rate/Year); // Time-discounting
    if (Plot=="On") {outputLog << "Time-discounted dForecastReal=" << dForecastReal << endl;};
    
    double dForecastVirtual=gsl_matrix_get(ReflexiveValues, j, t);
    if (Exploration==4) {
        Arma=0; Linear=0; Polynomial=0;
        MeanGen=0; MeanStep=0; y0=0; y1=0; Poly=1; DiscountSum=0; Count=0;
        vector<double> VecVirtual;
        if (Plot=="On") {outputLog << "We are in the case VirtualTool=" << VirtualTool;};
        if (t<(VirtualLag+1)*Future) { // If not enough steps backwards everyone (Arma, Linear, Polynomial) is given the Mean of past reference values
            for (int n=0; n<=t; n++) {MeanGen+=gsl_matrix_get(ReflexiveValues, j, n)/(t+1);};
            dForecastVirtual=MeanGen;
            if (Plot=="On") {outputLog << " and t=" << t << "<" << "(VirtualLag+1)*Future-1=" << (VirtualLag+1)*Future-1 << endl;};
        } // closes if
        else { // If enough steps backwards computing the mean over every (VirtualLag+1) past array
            dForecastVirtual=BinaryProjection (ReflexiveValues, t, VirtualTool, VirtualLag, Future);
        }; // closes else
        //if (dForecastVirtual<0) {dForecastVirtual=0.01;}; // Just in case because of Polynomial easily having exaggerated expectations
        if (dForecastVirtual<=0) {dForecastVirtual=gsl_matrix_get(ReflexiveValues, j, t);}; // Failsafe for Polynomial exaggerated expectations
        if (Plot=="On") {outputLog << "dForecastVirtual=" << dForecastVirtual << endl;};
        //dForecastVirtual*=1/(1-(Future)*Rate/Year); // Time-discounting
        if (Plot=="On") {outputLog << " Time-discounted dForecastVirtual=" << dForecastVirtual << endl;};
    }; // closes Exploration==4 condition
    
    //XXXExploration5
    if ((t-Future>=0) && (t%(2+Future/Month)==0)) {
        OffPolicy2=OffPolicy;
        if (Plot=="On") {outputLog << "We are in the case OffPolicy2=OffPolicy=" << OffPolicy << endl;};
    }
    else {
        OffPolicy2=0;
        if (Plot=="On") {outputLog << "We are in the case OffPolicy2=0" << endl;};
    };
    double dForecast5=gsl_matrix_get(ReflexiveValues, j, t); double MinimumCostFunction=9999999;
    int OptimalTool5=0; int OptimalLag5=0; int OptimalRWeight5=0;
    if (OffPolicy2==1) {
        for (int k1=0; k1<3; k1++) { // Tool
            for (int k2=0; k2<3; k2++) { // Lag
                for (int k3=0; k3<3; k3++) { // RWeight
                    dForecast5=BinaryProjection (ReflexiveValues, t-Future, k1, k2, Future);
                    double RWeight5=0;
                    if (Reflexivity<=0.5) {
                        if (k3==0) {RWeight5=0;}
                        else if (k3==1) {RWeight5=100*Reflexivity;}
                        else if (k3==2) {RWeight5=200*Reflexivity;};
                    }
                    else {
                        if (k3==0) {RWeight5=100*(2*Reflexivity-1);}
                        else if (k3==1) {RWeight5=100*Reflexivity;}
                        else if (k3==2) {RWeight5=100;};
                    };
                    double XBid = gsl_matrix_get (BiasedValues, j , t-Future);
                    double YBid = dForecast5;
                    double b = RWeight5;
                    double a = abs(100-b);
                    double ForecastReference5=(XBid*a + YBid*b)/(a+b);
                    double CostFunction=100*abs(gsl_matrix_get(ReflexiveValues, j, t)-ForecastReference5)/gsl_matrix_get(ReflexiveValues, j, t); // Comparison between what F() would have forecasted with k1, k2, k3 and present realization
                    if (MinimumCostFunction>=CostFunction) {MinimumCostFunction=CostFunction; OptimalTool5=k1; OptimalLag5=k2; OptimalRWeight5=k3;}
                    if (Plot=="On") {outputLog << "For Tool=" << k1 << ", Lag="<< k2 << ", RWeight="<< k3 << ": dForecast5=$" << dForecast5 << ", ForecastReference5=$" << ForecastReference5 << ", P(t)=$" << gsl_matrix_get(ReflexiveValues, j, t) << ", CostFunction=" << CostFunction << "%" << endl;};
                }; // closes k3
            }; // closes k2
        }; // closes k1
        if (Plot=="On") {outputLog << "MinimumCostFunction=" << MinimumCostFunction << " for OptimalTool5=" << OptimalTool5 << " (0: anti-trend, 1: mean-reversion, 2: trend), OptimalLag5=" << OptimalLag5 << ", OptimalRWeight5=" << OptimalRWeight5 << endl;};
    }; // closes OffPolicy2==1 condition
    
    // Computation of Mini and Maxi : minimum = (a<b) ? a : b;
    double RWeightReal=100*Reflexivity; double RWeightVirtual=100*Reflexivity;
    // Transforming Weight to real reflexivity coefficients RWeightReal
    if (Reflexivity<=0.5) {
        if (Weight==0) {RWeightReal=0;}
        else if (Weight==1) {RWeightReal=100*Reflexivity;}
        else if (Weight==2) {RWeightReal=200*Reflexivity;};
    }
    else {
        if (Weight==0) {RWeightReal=100*(2*Reflexivity-1);}
        else if (Weight==1) {RWeightReal=100*Reflexivity;}
        else if (Weight==2) {RWeightReal=100;};
    };
    //double XBid = gsl_matrix_get (BiasedValues, j , Time-1); if (t<Time-Future) {XBid = gsl_matrix_get (BiasedValues, j , t+Future);}; // Failsafe QQQ
    double XBid = gsl_matrix_get (BiasedValues, j , t); // Changed to present time (more rigorous approach)
    double YBid = dForecastReal;
    double b = RWeightReal;
    double a = abs(100-b);
    double ForecastReference=(XBid*a + YBid*b)/(a+b); // Now taking into account the dynamic reflexivity via Weight
    ForecastReal[j].push_back(ForecastReference); // Pushing the value into the vector
    if (int(ForecastReal[j].size()) - Future - 1 >0) {ForecastReal[j].erase (ForecastReal[j].begin());};
    if (Plot=="On") {
        outputLog << "Forecasted with ";
        if (Tool==0) {outputLog << "Anti-trend following";}
        else if (Tool==1) {outputLog << "Mean reversion";}
        else if (Tool==2) {outputLog << "Trend following";};
        outputLog << " ((Lag+1)*Future=" << (Lag+1)*Future << ", Future=" << Future << ")" << endl;
        outputLog << "ForecastReal[j][]={" ; for (int u=0; u<int(ForecastReal[j].size()); u++) {if (u<int(ForecastReal[j].size())-1) {outputLog << ForecastReal[j][u] << ", ";} else {outputLog << ForecastReal[j][u];};}; outputLog << "}" << endl;
    }; // closes Plot condition
    
    double ForecastReferenceVirtual=ForecastReference;
    if (Exploration==4) {
        if (Reflexivity<=0.5) {
            if (VirtualWeight==0) {RWeightVirtual=0;}
            else if (VirtualWeight==1) {RWeightVirtual=100*Reflexivity;}
            else if (VirtualWeight==2) {RWeightVirtual=200*Reflexivity;};
        }
        else {
            if (VirtualWeight==0) {RWeightVirtual=100*(2*Reflexivity-1);}
            else if (VirtualWeight==1) {RWeightVirtual=100*Reflexivity;}
            else if (VirtualWeight==2) {RWeightVirtual=100;};
        };
        //if (VirtualWeight==0) {RWeightVirtual=0*Reflexivity;} else if (VirtualWeight==1) {RWeightVirtual=50*Reflexivity;} else if (VirtualWeight==2) {RWeightVirtual=100*Reflexivity;}; // Transforming VirtualWeight to real reflexivity coefficients RWeightVirtual
        YBid = dForecastVirtual;
        b = RWeightVirtual;
        a = abs(100-b);
        ForecastReferenceVirtual=(XBid*a + YBid*b)/(a+b); // Now taking into account the dynamic reflexivity via Weight
        ForecastVirtual[j].push_back(ForecastReferenceVirtual); // Pushing the value into the vector
        if (int(ForecastVirtual[j].size()) - Future - 1 >0) {ForecastVirtual[j].erase (ForecastVirtual[j].begin());};
        if (Plot=="On") {outputLog << " to " << dForecastVirtual;};
        if (Plot=="On") {
            outputLog << ", forecasted with ";
            if (VirtualTool==0) {outputLog << "Anti-trend following";}
            else if (VirtualTool==1) {outputLog << "Mean reversion";}
            else if (VirtualTool==2) {outputLog << "Trend following";};
            outputLog << " ((VirtualLag+1)*Future=" << (VirtualLag+1)*Future << ", Future=" << Future << ")" << endl;
            outputLog << "ForecastVirtual[j][]={" ; for (int u=0; u<int(ForecastVirtual[j].size()); u++) {if (u<int(ForecastVirtual[j].size())-1) {outputLog << ForecastVirtual[j][u] << ", ";} else {outputLog << ForecastVirtual[j][u];};}; outputLog << "}" << endl;
        }; // closes Plot condition
    }; // closes Exploration==4 condition
    
    double ReflexiveVal_t = gsl_matrix_get (ReflexiveValues, j, t);
    //double ReflexiveVal2_t=ReflexiveVal_t*(1+Future*Rate/Year); // Time-discounting
    double ReflexiveVal2_t=ReflexiveVal_t;
    double Mini=ReflexiveVal2_t; double Maxi=ReflexiveVal2_t; int MiniTime=t; int MaxiTime=t+Future;
    if (ForecastReference<= ReflexiveVal2_t) {Mini=ForecastReference; Maxi=ReflexiveVal2_t; MiniTime=t+Future; MaxiTime=t;}
    else if (ForecastReference > ReflexiveVal2_t) {Maxi=ForecastReference; Mini=ReflexiveVal2_t; MaxiTime=t+Future; MiniTime=t;};
    Stocks[j].StockBidValue = Mini;
    Stocks[j].StockAskValue = Maxi;
    if (Trunk>-1) {
        Stocks[j].StockBidValue=DigitTrunk (Mini, Trunk, "Ceil");
        Stocks[j].StockAskValue=DigitTrunk (Maxi, Trunk, "Ceil");
    };
    
    if (Plot=="On") {
        outputLog << endl << "Reflexivity=" << Reflexivity << ", Weight=" << Weight << "=> RWeightReal=" << RWeightReal << "% as reflexivity coefficient, and Reflexive=" << Reflexive << endl;
        if (Exploration==4) {outputLog << "VirtualWeight=" << VirtualWeight << "=> RWeightVirtual=" << RWeightVirtual << "% as reflexivity coefficient..." << endl;};
        outputLog << endl << "Updating Bid and Ask:" << endl;
        outputLog <<"dForecastReal=" << dForecastReal << ", time-discounted P(t)=" << ReflexiveVal2_t << ", t=" << t << ", B(t)=" << gsl_matrix_get (BiasedValues, j , t) << ", ForecastReference=" << ForecastReference << endl;
        if (Exploration==4) {outputLog <<"dForecastVirtual=" << dForecastVirtual << ", time-discounted P(t)=" << ReflexiveVal2_t << ", t=" << t << ", ForecastReferenceVirtual=" << ForecastReferenceVirtual << endl;};
        outputLog << "Mini=Stocks[j].StockBidValue=" << Mini << ", MiniTime=" << MiniTime << endl;
        outputLog << "Maxi=Stocks[j].StockAskValue=" << Maxi << ", MaxiTime=" << MaxiTime << endl;
        
        
        
        /*************************************ABM I.7******************************************/
        outputLog << endl << "    ABM I.7: Calculating the return if t>Future as the percentage of difference between forecast and realization" << endl;
    }; // closes Plot condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.6 computed..." << endl;
    
    if (t==1) {
        dFResultVec[j].push_back(0);
        for (int i=0; i<FS*FA; i++) {FQ[j].push_back(0); FNumberA[j].push_back(0);};
        if (Plot=="On") {outputLog << "Initialization of FQ[], FNumberA[], by Agent " << AgentName << ", t=" << t << endl;};
    }; // closes if
    // We now calculate the return at time $t-T$ as the difference between the forecast and its realization
    if (t-Future<=0) {if (Plot=="On") {outputLog << "We are not in the case t>Future (t=" << t << ", Future=" << Future << ")" << endl;};}
    else if (t-Future>0) {
        if (Plot=="On") {
            outputLog << "We are in the case t=" << t << ">" << Future << "=Future (ForecastReal.size()=" << ForecastReal.size() << ")" << endl;
            outputLog << "At time t-Future=" << t-Future << ", real forecast was $" << ForecastReal[j][max(0,int(ForecastReal[j].size()) - Future - 1)] << " and at present time t=" << t << ", market price is $" << ReflexiveVal_t << endl;
            if (Exploration==4) {outputLog << " and virtual forecast was $" << ForecastVirtual[j][max(0,int(ForecastVirtual[j].size()) - Future - 1)] << endl;};
        }; // closes Plot condition
        double dFResultReal=0; double dFResultVirtual=0;
        //dFResultReal=-100*abs(ReflexiveVal_t - ForecastReal[int(ForecastReal.size()) - Future - 1])/ReflexiveVal_t; // minus the percentage of the absolute value of miss compared to market price.
        //dFResultVirtual=-100*abs(ReflexiveVal_t - ForecastVirtual[int(ForecastVirtual.size()) - Future - 1])/ReflexiveVal_t; // minus the percentage of the absolute value of miss compared to market price.
        dFResultReal=100*abs(gsl_matrix_get(ReflexiveValues, j, t) - ForecastReal[j][max(0,int(ForecastReal[j].size()) - Future - 1)])/gsl_matrix_get(ReflexiveValues, j, t); // minus the percentage of the absolute value of miss compared to market price.
        if (Exploration==4) {dFResultVirtual=100*abs(gsl_matrix_get(ReflexiveValues, j, t) - ForecastVirtual[j][max(0,int(ForecastVirtual[j].size()) - Future - 1)])/gsl_matrix_get(ReflexiveValues, j, t);}; // minus the percentage of the absolute value of miss compared to market price.
        //FResultReal.push_back(dFResultReal); // Recording result of forecast (ABM#6 for Forecast)
        //if (int(FResultReal.size()) - Future - 1 > 0) {FResultReal.erase (FResultReal.begin());}; // OOO
        //FResultVirtual.push_back(dFResultVirtual); // Recording result of forecast (ABM#6 for Forecast)
        //if (int(FResultVirtual.size()) - Future - 1 > 0) {FResultVirtual.erase (FResultVirtual.begin());}; // OOO
        
        // Finding the optimal action at time t-Future
        //if (OffPolicy2==1) // XXXExploration5
        
        
        
        if (Plot=="On") {
            outputLog << "So dFResultReal=" << dFResultReal << "%" << endl;
            if (Exploration==4) {outputLog << "dFResultVirtual=" << dFResultVirtual << "%" << endl;};
            //outputLog << "FResultReal[]={" ; for (int u=0; u<int(FResultReal.size()); u++) {if (u<int(FResultReal.size())-1) {outputLog << FResultReal[u] << ", ";} else {outputLog << FResultReal[u];};}; outputLog << "}" << endl;
            //outputLog << "FResultVirtual[]={" ; for (int u=0; u<int(FResultVirtual.size()); u++) {if (u<int(FResultVirtual.size())-1) {outputLog << FResultVirtual[u] << ", ";} else {outputLog << FResultVirtual[u];};}; outputLog << "}" << endl;
            
            
            /*************************************ABM I.8******************************************/
            outputLog << endl << "    ABM I.8: Returns above now discretized according to past returns (cf. reinforcement comparison & actor-critic methods)" << endl;
        }; // closes Plot condition
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.7 computed..." << endl;
        
        // We now discretise this result by giving it the following score -2, -1, +1, or +2 depending on where it stands in all past rewards

        double dFResultPercentile=0; int dFResultDisReal=0; int dFResultVecSize=int(dFResultVec[j].size());
        if (dFResultReal<0.001) {dFResultReal=0.001;}; // Failsafe
        for (int i=0; i<dFResultVecSize; i++) { // sorting the vector SmoothVec by ascending values
            if (dFResultReal<dFResultVec[j][i]) {dFResultVec[j].insert(dFResultVec[j].begin()+i, dFResultReal*(1 + 0.0001*VRan[25])); dFResultPercentile=i*1.0/dFResultVecSize; break;};
            if (i==dFResultVecSize-1) {dFResultVec[j].push_back(dFResultReal*(1 + 0.0001*VRan[25])); dFResultPercentile=1; break;};
        }; // closes for loop
        if (dFResultPercentile<0.05) {dFResultDisReal=+4;}
        else if ((dFResultPercentile>=0.05) && (dFResultPercentile<0.25)) {dFResultDisReal=2;}
        else if ((dFResultPercentile>=0.25) && (dFResultPercentile<0.5)) {dFResultDisReal=1;}
        else if ((dFResultPercentile>=0.5) && (dFResultPercentile<0.75)) {dFResultDisReal=-1;}
        else if ((dFResultPercentile>=0.75) && (dFResultPercentile<0.95)) {dFResultDisReal=-2;}
        else if (dFResultPercentile>=0.95) {dFResultDisReal=-4;};
        if (dFResultVecSize==1) {dFResultDisReal=1;}; // In the beginning it means nothing
        if ((dFResultReal<5) && (dFResultDisReal<0)) {dFResultDisReal=1;}; // Making sure a forecast smaller than 5% mismatch is not accounted as a negative reward XYZ333
        
        if (dFResultVecSize - Future - History >0) {dFResultVec[j].erase (dFResultVec[j].begin());}; // OOO
        FResultDisReal[j].push_back(dFResultDisReal); // LLL
        //if (int(FResultDisReal.size()) - Future - 1 >0) {FResultDisReal.erase (FResultDisReal.begin());}; // OOO
        

        if (Plot=="On") {
            //outputLog << "FRMean=" << FRMean << ", FRVariance=" << FRVariance << ", FRStdDev=" << FRStdDev << endl;
            //outputLog << "BBMinus=" << BBMinus << ", BBPlus=" << BBPlus << endl;
            outputLog << "dFResultDisReal=" << dFResultDisReal << " (dFResultReal=" << dFResultReal << ")" << endl;
            outputLog << "FResultDisReal[j][]={" ; for (int u=0; u<int(FResultDisReal[j].size()); u++) {if (u<int(FResultDisReal[j].size())-1) {outputLog << FResultDisReal[j][u] << ", ";} else {outputLog << FResultDisReal[j][u];};}; outputLog << "}" << endl;
        }; // closes Plot condition
        
        int dFResultDisVirtual=0;
        if (Exploration==4) {
            dFResultPercentile=0; dFResultVecSize=int(dFResultVec[j].size());
            if (dFResultVirtual<0.001) {dFResultVirtual=0.001;}; // Failsafe
            for (int i=0; i<dFResultVecSize; i++) { // sorting the vector by ascending values
                if (dFResultVirtual<dFResultVec[j][i]) {dFResultVec[j].insert(dFResultVec[j].begin()+i, dFResultVirtual*(1 + 0.0001*VRan[26])); dFResultPercentile=i*1.0/dFResultVecSize; break;};
                if (i==dFResultVecSize-1) {dFResultVec[j].push_back(dFResultVirtual*(1 + 0.0001*VRan[26])); dFResultPercentile=1; break;};
            }; // closes for loop
            if (dFResultPercentile<0.05) {dFResultDisVirtual=4;}
            else if ((dFResultPercentile>=0.05) && (dFResultPercentile<0.25)) {dFResultDisVirtual=2;}
            else if ((dFResultPercentile>=0.25) && (dFResultPercentile<0.5)) {dFResultDisVirtual=1;}
            else if ((dFResultPercentile>=0.5) && (dFResultPercentile<0.75)) {dFResultDisVirtual=-1;}
            else if ((dFResultPercentile>=0.75) && (dFResultPercentile<0.95)) {dFResultDisVirtual=-2;}
            else if (dFResultPercentile>=0.95) {dFResultDisVirtual=-4;};
            if (dFResultVecSize==1) {dFResultDisVirtual=1;}; // In the beginning it means nothing
            if (dFResultVecSize - Future - History >0) {dFResultVec[j].erase (dFResultVec[j].begin());}; // OOO
            //if (int(FResultDisVirtual.size()) - Future - 1 >0) {FResultDisVirtual.erase (FResultDisVirtual.begin());}; // OOO
            FResultDisVirtual[j].push_back(dFResultDisVirtual); // LLL
            if (Plot=="On") {
                //outputLog << "FRMeanVirtual=" << FRMeanVirtual << ", FRVarianceVirtual=" << FRVarianceVirtual << ", FRStdDevVirtual=" << FRStdDevVirtual << endl;
                //outputLog << "BBMinusVirtual=" << BBMinusVirtual << ", BBPlusVirtual=" << BBPlusVirtual << endl;
                outputLog << "dFResultDisVirtual=" << dFResultDisVirtual << " (dFResultVirtual=" << dFResultVirtual << ")" << endl;
                outputLog << "FResultDisVirtual[j][]={" ; for (int u=0; u<int(FResultDisVirtual[j].size()); u++) {if (u<int(FResultDisVirtual[j].size())-1) {outputLog << FResultDisVirtual[j][u] << ", ";} else {outputLog << FResultDisVirtual[j][u];};}; outputLog << "}" << endl;
            }; // closes Plot condition
        }; // closes if
        
        if ((NEBPositivity>VRan[17]) && (dFResultDisReal>0)) {
            double FirstdFResultDisReal=dFResultDisReal;
            dFResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
            if (Plot=="On") {outputLog << "NEBPositivity=" << NEBPositivity << " : positivity bias switched on => dFResultDisReal=" << FirstdFResultDisReal << "->" << dFResultDisReal << endl;};
        };
        
        if ((NEBNegativity>VRan[17]) && (dFResultDisReal<0)) {
            double FirstdFResultDisReal=dFResultDisReal;
            dFResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
            if (Plot=="On") {outputLog << "NEBNegativity=" << NEBNegativity << " : negativity bias switched on => dFResultDisReal=" << FirstdFResultDisReal << "->" << dFResultDisReal << endl;};
        };
        
        int dFResultDis5=4;
        if (OffPolicy2==1) {
            if (Plot=="On") {outputLog << "dFResultDis5=" << 4 << " (always) corresponding to a mismatch MinimumCostFunction=100*abs(P(t)-ˆP(t-T))/P(t)=" << MinimumCostFunction << "%" << endl;};
        };
        
        /*************************************ABM I.9******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM I.9: Updating Q(s,a) via SAM according to these discretized returns" << endl;};
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.8 computed..." << endl;
        
        // We now update the action-value function $Q(s,a)$ based on this return via the SAM or Sample Average Method (ABM#7)
        vector<int> OldStates = get_tensor_coord (FSIndex[j][max(0,int(FSIndex[j].size()) - Future - 1)], 3, FSDim); // Coord in s-tensor of past state $s_{t-T}$
        vector<int> OldRealActions = get_tensor_coord (FAIndexReal[j][max(0,int(FAIndexReal[j].size()) - Future - 1)], 3, FADim); //Coord in a-tensor of past action $a_{t-T}$
        vector<int> OldVirtualActions, OldActions5;
        int xxx=max(0, int(FResultDisReal[j].size()) - Future - 1);
        dFResultDisReal=FResultDisReal[j][xxx]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
        if (Exploration==4) {
            OldVirtualActions = get_tensor_coord (FAIndexVirtual[j][max(0,int(FAIndexVirtual[j].size()) - Future - 1)], 3, FADim); //Coord in a-tensor of past action $a_{t-T}$
            dFResultDisVirtual=FResultDisVirtual[j][max(0,int(FResultDisVirtual[j].size()) - Future - 1)]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
        };
        FAInd5[0]=OptimalTool5; FAInd5[1]=OptimalLag5; FAInd5[2]=OptimalRWeight5; int ATensorIndex5 = 0;
        if (OffPolicy2==1) {
            ATensorIndex5 = get_tensor_index (3, FAInd5, FADim);
            OldActions5 = get_tensor_coord (ATensorIndex5, 3, FADim);
        };
        
        if (Plot=="On") {
            outputLog << "OldStates[]={" ; for (int u=0; u<int(OldStates.size()); u++) {if (u<int(OldStates.size())-1) {outputLog << OldStates[u] << ", ";} else {outputLog << OldStates[u];};}; outputLog << "}" << endl;
            outputLog << "OldRealActions[]={" ; for (int u=0; u<int(OldRealActions.size()); u++) {if (u<int(OldRealActions.size())-1) {outputLog << OldRealActions[u] << ", ";} else {outputLog << OldRealActions[u];};}; outputLog << "}" << endl;
            if (Exploration==4) {outputLog << "OldVirtualActions[]={" ; for (int u=0; u<int(OldVirtualActions.size()); u++) {if (u<int(OldVirtualActions.size())-1) {outputLog << OldVirtualActions[u] << ", ";} else {outputLog << OldVirtualActions[u];};}; outputLog << "}" << endl;};
            if (OffPolicy2==1) {outputLog << "ATensorIndex5=" << ATensorIndex5 << endl;};
            outputLog << "FQ[], FNumberA[], initiated with size " << FNumberA.size() << " (=FS*FA)" << endl;
            outputLog << "We are at t=Future=" << t << endl;
        }; // closes Plot condition
        FQIndReal[0]=OldStates[0]; FQIndReal[1]=OldStates[1]; FQIndReal[2]=OldStates[2];
        FQIndReal[3]=OldRealActions[0]; FQIndReal[4]=OldRealActions[1]; FQIndReal[5]=OldRealActions[2];
        int QRealTensorIndex = 0;
        QRealTensorIndex = get_tensor_index (6, FQIndReal, FQDim); // Vector index of Q-tensor corresponding to $s_{t-T}$ x $a_{t-T}$
        int QVirtualTensorIndex=0;
        if (Exploration==4) {
            FQIndVirtual[0]=OldStates[0]; FQIndVirtual[1]=OldStates[1]; FQIndVirtual[2]=OldStates[2]; FQIndVirtual[3]=OldVirtualActions[0]; FQIndVirtual[4]=OldVirtualActions[1]; FQIndVirtual[5]=OldVirtualActions[2];
            QVirtualTensorIndex = get_tensor_index (6, FQIndVirtual, FQDim); // Vector index of Q-tensor corresponding to $s_{t-T}$ x $a_{t-T}$
        };
        int QTensorIndex5=0;
        if (OffPolicy2==1) {
            FQInd5[0]=OldStates[0]; FQInd5[1]=OldStates[1]; FQInd5[2]=OldStates[2]; FQInd5[3]=OldActions5[0]; FQInd5[4]=OldActions5[1]; FQInd5[5]=OldActions5[2];
            QTensorIndex5 = get_tensor_index (6, FQInd5, FQDim); // Vector index of Q-tensor corresponding to $s_{t-T}$ x $a_{t-T}$
        };
        if (Plot=="On") {
            outputLog << "FQIndReal[]={" ; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FQIndReal[u] << ", ";} else {outputLog << FQIndReal[u];};}; outputLog << "}" << endl;
            if (Exploration==4) {outputLog << "FQIndVirtual[]={" ; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FQIndVirtual[u] << ", ";} else {outputLog << FQIndVirtual[u];};}; outputLog << "}" << endl;};
            outputLog << "FQDim[]={" ; for (int u=0; u<6; u++) {if (u<6-1) {outputLog << FQDim[u] << ", ";} else {outputLog << FQDim[u];};}; outputLog << "}" << endl;
            outputLog << "QRealTensorIndex=" << QRealTensorIndex << endl;
            if (Exploration==4) {outputLog << "QVirtualTensorIndex=" << QVirtualTensorIndex << endl;};
            if (OffPolicy2==1) {outputLog << "QTensorIndex5=" << QTensorIndex5 << endl;};
        }; // closes Plot condition
        
        FNumberA[j][QRealTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
        if (FNumberA[j][QRealTensorIndex]==1) {FQ[j][QRealTensorIndex]=dFResultReal;} // No SAM if not yet populated (i.e =0)
        else {FQ[j][QRealTensorIndex]=(FQ[j][QRealTensorIndex]*(FNumberA[j][QRealTensorIndex]-1) + dFResultReal)/FNumberA[j][QRealTensorIndex];}; // SAM
        if (Plot=="On") {
            outputLog << "Number of times this real action was previously taken in that former state: FNumberA[j][QRealTensorIndex]-1=" << FNumberA[j][QRealTensorIndex]-1 << ", so FQ[j][QRealTensorIndex]=" << FQ[j][QRealTensorIndex] << endl;
        }; // closes plot condition
        
        if (Exploration==4) {
            FNumberA[j][QVirtualTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
            if (FNumberA[j][QVirtualTensorIndex]==1) {FQ[j][QVirtualTensorIndex]=dFResultVirtual;} // No SAM if not yet populated
            else {FQ[j][QVirtualTensorIndex]=(FQ[j][QVirtualTensorIndex]*(FNumberA[j][QVirtualTensorIndex]-1) + dFResultVirtual)/FNumberA[j][QVirtualTensorIndex];}; // SAM
            if (Plot=="On") {
                outputLog << "Number of times this virtual action was previously taken in that former state: FNumberA[j][QVirtualTensorIndex]-1=" << FNumberA[j][QVirtualTensorIndex]-1 << ", so FQ[j][QVirtualTensorIndex]=" << FQ[j][QVirtualTensorIndex] << endl;
            }; // closes plot condition
        }; // closes if
        
        if (OffPolicy2==1) {
            FNumberA[j][QTensorIndex5]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
            if (FNumberA[j][QTensorIndex5]==1) {FQ[j][QTensorIndex5]=MinimumCostFunction;} // No SAM if not yet populated
            else {FQ[j][QTensorIndex5]=(FQ[j][QTensorIndex5]*(FNumberA[j][QTensorIndex5]-1) + MinimumCostFunction)/FNumberA[j][QTensorIndex5];}; // SAM
            if (Plot=="On") {
                outputLog << "Number of times this Exploration5 action was previously taken in that former state: FNumberA[j][QTensorIndex5]-1=" << FNumberA[j][QTensorIndex5]-1 << ", so FQ[j][QTensorIndex5]=" << FQ[j][QTensorIndex5] << endl;
            }; // closes plot condition
        }; // closes if
        
        
        
        /*************************************ABM I.10******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM I.10: Updating pi(s,a) according to these discretized returns" << endl;};
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.9 computed..." << endl;
        
        // We now update existing policy based on this discretized results FResultDisReal and FResultDisVirtual (ABM#8)
        FSInd[0]=OldStates[0]; FSInd[1]=OldStates[1]; FSInd[2]=OldStates[2];
        int OldSTensorIndex = get_tensor_index (3, FSInd, FSDim); // Vector index of the s-tensor associated with $s_{t-T}$
        int OldSCount=0;
        for (int i=0; i<int(FSIndex[j].size())-Future; i++) {if (FSIndex[j][i]==OldSTensorIndex) {OldSCount+=1;};}; // Nb of times $s_{t-T}$ seen
        if (Plot=="On") {
            outputLog << "State at time t-Future: FSInd[0]=OldStates[0]=" << FSInd[0] << ", FSInd[1]=OldStates[1]=" << FSInd[1] << ", FSInd[2]=OldStates[2]=" << FSInd[2] << endl;
            outputLog << "OldSTensorIndex=" << OldSTensorIndex << endl;
            outputLog << "Number of times state at FSIndex.size()-Future encountered: OldSCount=" << OldSCount << endl;
        }; // closes Plot condition
        FQIndReal[0]=OldStates[0]; FQIndReal[1]=OldStates[1]; FQIndReal[2]=OldStates[2]; // Using results from above on Q
        FQIndReal[3]=0; FQIndReal[4]=0; FQIndReal[5]=0; // As start of the Pi vector index to span
        int PiSTensorIndex = get_tensor_index (6, FQIndReal, FQDim); // Vector index of the Q-tensor associated with $s_{t-T}$
        if (Plot=="On") {
            outputLog << "FQIndReal[0]=OldStates[0]=" << FQIndReal[0] << ", FQIndReal[1]=OldStates[1]=" << FQIndReal[1] << ", FQIndReal[2]=OldStates[2]=" << FQIndReal[2] << endl;
            outputLog << "FQIndReal[3]=" << FQIndReal[3] << ", FQIndReal[4]=" << FQIndReal[4] << ", FQIndReal[5]=" << FQIndReal[5] << ", which is PiSTensorIndex=" << PiSTensorIndex << " from which we span all actions to update their probability" << endl;
            //outputLog << "State s was encountered OldSCount=" << OldSCount << " times" << endl;
            outputLog << "Results was dFResultDisReal=" << dFResultDisReal << endl;
            if (Exploration==4) {outputLog << "And dFResultDisVirtual=" << dFResultDisVirtual << endl;};
            if (OffPolicy2==1) {outputLog << "And dFResultDis5=" << dFResultDis5 << endl;};
        }; // closes Plot condition
        double Sum=0;
        for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
            //if (Plot=="On") {outputLog << "Old FPi[j][" << i << "]=" << floor(100*FPi[j][i]) << "%" << endl;};
            Sum+=FPi[j][i];
        };
        if (Plot=="On") {
            outputLog << "Sum=" << 100*Sum << endl;
            if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorOriginalF[]" << endl;};
        };
        Sum=0;
        
        //if (1-NEB8>VRan[19]) {
        // First for real
        for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
            if (Plot=="On") {outputLog << "Old real FPi[j][" << i << "]=" << 100*FPi[j][i] << "% => ";};
            if ((i==QRealTensorIndex) && (dFResultDisReal>=0)) {
                for (int f=0; f<dFResultDisReal; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]) + NEBLearningRate/FNumberA[j][QRealTensorIndex];}; // DDD3
            } // closes if
            else if ((i!=QRealTensorIndex) && (dFResultDisReal>=0)) {
                for (int f=0; f<dFResultDisReal; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]);};
            } // closes else if
            else if ((i==QRealTensorIndex) && (dFResultDisReal<0)) {
                for (int f=0; f<abs(dFResultDisReal); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]);};
            } // closes else if
            else if ((i!=QRealTensorIndex) && (dFResultDisReal<0)) {
                for (int f=0; f<abs(dFResultDisReal); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QRealTensorIndex]) + NEBLearningRate/(FA-1)/FNumberA[j][QRealTensorIndex];};
            }; // closes else if
            if (Plot=="On") {outputLog << "New real FPi[j][" << i << "]=" << 100*FPi[j][i] << "%" << endl;};
            Sum+=FPi[j][i];
        }; // closes i loop
        if (Plot=="On") {
            outputLog << "(QRealTensorIndex=" << QRealTensorIndex << " was the real action update), FNumberA[j][QRealTensorIndex]=" << FNumberA[j][QRealTensorIndex] << endl;
            outputLog << "Sum=" << 100*Sum << endl;
            if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorRealF[]: For OldSCount=" << OldSCount << ", dFResultDisReal=" << dFResultDisReal << endl;};
        }; // closes Plot condition
        Sum=0;
        //}; // closes if
        // Second for virtual
        //if ((Exploration==4) && (1-NEB8>VRan[19])) { // This is for off-policy Watkin's Q-learning
        if (Exploration==4) { // This is for off-policy Watkin's Q-learning
            for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
                if (Plot=="On") {outputLog << "Old virtual FPi[j][" << i << "]=" << 100*FPi[j][i] << "% => ";};
                if ((i==QVirtualTensorIndex) && (dFResultDisVirtual>=0)) {
                    for (int f=0; f<dFResultDisVirtual; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/FNumberA[j][QVirtualTensorIndex];};
                } // closes if
                else if ((i!=QVirtualTensorIndex) && (dFResultDisVirtual>=0)) {
                    for (int f=0; f<dFResultDisVirtual; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]);};
                } // closes else if
                else if ((i==QVirtualTensorIndex) && (dFResultDisVirtual<0)) {
                    for (int f=0; f<abs(dFResultDisVirtual); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]);};
                } // closes else if
                else if ((i!=QVirtualTensorIndex) && (dFResultDisVirtual<0)) {
                    for (int f=0; f<abs(dFResultDisVirtual); f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/(FA-1)/FNumberA[j][QVirtualTensorIndex];};
                }; // closes else if
                if (Plot=="On") {outputLog << "New virtual FPi[j][" << i << "]=" << 100*FPi[j][i] << "%" << endl;};
                Sum+=FPi[j][i];
            }; // closes i loop
            if (Plot=="On") {
                outputLog << "(QVirtualTensorIndex=" << QVirtualTensorIndex << " was the real action update), FNumberA[j][QVirtualTensorIndex]=" << FNumberA[j][QVirtualTensorIndex] << endl;
                outputLog << "Sum=" << 100*Sum << endl; if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorVirtualF[]: For OldSCount=" << OldSCount << ", dFResultDisVirtual=" << dFResultDisVirtual << endl;};
            }; // closes Plot condition
            Sum=0;
        }; // closes if Exploration
        // Second for Exploration5
        if (OffPolicy2==1) { // This is for off-policy Watkin's Q-learning
            for (int i=PiSTensorIndex; i<PiSTensorIndex+FA; i++) {
                if (Plot=="On") {outputLog << "Old Exploration5 FPi[j][" << i << "]=" << 100*FPi[j][i] << "% => ";};
                if (i==QTensorIndex5) {
                    for (int f=0; f<dFResultDis5; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QTensorIndex5]) + NEBLearningRate/FNumberA[j][QTensorIndex5];};
                } // closes if
                else if (i!=QTensorIndex5) {
                    for (int f=0; f<dFResultDis5; f++) {FPi[j][i]=FPi[j][i]*(1-NEBLearningRate/FNumberA[j][QTensorIndex5]);};
                } // closes else if
                if (Plot=="On") {outputLog << "New Exploration5 FPi[j][" << i << "]=" << 100*FPi[j][i] << "%" << endl;};
                Sum+=FPi[j][i];
            }; // closes i loop
            if (Plot=="On") {
                outputLog << "(QTensorIndex5=" << QTensorIndex5 << " was the real action update), FNumberA[j][QTensorIndex5]=" << FNumberA[j][QTensorIndex5] << endl;
                outputLog << "Sum=" << 100*Sum << endl; if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "Error5F[]: For OldSCount=" << OldSCount << ", dFResultDis5=" << dFResultDis5 << endl;};
            }; // closes Plot condition
            Sum=0;
        }; // closes if Exploration
    }; // closes if (t-Future>0)
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 1.10 computed..." << endl;
    
    
    
    
    
    
    /*******************************************************************************/
    /*******************************************************************************/
    /*******************************************************************************/
    /***********************************RL TRADING**********************************/
    /*******************************************************************************/
    /*******************************************************************************/
    /*******************************************************************************/
    
    /*************************************ABM II.1**************************************************/
    if (Plot=="On") {outputLog << endl << endl << "    ABM 2.1: Defining present state s" << endl;};
    // Griding growth forecast results for RL states s
    double dMu=100*(ForecastReference-gsl_matrix_get(ReflexiveValues, j, t))/(gsl_matrix_get(ReflexiveValues, j, t));
    int dMuPosVecSize=int(dMuPosVec[j].size()); int dMuNegVecSize=int(dMuNegVec[j].size());
    double dMuPosPercentile=0; double dMuNegPercentile=0;
    if (t==1) {dMuPosVec[j].push_back(0); dMuPosPercentile=0.5; dMuNegVec[j].push_back(0); dMuNegPercentile=0.5;}
    else if ((t>1) && (dMu>=0)) {
        for (int i=0; i<dMuPosVecSize; i++) {
            if (dMu<dMuPosVec[j][i]) {dMuPosVec[j].insert(dMuPosVec[j].begin()+i, dMu*(1 + 0.0001*VRan[30])); dMuPosPercentile=i*1.0/dMuPosVecSize; break;};
            if (i==dMuPosVecSize-1) {dMuPosVec[j].push_back(dMu*(1 + 0.0001*VRan[30])); dMuPosPercentile=1; break;};
        }; // closes for loop
        
        /*
            if (dMuPosPercentile<0.05) {Mu=3;}
            else if ((dMuPosPercentile>=0.05) && (dMuPosPercentile<0.33)) {Mu=4;}
            else if ((dMuPosPercentile>=0.33) && (dMuPosPercentile<0.67)) {Mu=5;}
            else if (dMuPosPercentile>=0.67) {Mu=6;};
            */
        
        if (dMuPosPercentile<0.05) {Mu=1;}
        else if (dMuPosPercentile>=0.05) {Mu=2;};
        // NEB4 implementation
        //if ((NEB4p>VRan[16]) && (NEB4==0)) {if (dMuPosPercentile>=0.025) {Mu=2;};}; // Regressive conservatism: overestimating high values & likelihoods while underestimating low values & likelihoods
        //if ((NEB4p>VRan[16]) && (NEB4==1)) {if (dMuPosPercentile<=0.075) {Mu=1;};}; // Regressive conservatism: underestimating high values & likelihoods while overestimating low values & likelihoods
        
    } // closes else if
    else if ((t>1) && (dMu<0)) {
        for (int i=0; i<dMuNegVecSize; i++) {
            if (dMu<dMuNegVec[j][i]) {dMuNegVec[j].insert(dMuNegVec[j].begin()+i, dMu*(1 + 0.0001*VRan[30])); dMuNegPercentile=i*1.0/dMuNegVecSize; break;};
            if (i==dMuNegVecSize-1) {dMuNegVec[j].push_back(dMu*(1 + 0.0001*VRan[30])); dMuNegPercentile=1; break;};
        }; // closes for loop
        
        /*
            if (dMuNegPercentile<0.33) {Mu=0;}
            else if ((dMuNegPercentile>=0.33) && (dMuNegPercentile<0.67)) {Mu=1;}
            else if ((dMuNegPercentile>=0.67) && (dMuNegPercentile<0.95)) {Mu=2;}
            else if (dMuNegPercentile>=0.95) {Mu=3;};
            */
        
        if (dMuNegPercentile<=0.95) {Mu=0;}
        else if (dMuNegPercentile>0.95) {Mu=1;};
        // NEB4 implementation
        //if ((NEB4p>VRan[16]) && (NEB4==0)) {if (dMuNegPercentile<=0.975) {Mu=0;}}; // Regressive conservatism: overestimating high values & likelihoods while underestimating low values & likelihoods
        //if ((NEB4p>VRan[16]) && (NEB4==1)) {if (dMuNegPercentile>=0.925) {Mu=1;};}; // Regressive conservatism: underestimating high values & likelihoods while overestimating low values & likelihoods
        
    }; // closes else if
    if (int(dMuPosVec[j].size()) - Future - History >0) {dMuPosVec[j].erase (dMuPosVec[j].begin());}; // OOO
    if (int(dMuNegVec[j].size()) - Future - History >0) {dMuNegVec[j].erase (dMuNegVec[j].begin());}; // OOO
    
    // NEB10 implementation
    //if ((NEB10>VRan[18]) && (Mu<2)) {Mu+=1;
    //    if (Plot=="On") {outputLog << "NEB10=" << NEB10 << " => DelayDiscounting : Mu=" << Mu-1 << "->" << Mu << endl;};
    //};
    
    
    // Now griding volatility forecast results for RL states s
    //double dSig=100*Lvol/(gsl_matrix_get(ReferenceValues, j, t));
    double dSig=Lvol;
    int dSigVecSize=int(dSigVec[j].size()); double dSigPercentile=0;
    if (t==1) {dSigVec[j].push_back(dSig); dSigPercentile=0.5;}
    else {
        for (int i=0; i<dSigVecSize; i++) {
            if (dSig<dSigVec[j][i]) {dSigVec[j].insert(dSigVec[j].begin()+i, dSig*(1 + 0.0001*VRan[31])); dSigPercentile=i*1.0/dSigVecSize; break;};
            if (i==dSigVecSize-1) {dSigVec[j].push_back(dSig*(1 + 0.0001*VRan[31])); dSigPercentile=1; break;};
        }; // closes for loop
        if (dSigPercentile<0.33) {Sig=0;}
        else if ((dSigPercentile>=0.33) && (dSigPercentile<0.67)) {Sig=1;}
        else if (dSigPercentile>=0.67) {Sig=2;};
    } // closes else
    if (int(dSigVec[j].size()) - Future - History >0) {dSigVec[j].erase (dSigVec[j].begin());}; // OOO
    
    
    
    
    if (t==1) {RFAFirst=RFA;}; // To make sure that it is the market price at t=1 that acts as standard
    if ((RFAFirst==0) && (t>=2*Week)) {RFAFirst=RFA;}; // Failsafe if agent has RFAFirst=0
    /*
        if ((RFA>=0) && (RFA<0.25*RFAFirst)) {RFALevel=0;}
        else if ((RFA>=0.25*RFAFirst) && (RFA<0.6*RFAFirst)) {RFALevel=1;}
        else if (RFA>=0.6*RFAFirst) {RFALevel=2;};
        */
    if ((RFA>=0) && (RFA<0.6*RFAFirst)) {RFALevel=0;}
    else if (RFA>=0.6*RFAFirst) {RFALevel=1;};
    
    // NEB10 implementation
    //if (NEB10>VRan[18]) {RFALevel=1;
    //if (Plot=="On") {outputLog << "NEB10=" << NEB10 << " => DelayDiscounting : RFALevel set at 1" << endl;};
    //};
    
    
    if (Plot=="On") {RFALog[j].push_back(RFALevel);}; // LLL
    // Now griding RBA since it is a parameter for shorting stocks
    RBA=StockHoldings();
    if (t==1) {RBAFirst=RBA;}; // To make sure that it is the market price at t=1 that acts as standard
    if ((RBAFirst==0) && (t>=2*Week)) {RBAFirst=RBA;}; // Failsafe if agent has Q=0 stock holdings and hence RBAFirst=0
    /*
        if ((RBA>=0) && (RBA<0.25*RBAFirst)) {RBALevel=0;}
        else if ((RBA>=0.25*RBAFirst) && (RBA<0.6*RBAFirst)) {RBALevel=1;}
        else if (RBA>=0.6*RBAFirst) {RBALevel=2;};
        */
    if ((RBA>=0) && (RBA<0.6*RBAFirst)) {RBALevel=0;}
    else if (RBA>=0.6*RBAFirst) {RBALevel=1;};
    
    // NEB10 implementation
    //if (NEB10>VRan[18]) {RBALevel=1;
    //if (Plot=="On") {outputLog << "NEB10=" << NEB10 << " => DelayDiscounting : RBALevel set at 1" << endl;};
    //};
    
    if (Plot=="On") {RBALog[j].push_back(RBALevel);}; // LLL
    // Now griding Liquid as it is a parameter for the Law of Supply and Offer and Gesture parameter
    int LiquidPercentVecSize=int(LiquidPercentVec[j].size()); double LiquidPercentPercentile=0;
    if (t==1) {LiquidPercentVec[j].push_back(LiquidPercent); LiquidPercentPercentile=0.5;}
    else {
        for (int i=0; i<LiquidPercentVecSize; i++) { // Important: LiquidPercentVec built out of non-zero values of LiquidPercent!
            if ((LiquidPercent<LiquidPercentVec[j][i]) && (LiquidPercent!=0)) {LiquidPercentVec[j].insert(LiquidPercentVec[j].begin()+i, LiquidPercent*(1 + 0.0001*VRan[32])); LiquidPercentPercentile=i*1.0/LiquidPercentVecSize; break;};
            if ((i==LiquidPercentVecSize-1) && (LiquidPercent!=0)) {LiquidPercentVec[j].push_back(LiquidPercent*(1 + 0.0001*VRan[32])); LiquidPercentPercentile=1; break;};
        }; // closes for loop
        /*
            if (LiquidPercentPercentile<0.33) {Liquid=1;}
            else if ((LiquidPercentPercentile>=0.33) && (LiquidPercentPercentile<0.67)) {Liquid=2;}
            else if (LiquidPercentPercentile>=0.67) {Liquid=3;};
            if (LiquidPercent==0) {Liquid=0;};
            */
        if (LiquidPercentPercentile<0.33) {Liquid=1;}
        else if (LiquidPercentPercentile>=0.33) {Liquid=2;};
        if (LiquidPercent==0) {Liquid=0;};
    } // closes else
    if (int(LiquidPercentVec[j].size()) - Future - History >0) {LiquidPercentVec[j].erase (LiquidPercentVec[j].begin());}; // OOO

    
    if (t==1) {Mu=1; Sig=1; Liquid=0;};
    TSInd[0]=Mu; TSInd[1]=Sig; TSInd[2]=RFALevel; TSInd[3]=RBALevel; TSInd[4]=Liquid;
    int TSid=get_tensor_index (5, TSInd, TSDim);// This gives the state index tensor for our second RL algorithm
    TSIndex[j].push_back(TSid); // Back up
    if (int(TSIndex[j].size()) - Future - 1 >0) {TSIndex[j].erase (TSIndex[j].begin());}; // OOO
    if (Plot=="On") {
        outputLog << "dMu=" << dMu << ", dSig=" << dSig << ", RFA=" << RFA << " (RFAFirst=" << RFAFirst << "), RBA=" << RBA << " (RBAFirst=" << RBAFirst << "), LiquidPercent=" << LiquidPercent << endl;
        outputLog << "Mu=" << Mu << ", Sig=" << Sig << ", RFALevel=" << RFALevel << ", RBALevel=" << RBALevel << ", Liquid=" << Liquid << endl;
        outputLog << "TSDim[]={" ; for (int u=0; u<5; u++) {if (u<5-1) {outputLog << TSDim[u] << ", ";} else {outputLog << TSDim[u];};}; outputLog << "}" << endl;
        outputLog << "TSInd[]={" ; for (int u=0; u<5; u++) {if (u<5-1) {outputLog << TSInd[u] << ", ";} else {outputLog << TSInd[u];};}; outputLog << "}" << endl;
        outputLog << "TSid=" << TSid << endl;
        outputLog << "TSIndex[j][]={" ; for (int u=0; u<int(TSIndex[j].size()); u++) {if (u<int(TSIndex[j].size())-1) {outputLog << TSIndex[j][u] << ", ";} else {outputLog << TSIndex[j][u];};}; outputLog << "}" << endl;
        MuLog[j].push_back(Mu); // LLL
        SigLog[j].push_back(Sig); // LLL
        LiquidityLog[j].push_back(Liquid); // LLL
        
        
        /*************************************ABM II.2******************************************/
        outputLog << endl << "    ABM 2.2: Initializing the start policy pi_0" << endl;
    }; // closes Plot condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.1 computed..." << endl;
    
    if (t==1) {
        TPiInd[0]=1; TPiInd[1]=0; TPiInd[2]=1; TPiInd[3]=1; TPiInd[4]=1; // State s
        TPiInd[5]=1; TPiInd[6]=1; // Actions a
        for (int i=0; i<TS*TA; i++) {TPi[j].push_back(1.0/TA);}; // HHH : All actions have equiprobability
        if (Plot=="On") {outputLog << "Policy initialized... all " << TA << " possible actions are given equiprobability 1/" << TA << "=" << 1.0/TA << endl;};
    }
    else if (t>1) {
        if (Plot=="On") {outputLog << "Policy already initialized..." << endl;};
    };
    if (Plot=="On") {
        outputLog << "(To reduce output size, the PiTensorIndex vector indices and their corresponding TPiInd[] coordinates in pi-tensor space are not displayed)" << endl;
        
        
        /*************************************ABM II.3******************************************/
        outputLog << endl << "    ABM 2.3: Using current policy pi to select real action a=(Quant,Pinch) in present state s" << endl;
    }; // closes Plot condition
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.2 computed...";
    
    // We now use the current policy to perform an action $a$ in our state $s$. Note: There is no need to use \verb?PiStack?, all we need is when we have a vector of action each with a probability in $[0,1]$, we generate a random number from a uniform distribution: if that number is larger than the first probability it is subtracted from it, and if that new number is still larger than the second probability, we do likewise until we arrive to the probability corresponding to action $a$ to be picked up. This is because the sum of all these probabilities is $1$.
    TTenPiId[0]=Mu; TTenPiId[1]=Sig; TTenPiId[2]=RFALevel; TTenPiId[3]=RBALevel; TTenPiId[4]=Liquid; TTenPiId[5]=0; TTenPiId[6]=0;
    int Tid=0;
    Tid = get_tensor_index (7, TTenPiId, TPiDim); // Vector index of pi-tensor representation corresponding to present s and a=(0,0)
    if (Plot=="On") {
        outputLog << "TTenPiId[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TTenPiId[u] << ", ";} else {outputLog << TTenPiId[u];};}; outputLog << "}" << endl;
        outputLog << "Tid=" << Tid << endl;
    }; // closes Plot condition
    
    double Tix=VRan[6];
    double Tix2=VRan[7];
    Act=0; Pos=Tid;
    if (Plot=="On") {outputLog << "Tix=" << Tix << ", Tix2=" << Tix2 << endl;};
    // outputDebug << ": a, ";
    if (Plot=="On") {outputLog << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.2" << endl;};
    
    if (Exploration==1) { // Epsilon-greedy: Choosing always the action with highest probability, but once in a while (with a small probability epsilon) choosing a (uniform) random one (which is a problem: the worst potential action can be taken as likely as the best, with in some cases devastating returns).
        // outputDebug << "b1, ";
        if (Plot=="On") {outputLog << "Epsilon-greedy method selected" << endl;};
        if (Tix<Epsilon) { // Exploration
            Pos=Tid + int(floor(TA*Tix2));
            Act=TPi[j][Pos];
            if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }
        else { // Exploitation
            for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                if (TPi[j][n]>=Act) {Pos=n; Act=TPi[j][Pos];}; // We try and find the action with greatest probability (greedy move) PPP
            }; // closes for
            vector<int> PosVec; // Indices of equiprobable actions
            for (int n=Tid; n<Tid+TA; n++) {if (TPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
            Pos=PosVec[int(floor(VRan[24]*PosVec.size()))]; // We choose a random indice among these!
            PosVec.erase(PosVec.begin(),PosVec.end());
            if (Plot=="On") {outputLog << "Exploitation epsilon-greedy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else
    } // closes if
    
    else if (Exploration==2) { // Softmax: Same as above except the random choice in now not (uniform) equiprobable for all actions, but graded according to their estimated value. One method (using a Gibbs/Boltzmann distribution $e^{Q_t(a) / \tau} / \sum_{b=1}^n e^{Q_t(b) / \tau}$, with $\tau$ the temperature) allows to even shift this grading all the way back to epsilon greedy methods (when $\tau \rightarrow 0$).
        // outputDebug << "b2, ";
        if (Plot=="On") {outputLog << "Softmax method selected" << endl;};
        if (Tix<Epsilon) { // Exploration
            Pos=Tid + int(floor(TA*Tix2));
            Act=TPi[j][Pos];
            if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }
        else { // Exploitation
            vector<int> N = Shuffle(TA);
            if (Plot=="On") {outputLog << "From n=" << Tid << " to " << Tid+TA-1 << ": Tix2=" << Tix2 << endl;};
            for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                int Newn = N[n-Tid]+Tid;
                if (Plot=="On") {outputLog << "For n=" << Newn << ", TPi[j][Newn]=" << floor(100*TPi[j][Newn]) << "%: Tix2-=TPi[j][Newn]=" << Tix2 << endl;};
                Tix2-=TPi[j][Newn];
                if (Tix2<=0) {Pos=Newn; Act=TPi[j][Pos]; break;}; // We try and find a new action according to its graded probability
            }; // closes for
            if (Plot=="On") {outputLog << "Exploitation softmax: Chose action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else
    } // closes else if
    
    else if (Exploration==3) { // Pursuit: The probability of selecting the greedy action $a_{t+1}=a_{t+1}^{\ast}$ at next time step is incremented a fraction $\beta$ of the way toward $1$, while the probabilities of selecting the other actions are decremented toward $0$, respectively shown by the following two expressions: $\pi_{t+1}(a_{t+1}^{\ast}) = \pi_{t}(a_{t+1}^{\ast}) + \beta \left[ 1- \pi_{t}(a_{t+1}^{\ast}) \right]$ and $\pi_{t+1}(a) = \pi_{t}(a) + \beta \left[ 0- \pi_{t}(a) \right] , \hspace{5mm} \forall a\neq a_{t+1}^{\ast}$.
        // outputDebug << "b3, ";
        if (Plot=="On") {outputLog << "Pursuit method selected" << endl;};
        if (Tix<Epsilon*(1-(double(t)/Time))) { // Exploration
            Pos=Tid + int(floor(TA*Tix2));
            Act=TPi[j][Pos];
            if (Plot=="On") {outputLog << "Exploration: Chose random action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }
        else { // Exploitation
            for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
                if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
                if (TPi[j][n]>=Act) {Pos=n; Act=TPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
            }; // closes for
            vector<int> PosVec; // Indices of equiprobable actions
            for (int n=Tid; n<Tid+TA; n++) {if (TPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
            Pos=PosVec[int(floor(VRan[24]*PosVec.size()))]; // We choose a random indice among these!
            PosVec.erase(PosVec.begin(),PosVec.end());
            if (Plot=="On") {outputLog << "Exploitation pursuit: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
        }; // closes else
    } // closes else if
    
    else if (Exploration==4) { // This is for off-policy Watkin's Q-learning
        // outputDebug << "b4, ";
        if (Plot=="On") {outputLog << "Off-policy Watkin's Q-learning method selected" << endl;};
        for (int n=Tid; n<Tid+TA; n++) { // Running over the actions a associated with that state s
            if (Plot=="On") {outputLog << "For n=" << n << ", TPi[j][n]=" << floor(100*TPi[j][n]) << "%" << endl;};
            if (TPi[j][n]>=Act) {Pos=n; Act=TPi[j][Pos];}; // We try and find the action with greatest probability (greedy move)
        }; // closes for
        vector<int> PosVec; // Indices of equiprobable actions
        for (int n=Tid; n<Tid+TA; n++) {if (TPi[j][n]==Act) {PosVec.push_back(n);};}; // We record all the indices of equiprobable actions
        Pos=PosVec[int(floor(VRan[24]*PosVec.size()))]; // We choose a random indice among these!
        PosVec.erase(PosVec.begin(),PosVec.end());
        if (Plot=="On") {outputLog << "Exploitation Off-policy: Chose greediest action at n=Pos=" << Pos << " with probability Act=" << Act << endl;};
    }; // closes else if
    // outputDebug << "b, ";
    
    
    int Greenlight=1;
    if (Plot=="On") {outputLog << "Model-based RL switches Greenlight=1 if current state s is good wrt. previous states representations according to arg max_a Q(s,a)" << endl;};
    if (Plot=="On") {outputLog << "TradingWindowClock[j]/TradingWindow=" << TradingWindowClock[j] << "/" << TradingWindow << endl << "Greenlight=" << Greenlight << endl;};
    if ((TradingFrequencyCond=="On") && (t-Future>0)) {
        Greenlight=0; TradingWindowClock[j]+=1;
        if (TradingWindowClock[j]==1) {
            Qs[j].erase(Qs[j].begin(),Qs[j].end()); // Completely updating the model-based RL framework at every new trading window
            double Max=-999999;
            for (int u=0; u<TS*TA; u++) {
                if (u%TA==0) {Max=-999999;};
                if (TQ[j][u]>=Max) {Max=TQ[j][u];};
                if ((u%TA==TA-1) && (Max!=0)) {Qs[j].push_back(Max);};
                //if (u%TA==TA-1) {Qs[j].push_back(Max);};
            }; // closes u loop
            //string A = "TQ_"; string B=to_string(t); string C=".txt"; A+=B+C; if (j==0) {PlotSTL(Qs[j], A.c_str(), "0");};
        }; // closes if
        double TSPercentile=0; double QsPresent=-999999; // QsPresent=arg max_a Q(s,a) for present state s
        for (int u=Tid; u<Tid+TA; u++) {if (TQ[j][u]>=QsPresent) {QsPresent=TQ[j][u];};}; // Defining QsPresent
        for (int u=0; u<int(Qs[j].size()); u++) {if (QsPresent>=Qs[j][u]) {TSPercentile+=1;};};
        TSPercentile*=100.0/int(Qs[j].size());
        if (Qs[j].size()==0) {TSPercentile=0;}; // Failsafe
        if (((TradingWindowClock[j]*100.0/TradingWindow)>=TSPercentile) || (Qs[j].size()==0) || (TradingWindowClock[j]==TradingWindow)) {Greenlight=1; TradingWindowClock[j]=0;};
        if (Plot=="On") {outputLog << "Qs[j].size()=" << Qs[j].size() << ", QsPresent=" << QsPresent << endl << "100*t/T=" << (TradingWindowClock[j]*100.0/TradingWindow) << " >= TSPercentile=" << TSPercentile << endl << "TradingWindowClock[j]/TradingWindow=" << TradingWindowClock[j] << "/" << TradingWindow << endl << "Greenlight=" << Greenlight << endl;};
    }; // closes if
    // outputDebug << "c, ";
    
    vector<int> TenPiCoo2 = get_tensor_coord (Pos, 7, TPiDim); // Getting the tensor coordinates from that index
    Quant=TenPiCoo2[5]; // Selection of first action from policy
    Pinch=TenPiCoo2[6]; // Selection of second action from policy
    if (Greenlight==0) {
        Quant=1; TenPiCoo2[5]=Quant;
        if (Plot=="On") {outputLog << "Greenlight is OFF..." << endl;}; // Holding position if the current state is not great enough according to BestStates
    }
    else {if (Plot=="On") {
        outputLog << "Greenlight is ON..." << endl;};
        ExitHorizons.push_back(t+Future); // JJJ4
    };

    
    if (NEBLossAversion>VRan[13]) {
        int FirstPinch=Pinch;
        if (Pinch==0) {Pinch=1;} // Loss aversion asking more to sell than to buy
        else {Pinch=2;};
        if (Plot=="On") {outputLog << "NEBLossAversion=" << NEBLossAversion << " : loss aversion => Pinch=" << FirstPinch << "->" << Pinch << endl;};
    }; // closes if
    
    if (NEBFear>VRan[13]) {
        int FirstQuant=Quant;
        if ((Sig==2) || (RFALevel==0) || (Liquid==0)) {Quant=0;}
        if (Plot=="On") {outputLog << "NEBFear=" << NEBFear << " : fear => Quant=" << FirstQuant << "->" << Quant << endl;};
    }; // closes if
    
    if (NEBGreed>VRan[13]) {
        int FirstQuant=Quant;
        if (Mu==2) {Quant=2;}
        if (Plot=="On") {outputLog << "NEBGreed=" << NEBGreed << " : fear => Quant=" << FirstQuant << "->" << Quant << endl;};
    }; // closes if
    
    
    // Cluster trading
    if ((AgentName<ClusterLimit) && (AgentName!=LeaderAgent) && (t>=LearningPhase)) {// WWW
        Quant=LeaderQuant;
        Pinch=LeaderPinch;
        cout << AgentName << ", ";
        if (Plot=="On") {outputLog << "i=" << AgentName << ", LeaderAgent=" << LeaderAgent << ", ClusterLimit=" << ClusterLimit << " => Quant=LeaderQuant=" << LeaderQuant << ", Pinch=LeaderPinch=" << LeaderPinch << " from LeaderAgent=" << LeaderAgent << endl;};
    };
    
    
    
    if (Plot=="On") {
        outputLog << "TenPiCoo2[]={" ; for (int u=0; u<int(TenPiCoo2.size()); u++) {if (u<int(TenPiCoo2.size())-1) {outputLog << TenPiCoo2[u] << ", ";} else {outputLog << TenPiCoo2[u];};}; outputLog << "}, ";
        outputLog << "Quant=TenPiCoo2[5]=" << TenPiCoo2[5] << ", Pinch=TenPiCoo2[6]=" << TenPiCoo2[6] << endl;
    }; // closes Plot condition
    // Defining the real action as tensor index to take according to above policy which selected Quant and Pinch
    if (t==1) {TAIndReal[0]=1; TAIndReal[1]=0;}
    else {TAIndReal[0]=Quant; TAIndReal[1]=Pinch;}; // Chosen in the policy above
    if (Plot=="On") {
        outputLog << "TADim[]={" ; for (int u=0; u<2; u++) {if (u<2-1) {outputLog << TADim[u] << ", ";} else {outputLog << TADim[u];};}; outputLog << "}" << endl;
        outputLog << "TAIndReal[]={" ; for (int u=0; u<2; u++) {if (u<2-1) {outputLog << TAIndReal[u] << ", ";} else {outputLog << TAIndReal[u];};}; outputLog << "}" << endl;
    }; // closes Plot condition
    // outputDebug << "d, ";
    
    ARealTensorIndex = get_tensor_index (2, TAIndReal, TADim); // Vector index in the tensor of all possible RL actions
    TAIndexReal[j].push_back(ARealTensorIndex); // Recording real action a_t (ABM#2 for Forecast)
    if (int(TAIndexReal[j].size()) - Future - 1 >0) {TAIndexReal[j].erase (TAIndexReal[j].begin());}; // OOO
    if (Plot=="On") {
        outputLog << "ARealTensorIndex=" << ARealTensorIndex << endl;
        outputLog << "TAIndexReal[j][]={" ; for (int u=0; u<int(TAIndexReal[j].size()); u++) {if (u<int(TAIndexReal[j].size())-1)  {outputLog << TAIndexReal[j][u] << ", ";} else {outputLog << TAIndexReal[j][u];};}; outputLog << "}, ";
        // Update the Bid and Ask that will be filed in the OB to take into account the Pinch
        outputLog << endl << "Updating the Bid and Ask that will be filed in the OB to take into account the Pinch" << endl;
    }; // closes Plot condition
    //double Spread=abs(Stocks[j].StockAskValue - Stocks[j].StockBidValue); // GGG
    double Spread=VSpread; // VSpread is the absolute value of average top bid minus average top ask
    if (Plot=="On") {
        outputLog << "ARealTensorIndex=" << ARealTensorIndex << endl;
        outputLog << "Stocks[j].StockBidValue=Mini=" << Stocks[j].StockBidValue << endl;
        outputLog << "Stocks[j].StockAskValue=Maxi=" << Stocks[j].StockAskValue << endl;
    }; // closes Plot condition
    // outputDebug << "e, ";
    
    // Gesture=0.2+0.8r
    // Stocks[j].StockBidValue = min(Pt,ˆPt)±Gesture*Spread
    // Stocks[j].StockAskValue = max(Pt,ˆPt)±Gesture*Spread
    
    double TransactionVirtualBid=ReflexiveVal_t;
    double TransactionVirtualAsk=ReflexiveVal_t;
    if (Exploration==4) {TransactionVirtualBid=Stocks[j].StockBidValue; TransactionVirtualAsk=Stocks[j].StockAskValue;}; // Backup to compute the virtual return (see below)
    if (Plot=="On") {outputLog << "Since Spread=" << Spread << "$, Pinch=" << Pinch << ", and Gesture=" << floor(Gesture*100-100) << "% of the Spread, we have:" << endl;};
    
    
    if (Pinch==0) {Stocks[j].StockBidValue+=Gesture*Spread; Stocks[j].StockAskValue-=Gesture*Spread;} // This corresponds to Pinch=-Gesture (more desperate)
    else if (Pinch==2) {Stocks[j].StockBidValue-=Gesture*Spread; Stocks[j].StockAskValue+=Gesture*Spread;}; // This corresponds to Pinch=+Gesture (less desperate)
    
    if (Trunk>-1) {
        if (Pinch==0) { // This corresponds to Pinch=-Gesture (more desperate)
            Stocks[j].StockBidValue=DigitTrunk (Mini+Gesture*Spread, Trunk, "Ceil");
            Stocks[j].StockAskValue=DigitTrunk (Maxi-Gesture*Spread, Trunk, "Ceil");
        }
        else if (Pinch==2) { // This corresponds to Pinch=+Gesture (less desperate)
            Stocks[j].StockBidValue=DigitTrunk (Mini-Gesture*Spread, Trunk, "Ceil");
            Stocks[j].StockAskValue=DigitTrunk (Maxi+Gesture*Spread, Trunk, "Ceil");
        };
        if (abs(Stocks[j].StockBidValue-ceil(Stocks[j].StockBidValue))>0) {cout << "ISSUE : Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << endl;};
        if (abs(Stocks[j].StockAskValue-ceil(Stocks[j].StockAskValue))>0) {cout << "ISSUE : Stocks[j].StockAskValue=" << Stocks[j].StockAskValue << endl;};
    };
    
    if (Stocks[j].StockBidValue<0) {Stocks[j].StockBidValue=0;}; // Failsafe for the Pinch effect
    if (Stocks[j].StockAskValue<0) {Stocks[j].StockAskValue=0;}; //Failsafe for the Pinch effect
    if (Plot=="On") {
        outputLog << "Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << endl;
        outputLog << "Stocks[j].StockAskValue=" << Stocks[j].StockAskValue << endl;
    }; // closes Plot condition
    double RealTransac=ReflexiveVal_t;
    int QuantityFraction=1;
    double AskVal=Stocks[j].StockQuantity/QuantityFraction;
    double BidVal=RFA/((NumberOfStocks+1)*Stocks[j].StockAskValue)/QuantityFraction; // Cannot buy more than its RFA divided by number of stocks
    //double BidVal=RFA/((NumberOfStocks)*Stocks[j].StockAskValue); // Cannot buy more than its RFA divided by number of stocks // GGG
    if ((RFA<0) || (Stocks[j].StockAskValue==0)) {BidVal=0;}; // Failsafe
    if (Quant==0) {QTradeReal=-int(AskVal); RealTransac=Stocks[j].StockAskValue;}
    else if (Quant==1) {QTradeReal=0; RealTransac=ReflexiveVal_t;}
    else if (Quant==2) {QTradeReal=int(BidVal); RealTransac=Stocks[j].StockBidValue;}
    QuantitiesReal[j].push_back(QTradeReal); // By default, but if there is an actual transaction then it is in OB's trading
    if (int(QuantitiesReal[j].size()) - Future - 1 >0) {QuantitiesReal[j].erase (QuantitiesReal[j].begin());}; // OOO
    TransactionPriceReal[j].push_back(RealTransac); // By default, but if there is an actual transaction then it is in OB's trading
    if (int(TransactionPriceReal[j].size()) - Future - 1 >0) {TransactionPriceReal[j].erase (TransactionPriceReal[j].begin());}; // OOO
    // outputDebug << "f, ";
    
    
    //ofstream outputTest("/Users/admin/Documents/GNT/SYMBA/ZULU.txt", ofstream::app); // ZZZ
    //if ((QTradeReal>0) && (QTradeReal>RFA/Stocks[j].StockBidValue)) {outputTest << "ErrorX Bid at i=" << AgentName << ", j=" << j << ", t=" << t << ": QTradeReal=" << QTradeReal << ", RFA/Stocks[j].StockBidValue=" << RFA << "/" << Stocks[j].StockBidValue << "=" << RFA/Stocks[j].StockBidValue << endl;};
    //if ((QTradeReal<0) && (abs(QTradeReal)>Stocks[j].StockQuantity)) {outputTest << "ErrorY Ask at i=" << AgentName << ", j=" << j << ", t=" << t << ": QTradeReal=" << QTradeReal << ", Stocks[j].StockQuantity=" << Stocks[j].StockQuantity << endl;};
    if (Plot=="On") {
        outputLog << "Since Quant=" << Quant << ", QTradeReal=" << QTradeReal << " (Stocks[j].StockQuantity=" << Stocks[j].StockQuantity << ", RFA=" << RFA << ", and Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << ")";
        outputLog << endl << "QuantitiesReal[j][]={" ; for (int u=0; u<int(QuantitiesReal[j].size()); u++) {if (u<int(QuantitiesReal[j].size())-1)  {outputLog << QuantitiesReal[j][u] << ", ";} else {outputLog << QuantitiesReal[j][u];};}; outputLog << "}" << endl;
        outputLog << "RealTransac=" << RealTransac;
        outputLog << endl << "TransactionPriceReal[j][]={" ; for (int u=0; u<int(TransactionPriceReal[j].size()); u++) {if (u<int(TransactionPriceReal[j].size())-1)  {outputLog << TransactionPriceReal[j][u] << ", ";} else {outputLog << TransactionPriceReal[j][u];};}; outputLog << "}" << endl;
        QuantRealLog[j].push_back(Quant); // LLL
        PinchRealLog[j].push_back(Pinch); // LLL
    }; // closes Plot condition
    // outputDebug << "g" << endl;
    
    
    
    // Exploration5 retro-infering the right actions OptimalQuant and OptimalPinch at time t-Future by off-policy learning
    int OptimalQuant=0; int OptimalPinch=0; int TAid5 = 0;
    if (OffPolicy2==1) {
        // Finding OptimalQuant
        if (gsl_matrix_get (ReflexiveValues, j, t)>gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalQuant=2; // Buy
            if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)>gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=0;} // Liberal gesture
            else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)==gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=1;} // Neutral gesture
            else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)<gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=2;}; // Tough gesture
        }
        
        else if (gsl_matrix_get (ReflexiveValues, j, t)==gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalQuant=1; OptimalPinch=1;}// Hold, neutral gesture
        
        else if (gsl_matrix_get (ReflexiveValues, j, t)<gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalQuant=0; // Sell
            if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)>gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=2;} // Tough gesture
            else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)==gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=1;} // Neutral gesture
            else if (gsl_matrix_get (ReflexiveValues, j, t-Future+1)<gsl_matrix_get (ReflexiveValues, j, t-Future)) {OptimalPinch=0;}; // Liberal gesture
        };
        // Finding OptimalPinch (issue with hold because all Pinch should technically be updated equiprobably, but it's ok if we select one only)
        if (Plot=="On") {outputLog << "Since P(t)=" << gsl_matrix_get (ReflexiveValues, j, t) << ", P(t-T)=" << gsl_matrix_get (ReflexiveValues, j, t-Future) << ", P(t-T+1)=" << gsl_matrix_get (ReflexiveValues, j, t-Future+1) << ": OptimalQuant=" << OptimalQuant << " (0: short, 1: hold, 2: long) and OptimalPinch=" << OptimalPinch << " (0: liberal, 1: neutral, 2: tough gesture)" << endl;};
        TAInd5[0]=OptimalQuant; TAInd5[1]=OptimalPinch;
        if (Plot=="On") {outputLog << "TAInd5[0]=OptimalQuant=" << TAInd5[0] << ", TAInd5[1]=OptimalPinch=" << TAInd5[1] << endl;};
        TAid5 = get_tensor_index (2, TAInd5, TADim); // Vector index in the tensor of all possible RL actions
        //TAIndex5[j].push_back(TAid5); // Recording real action a_t (ABM#2 for Trading)
        //if (int(TAIndex5[j].size()) - Future - 1 >0) {TAIndex5[j].erase (TAIndex5[j].begin());};
    }; // closes if
    
    
    
    
    /*************************************ABM II.4******************************************/
    if (Plot=="On") {outputLog << endl << "    ABM 2.4: Randomly selecting a virtual action a=(VirtualQuant,VirtualPinch)" << endl;};
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.3 computed..." << endl;
    
    // Defining the virtual action as tensor index to take (not based on any exploration policy)
    int VirtualQuant, VirtualPinch;
    if (Exploration==4) {
        VirtualQuant=int(1000*VRan[8])%3; // 0, 1, 2, 3, 4, 5, or 6
        VirtualPinch=int(1000*VRan[9])%3; // 0, 1, or 2
        for (int i=0; i<100; i++) {int VQ=int(1000*VRan[10])%3; if (VQ!=Quant) {VirtualQuant=VQ; break;};}; // Making sure virtual action not like real one
        for (int i=0; i<100; i++) {int VP=int(1000*VRan[11])%3; if (VP!=Lag) {VirtualPinch=VP; break;};}; // Making sure virtual action not like real one
        if (Plot=="On") {outputLog << "VirtualQuant=" << VirtualQuant << ", VirtualPinch=" << VirtualPinch << endl;};
        double SpreadVirtual=abs(TransactionVirtualAsk - TransactionVirtualBid);
        if (VirtualPinch==0) {TransactionVirtualBid+=Gesture*SpreadVirtual; TransactionVirtualAsk-=Gesture*SpreadVirtual;} // This corresponds to VirtualPinch=-Gesture
        else if (VirtualPinch==2) {TransactionVirtualBid-=Gesture*SpreadVirtual; TransactionVirtualAsk+=Gesture*SpreadVirtual;}; // This corresponds to VirtualPinch=+Gesture
        if (TransactionVirtualBid<0) {TransactionVirtualBid=0;}; // Failsafe for the VirtualPinch effect
        if (TransactionVirtualAsk<0) {TransactionVirtualAsk=0;}; // Failsafe for the VirtualPinch effect
        double VirtualTransac=ReflexiveVal_t;
        if (VirtualQuant==0) {QTradeVirtual=-Stocks[j].StockQuantity; VirtualTransac=TransactionVirtualAsk;}
        else if (VirtualQuant==1) {QTradeVirtual=0; VirtualTransac=ReflexiveVal_t;}
        else if (VirtualQuant==2) {QTradeVirtual=int(floor(1*RFA/Stocks[j].StockBidValue)); VirtualTransac=TransactionVirtualBid;};
        if ((Stocks[j].StockBidValue==0) && (VirtualQuant==2)) {QTradeVirtual=0;} // For the case RFA/Stocks[j].StockBidValue=#IND
        if (Plot=="On") {outputLog << "Since VirtualQuant=" << VirtualQuant << ", QTradeVirtual=" << QTradeVirtual << " (Stocks[j].StockQuantity=" << Stocks[j].StockQuantity << ", RFA=" << RFA << ", and Stocks[j].StockBidValue=" << Stocks[j].StockBidValue << ")" << endl;};
        //if (t==1) {TAIndVirtual[0]=3; TAIndVirtual[1]=0;}
        TAIndVirtual[0]=VirtualQuant; TAIndVirtual[1]=VirtualPinch; // Chosen in the policy above
        if (Plot=="On") {outputLog << "TAIndVirtual[0]=VirtualQuant=" << TAIndVirtual[0] << ", TAIndVirtual[1]=VirtualPinch=" << TAIndVirtual[1] << endl;};
        int TAid = get_tensor_index (2, TAIndVirtual, TADim); // Vector index in the tensor of all possible RL actions
        TAIndexVirtual[j].push_back(TAid); // Recording real action a_t (ABM#2 for Trading)
        if (int(TAIndexVirtual[j].size()) - Future - 1 >0) {TAIndexVirtual[j].erase (TAIndexVirtual[j].begin());}; // OOO
        QuantitiesVirtual[j].push_back(QTradeVirtual); // This doest need to be in OB's trading
        if (int(QuantitiesVirtual[j].size()) - Future - 1 >0) {QuantitiesVirtual[j].erase (QuantitiesVirtual[j].begin());}; // OOO
        TransactionPriceVirtual[j].push_back(VirtualTransac); // By default, but if there is an actual transaction then it is in OB's trading
        if (int(TransactionPriceVirtual[j].size()) - Future - 1 >0) {TransactionPriceVirtual[j].erase (TransactionPriceVirtual[j].begin());}; // OOO
        if (Plot=="On") {
            outputLog << "TAid=" << TAid << endl;
            outputLog << "TAIndexVirtual[j][]={" ; for (int u=0; u<int(TAIndexVirtual[j].size()); u++) {if (u<int(TAIndexVirtual[j].size())-1)  {outputLog << TAIndexVirtual[j][u] << ", ";} else {outputLog << TAIndexVirtual[j][u];};}; outputLog << "}" << endl;
            outputLog << "QTradeVirtual=" << QTradeVirtual;
            outputLog << endl << "QuantitiesVirtual[j][]={" ; for (int u=0; u<int(QuantitiesVirtual[j].size()); u++) {if (u<int(QuantitiesVirtual[j].size())-1)  {outputLog << QuantitiesVirtual[j][u] << ", ";} else {outputLog << QuantitiesVirtual[j][u];};}; outputLog << "}" << endl;
            outputLog << "VirtualTransac=" << VirtualTransac << " (=ReflexiveVal_t)" << endl;
            outputLog << "TransactionPriceVirtual[j][]={" ; for (int u=0; u<int(TransactionPriceVirtual[j].size()); u++) {if (u<int(TransactionPriceVirtual[j].size())-1)  {outputLog << TransactionPriceVirtual[j][u] << ", ";} else {outputLog << TransactionPriceVirtual[j][u];};}; outputLog << "}" << endl;
            QuantVirtualLog[j].push_back(VirtualQuant); // LLL
            PinchVirtualLog[j].push_back(VirtualPinch); // LLL
        }; // closes Plot condition
    }; // closes Exploration==4 condition
    
    /*************************************ABM II.5******************************************/
    if (Plot=="On") {outputLog << endl << "    ABM 2.5: Calculating the return if t>Future as the cashflow difference if the action had been taken and if it had not" << endl;};
    
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.4 computed..." << endl;
    double DividendYield=2.276*0.01/Year;
    double BrokerFees=0.01*0.01; // 0.01% on the LSE
    if (t==1) {
        dTResultVec[j].push_back(0);
        for (int i=0; i<TS*TA; i++) {TQ[j].push_back(0); TNumberA[j].push_back(0);};
        if (Plot=="On") {outputLog << "Initialization of TQ[j][], TNumberA[j][], by Agent " << AgentName << ", t=" << t << endl;};
    };
    // We now calculate the return at time $t-T$. Let TransacQ(t-T) be the quantity of transaction at time $t-T$, TransacRFA(t-T) be the cash flow from RFA due to transaction at time $t-T$. Then the return is given by NAV(t) - NAV'(t), where NAV'(t) is the NAV(t) if the transaction at time $t-T$ had not taken place (Q(t-T) [P(t) - T(t-T)]), given by RFA(t) - TransacRFA(t-T) + (Q(t) - TransacQ(t-T))*MarketPrice(t). So the return is simply TransacRFA(t-T) + TransacQ(t-T)*MarketPrice(t)
    if (t-Future<=0) {if (Plot=="On") {outputLog << "We are not in the case t>Future (t=" << t << ", Future=" << Future << ")" << endl;}}
    else if (t-Future>0) {
        if (Plot=="On") {outputLog << "We are in the case t=" << t << ">" << Future << "=Future" << endl;};
        double dTResultReal=0; double dTResultVirtual=0;
        int Q=QuantitiesReal[j][max(0,int(QuantitiesReal[j].size()) - Future - 1)];
        double PastTransacPt=TransactionPriceReal[j][max(0,int(TransactionPriceReal[j].size()) - Future - 1)];
        dTResultReal=Q*(ReflexiveVal_t-PastTransacPt); // TransacRFA(t-T) + TransacQ(t-T)*MarketPrice(t). Note: the transaction price should be the one realized at OB (may be cancelled!) // JJJ
        dTResultReal+=Q*PastTransacPt*DividendYield*Future; // Dividends
        dTResultReal-=2*abs(Q)*PastTransacPt*BrokerFees; // Broker fees, one for in and one for out of position
        if (Exploration==4) {dTResultVirtual=QuantitiesVirtual[j][max(0,int(QuantitiesVirtual[j].size()) - Future - 1)] * (ReflexiveVal_t - TransactionPriceVirtual[j][max(0,int(TransactionPriceVirtual[j].size()) - Future - 1)]);};
        //TResultReal.push_back(dTResultReal);
        //if (int(TResultReal.size()) - Future - 1 >0) {TResultReal.erase (TResultReal.begin());}; // OOO
        //TResultVirtual.push_back(dTResultVirtual);
        //if (int(TResultVirtual.size()) - Future - 1 >0) {TResultVirtual.erase (TResultVirtual.begin());}; // OOO
        double dYield=100.0*(dTResultReal*NumberOfStocks*Year/Future)/(RFA+RBA); // Returns made by this action in % of agent NAV extended as adjusted return to all year
        double dMarketYield=(MarketPerformance-100)*Year/TimeSinceJan1st; // Or should we use Reflexive_t instead?
        double dTResult5=0;
        // Not really necessary since dTdisResult5 is always set at 4, regardless of dTResult5
        if (OffPolicy2==1) {
            int PseudoQuantity=ceil(RFAFirst/((NumberOfStocks+1)*ReflexiveVal_t)/QuantityFraction);
            dTResult5=PseudoQuantity * (ReflexiveVal_t - gsl_matrix_get (ReflexiveValues, j, t-Future)); // Approximation of Exploration5 by QuantitiesReal // JJJ
            dTResult5+=PseudoQuantity * gsl_matrix_get (ReflexiveValues, j, t-Future) * DividendYield * Future;
            dTResult5-=2*abs(PseudoQuantity) * gsl_matrix_get (ReflexiveValues, j, t-Future) * BrokerFees; // One for in and one for out of position
            if (Plot=="On") {outputLog << "PseudoQuantity=" << PseudoQuantity << ", P(t)=" << ReflexiveVal_t << ", P(t-T)=" << gsl_matrix_get (ReflexiveValues, j, t-Future) << " => dTResult5=" << dTResult5 << endl;};
        };
        if (Plot=="On") {
            outputLog << "At time t=" << t-Future << ", real placed order was $" << TransactionPriceReal[j][max(0,int(TransactionPriceReal[j].size()) - Future - 1)] << " (for a quantity " << QuantitiesReal[j][max(0,int(QuantitiesReal[j].size()) - Future - 1)] << " thereof)" << endl;
            if (Exploration==4) {outputLog << "At time t=" << t-Future << ", virtual placed order was $" << TransactionPriceVirtual[j][max(0,int(TransactionPriceVirtual[j].size()) - Future - 1)] << " (for a quantity " << QuantitiesVirtual[j][max(0,int(QuantitiesVirtual[j].size()) - Future - 1)] << " thereof)" << endl;};
            outputLog << "Finally today it is $" << ReflexiveVal_t << endl;
            outputLog << "So we get dTResultReal=$" << dTResultReal << ", dYield=" << dYield << "%, dMarketYield=" << dMarketYield << "%" << endl;
            if (Exploration==4) {outputLog << "And dTResultVirtual=$" << dTResultVirtual << endl;};
            //outputLog << "TResultReal[]={" ; for (int u=0; u<int(TResultReal.size()); u++) {if (u<int(TResultReal.size())-1)  {outputLog << TResultReal[u] << ", ";} else {outputLog << TResultReal[u];};}; outputLog << "}" << endl;
            //outputLog << "TResultVirtual[]={" ; for (int u=0; u<int(TResultVirtual.size()); u++) {if (u<int(TResultVirtual.size())-1)  {outputLog << TResultVirtual[u] << ", ";} else {outputLog << TResultVirtual[u];};}; outputLog << "}" << endl;
            
            
            
            
            /*************************************ABM II.6******************************************/
            outputLog << endl << "    ABM 2.6: Returns above now discretized according to past returns (cf. reinforcement comparison & actor-critic methods)" << endl;
        }; // closes Plot condition

        
        double dTResultPercentile=0; int dTResultVecSize=int(dTResultVec[j].size()); int dTResultDisReal=0;
        bool CondRecord=0;
        for (int i=0; i<int(ExitHorizons.size()); i++) {
            if (t==ExitHorizons[i]) {CondRecord=1;};
        };
        if (CondRecord==1) { // JJJ3
            for (int i=0; i<dTResultVecSize; i++) { // sorting the vector SmoothVec by ascending values
                if (dTResultReal<dTResultVec[j][i]) {
                    dTResultVec[j].insert(dTResultVec[j].begin()+i, dTResultReal*(1 + 0.0001*VRan[27]));
                    dTResultPercentile=i*1.0/dTResultVecSize;
                    break;
                };
                if (i==dTResultVecSize-1) {
                    dTResultVec[j].push_back(dTResultReal*(1 + 0.0001*VRan[27]));
                    dTResultPercentile=1;
                    break;
                };
            }; // closes for loop
            if ((dTResultReal>0) && (dTResultDisReal<0)) {dTResultDisReal=1;}; // Making sure a net profit is not accounted as a negative reward XYZ333
        }; // closes CondRecord condition
        
        if (dTResultPercentile<0.05) {dTResultDisReal=-4;}
        else if ((dTResultPercentile>=0.05) && (dTResultPercentile<0.25)) {dTResultDisReal=-2;}
        else if ((dTResultPercentile>=0.25) && (dTResultPercentile<0.5)) {dTResultDisReal=-1;}
        else if ((dTResultPercentile>=0.5) && (dTResultPercentile<0.75)) {dTResultDisReal=1;}
        else if ((dTResultPercentile>=0.75) && (dTResultPercentile<0.95)) {dTResultDisReal=2;}
        else if (dTResultPercentile>=0.95) {dTResultDisReal=4;};
        if (dTResultVecSize==1) {dTResultDisReal=1;}; // In the beginning it means nothing
        if (dTResultVecSize - Future - History >0) {dTResultVec[j].erase (dTResultVec[j].begin());}; // OOO
        

        if ((NEBPositivity>VRan[17]) && (dTResultDisReal>0)) {
            double FirstdTResultDisReal=dTResultDisReal;
            dTResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
            if (Plot=="On") {outputLog << "NEBPositivity=" << NEBPositivity << " : positivity bias switched on => dTResultDisReal=" << FirstdTResultDisReal << "->" << dTResultDisReal << endl;};
        };
        
        if ((NEBNegativity>VRan[17]) && (dTResultDisReal<0)) {
            double FirstdTResultDisReal=dTResultDisReal;
            dTResultDisReal+=1; // Positivity bias: better recalling pleasant memories than unpleasant ones
            if (Plot=="On") {outputLog << "NEBNegativity=" << NEBNegativity << " : negativity bias switched on => dTResultDisReal=" << FirstdTResultDisReal << "->" << dTResultDisReal << endl;};
        };
        
        if ((dTResultReal>0) && (dTResultDisReal<0)) {
            dTResultDisReal=1; // Making sure a net profit is not accounted as a negative reward
            if (Plot=="On") {outputLog << "Making sure a net profit is not accounted as a negative reward => dTResultDisReal=1" << endl;};
        };
        TResultDisReal[j].push_back(dTResultDisReal); // LLL
        
        //if (int(TResultDisReal.size()) - Future - 1 >0) {TResultDisReal.erase (TResultDisReal.begin());}; // OOO
        /*
            ofstream outputCheck("/Users/admin/Documents/GNT/SYMBA/Check.txt", ofstream::app);
            if (AgentName==5) {
            outputCheck << "At t=" << t << ": dTResultReal=" << dTResultReal*(1 + 0.0001*VRan[27]) << ", dTResultRealPercentile=" << ceil(100*dTResultRealPercentile) << "% => dTResultDisReal=" << dTResultDisReal << endl;
            outputCheck << "dTResultVec[]={" ; for (int u=0; u<int(dTResultVec.size()); u++) {if (u<int(dTResultVec.size())-1) {outputCheck << dTResultVec[u] << ", ";} else {outputCheck << dTResultVec[u];};}; outputCheck << "}" << endl;
            };
            */
        if (Plot=="On") {
            //outputLog << "TRMean=" << TRMean << ", TRVariance=" << TRVariance << ", TRStdDev=" << TRStdDev << endl;
            //outputLog << "BBMinus=" << BBMinus << ", BBPlus=" << BBPlus << endl;
            outputLog << "dTResultDisReal=" << dTResultDisReal << " (dTResultReal=" << dTResultReal << ")" << endl;
            outputLog << "TResultDisReal[j][]={" ; for (int u=0; u<int(TResultDisReal[j].size()); u++) {if (u<int(TResultDisReal[j].size())-1) {outputLog << TResultDisReal[j][u] << ", ";} else {outputLog << TResultDisReal[j][u];};}; outputLog << "}" << endl;
        }; // closes Plot condition
        // We do the same for the virtual results
        
        int dTResultDisVirtual=0;
        if (Exploration==4) {
            double dTResultPercentile=0; int dTResultVecSize=int(dTResultVec[j].size());
            for (int i=0; i<dTResultVecSize; i++) { // sorting the vector SmoothVec by ascending values
                if (dTResultVirtual<dTResultVec[j][i]) {dTResultVec[j].insert(dTResultVec[j].begin()+i, dTResultVirtual*(1 + 0.0001*VRan[28])); dTResultPercentile=i*1.0/dTResultVecSize; break;};
                if (i==dTResultVecSize-1) {dTResultVec[j].push_back(dTResultVirtual*(1 + 0.0001*VRan[28])); dTResultPercentile=1; break;};
            }; // closes for loop
            if (dTResultPercentile<0.05) {dTResultDisVirtual=-4;}
            else if ((dTResultPercentile>=0.05) && (dTResultPercentile<0.25)) {dTResultDisVirtual=-2;}
            else if ((dTResultPercentile>=0.25) && (dTResultPercentile<0.5)) {dTResultDisVirtual=-1;}
            else if ((dTResultPercentile>=0.5) && (dTResultPercentile<0.75)) {dTResultDisVirtual=1;}
            else if ((dTResultPercentile>=0.75) && (dTResultPercentile<0.95)) {dTResultDisVirtual=2;}
            else if (dTResultPercentile>=0.95) {dTResultDisVirtual=4;};
            if (dTResultVecSize==1) {dTResultDisVirtual=1;}; // In the beginning it means nothing
            if (dTResultVecSize - Future - History >0) {dTResultVec[j].erase (dTResultVec[j].begin());}; // OOO
            TResultDisVirtual[j].push_back(dTResultDisVirtual); // LLL
            if (Plot=="On") {
                outputLog << "dTResultDisVirtual=" << dTResultDisVirtual << " (dTResultVirtual=" << dTResultVirtual << ")" << endl;
                outputLog << "TResultDisVirtual[j][]={" ; for (int u=0; u<int(TResultDisVirtual[j].size()); u++) {if (u<int(TResultDisVirtual[j].size())-1) {outputLog << TResultDisVirtual[j][u] << ", ";} else {outputLog << TResultDisVirtual[j][u];};}; outputLog << "}" << endl;
            }; // closes plot condition
        }; // closes if
        
        double dTResultdis5=4;
        if (OffPolicy2==1) {
            if (Plot=="On") {outputLog << "dTResultdis5=" << dTResultdis5 << " (always)" << endl;};
        };
        
        
        /*************************************ABM II.7******************************************/
        if (Plot=="On") {outputLog << endl << "    ABM 2.7: Updating Q(s,a) based on these returns via SAM" << endl;};
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.6 computed..." << endl;
        
        // We now update the action-value function $Q(s,a)$ based on this return via the SAM or Sample Average Method
        vector<int> OldStates = get_tensor_coord (TSIndex[j][max(0,int(TSIndex[j].size()) - Future - 1)], 5, TSDim); // Coord in s-tensor of past state $s_{t-T}$
        vector<int> OldRealActions = get_tensor_coord (TAIndexReal[j][max(0,int(TAIndexReal[j].size()) - Future - 1)], 2, TADim); //Coord in a-tensor of past action $a_{t-T}$
        vector<int> OldVirtualActions, OldActions5;
        dTResultDisReal=TResultDisReal[j][max(0,int(TResultDisReal[j].size()) - Future - 1)]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
        
        if (Exploration==4) {
            OldVirtualActions = get_tensor_coord (TAIndexVirtual[j][max(0,int(TAIndexVirtual[j].size()) - Future - 1)], 2, TADim); //Coord in a-tensor of past action
            dTResultDisVirtual=TResultDisVirtual[j][max(0,int(TResultDisVirtual[j].size()) - Future - 1)]; // Past value t-Future time steps ago that will be used to update both Q() and pi() // TTT1
        };
        
        if (OffPolicy2==1) {OldActions5=get_tensor_coord (TAid5, 2, TADim);};
        
        if (Plot=="On") {
            outputLog << "OldStates[]={" ; for (int u=0; u<int(OldStates.size()); u++) {if (u<int(OldStates.size())-1) {outputLog << OldStates[u] << ", ";} else {outputLog << OldStates[u];};}; outputLog << "}" << endl;
            outputLog << "OldRealActions[]={" ; for (int u=0; u<int(OldRealActions.size()); u++) {if (u<int(OldRealActions.size())-1) {outputLog << OldRealActions[u] << ", ";} else {outputLog << OldRealActions[u];};}; outputLog << "}" << endl;
            if (Exploration==4) {outputLog << "OldVirtualActions[]={" ; for (int u=0; u<int(OldVirtualActions.size()); u++) {if (u<int(OldVirtualActions.size())-1) {outputLog << OldVirtualActions[u] << ", ";} else {outputLog << OldVirtualActions[u];};}; outputLog << "}" << endl;};
            outputLog << "We are at t=Future=" << t << endl;
        }; // closes Plot condition
        TQIndReal[0]=OldStates[0]; TQIndReal[1]=OldStates[1]; TQIndReal[2]=OldStates[2]; TQIndReal[3]=OldStates[3]; TQIndReal[4]=OldStates[4];//s
        TQIndReal[5]=OldRealActions[0]; TQIndReal[6]=OldRealActions[1]; // a
        if (Exploration==4) {TQIndVirtual[0]=OldStates[0]; TQIndVirtual[1]=OldStates[1]; TQIndVirtual[2]=OldStates[2]; TQIndVirtual[3]=OldStates[3]; TQIndVirtual[4]=OldStates[4]; // s
            TQIndVirtual[5]=OldVirtualActions[0]; TQIndVirtual[6]=OldVirtualActions[1];}; // a
        if (OffPolicy2==1) {TQInd5[0]=OldStates[0]; TQInd5[1]=OldStates[1]; TQInd5[2]=OldStates[2]; TQInd5[3]=OldStates[3]; TQInd5[4]=OldStates[4]; // s
            TQInd5[5]=OldActions5[0]; TQInd5[6]=OldActions5[1];}; // a
        int QRealTensorIndex = get_tensor_index (7, TQIndReal, TQDim); // Vector index of Q-tensor <=> $s_{t-T}$ x $a_{t-T}$
        int QVirtualTensorIndex=0; int QTensorIndex5=0;
        if (Exploration==4) {QVirtualTensorIndex = get_tensor_index (7, TQIndVirtual, TQDim);}; // Vector index of Q-tensor <=> $s_{t-T}$ x $a_{t-T}$
        if (OffPolicy2==1) {QTensorIndex5=get_tensor_index (7, TQInd5, TQDim);};
        if (Plot=="On") {
            outputLog << "TQIndReal[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQIndReal[u] << ", ";} else {outputLog << TQIndReal[u];};}; outputLog << "}" << endl;
            outputLog << "TQDim[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQDim[u] << ", ";} else {outputLog << TQDim[u];};}; outputLog << "}" << endl;
            outputLog << "QRealTensorIndex=" << QRealTensorIndex << endl;
            if (Exploration==4) {outputLog << "TQIndVirtual[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQIndVirtual[u] << ", ";} else {outputLog << TQIndVirtual[u];};}; outputLog << "}" << endl;
                outputLog << "QVirtualTensorIndex=" << QVirtualTensorIndex << endl;};
            if (OffPolicy2==1) {outputLog << "TQInd5[]={" ; for (int u=0; u<7; u++) {if (u<7-1) {outputLog << TQInd5[u] << ", ";} else {outputLog << TQInd5[u];};}; outputLog << "}" << endl;
                outputLog << "QTensorIndex5=" << QTensorIndex5 << endl;};
        }; // closes Plot condition
        
        if (CondRecord==1) {
            TNumberA[j][QRealTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
            if (TNumberA[j][QRealTensorIndex]<=1) {TQ[j][QRealTensorIndex]=dTResultReal;} // No SAM if not yet populated (i.e =0)
            else {TQ[j][QRealTensorIndex]=(TQ[j][QRealTensorIndex]*(TNumberA[j][QRealTensorIndex]-1) + dTResultReal)/TNumberA[j][QRealTensorIndex];}; // SAM
            if (Plot=="On") {
                outputLog << "Number of times this real action was previously taken in that former state: TNumberA[j][QRealTensorIndex]-1=" << TNumberA[j][QRealTensorIndex]-1 << ", so TQ[j][QRealTensorIndex]=" << TQ[j][QRealTensorIndex] << endl;
            }; // closes Plot condition
        };
        
        if (Exploration==4) {
            TNumberA[j][QVirtualTensorIndex]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
            if (TNumberA[j][QVirtualTensorIndex]<=1) {TQ[j][QVirtualTensorIndex]=dTResultVirtual;} // No SAM if not yet populated
            else {TQ[j][QVirtualTensorIndex]=(TQ[j][QVirtualTensorIndex]*(TNumberA[j][QVirtualTensorIndex]-1) + dTResultVirtual)/TNumberA[j][QVirtualTensorIndex];}; // SAM
            if (Plot=="On") {
                outputLog << "Number of times this virtual action was previously taken in that former state: TNumberA[j][QVirtualTensorIndex]-1=" << TNumberA[j][QVirtualTensorIndex]-1 << ", so TQ[j][QVirtualTensorIndex]=" << TQ[j][QVirtualTensorIndex] << endl;
            }; // closes Plot condition
        }; // closes if
        
        if (OffPolicy2==1) {
            TNumberA[j][QTensorIndex5]+=1; // Adding this action as taken in the ActionNumber tensor so as to count for SAM
            if (TNumberA[j][QTensorIndex5]<=1) {TQ[j][QTensorIndex5]=dTResult5;} // No SAM if not yet populated
            else {TQ[j][QTensorIndex5]=(TQ[j][QTensorIndex5]*(TNumberA[j][QTensorIndex5]-1) + dTResult5)/TNumberA[j][QTensorIndex5];}; // SAM
            if (Plot=="On") {
                outputLog << "Number of times this off-policy action was previously taken in that former state: TNumberA[j][QTensorIndex5]-1=" << TNumberA[j][QTensorIndex5]-1 << ", so TQ[j][QTensorIndex5]=" << TQ[j][QTensorIndex5] << endl;
            }; // closes Plot condition
        };
        
        
        
        
        /*************************************ABM II.8******************************************/
        if (Plot=="On") {
            outputLog << endl << "    ABM 2.8: Updating existing policy pi(s,a) based on these discretized results" << endl;
            outputLog << "TA=" << TA << ", TS=" << TS << endl;
        };
        // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.7 computed..." << endl;
        
        // We now update existing policy pi(s,a) based on these discretized results
        TSInd[0]=OldStates[0]; TSInd[1]=OldStates[1]; TSInd[2]=OldStates[2]; TSInd[3]=OldStates[3]; TSInd[4]=OldStates[4];
        int OldSTensorIndex = get_tensor_index (5, TSInd, TSDim); // Vector index of the s-tensor associated with $s_{t-T}$
        if (Plot=="On") {
            outputLog << "State at time t-Future: TSInd[0]=OldStates[0]=" << TSInd[0] << ", TSInd[1]=OldStates[1]=" << TSInd[1] << ", TSInd[2]=OldStates[2]=" << TSInd[2] << ", TSInd[3]=OldStates[3]=" << TSInd[3] << ", TSInd[4]=OldStates[4]=" << TSInd[4] << endl;
            outputLog << "OldSTensorIndex=" << OldSTensorIndex << endl;
        }; // closes Plot condition
        int OldSCount=0;
        for (int i=0; i<int(TSIndex[j].size())-Future; i++) {if (TSIndex[j][i]==OldSTensorIndex) {OldSCount+=1;};}; // Nb of times $s_{t-T}$ seen
        TQIndReal[0]=OldStates[0]; TQIndReal[1]=OldStates[1]; TQIndReal[2]=OldStates[2]; TQIndReal[3]=OldStates[3]; TQIndReal[4]=OldStates[4]; TQIndReal[5]=0; TQIndReal[6]=0; // Using results from above on Q for s but with zero action
        if (Plot=="On") {
            outputLog << "Number of times state at t-Future encountered: OldSCount=" << OldSCount << endl;
            outputLog << "TQIndReal[0]=OldStates[0]=" << TQIndReal[0] << ", TQIndReal[1]=OldStates[1]=" << TQIndReal[1] << ", TQIndReal[2]=OldStates[2]=" << TQIndReal[2] << ", TQIndReal[3]=OldStates[3]=" << TQIndReal[3] << ", TQIndReal[4]=OldStates[4]=" << TQIndReal[4] << ", TQIndReal[5]=" << 0 << ", TQIndReal[6]=" << 0 << endl;
        }; // closes Plot condition
        int PiSTensorIndex = get_tensor_index (7, TQIndReal, TQDim); // Vector index of the Q-tensor associated with $s_{t-T}$
        if (Plot=="On") {
            outputLog << "PiSTensorIndex=" << PiSTensorIndex << " from which we span all actions to update their probability" << endl;
            //outputLog << "State s was encountered OldSCount=" << OldSCount << " times" << endl;
            outputLog << "Results where dTResultDisReal=" << dTResultDisReal << ", but with CondRecord=" << CondRecord << " => 0: not taken into account, 1: taken into account" << endl;
            if (Exploration==4) {outputLog << "And dTResultDisVirtual=" << dTResultDisVirtual << endl;};
            if (OffPolicy2==1) {outputLog << "And dTResultdis5=" << dTResultdis5 << endl;};
        }; // closes Plot condition
        double Sum=0;
        for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
            //if (Plot=="On") {outputLog << "Old TPi[j][" << i << "]=" << floor(100*TPi[j][i]) << "%" << endl;};
            Sum+=TPi[j][i];
        };
        if (Plot=="On") {
            outputLog << "Sum=" << 100*Sum << endl;
            if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorOriginalT[]" << endl;};
        };
        Sum=0;
        
        if (CondRecord==1) { // JJJ3
            // First for real
            for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                if (Plot=="On") {outputLog << "Old real TPi[j][" << i << "]=" << 100*TPi[j][i] << "% => ";};
                if ((i==QRealTensorIndex) && (dTResultDisReal>=0)) {
                    for (int f=0; f<dTResultDisReal; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]) + NEBLearningRate/TNumberA[j][QRealTensorIndex];};
                } // closes if
                else if ((i!=QRealTensorIndex) && (dTResultDisReal>=0)) {
                    for (int f=0; f<dTResultDisReal; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]);};
                } // closes else if
                else if ((i==QRealTensorIndex) && (dTResultDisReal<0)) {
                    for (int f=0; f<abs(dTResultDisReal); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]);};
                } // closes else if
                else if ((i!=QRealTensorIndex) && (dTResultDisReal<0)) {
                    for (int f=0; f<abs(dTResultDisReal); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QRealTensorIndex]) + NEBLearningRate/(TA-1)/TNumberA[j][QRealTensorIndex];};
                }; // closes else if
                if (Plot=="On") {outputLog << "New real TPi[j][" << i << "]=" << 100*TPi[j][i] << "%" << endl;};
                Sum+=TPi[j][i];
            }; // closes i loop
            if (Plot=="On") {
                outputLog << "(QRealTensorIndex=" << QRealTensorIndex << " was the real action update), TNumberA[j][QRealTensorIndex]=" << TNumberA[j][QRealTensorIndex] << endl;
                outputLog << "Sum=" << 100*Sum << endl;
                if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorRealT[]: For OldSCount=" << OldSCount << ", dTResultDisReal=" << dTResultDisReal << endl;};
            }; // closes Plot condition
            Sum=0;
        }; // closes if
        // Second for virtual
        //if ((Exploration==4) && (1-NEB8>VRan[19])) { // This is for off-policy Watkin's Q-learning
        if (Exploration==4) { // This is for off-policy Watkin's Q-learning
            for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                if (Plot=="On") {outputLog << "Old virtual TPi[j][" << i << "]=" << 100*TPi[j][i] << "% => ";};
                if ((i==QVirtualTensorIndex) && (dTResultDisVirtual>=0)) {
                    for (int f=0; f<dTResultDisVirtual; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/TNumberA[j][QVirtualTensorIndex];};
                } // closes if
                else if ((i!=QVirtualTensorIndex) && (dTResultDisVirtual>=0)) {
                    for (int f=0; f<dTResultDisVirtual; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]);};
                } // closes else if
                else if ((i==QVirtualTensorIndex) && (dTResultDisVirtual<0)) {
                    for (int f=0; f<abs(dTResultDisVirtual); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]);};
                } // closes else if
                else if ((i!=QVirtualTensorIndex) && (dTResultDisVirtual<0)) {
                    for (int f=0; f<abs(dTResultDisVirtual); f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QVirtualTensorIndex]) + NEBLearningRate/(TA-1)/TNumberA[j][QVirtualTensorIndex];};
                }; // closes else if
                if (Plot=="On") {outputLog << "New virtual TPi[j][" << i << "]=" << 100*TPi[j][i] << "%" << endl;};
                //if (TPi[i]<0) {outputLog << "ErrorVirtual2a: at i=" << i << ", TPi[i]=" << TPi[i] << endl;};
                Sum+=TPi[j][i];
            }; // closes i loop
            if (Plot=="On") {
                outputLog << "(QVirtualTensorIndex=" << QVirtualTensorIndex << " was the virtual action update), TNumberA[j][QVirtualTensorIndex]=" << TNumberA[j][QVirtualTensorIndex] << endl;
                outputLog << "Sum=" << 100*Sum << endl;
                if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorVirtualT[]: For OldSCount=" << OldSCount << ", dTResultDisVirtual=" << dTResultDisVirtual << endl;};
            }; // closes Plot condition
            Sum=0;
        }; // closes if Exploration
        
        
        // Second for Exploration5
        if (OffPolicy2==1) {
            for (int i=PiSTensorIndex; i<PiSTensorIndex+TA; i++) {
                if (Plot=="On") {outputLog << "Old Exploration5 TPi[j][" << i << "]=" << 100*TPi[j][i] << "% => ";};
                if (i==QTensorIndex5) {
                    for (int f=0; f<dTResultdis5; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QTensorIndex5]) + NEBLearningRate/TNumberA[j][QTensorIndex5];};
                } // closes if
                else if (i!=QTensorIndex5) {
                    for (int f=0; f<dTResultdis5; f++) {TPi[j][i]=TPi[j][i]*(1-NEBLearningRate/TNumberA[j][QTensorIndex5]);};
                } // closes else if
                if (Plot=="On") {outputLog << "New Exploration5 TPi[j][" << i << "]=" << 100*TPi[j][i] << "%" << endl;};
                //if (TPi[i]<0) {outputLog << "ErrorVirtual2a: at i=" << i << ", TPi[i]=" << TPi[i] << endl;};
                Sum+=TPi[j][i];
            }; // closes i loop
            if (Plot=="On") {
                outputLog << "(QTensorIndex5=" << QTensorIndex5 << " was the off-policy action update), TNumberA[j][QTensorIndex5]=" << TNumberA[j][QTensorIndex5] << endl;
                outputLog << "Sum=" << 100*Sum << endl;
                if ((int(100*Sum)<99) || (int(100*Sum)>101)) {outputLog << "ErrorOffPolicyT[]: For OldSCount=" << OldSCount << ", dTResultdis5=" << dTResultdis5 << endl;};
            }; // closes Plot condition
            Sum=0;
        }; // closes if Exploration
        
        
        
        
    }; // closes if (t-Future>0)
    
    if (t==Future-1) { // Cancel all previous learnings (both Q(s,a) and Pi(s,a)) as they cannot be considered valid (either in state s which is represented according to Past*Future, or forecast which are for Lag*Future with Lag=1,2,3)
        for (int i=0; i<FS*FA; i++) {FQ[j][i]=0; FNumberA[j][i]=0; FPi[j][i]=1.0/FA;};
        for (int i=0; i<TS*TA; i++) {TQ[j][i]=0; TNumberA[j][i]=0; TPi[j][i]=1.0/TA;};
    };
    if (t<Future-1) {QTradeReal=0;}; // So agent doesn't trade before time
    if (Plot=="On") {outputLog << endl << endl;};
    // outputDebug << "t=" << t << ", i=" << AgentName << ", j=" << j << ": ABM 2.8 computed..." << endl;
    
    
    
}; // closes method RL()