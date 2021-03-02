#pragma once
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include "stock.hpp"

using namespace std;


class Agent {
public:
    // GENERAL VARIABLES
    int AgentName; // Name of agent
    double RFA, RBA, RFAFirst, RBAFirst, NAVOnJanFirst; // Risk-free asset (Zero Coupon Bond, ZCB), Risk-based assets (total stock holdings)
    int Bankruptcy, LiquidationFloor, TradingWindow; // Bankruptcy parameters and trading window
    vector<Stock> Stocks; // Collection of stocks
    gsl_matrix* BiasedValues; // J*T matrix built via CointegratedWalk()
    vector<int> TradingWindowClock; // Clock counter reinitialized at each trading publication
    vector<int> ExitHorizons; // Times at which the agent must exit its position
    vector<vector<double>> Qs; // Qs[j][i] is the arg max_a Q(s,a) of state s of 5-index i for stock j for model-based RL (updated every TradingWindow)
    
    // COINTEGRATION VARIABLES
    vector<double> Leash, LeashVol, Accuracy; // One for each j stock. Accuracy is the prob. of bubble burst, on average 0.2,0.75,1,2 per year for Accuracy=13,11,10,9 resp.
    
    // BIAS VARIABLES
    double Gesture; // Commercial gesture in percent added (or subtracted!) to the Bid, subtracted (or added!) to the Ask
    double Versatility; // How often (number of time steps) on average does the agent change his weight
    double Reflexivity; // Reflexivity coefficient
    
    // RL VARIABLE
    //double DiscountFactor; // For the ARMA forecast
    int Future, History; // History is the past lag over which the agent performs comparison methods for FSIndex, FAIndexReal, FAIndexVirtual, ForecastReal, ForecastVirtual, FResultDisReal, FResultDisVirtual, AND TSIndex, TAIndexReal, QuantitiesReal, TransactionPriceReal, TAIndexVirtual, QuantitiesVirtual, TransactionPriceVirtual, TResultDisReal, TResultDisVirtual and also FResultReal, FResultVirtual, TResultReal, TResultVirtual (formerly HistoryRes).
    int Rough, Smooth, Reflexive; // Forecast states
    int Tool, Lag, Weight;// Forecast actions
    int Mu, Sig, RFALevel, RBALevel, Liquid; // Trading states
    int Quant, Pinch; // Trading actions
    int QTradeReal, QTradeVirtual; // Trading quantity to file in the OB
    int Exploration; // Epsilon-Greedy (for Exploration=1), Softmax (for Exploration=2), Pursuit (for Exploration=3), Off-Policy Watkin's Q-learning (for Exploration=4), SAM (none)
    int OffPolicy, OffPolicy2; // 0: no off-policy, 1: off-policy
    double Epsilon; // Epsilon-greedy method parameter (probability value below which exploration is performed)
    //int Update; // Action-value functions Q(s,a) and policies Pi(s,a) updating methods
    
    // Neuroeconomics biases
    int NEB; // Given the value {0,1,2,3,4,5,6} if belonging to class {no NEB, NEB124, NEB356, NEB89, NEB123456, NEB12489, NEB35689}.
    double Human, NEBLossAversion, NEBPositivity, NEBNegativity, NEBFear, NEBGreed, NEBLearningRate;
    
    // RL INDICES // DDD
    int FS, FA, TS, TA; // Number of states S and actions A
    int FSDim[3], FADim[3], FQDim[6], FPiDim[6], TSDim[5], TADim[2], TQDim[7], TPiDim[7]; // Dimensions specifications
    //int FSDim[2], FADim[2], FQDim[4], FPiDim[4], TSDim[5], TADim[2], TQDim[7], TPiDim[7]; // Dimensions specifications // DDD2
    int FSInd[3], FAIndReal[3], FAIndVirtual[3], FAInd5[3], FQIndReal[6], FQIndVirtual[6], FQInd5[6], FPiInd[6];
    int TSInd[5], TAIndReal[2], TAIndVirtual[2], TAInd5[2], TQIndReal[7], TQIndVirtual[7], TQInd5[7], TPiInd[7];
    int FTenPiId[6], TTenPiId[7]; // To select an exploratory policy
    // int FTenPiId[4], TTenPiId[7]; // To select an exploratory policy // DDD2
    
    // MMM
    vector<vector<double>> dRoughVec, dSmoothVec, dMuPosVec, dMuNegVec, dSigVec, LiquidPercentVec, dFResultVec, dTResultVec; // Past values to determine median & percentiles and discrete score via History // MMM
    vector<vector<int>> FSIndex, TSIndex, FAIndexReal, FAIndexVirtual, TAIndexReal, TAIndexVirtual; // Real and virtual actions taken (identified by its vector index) at time t // MMM
    vector<vector<double>> FPi, TPi, FQ, TQ; // Action-value function Q(s)=E[R|s,a] as a SxA matrix // MMM
    vector<vector<double>> ForecastReal, ForecastVirtual, Forecast5; // MMM
    vector<vector<int>> FNumberA, TNumberA; // The number of times a was performed in s for SAM // MMM
    vector<vector<int>> QuantitiesReal, QuantitiesVirtual; // Quantities of stock at each time step that are willed to be transacted // MMM
    vector<vector<double>> TransactionPriceReal, TransactionPriceVirtual; // Transaction prices // MMM
    
    // LOG // MMM
    vector<vector<double>> LvolLog, SvolLog; // MMM
    vector<vector<int>> RoughLog, SmoothLog, ReflexiveLog, ToolRealLog, LagRealLog, WeightRealLog, ToolVirtualLog, LagVirtualLog, WeightVirtualLog, MuLog, SigLog, RFALog, RBALog, LiquidityLog, PinchRealLog, QuantRealLog, PinchVirtualLog, QuantVirtualLog, FResultDisReal, FResultDisVirtual, TResultDisReal, TResultDisVirtual; // MMM
    
    
    
    // Returns the index of a vector<double> Tensor (decomposed from a tensor of dimensions vector<int> Dimensions) corresponding to the vector<int> Indices. This can be used to update the corresponding element in Tensor, or to simply access it
    int get_tensor_index (const int IndSize, int Indices[], int Dimensions[]);
    
    // Returns the D-dim tensor coordinates corresponding to the Index of the STL vector representation of that tensor
    vector<int> get_tensor_coord (int Index, int DimSize, int Dimensions[]);
    
    double Capital();
    
    double StockHoldings();
    
    void Liquidation(int i, vector<Agent> &Market);

    // RL algorithm for forecasting and placing order in the OB (SS6)
    void RL(
        int j, int t, double Rate, gsl_matrix* ReflexiveValues, double VSpread, double LiquidPercent, int Time,
        int NumberOfStocks, string TradingFrequencyCond, string Plot, string VersatilityCondition,
        double MarketPerformance, int TimeSinceJan1st, int LearningPhase, int LeaderAgent, int LeaderQuant,
        int LeaderPinch, int ClusterLimit, int Trunk
    );
}; // closes Agent class