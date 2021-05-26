#pragma once
#include <vector>
#include <filesystem>
#include "agent.hpp"

using namespace std;


/// A line in the order book.
struct Order {
    bool Nature;  ///< 0 for Bid and 1 for Ask
    int Agent;    ///< Agent placing the order
    int Stock;    ///< Stock to trade
    int Q;        ///< Quantity of the stock to trade
    double Price; ///< Price or less for bids, Price or more for asks
};


/// Order book for a given stock
class StockOrderBook {
public:
    vector<Order> Bids;     ///< Collection of all bids
    vector<Order> Asks;     ///< Collection of all asks
    int MetaorderLastAgent; ///< Agent placing the order
    int MetaorderNature;    ///< 0 for Bid and 1 for Ask
    double Credit;          ///< Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
    
    void Sort();
    vector<double> Clear(vector<Agent>& Market, int t, int j, string Plot, const filesystem::path& output_dir);
    int MetaorderInjection(vector<Agent> &Market, int SharesOutstanding, int MetaorderImpact, int OBLevelSize,
        const filesystem::path& output_dir);
    
}; // closes class StockOrderBook
