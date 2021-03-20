#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "global.hpp"
#include "order_book.hpp"
#include "agent.hpp"
#include "utility.hpp"

using namespace std;


/// Sorting the order book for the given stock, with Bids of decreasing prices and Asks of increasing prices (SS9)
void StockOrderBook::Sort() {
    vector<Order> SortedBids, SortedAsks, SortedBidsNoId, SortedAsksNoId;
    vector<double> SortedListBids, SortedListAsks;
    // Saving the prices only of these orders Bids and Asks into two separate vectors SortedListBids and SortedListAsks
    for (int k=0; k<int(Bids.size()); k++) {SortedListBids.push_back(Bids[k].Price);}; // closes k loop
    for (int k=0; k<int(Asks.size()); k++) {SortedListAsks.push_back(Asks[k].Price);}; // closes k loop
    // Sorting them
    sort(SortedListBids.begin(), SortedListBids.begin() + Bids.size()); // sort() is by ascending order, bad for Bids
    reverse(SortedListBids.begin(), SortedListBids.begin() + Bids.size()); // reverse() puts it by descending order, good for Bids
    sort(SortedListAsks.begin(), SortedListAsks.begin() + Asks.size()); // sort() is by ascending order, good for Asks
    // Sorting the orders Bids and Asks according to these two sorted lists of prices SortedListBids and SortedListAsks, and saving them in new sorted orders SortedBids and SortedAsks
    for (int k=0; k<int(SortedListBids.size()); k++) {
        for (int m=0; m<int(Bids.size()); m++) {
            if (Bids[m].Price==SortedListBids[k]) {
                SortedBids.push_back(Bids[m]);
                Bids.erase(Bids.begin()+m); // To ensure this order is not taken again in the next count if another bid has the same price
                m=int(SortedListBids.size());
            }; // closes if
        }; // closes m loop
    }; // closes k loop
    for (int k=0; k<int(SortedListAsks.size()); k++) {
        for (int m=0; m<int(Asks.size()); m++) {
            if (Asks[m].Price==SortedListAsks[k]) {
                SortedAsks.push_back(Asks[m]);
                Asks.erase(Asks.begin()+m); // To ensure this order is not taken again in the next count if another ask has the same price
                m=int(SortedListAsks.size());
            }; // closes if
        }; // closes m loop
    }; // closes k loop
    Bids=SortedBids;
    Asks=SortedAsks;
}; // closes method Sort()


vector<double> StockOrderBook::Clear(vector<Agent>& Market, int t, int j, string Plot) {
    // Passing by ref to modify Market (SS11)
    double TempClearingPrice=0;
    double TotalStockQuantityTraded=0;
    int Count=0;
    double TempClearingPriceLast=0;
    double TotalStockQuantityTradedLast=0;
    int CountLast=0;
    int LevelCount=0;
    double BrokerFees=0.01*0.01; // 0.01% on the LSE
    //double AverageBids=0;
    //double AverageAsks=0;
    //ofstream outputTest("/Users/admin/Documents/GNT/SYMBA/ZULU.txt", ofstream::app); int NumberOfAgents=30; int NumberOfStocks=2;
    ofstream outputLog(Machine+"SimLog.txt", ofstream::app);
    if (Plot=="On") {outputLog << endl << endl << "***************CLEAR() AT TIME t=" << t << "**************" << endl;};
    for (int k=-1; k<int(min(Bids.size(), Asks.size())); k++) {
        if (int(min(Bids.size(), Asks.size()))==0) {/*outputLog << "Function Clear() not activated" << endl;*/ break;}
        if (k==-1) {k=0;};
        //if (Plot=="On") {outputLog << "New loop: k=" << k << ", Bids.size()=" << Bids.size() << ", Asks.size()=" << Asks.size() << ", min(Bids.size(), Asks.size())=" << min(Bids.size(), Asks.size()) << endl;};
        int TempQuantity=int(min (Bids[k].Q, Asks[k].Q)); // KKK
        //int TempQuantity=max(0, int(min (Bids[k].Q, Asks[k].Q))); // KKK
        TempClearingPrice = (Bids[k].Price + Asks[k].Price)/2; // Transaction taken at midprice, a neutral and common solution
        bool C1 = Bids[k].Price>=Asks[k].Price; // Bid is larger than Ask at a given level k of the order book
        bool C2 = Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity > 0; // Seller has non-zero quantity of this stock
        bool C3 = Market[Bids[k].Agent].RFA >= TempClearingPrice * TempQuantity; // Buyer has enough RFA to buy
        bool C4 = Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity >= TempQuantity; // Seller has enough quantities to sell
        if (C1 && C2 && C3 && C4) {
            //int QBidsBefore=Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity;
            //int QAsksBefore=Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity;
            // This is the OB clearing per se
            Market[Bids[k].Agent].RFA -= TempClearingPrice * TempQuantity;
            Market[Asks[k].Agent].RFA += TempClearingPrice * TempQuantity;
            Market[Bids[k].Agent].RFA*=1-BrokerFees; // Broker fees
            Market[Asks[k].Agent].RFA*=1-BrokerFees; // Broker fees
            Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity += TempQuantity;
            Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity -= TempQuantity;
            Market[Bids[k].Agent].TradingWindowClock.push_back(Market[Bids[k].Agent].Future+t); // JJJ3
            Market[Asks[k].Agent].TradingWindowClock.push_back(Market[Asks[k].Agent].Future+t); // JJJ3
            TotalStockQuantityTraded+=TempQuantity;
            Count+=1;
            //for (int i=0; i<NumberOfAgents ; i++) {for (int j=0; j<NumberOfStocks; j++) {if (Market[i].Stocks[j].StockQuantity <0) {outputTest << "Error Z Clear(): TempQuantity=min (Bids[k].Q, Asks[k].Q)=" << TempQuantity << ", Bids[k].Q=" << Bids[k].Q << ", Asks[k].Q=" << Asks[k].Q << ", Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity=" << QBidsBefore << "->" << Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity << ", Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity=" << QAsksBefore << "->" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity << endl;};};};
            // This is for method RL()
            int TransactionPriceRealSize = int(Market[Bids[k].Agent].TransactionPriceReal[j].size());
            int TransactionPriceRealSize2 = int(Market[Asks[k].Agent].TransactionPriceReal[j].size());
            int QuantitiesRealSize = int(Market[Bids[k].Agent].QuantitiesReal[j].size());
            int QuantitiesRealSize2 = int(Market[Asks[k].Agent].QuantitiesReal[j].size());
            Market[Bids[k].Agent].TransactionPriceReal[j][TransactionPriceRealSize-1]=TempClearingPrice;
            Market[Asks[k].Agent].TransactionPriceReal[j][TransactionPriceRealSize2-1]=TempClearingPrice;
            Market[Bids[k].Agent].QuantitiesReal[j][QuantitiesRealSize-1]=TempQuantity;
            Market[Asks[k].Agent].QuantitiesReal[j][QuantitiesRealSize2-1]=-TempQuantity;
            
            // This computes the average of all cleared Bids and all cleared Asks
            //AverageBids+=Bids[k].Price;
            //AverageAsks+=Asks[k].Price;
            
            TempClearingPriceLast=TempClearingPrice;
            TotalStockQuantityTradedLast=TotalStockQuantityTraded;
            CountLast=Count;
            
            // Output
            if (Plot=="On") {
                outputLog << "Level " << LevelCount << ": Agent " << Bids[k].Agent << " with Q=" << Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity - TempQuantity << "->" << Market[Bids[k].Agent].Stocks[Bids[k].Stock].StockQuantity << " longs " << TempQuantity << " of Stock " << Bids[k].Stock << " from Agent " << Asks[k].Agent << " with Q=" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity+TempQuantity << "->" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity << " for $" << TempClearingPrice << " (midprice from " << Asks[k].Price << "~" << Bids[k].Price << ")" << endl;
            }; // closes if
            
            // Proper draft of all orders
            if (Bids[k].Q>TempQuantity) {Bids[k].Q-=TempQuantity; Asks.erase(Asks.begin()+k); k-=1;} // Buyer stays in business if still eager to buy more from next seller
            else if (Asks[k].Q>TempQuantity) {Asks[k].Q-=TempQuantity; Bids.erase(Bids.begin()+k); k-=1;}; // Seller stays in business if still eager to sell more to next buyer
            //if (Plot=="On") {outputLog << "Now k=" << k << ", Bids.size()=" << Bids.size() << ", Asks.size()=" << Asks.size() << ", min(Bids.size(), Asks.size())=" << min(Bids.size(), Asks.size()) << endl;};
        } // closes if
        
        else if (!C1) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since bid=" << Bids[k].Price << " < ask=" << Asks[k].Price << endl;};}
        else if (!C2) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since seller out of stock" << endl;};}
        else if (!C3) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since buyer hasn't enough RFA to buy ($" << Market[Bids[k].Agent].RFA << ")" << endl;};}
        else if (!C4) {if (Plot=="On") {outputLog << "Level " << LevelCount << ": not cleared since seller hasn't enough quantities to sell (" << Market[Asks[k].Agent].Stocks[Asks[k].Stock].StockQuantity << ")" << endl;};};
        
        /*
            if ((k<0) || (OBSize<=0)) {
            //if (Plot=="On") {outputLog << "We are in a case (k<0) || (OBSize<=0)" << endl;};
            k=0; OBSize=0;
            //if (Plot=="On") {outputLog << "So we change k=" << k << ", OBSize=" << OBSize << endl;};
            }; // closes if
            */
        LevelCount+=1;
    }; // closes k loop
    
    // Outputting clearing price
    //if ((Count!=0) && ((Plot=="On") || (Plot=="Mini"))) {outputLog << endl << "CLEARING PRICE FOR STOCK " << OBStock << " AT TIME t=" << t << ": $" << TempClearingPrice << endl<< endl << endl;};
    
    // Returning result
    //if (Count>0) {AverageBids/=Count; AverageAsks/=Count;};
    vector<double> Res;
    Res.push_back(TempClearingPriceLast); // Market price given by the latest (lowest) clearing price
    Res.push_back(TotalStockQuantityTradedLast); // Total stock quantity that was traded at this time step t
    Res.push_back(CountLast); // Indicator of emptyness of OB
    //ofstream outputOBLevels(Machine+"OBlevels.csv", ofstream::app);
    //outputOBLevels << CountLast << endl;
    //Res.push_back(AverageBids); // Average of all Bids
    //Res.push_back(AverageAsks); // Average of all Asks
    if (Plot=="On") {outputLog << endl << "TempClearingPriceLast=" << Res[0] << ", TotalStockQuantityTradedLast=" << Res[1] << ", CountLast=" << Res[2] << endl;};
    outputLog.close();
    //outputOBLevels.close();
    return Res;
}; // closes method Clear()


// Selecting a random level in the OB and shifting it to a metaorder (if at least one business order in the OB)
int StockOrderBook::MetaorderInjection(vector<Agent> &Market, int SharesOutstanding, int MetaorderImpact, int OBLevelSize) {
    ofstream outputMetalog(Machine+"MetaorderLog.txt", ofstream::app);
    int Res=-1;
    if (OBLevelSize>-1) {
        vector<int> J = Shuffle(OBLevelSize+1);
        int MetaorderLevel=J[0]; // OB level where the metaorder will be injected (whether in the ask or the bid order)
        int BidOrAskMetaorder=int(time(0))%2; // If 0 then bid, if 1 then ask
        if ((BidOrAskMetaorder==0) && (OBLevelSize>-1)) {
            Bids[MetaorderLevel].Q=int(MetaorderImpact*SharesOutstanding/100);
            Market[Bids[MetaorderLevel].Agent].RFA+=ceil(Bids[MetaorderLevel].Q*(Bids[MetaorderLevel].Price+Asks[MetaorderLevel].Price)/2); // This agent is made richer in order to buy
            Res=MetaorderLevel;
            outputMetalog << "Bid metaorder injected at OB level " << MetaorderLevel << ": agent " << Bids[MetaorderLevel].Agent << " buys " << Bids[MetaorderLevel].Q << " for £" << Bids[MetaorderLevel].Price << endl;
            MetaorderLastAgent=Bids[MetaorderLevel].Agent; // Agent placing the order
            MetaorderNature=0; // 0 for Bid and 1 for Ask
            Credit=ceil(Bids[MetaorderLevel].Q*(Bids[MetaorderLevel].Price+Asks[MetaorderLevel].Price)/2); // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
        }
        else if ((BidOrAskMetaorder==1) && (OBLevelSize>-1)) {
            Asks[MetaorderLevel].Q=int(MetaorderImpact*SharesOutstanding/100);
            Market[Asks[MetaorderLevel].Agent].Stocks[0].StockQuantity+=Asks[MetaorderLevel].Q+1; // This agent is made richer in order to sell
            Res=MetaorderLevel;
            outputMetalog << "Ask metaorder injected at OB level " << MetaorderLevel << ": agent " << Asks[MetaorderLevel].Agent << " sells " << Asks[MetaorderLevel].Q << " for £" << Asks[MetaorderLevel].Price << endl;
            MetaorderLastAgent=Asks[MetaorderLevel].Agent; // Agent placing the order
            MetaorderNature=1; // 0 for Bid and 1 for Ask
            Credit=double(Asks[MetaorderLevel].Q+1); // Extra RFA (for a buy) or stock quantity (for a sell) temporarily credited to the agent
        };
    };
    outputMetalog.close();
    return Res;
}; // closes method MetaorderInjection()