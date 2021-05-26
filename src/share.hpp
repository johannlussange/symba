#pragma once
#include <string>
#include <filesystem>
#include <gsl/gsl_matrix.h>

using namespace std;


class Share {
public:
    string Symbol, Title, Exchange, File, Country, Currency, Sector, InPF;
    gsl_matrix* Data; /// Prices & cashflows
    
    void Gen(const filesystem::path& Name, gsl_matrix* Macroeconomics, int FirstDate);    
    
    // Broker fees of IG are downloaded at https://www.ig.com/fr/actions/conditions-actions
    // This outputs files LSEnew.txt, NYSEnew.txt, NASDAQnew.txt which are the stocks offered by IG. These files must be formatted to Legacy MAC OS (CR) text format, and can in turn be used as direct data feed
    void ProcessBrokerData(const filesystem::path& output_dir);

};