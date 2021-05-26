#include <vector>
#include <fstream>
#include <sstream>
#include "share.hpp"
#include "global.hpp"
#include "cli.hpp"


void Share::Gen(const filesystem::path& Name, gsl_matrix* Macroeconomics, int FirstDate) {
    //int Numeraire=5;
    //if (Currency=="USD") {Numeraire=3;} else if (Currency=="GBP") {Numeraire=4;}; // This is switched off so as to study local market microstructure apart from FX perturbations // FXX
    vector<double> LocalDates, LocalOpen, LocalClose, LocalHigh, LocalLow, LocalVolumes;
    ifstream FileInput(Name); string Line; int LineNb=0;
    while (getline (FileInput, Line)) {LineNb++;
        istringstream LineStream(Line); string Item; int ItemNb=0;
        while (getline (LineStream, Item, ',')) {ItemNb++;
            if ((LineNb>1) && (ItemNb==3)) {LocalDates.push_back(atof(Item.c_str()));};
            if ((LineNb>1) && (ItemNb==5)) {LocalOpen.push_back(atof(Item.c_str()));};
            if ((LineNb>1) && (ItemNb==8)) {LocalClose.push_back(atof(Item.c_str()));};
            if ((LineNb>1) && (ItemNb==9)) {LocalVolumes.push_back(atof(Item.c_str()));};
            if ((LineNb>1) && (ItemNb==6)) {LocalHigh.push_back(atof(Item.c_str()));};
            if ((LineNb>1) && (ItemNb==7)) {LocalLow.push_back(atof(Item.c_str()));};
        };
    };
    FileInput.close();
    // Computing the data size from FirstDate to today
    int Size=0; bool Cond=0;
    for (int i=0; i<int(Macroeconomics->size2); i++) {
        if (gsl_matrix_get (Macroeconomics, 0, i)==FirstDate) {Cond=1;};
        if (Cond==1) {Size+=1;};
    };
    // Populating gsl_matrix* Data
    int PastLag=120+Year; // Steps before FirstDate necessary for BinaryProjector() to start at FirstDate
    int FirstStep=0; // STL vector index where data matches FirstDate
    int VectorSize=int(LocalDates.size());
    Data = gsl_matrix_calloc (6, Size+PastLag);
    Cond=0;
    for (int i=0; i<Size; i++) {gsl_matrix_set(Data, 0, i+PastLag, gsl_matrix_get (Macroeconomics, 0, i));}; //Dates
    for (int i=0; i<Size; i++) {
        for (int k=0; k<VectorSize; k++) {
            if ((LocalDates[k]==gsl_matrix_get (Macroeconomics, 0, i)) && Cond==0) {FirstStep=k; Cond=1;};
            if (LocalDates[k]==gsl_matrix_get (Macroeconomics, 0, i)) {
                gsl_matrix_set(Data, 1, i+PastLag, LocalOpen[k]); //Open in €
                gsl_matrix_set(Data, 2, i+PastLag, LocalClose[k]); //Close in €
                gsl_matrix_set(Data, 3, i+PastLag, LocalVolumes[k]); // Daily volume
                gsl_matrix_set(Data, 4, i+PastLag, LocalHigh[k]); //High in €
                gsl_matrix_set(Data, 5, i+PastLag, LocalLow[k]); //Low in €
                break;
            };
        }; // closes k loop
    }; // closes i loop
    // Stiching so as to match Macroeconomics and other exchanges data
    for (int i=1; i<Size; i++) {
        if (gsl_matrix_get(Data, 1, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 1, i+PastLag, gsl_matrix_get(Data, 1, i+PastLag-1));};
        if (gsl_matrix_get(Data, 2, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 2, i+PastLag, gsl_matrix_get(Data, 2, i+PastLag-1));};
        if (gsl_matrix_get(Data, 3, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 3, i+PastLag, gsl_matrix_get(Data, 3, i+PastLag-1));};
        if (gsl_matrix_get(Data, 4, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 4, i+PastLag, gsl_matrix_get(Data, 4, i+PastLag-1));};
        if (gsl_matrix_get(Data, 5, i+PastLag)<0.0000001) {gsl_matrix_set(Data, 5, i+PastLag, gsl_matrix_get(Data, 5, i+PastLag-1));};
    };
    // Populating the region 0<t<PastLag for BinaryProjector
    if (FirstStep>PastLag) {
        for (int k=FirstStep-PastLag; k<FirstStep; k++) {
            gsl_matrix_set(Data, 0, k-FirstStep+PastLag, LocalDates[k]);
            gsl_matrix_set(Data, 1, k-FirstStep+PastLag, LocalOpen[k]);
            gsl_matrix_set(Data, 2, k-FirstStep+PastLag, LocalClose[k]);
            gsl_matrix_set(Data, 3, k-FirstStep+PastLag, LocalVolumes[k]);
            gsl_matrix_set(Data, 4, k-FirstStep+PastLag, LocalHigh[k]);
            gsl_matrix_set(Data, 5, k-FirstStep+PastLag, LocalLow[k]);
        };
    };
    
    // Stitching stock splits
    for (int t=1; t<Size+PastLag; t++) {
        double Split=1;
        if (gsl_matrix_get(Data, 2, t-1)>0.00001) {Split=gsl_matrix_get(Data, 2, t)/gsl_matrix_get(Data, 2, t-1);};
        if (((Split>=4) || (Split<=0.25)) && (Split>0)) {
            for (int i=t; i<Size+PastLag; i++) {gsl_matrix_set(Data, 2, i, gsl_matrix_get(Data, 2, i)/Split);};
        }; // if
    }; // t loop
}


// Broker fees of IG are downloaded at https://www.ig.com/fr/actions/conditions-actions
// This outputs files LSEnew.txt, NYSEnew.txt, NASDAQnew.txt which are the stocks offered by IG. These files must be formatted to Legacy MAC OS (CR) text format, and can in turn be used as direct data feed
void Share::ProcessBrokerData () {
    InPF="Out";
    vector<vector<string>> W;
    auto Path = cli::args::output_dir / "Symba/CSV/Tiered Margin_cfd.csv";
    ifstream FileInput(Path); string Line; int LineNb=0;
    vector<string> V1, V2, V3, V4;
    while (getline (FileInput, Line)) {LineNb++;
        istringstream LineStream(Line); string Item; int ItemNb=0;
        while (getline (LineStream, Item, '\t')) {ItemNb++;
            if ((ItemNb==1) && (LineNb>5)) {V1.push_back(Item.c_str());}; // Symbol
            if ((ItemNb==2) && (LineNb>5)) {V2.push_back(Item.c_str());}; // Title
            if ((ItemNb==3) && (LineNb>5)) {V3.push_back(Item.c_str());}; // Country
            if ((ItemNb==5) && (LineNb>5)) {V4.push_back(Item.c_str());}; // Can go short?
        };
    };
    // Reuters RIC: .L (London Stock Exchange), .O (NASDAQ), .N (NYSE), .P (Nyse ARCA), .PK (OTC Market Group), .OB (?), .K (New York Consolidated), .A (American Stock Exchange)
    auto Name = cli::args::output_dir / "Symba/CSV" / Exchange / File;
    string Suffix=".L";
    if (Exchange=="NYSE") {Suffix=".N";} // Some stocks will taken as .K on the New York Consolidated!
    else if (Exchange=="NASDAQ") {Suffix=".O";};
    W.push_back(V1); W.push_back(V2); W.push_back(V3); W.push_back(V4);
    ifstream FileInput2(Name); string Line2; LineNb=0;
    while (getline (FileInput2, Line2)) {LineNb++;
        istringstream LineStream(Line2); string Item; int ItemNb=0;
        while (getline (LineStream, Item, ',')) {ItemNb++;
            if ((LineNb==2) && (ItemNb==1)) {
                for (int i=0; i<int(W[0].size()); i++) {
                    if (((Item.c_str()==W[0][i]) || (Item.c_str()+Suffix==W[0][i])) && (W[2][i]==Country) && (W[3][i]=="Yes")) {
                        Symbol=W[0][i]; Title=W[1][i]; InPF="In";
                    };
                }; // closes i loop
            } // closes if
        }; // closes while
    }; // closes while
    if (InPF=="In") {
        ofstream outputLSE(cli::args::output_dir / "Symba/CSV/LSEnew.txt", ofstream::app);
        ofstream outputNYSE(cli::args::output_dir / "Symba/CSV/NYSEnew.txt", ofstream::app);
        ofstream outputNASDAQ(cli::args::output_dir / "Symba/CSV/NASDAQnew.txt", ofstream::app);
        if (Exchange=="LSE") {outputLSE << File << ',' << Symbol << ',' << Title << ',' << Exchange << ',' << Country << ',' << Currency << endl;}
        if (Exchange=="NYSE") {outputNYSE << File << ',' << Symbol << ',' << Title << ',' << Exchange << ',' << Country << ',' << Currency << endl;}
        if (Exchange=="NASDAQ") {outputNASDAQ << File << ',' << Symbol << ',' << Title << ',' << Exchange << ',' << Country << ',' << Currency << endl;}
    };
}