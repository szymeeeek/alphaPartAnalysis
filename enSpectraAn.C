#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>
using namespace std;

Bool_t Widmo(string filename) {
    //energy calib and histo definition
    string outFile = "fitParams_absCounts.txt";
    string params[3] = {"amplitude", "mean", "sigma"};
    const Int_t nBins = 1024;

    double minChannel = 0, maxChannel = nBins-1;

    double m = 7.917;
    double b = -6e-14;

    double energyBins[nBins + 1];
    for (int i = 0; i <= nBins; i++) {
        energyBins[i] = m * (minChannel + i) + b;
    }

    TH1F *hist = new TH1F(Form("Widmo%s", filename.c_str()), Form("Widmo %s", filename.c_str()), nBins, energyBins);
    TF1 *gaus = new TF1("gaus", "gaus", m*minChannel+b, m*maxChannel+b);
    gaus->SetParameters(500, 3000, 100);
    TFile *file = new TFile("integralsTest.root", "UPDATE");

    //file opening and reading it
    fstream myfile(filename, ios::in);
    if (!myfile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return kFALSE;
    }

    Int_t counts[nBins] = {};
    Double_t errors[nBins] = {};

    Int_t i = 0;
    Double_t dummy = 0;
    while (myfile >> counts[i]) {
        errors[i] = sqrt(counts[i]);
        i++;
    }
    myfile.close();

    //histo filling, fitting and formatting
    for (Int_t i = 0; i<nBins; i++) {
        hist->SetBinContent(i, counts[i]);   
        hist->SetBinError(i, errors[i]);
    }

    fstream myfile1;
    myfile1.open(outFile, ios::out | ios::app);
    TFitResultPtr results = hist->Fit(gaus, "");
    //     cout << "Chi2 = " << results->Chi2() << endl;
    //     cout << "NDF = " << results->Ndf() << endl;
    //     cout << "EDM = " << results->Edm() << endl;

    Double_t absCounts = hist->Integral(0, 7500);

    myfile1<<filename<<" INT: "<<absCounts<<endl;
    cout<<"-----------------------"<<absCounts<<endl;

    for (int j=0;j<gaus->GetNpar();j++) {
        Float_t value = gaus->GetParameter(j);
        Float_t valEr = gaus->GetParError(j);
        myfile1<<filename<<" "<<params[j]<<" = "<<value<<" "<<params[j]<<" error = "<<valEr<<endl;
    }

    TCanvas *c1 = new TCanvas(Form("%s",filename.c_str()), "Widmo", 800, 600);
    c1->cd();
    hist->Draw();
    hist->GetXaxis()->SetTitle("energy [keV]");
    hist->GetYaxis()->SetTitle("counts [a.u.]");
    hist->SetStats(kFALSE);

    hist->Write();
    c1->Write();
    file->Write();
    file->Close();
    myfile1.close();

    return kTRUE;
}


Bool_t enSpectraAn(){
    vector<string> x = {"15.dat", "45.dat", "75.dat", "105.dat", "135.dat", "165.dat", "195.dat", "225.dat"};
    vector<string> p = {"0022.dat", "00115.dat", "0071.dat", "0091.dat", "003.dat", "0037.dat", "0098.dat", "006.dat", "0048.dat"};

    for(int i = 0; i<x.size(); i++){
        Widmo(x.at(i));
    }

    for(int i = 0; i<p.size(); i++){
        Widmo(p.at(i));
    }

    return kTRUE;
}
    
