#include <string>
#include <vector>
#include "betheSim.C"

Bool_t graphs(){
    //initial values
    const Int_t nGraphs = 2;
    string parameters[nGraphs] = {"x", "p"};

    //graphs and canvas defining
    TFile *file = new TFile("graphs.root", "UPDATE");
    TGraphErrors *enDep[nGraphs] = {nullptr};
    TCanvas *c1 = new TCanvas();
    //c1->Divide(1,2);
    TMultiGraph *mg = new TMultiGraph();

    //data reading
    Double_t x[8] = {22.5, 19.5, 16.5, 13.5, 10.5, 7.5, 4.5, 1.5}; //mm
    Double_t xEn[8] = {3407.66, 3612.94, 3793.33, 3959.86, 4120.72, 4286.72, 4449.63, 4589.51};
    Double_t p[8] = {-0.098, -0.091, -0.071, -0.06, -0.048, -0.037, -0.03, -0.022}; //MPa
    Double_t pEn[8] = {5380.97, 5028.96, 4099.4, 3653.2, 3025.82, 2282.05, 1701.02, 971.534};
    Double_t pInt[8] = {15875 / 28.37, 17972 / 32.2, 17257 / 31.27, 14926 / 27.27, 17563 / 32.68, 19874 / 37.78, 23531 / 46.63, 22025 / 49.3};
    Double_t dp[8];
    Double_t dD[8];
    Double_t pEnx[8];
    Double_t dx[8];
    Double_t xEnx[8];

    Double_t p0 = 0.101325;
    Double_t t1 = (40/p0*0.01);
    
    for(int s = 0; s<8; s++){
        p[s] = ((101325+p[s]*1e6)*40)/101325;
        dp[s] = TMath::Sqrt((t1*t1)+(((p0+p[s])/p0)*0.001)*(((p0+p[s])/p0)*0.1));
        cout<<p[s]<<"+/-"<<dp[s]<<endl;
    }
    for(int s = 0; s<7; s++){
        dD[s] = p[s+1]-p[s];
        cout<<dD[s]<<endl;
    }
    dD[7] = 1;
    

    //filling graphs
    
    enDep[0] = new TGraphErrors(8, p, pInt, dp, 0);
    c1->cd();
    // gPad->SetLogx();
    // gPad->SetLogy();
    // enDep[0]->SetTitle(" ; d (mm); N (1/s)");
    // enDep[0]->SetMarkerStyle(8);
    //enDep[0]->Draw("ap");
    enDep[0]->Write();
    for(int s = 0; s<8; s++){
        xEnx[s] = -1*xEn[s]/3+1600;
        cout<<xEnx[s]<<endl;
        dx[s] = 0.01;
    }
    enDep[1] = new TGraphErrors(8, x, xEn, dx, 0);
    c1->cd();
    gPad->SetLogx();
    gPad->SetLogy();
    enDep[1]->SetTitle(" ; d (mm); energy (keV)");
    enDep[1]->SetMarkerStyle(8);
    enDep[1]->Draw("ap");
    enDep[1]->Write();

    for(int i = 0; i<8; i++){
        xEn[i] = 5486-xEn[i];
        cout<<pEn[i]<<endl;
    }
    cout<<pEn[7]<<endl;
    pEnx[7] = pEn[7]/dD[7];
    

    // TCanvas *c2 = new TCanvas();
    // TGraphErrors *fibiabsvdf = new TGraphErrors(8, x, xEn, dx, 0);
    // c2->cd();
    // gPad->SetLogx();
    // fibiabsvdf->SetTitle(" ; d (mm); dE/dx (keV/mm)");
    // fibiabsvdf->Draw("ap");
    // fibiabsvdf->SetMarkerStyle(8);
    // fibiabsvdf->Write();

    c1->Write();
    file->Close();
    return kTRUE;
}