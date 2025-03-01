#include <string>
#include <vector>
#include "betheSim.C"

Bool_t graphs(){
    TFile *file = new TFile("graphs3.root", "RECREATE");
    TCanvas *c1 = new TCanvas("c1", "Graphs", 800, 600);
    
    // Dane eksperymentalne
    Double_t x[8] = {22.5, 19.5, 16.5, 13.5, 10.5, 7.5, 4.5, 1.5}; // mm
    Double_t xEn[8] = {3407.66, 3612.94, 3793.33, 3959.86, 4120.72, 4286.72, 4449.63, 4589.51}; // keV
    Double_t dx[8] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    
    Double_t p[8] = {-0.098, -0.091, -0.071, -0.06, -0.048, -0.037, -0.03, -0.022}; // MPa
    Double_t pEn[8] = {5380.97, 5028.96, 4099.4, 3653.2, 3025.82, 2282.05, 1701.02, 971.534}; // keV
    Double_t dp[8] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    Double_t pInt[8] = {15875 / 28.37, 17972 / 32.2, 17257 / 31.27, 14926 / 27.27, 17563 / 32.68, 19874 / 37.78, 23531 / 46.63, 22025 / 49.3};
    
    Double_t p0 = 0.101325;
    Double_t t1 = (40/p0*0.01);

    for(int s = 0; s<8; s++){
        p[s] = ((101325+p[s]*1e6)*40)/101325;
        dp[s] = TMath::Sqrt((t1*t1)+(((p0+p[s])/p0)*0.001)*(((p0+p[s])/p0)*0.1));
        cout<<p[s]<<"+/-"<<dp[s]<<endl;
    }

    // Wykres N(x) - liczba zliczeń w funkcji grubości warstwy powietrza
    TCanvas *cN = new TCanvas("cN", "N(x)", 800, 600);
    TGraphErrors *gN = new TGraphErrors(8, p, pInt, dp, nullptr);
    gN->SetTitle("Liczba zliczeń w funkcji efektywnej grubości warstwy powietrza; d (mm); N (1/s)");
    gN->SetMarkerStyle(20);
    gN->SetMarkerColor(kBlue);
    gN->SetLineColor(kBlue);
    gN->Draw("AP");
    gN->Write();

    // Obliczenie strat energii dE/dx jako różnica energii w kolejnych punktach podzielona przez odległość
    Double_t dEdx[8];
    Double_t dd[8];
    for (int i = 0; i < 7; i++) {
        dEdx[i] = (xEn[i+1] - xEn[i]) / 3;
        dd[i] = (x[i] + x[i+1]); // Średnia pozycja
    }
    

    // Wykres dE/dx(d) - straty energii w funkcji grubości warstwy powietrza
    TCanvas *cDE = new TCanvas("cDE", "dE/dx(d)", 800, 600);
    TGraphErrors *gdE = new TGraphErrors(8, x, dEdx, dx, 0);
    gdE->SetTitle("Straty energii w funkcji grubości warstwy powietrza; d (mm); dE/dx (keV/mm)");
    gdE->SetMarkerStyle(21);
    gdE->SetMarkerColor(kRed);
    gdE->SetLineColor(kRed);
    gdE->Draw("AP");
    gdE->Write();

    // Wykres E(x) - energia cząstek w funkcji grubości warstwy powietrza
    TCanvas *cE = new TCanvas("cE", "E(x)", 800, 600);
    TGraphErrors *gE = new TGraphErrors(8, x, xEn, dx, 0);
    gE->SetTitle("Energia cząstek w funkcji grubości warstwy powietrza; d (mm); E (keV)");
    gE->SetMarkerStyle(22);
    gE->SetMarkerColor(kGreen);
    gE->SetLineColor(kGreen);
    gE->Draw("AP");
    gE->Write();

    // Wykresy dla zależności od odległości x
    // dE/dx(d)
    TCanvas *cDX = new TCanvas("cDX", "dE/dx(x)", 800, 600);
    TGraphErrors *gdX = new TGraphErrors(8, x, xEn, dx, 0);
    gdX->SetTitle("Straty energii w funkcji odległości źródła; x (mm); dE/dx (keV/mm)");
    gdX->SetMarkerStyle(23);
    gdX->SetMarkerColor(kMagenta);
    gdX->SetLineColor(kMagenta);
    gdX->Draw("AP");
    gdX->Write();

    // E(x)
    TCanvas *cEX = new TCanvas("cEX", "E(x)", 800, 600);
    TGraphErrors *gEX = new TGraphErrors(8, x, xEn, dx, 0);
    gEX->SetTitle("Energia cząstek w funkcji odległości; x (mm); E (keV)");
    gEX->SetMarkerStyle(24);
    gEX->SetMarkerColor(kOrange);
    gEX->SetLineColor(kOrange);
    gEX->Draw("AP");
    gEX->Write();

    for(int i = 0; i<8; i++){
        Double_t temp = xEn[i];
        xEn[i] = 5486-temp;
    }

    // Wykres piku Bragga - pokazujący maksymalne straty energii
    TCanvas *cBragg = new TCanvas("cBragg", "Bragg Peak", 800, 600);
    TGraph *gBragg = new TGraph(7, dd, xEn);
    gBragg->SetTitle("Pik Bragga - straty energii w funkcji głębokości; d (mm); dE/dx (keV/mm)");
    gBragg->SetMarkerStyle(23);
    gBragg->SetMarkerColor(kMagenta);
    gBragg->SetLineColor(kMagenta);
    gBragg->Draw("AP");
    gBragg->Write();

    file->Close();
    return kTRUE;
}
