#include <TH2F.h>
#include <cmath>
#include <iostream>
using namespace std;

//energy unit [MeV]
//distance unit [m]

Double_t velo(Double_t T, Double_t m){
    Double_t E = T + m;
    return sqrt(E*E - m*m)/E;
}

Double_t betheBloch(Double_t mean, Int_t Z, Double_t I, Double_t v){
    const Double_t m_e = 0.511;
    const Double_t Na = 6.02214076e23;
    const Double_t re = 2.8179403262e-15;
    const Double_t rho = 1e6; //kg/m3

    Double_t K = 4 * TMath::Pi() * Na * re * re * m_e;  
    Double_t beta = v;  
    Double_t gamma = 1 / sqrt(1 - beta*beta);  

    Double_t T_max = 2 * m_e * (gamma * beta)*(gamma * beta);

    Double_t t1 = K * Z * Z * mean * (rho / (beta * beta));
    Double_t temp = (T_max / I)*(T_max / I);
    Double_t t2 = 0.5 * log(temp) - beta*beta;

    return t1 * t2;
}


Bool_t betheSim(Double_t E0 = 5486, Long_t steps = 1e7){
    TFile *simulation = new TFile("graphs.root", "UPDATE");
    Double_t zN = 7;
    Double_t zO = 8;
    Double_t mN = 14;
    Double_t mO = 16;
    Double_t ZAO = zO/mO;
    Double_t ZAN = zN/mN;

    Double_t IN = 93.1;
    Double_t IO = 104.8;
    Double_t mA = 3727; //MeV

    Double_t wN = (0.78*mN)/(0.78*mN+0.22*mO);
    Double_t wO = (0.22*mO)/(0.78*mN+0.22*mO);

    vector <Double_t> E;
    E.push_back(E0);

    Double_t x0 = 0;
    vector <Double_t> x;
    x.push_back(x0);

    TH2D *sim = new TH2D("sim", "Bragg Peak", 1000, 0, steps * 0.00001, 1000, 0, E0+1000);
    sim->Fill(x0, wO*betheBloch(ZAO, 2, IO, velo(E0, mA))+wN*betheBloch(ZAN, 2, IN, velo(E0, mA)));
    std::cout<<sim->GetBinContent(0)<<std::endl;
    
    cout<<wO*betheBloch(ZAO, 2, IO, velo(E0, mA))+wN*betheBloch(ZAN, 2, IN, velo(E0, mA))<<endl;

    // for(Int_t i = 0; i<steps-1; i++){
    //     x.push_back(x.at(i)+0.001);
    //     Double_t temp = (wO*betheBloch(ZAO, 2, IO, velo(E.at(i), mA))+wN*betheBloch(ZAN, 2, IN, velo(E.at(i), mA)));
    //     E.push_back(E.at(i)+temp*0.001);

    //     sim->Fill(x.at(i), (wO*betheBloch(ZAO, 2, IO, velo(E.at(i+1), mA))+wN*betheBloch(ZAN, 2, IN, velo(E.at(i+1), mA))));
    //     //SetBinContent(x.at(i+1), (wO*betheBloch(ZAO, 2, IO, velo(E.at(i+1), mA))+wN*betheBloch(ZAN, 2, IN, velo(E.at(i+1), mA))));

    //     if(i % 100000 == 0){
    //         std::cout<<"iter: "<<i<<std::endl;
    //         std::cout<<"Last bin: "<<sim->GetBinContent(i)<<std::endl;
    //         std::cout<<E.at(i)<<std::endl;
    //     }
    //     if (E.at(i) < 1) break;

    // }
    
    for(Int_t i = 0; i < steps-1; i++){
        x.push_back(x.at(i) + 0.01);
        
        Double_t dEdx = wO*betheBloch(ZAO, 2, IO, velo(E.at(i), mA)) + 
                        wN*betheBloch(ZAN, 2, IN, velo(E.at(i), mA));

        E.push_back(E.at(i) + dEdx * 0.001);  

        if (E.at(i) <= 0) break; 

        sim->Fill(x.at(i), (E0-E.at(i)));

        if(i % 100000 == 0){
            std::cout << "iter: " << i << ", E: " << E.at(i) << ", dE/dx: " << dEdx << ", x"<<x.at(i)<< std::endl;
        }
    }

    cout<<sim->Integral()<<endl;

    TCanvas *c1 = new TCanvas();
    c1->cd();
    gPad->SetLogx();
    sim->SetTitle(" ; d [mm]; dE/dx [keV/mm]");
    sim->Draw("COLZ");
    sim->Write();
    simulation->Close();

    return kTRUE;
}

