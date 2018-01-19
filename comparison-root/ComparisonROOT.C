
// ROOT libs
#include "TTree.h"
#include "ROOT/TDataFrame.hxx"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"

// Std libs
#include <vector>
#include <ctime>
#include <iostream>

using namespace std;

float get_mindR(const vector<float> &eta1, const vector<float> &eta2, const vector<float> &phi1, const vector<float> &phi2) {

  int n1 = eta1.size();
  int n2 = eta2.size();
  bool sameCollection = (eta1==eta2);
  float dRmin=100;

  for (auto i=0 ; i<n1 ; i++){
    if (sameCollection) n2=i;
    for (auto j=0; j<n2; j++){
      float deta2 = TMath::Power( eta1[i]-eta2[j], 2);
      float dphi2 = TMath::Power( phi1[i]-phi2[j], 2);
      float dR    = TMath::Sqrt(deta2+dphi2);
      if (dR<dRmin) dRmin=dR;
    }
  }

  if (dRmin==100) return -1.0;
  else            return dRmin;
  
};


void ComparisonROOT(bool isPara){

  clock_t t0 = clock();
  
  if (isPara) ROOT::EnableImplicitMT();
  
  // Load input file
  auto fileName = "../VectorNtuple_4topSM.root";
  auto treeName = "nominal_Loose";
  ROOT::Experimental::TDataFrame d(treeName, fileName);
  clock_t t_load = clock();

  
  // Apply mindR function and add the branch
  auto d_with_dR = d.
    Define("mindR_mj", get_mindR, {"mu_eta" ,"jet_eta","mu_phi" ,"jet_phi"} ).
    Define("mindR_jj", get_mindR, {"jet_eta","jet_eta","jet_phi","jet_phi"} ).
    Define("mindR_ej", get_mindR, {"el_eta" ,"jet_eta","el_phi" ,"jet_phi"} ).    
    Define("mindR_ee", get_mindR, {"el_eta" ,"el_eta" ,"el_phi" ,"el_phi"});    
  clock_t t_branch = clock();

  // Try some plotting with new computed variables
  auto hdRmj = d_with_dR.Filter("mindR_mj>0").Histo1D({"h1","h1",100,0,6}, "mindR_mj");
  auto hdRjj = d_with_dR.Filter("mindR_jj>0").Histo1D({"h2","h2",100,0,6}, "mindR_jj");
  auto hdRej = d_with_dR.Filter("mindR_ej>0").Histo1D({"h3","h3",100,0,6}, "mindR_ej");
  auto hdRee = d_with_dR.Filter("mindR_ee>0").Histo1D({"h4","h4",100,0,6}, "mindR_ee");
  clock_t t_hist = clock();

  TCanvas *c2 = new TCanvas();
  c2->cd();
  hdRmj->DrawNormalized();
  hdRjj->DrawNormalized("same");
  hdRej->DrawNormalized("same");
  hdRee->DrawNormalized("same");
  c2->Print("c2.png");
  clock_t t_draw = clock();

  // Timing printing
  double elapsed_tot  = double(t_draw - t0    ) / CLOCKS_PER_SEC;
  double elapsed_load = double(t_load - t0    ) / CLOCKS_PER_SEC;
  double elapsed_hist = double(t_hist - t_load) / CLOCKS_PER_SEC;
  double elapsed_draw = double(t_draw - t_hist) / CLOCKS_PER_SEC;
  cout << "total delta: "  << elapsed_tot  << endl;
  cout << "  -> loading: " << elapsed_load << endl;
  cout << "  -> histo  : " << elapsed_hist << endl;
  cout << "  -> drawing: " << elapsed_draw << endl;
  
  return;

}
