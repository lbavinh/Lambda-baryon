#define Vinh_cxx
#include "Vinh_lambda.h"
#include "TProfile.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TPave.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

static const double SoL=299792458;
static const double m2pion=pow(139.57061e-3,2);
static const double m2kaon=pow(493.677e-3,2);
static const double mkaon=493.677e-3;
static const double m2proton=0.880354511;
static const double pTmax=3.8;
static const double pTmin=0.4;
static const double pi=3.141592654;
static const float xmin=0.95;
static const float xmax=1.25;


TH1F *hinvLambdaE1 = new TH1F("hinvLambdaE1","East Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambdaW1 = new TH1F("hinvLambdaW1","West Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambda1 = new TH1F("hinvLambda1","#Lambda-baryon",300,0.98,1.28);

TH1F *hinvLambdaE2 = new TH1F("hinvLambdaE2","East Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambdaW2 = new TH1F("hinvLambdaW2","West Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambda2 = new TH1F("hinvLambda2","#Lambda-baryon",300,0.98,1.28);

TH1F *hinvLambdaE3 = new TH1F("hinvLambdaE3","East Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambdaW3 = new TH1F("hinvLambdaW3","West Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambda3 = new TH1F("hinvLambda3","#Lambda-baryon",300,0.98,1.28);

TH1F *hinvLambdaE4 = new TH1F("hinvLambdaE4","East Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambdaW4 = new TH1F("hinvLambdaW4","West Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambda4 = new TH1F("hinvLambda4","#Lambda-baryon",300,0.98,1.28);

TH1F *hinvLambdaE5 = new TH1F("hinvLambdaE5","East Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambdaW5 = new TH1F("hinvLambdaW5","West Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambda5 = new TH1F("hinvLambda5","#Lambda-baryon",300,0.98,1.28);

TH1F *hinvLambdaE6 = new TH1F("hinvLambdaE6","East arm #Lambda^{0}, #sqrt{s_{NN}}=62 GeV, cent: 10-60%, p_{T}^{#Lambda}>1.8 GeV/c",150,0.98,1.28);
TH1F *hinvLambdaW6 = new TH1F("hinvLambdaW6","West arm #Lambda^{0}, #sqrt{s_{NN}}=62 GeV, cent: 10-60%, p_{T}^{#Lambda}>1.8 GeV/c",150,0.98,1.28);
TH1F *hinvLambda6 = new TH1F("hinvLambda6","Both arm #Lambda^{0}, #sqrt{s_{NN}}=62 GeV, cent: 10-60%, p_{T}^{#Lambda}>1.8 GeV/c",150,0.98,1.28);

TH1F *hinvSTE1 = new TH1F("hinvSTE1","East Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvSTW1 = new TH1F("hinvSTW1","West Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvST1 = new TH1F("hinvST1","#ST-baryon",300,0.98,1.28);

TH1F *hinvSTE2 = new TH1F("hinvSTE2","East Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvSTW2 = new TH1F("hinvSTW2","West Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvST2 = new TH1F("hinvST2","#ST-baryon",300,0.98,1.28);

TH1F *hinvSTE3 = new TH1F("hinvSTE3","East Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvSTW3 = new TH1F("hinvSTW3","West Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvST3 = new TH1F("hinvST3","#ST-baryon",300,0.98,1.28);

TH1F *hinvSTE4 = new TH1F("hinvSTE4","East Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvSTW4 = new TH1F("hinvSTW4","West Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvST4 = new TH1F("hinvST4","#ST-baryon",300,0.98,1.28);

TH1F *hinvSTE5 = new TH1F("hinvSTE5","East Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvSTW5 = new TH1F("hinvSTW5","West Arm #ST-baryon",300,0.98,1.28);
TH1F *hinvST5 = new TH1F("hinvST5","#ST-baryon",300,0.98,1.28);

TH1F *hinvSTE6 = new TH1F("hinvSTE6","East arm neg.charged particle-proton pairs, #sqrt{s_{NN}}=62 GeV, cent: 10-60%, p_{T}^{#Lambda}>1.8 GeV/c",150,0.98,1.28);
TH1F *hinvSTW6 = new TH1F("hinvSTW6","West arm neg.charged particle-proton pairs, #sqrt{s_{NN}}=62 GeV, cent: 10-60%, p_{T}^{#Lambda}>1.8 GeV/c",150,0.98,1.28);
TH1F *hinvST6 = new TH1F("hinvST6","Both arm neg.charged particle-proton pairs, #sqrt{s_{NN}}=62 GeV, cent: 10-60%, p_{T}^{#Lambda}>1.8 GeV/c",150,0.98,1.28);

TH1F *hpE = new TH1F("hpE","Proton East",200,-0.2,1.5);
TH1F *hpW = new TH1F("hpW","Proton West",200,-0.2,1.5);

TH1F *hm2E = new TH1F("hm2E","mass^{2} for Tof.East",200,-0.2,1.5);
TH1F *hm2W = new TH1F("hm2W","mass^{2} for Tof.West",200,-0.2,1.5);

TH1F *hpiW = new TH1F("hpiW","West Arm #pi^{-}",200,-0.2,1.5);
TH1F *hpiE = new TH1F("hpiE","East Arm #pi^{-}",200,-0.2,1.5);

TH1F *hnegW = new TH1F("hnegW","West Arm negative particles",200,-0.2,1.5);
TH1F *hnegE = new TH1F("hnegE","East Arm negative particles",200,-0.2,1.5);

struct trk{
  trk(float Phi0=0, float Theta0=0, float Pt=0):phi0(Phi0), the0(Theta0), pt(Pt) {}
  float phi0;
  float the0;
  float pt;
};
void Vinh::loop_a_file(char *file) {
  TFile *treefile = TFile::Open(file);
  TTree *tree = (TTree*)treefile->Get("mtree");
  
  if(tree == 0) {
    cout <<"htree is not found in "<<file<<endl;
    treefile->Close();
    return;
  }
  cout << file <<" is opened"<<endl;
  Init(tree);
  Loop();
  treefile->Close();
  cout << file  <<" file processed"<<endl;
}
void Vinh::Loop(){
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      ana_event(jentry,ientry);
   }
}
float Vinh::IsPionE(float m2tof, float pt){
  float ispion=-9999.0;
  float sigma=1.0;
  float m2exp=0.01948816;
  sigma=-0.0108109 + 0.0310925*pt - 0.00757494*pt*pt + 0.00519401*pt*pt*pt;
  ispion=(m2tof-m2exp)/sigma;
  return  ispion;
}
float Vinh::IsPionW(float m2tof, float pt){
  float sigma=1.0;
  float m2exp=0.01948816;
  float ispion=-9999.0;
  sigma = 0.00336212 - 0.00674774*pt + 0.0175675*pt*pt - 0.00208485*pt*pt*pt;
  ispion=(m2tof-m2exp)/sigma;
  return  ispion;
}
float Vinh::IsProtonE(float m2tof, float pt){
  float sigma=1.0;
  float m2exp=0.880;
  float isprot=-9999.0;
  sigma = -0.0108109 + 0.0310925*pt - 0.00757494*pt*pt + 0.00519401*pt*pt*pt;
  isprot=(m2tof-m2exp)/sigma;
  return  isprot;
}
float Vinh::IsProtonW(float m2tof, float pt){
  float isprot=-9999.0;
  float sigma=1.0;
  float m2exp=0.880;
  sigma = 0.00336212 - 0.00674774*pt + 0.0175675*pt*pt - 0.00208485*pt*pt*pt;
  isprot=(m2tof-m2exp)/sigma;
  return  isprot;
}
float Vinh::IsKaonW(float m2tof, float pt){
  float sigma=1.0;
  float m2exp=0.24371698;
  float iskaon=-9999.0;  
  sigma=0.0502436-0.104543*pt+0.0950588*pt*pt-0.0194229*pt*pt*pt;
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}
float Vinh::IsKaonE(float m2tof, float pt){
  float sigma=1.0;
  float m2exp=0.24371698;
  float iskaon=-9999.0;
  sigma = 0.0277549 - 0.0394429*pt + 0.0379389*pt*pt - 0.00184303*pt*pt*pt;
  iskaon =(m2tof-m2exp)/sigma;
  return  iskaon;
}
void Vinh::ana_event(int jentry, int ientry) {
    // centrality cut and vertex +/- 30 cm cut
    if(fabs(bbcz)<30 && cent>=10 && cent<=60){//event selection, centrality selection
        if(ientry%100000==0) cout << ientry << endl;
        vector<trk> pW,pE,piW,piE,negW,negE;
        for(int i=0;i<mh;i++) {//track loop
          float pT = p[i]*sin(the0[i]);
          if(pT>pTmin && pT<pTmax){//track's momentum selection
            float m2 = pow(pT,2)*(pow(ttof[i]*(1e-7)*SoL/pltof[i],2)-1);
            if(dcarm[i]==0 && pltof[i]>0 && etof[i]>0.002 && sigtof[i]<3){//Tof East selection
              hm2E -> Fill(m2);
              if(fabs(IsProtonE(m2,pT))<2.0  && fabs(IsKaonE(m2,pT))>3.0  && m2>0.6 && m2<1.4 && pT>=0.5 && charge[i]>0){//Proton selection
                pE.push_back(trk(phi0[i],the0[i],pT));
                hpE -> Fill(m2);
              }
              if(fabs(IsPionE(m2,pT))<2.0 && fabs(IsKaonE(m2,pT))>3.0  && m2<0.13 && m2>-0.25 && charge[i]<0){//Neg Pion selection
                piE.push_back(trk(phi0[i],the0[i],pT));
                hpiE -> Fill(m2);
              }
              if(charge[i]<0){//Neg selection
                negE.push_back(trk(phi0[i],the0[i],pT));
                hnegE -> Fill(m2);
              }
            }//end of Tof East selection
            if(dcarm[i]==1 && pltof[i]>0 && etof[i]>60 && etof[i]<600 && sigtof[i]<3){// Tof West selection                
              hm2W -> Fill(m2);
              if(fabs(IsProtonW(m2,pT))<2.0  && fabs(IsKaonW(m2,pT))>3.0  && m2>0.6 && m2<1.4 && pT>=0.5 && charge[i]>0){//Proton selection
                pW.push_back(trk(phi0[i],the0[i],pT));
                hpW -> Fill(m2);
              }
              if(fabs(IsPionW(m2,pT))<2.0 && fabs(IsKaonW(m2,pT))>3.0  && m2<0.13 && m2>-0.25 && charge[i]<0){//Neg Pion selection
                piW.push_back(trk(phi0[i],the0[i],pT));
                hpiW -> Fill(m2);
              }
              if(charge[i]<0){//Neg selection
                negW.push_back(trk(phi0[i],the0[i],pT));
                hnegW -> Fill(m2);
              }
            }// Tof West selection   
          }//track's momentum selection
        }// end of track loop

        if((piE.size()>0 && pE.size()>0) || (negE.size()>0 && pE.size()>0)) {//loop for east lamda
          for(unsigned int p = 0; p<pE.size(); p++)
	          for(unsigned int pi = 0; pi<piE.size(); pi++){
	            float phi0p = pE[p].phi0;
              float phi0pi = piE[pi].phi0;
              float the0p = pE[p].the0;
              float the0pi = piE[pi].the0;
              float pTp = pE[p].pt;
              float pTpi = piE[pi].pt;
              float pxp = pTp*cos(phi0p);
              float pyp = pTp*sin(phi0p);
              float pzp = pTp/tan(the0p);
              float pxpi = pTpi*cos(phi0pi); 
              float pypi = pTpi*sin(phi0pi);  
              float pzpi = pTpi/tan(the0pi);
              float pP  = sqrt(pTp*pTp + pzp*pzp);
              float EP  = sqrt(pP*pP + m2proton);
              float pPi  = sqrt(pTpi*pTpi + pzpi*pzpi);
              float EPi  = sqrt(pPi*pPi + m2pion);
              float E  = EP + EPi;
              float px = pxp + pxpi;
              float py = pyp + pypi;
              float pz = pzp + pzpi;
              float pt = sqrt(px*px + py*py);
              float Minv = sqrt(E*E-px*px-py*py-pz*pz);

              float phi0neg = negE[p].phi0;
              float the0neg = negE[p].the0;
              float pTneg = negE[p].pt;
              float pxneg = pTneg*cos(phi0neg);
              float pyneg = pTneg*sin(phi0neg);
              float pzneg = pTneg/tan(the0neg);
              float pneg  = sqrt(pTneg*pTneg + pzneg*pzneg);
              float Eneg  = sqrt(pneg*pneg + m2pion);
              float EST  = EP + Eneg;
              float pxST = pxp + pxneg;
              float pyST = pyp + pyneg;
              float pzST = pzp + pzneg;
              float ptST = sqrt(pxST*pxST + pyST*pyST);
              float MinvST = sqrt(EST*EST-pxST*pxST-pyST*pyST-pzST*pzST);                            

              if(pt>0.4 && pt<2. && 3*pTpi<pTp) {hinvLambdaE1->Fill(Minv);hinvLambda1->Fill(Minv);}
              if(pt>0.4 && pt<2.4 && 3*pTpi<pTp) {hinvLambdaE2->Fill(Minv);hinvLambda2->Fill(Minv);}
              if(pt>1.5 && pt<2.5 && 3*pTpi<pTp) {hinvLambdaE3->Fill(Minv);hinvLambda3->Fill(Minv);}
              if(pt>1.2 && pt<2.8 && 3*pTpi<pTp) {hinvLambdaE4->Fill(Minv);hinvLambda4->Fill(Minv);}
              if(pt>1.8 && pt<3.2 && 3*pTpi<pTp) {hinvLambdaE5->Fill(Minv);hinvLambda5->Fill(Minv);}
              if(pt>1.8 && 3*pTpi<pTp) {hinvLambdaE6->Fill(Minv);hinvLambda6->Fill(Minv);}

              if(pt>0.4 && pt<2. && 3*pTneg<pTp) {hinvSTE1->Fill(MinvST);hinvST1->Fill(MinvST);}
              if(pt>0.4 && pt<2.4 && 3*pTneg<pTp) {hinvSTE2->Fill(MinvST);hinvST2->Fill(MinvST);}
              if(pt>1.5 && pt<2.5 && 3*pTneg<pTp) {hinvSTE3->Fill(MinvST);hinvST3->Fill(MinvST);}
              if(pt>1.2 && pt<2.8 && 3*pTneg<pTp) {hinvSTE4->Fill(MinvST);hinvST4->Fill(MinvST);}
              if(pt>1.8 && pt<3.2 && 3*pTneg<pTp) {hinvSTE5->Fill(MinvST);hinvST5->Fill(MinvST);}
              if(pt>1.8 && 3*pTneg<pTp) {hinvSTE6->Fill(MinvST);hinvST6->Fill(MinvST);}

          }
        }//end of loop for east lamda
        if((piW.size()>0 && pW.size()>0) || (negW.size()>0 && pW.size()>0)) {//loop for west lamda
          for(unsigned int p = 0; p<pW.size(); p++)
	          for(unsigned int pi = 0; pi<piW.size(); pi++){
	            float phi0p = pW[p].phi0;float phi0pi = piW[pi].phi0;
              float the0p = pW[p].the0;float the0pi = piW[pi].the0;
              float pTp = pW[p].pt;float pTpi = piW[pi].pt;
              float pxp = pTp*cos(phi0p);float pyp = pTp*sin(phi0p);float pzp = pTp/tan(the0p);
              float pxpi = pTpi*cos(phi0pi);float pypi = pTpi*sin(phi0pi);float pzpi = pTpi/tan(the0pi);
              float pP  = sqrt(pTp*pTp + pzp*pzp);float EP  = sqrt(pP*pP + m2proton);
              float pPi  = sqrt(pTpi*pTpi + pzpi*pzpi);float EPi  = sqrt(pPi*pPi + m2pion);
              float E  = EP + EPi;
              float px = pxp + pxpi;float py = pyp + pypi;float pz = pzp + pzpi;
              float pt = sqrt(px*px + py*py);
              float Minv = sqrt(E*E-px*px-py*py-pz*pz);

              float phi0neg = negW[p].phi0;
              float the0neg = negW[p].the0;
              float pTneg = negW[p].pt;
              float pxneg = pTneg*cos(phi0neg);
              float pyneg = pTneg*sin(phi0neg);
              float pzneg = pTneg/tan(the0neg);
              float pneg  = sqrt(pTneg*pTneg + pzneg*pzneg);
              float Eneg  = sqrt(pneg*pneg + m2pion);
              float EST  = EP + Eneg;
              float pxST = pxp + pxneg;
              float pyST = pyp + pyneg;
              float pzST = pzp + pzneg;
              float ptST = sqrt(pxST*pxST + pyST*pyST);
              float MinvST = sqrt(EST*EST-pxST*pxST-pyST*pyST-pzST*pzST);    

              if(pt>0.4 && pt<2. && 3*pTpi<pTp) {hinvLambdaW1->Fill(Minv);hinvLambda1->Fill(Minv);}
              if(pt>0.4 && pt<2.4 && 3*pTpi<pTp) {hinvLambdaW2->Fill(Minv);hinvLambda2->Fill(Minv);}
              if(pt>1.5 && pt<2.5 && 3*pTpi<pTp) {hinvLambdaW3->Fill(Minv);hinvLambda3->Fill(Minv);}
              if(pt>1.2 && pt<2.8 && 3*pTpi<pTp) {hinvLambdaW4->Fill(Minv);hinvLambda4->Fill(Minv);}
              if(pt>1.8 && pt<3.2 && 3*pTpi<pTp) {hinvLambdaW5->Fill(Minv);hinvLambda5->Fill(Minv);}
              if(pt>1.8 && 3*pTpi<pTp) {hinvLambdaW6->Fill(Minv);hinvLambda6->Fill(Minv);}

              if(pt>0.4 && pt<2. && 3*pTneg<pTp) {hinvSTW1->Fill(MinvST);hinvST1->Fill(MinvST);}
              if(pt>0.4 && pt<2.4 && 3*pTneg<pTp) {hinvSTW2->Fill(MinvST);hinvST2->Fill(MinvST);}
              if(pt>1.5 && pt<2.5 && 3*pTneg<pTp) {hinvSTW3->Fill(MinvST);hinvST3->Fill(MinvST);}
              if(pt>1.2 && pt<2.8 && 3*pTneg<pTp) {hinvSTW4->Fill(MinvST);hinvST4->Fill(MinvST);}
              if(pt>1.8 && pt<3.2 && 3*pTneg<pTp) {hinvSTW5->Fill(MinvST);hinvST5->Fill(MinvST);}
              if(pt>1.8 && 3*pTneg<pTp) {hinvSTW6->Fill(MinvST);hinvST6->Fill(MinvST);}                              

          }
        }//end of loop for east lamda           
    }// end of event selection, end of centrality selection    
}
void Vinh::ana_end(TString outFile) {

  TCanvas*c7 = new TCanvas ("c7","Test canvas2",200,10,1600,1200);
  c7 -> Divide(4,2);
  c7 -> cd(1);
  hpW -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpW -> GetYaxis() -> SetTitle ("Entries");
  hpW -> Draw();  

  c7 -> cd(5);
  hpE -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpE -> GetYaxis() -> SetTitle ("Entries");
  hpE -> Draw();

  c7 -> cd(2);
  hpiW -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpiW -> GetYaxis() -> SetTitle ("Entries");
  hpiW -> Draw();  

  c7 -> cd(6);
  hpiE -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpiE -> GetYaxis() -> SetTitle ("Entries");
  hpiE -> Draw();

  c7 -> cd(3);
  hnegW -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hnegW -> GetYaxis() -> SetTitle ("Entries");
  hnegW -> Draw();  

  c7 -> cd(7);
  hnegE -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hnegE -> GetYaxis() -> SetTitle ("Entries");
  hnegE -> Draw();

  c7->cd(4);
  hm2E -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hm2E -> GetYaxis() -> SetTitle ("Entries");
  hm2E -> Draw();  

  c7->cd(8);
  hm2W -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hm2W -> GetYaxis() -> SetTitle ("Entries");
  hm2W -> Draw();

  TFile *d_outfile = new TFile(outFile.Data(),"recreate");
  hpiE-> Write();
  hpiW-> Write();

  hpE -> Write();
  hpW -> Write();

  hnegE -> Write();
  hnegW -> Write();

  hm2E -> Write();
  hm2W -> Write();

  hinvLambdaE1 -> Write();  hinvLambdaW1 -> Write();  hinvLambda1 -> Write();
  hinvLambdaE2 -> Write();  hinvLambdaW2 -> Write();  hinvLambda2 -> Write();
  hinvLambdaE3 -> Write();  hinvLambdaW3 -> Write();  hinvLambda3 -> Write();
  hinvLambdaE4 -> Write();  hinvLambdaW4 -> Write();  hinvLambda4 -> Write();
  hinvLambdaE5 -> Write();  hinvLambdaW5 -> Write();  hinvLambda5 -> Write();
  hinvLambdaE6 -> Write();  hinvLambdaW6 -> Write();  hinvLambda6 -> Write();

  hinvSTE1 -> Write();  hinvSTW1 -> Write();  hinvST6 -> Write();
  hinvSTE2 -> Write();  hinvSTW2 -> Write();  hinvST2 -> Write();
  hinvSTE3 -> Write();  hinvSTW3 -> Write();  hinvST3 -> Write();
  hinvSTE4 -> Write();  hinvSTW4 -> Write();  hinvST4 -> Write();
  hinvSTE5 -> Write();  hinvSTW5 -> Write();  hinvST5 -> Write();
  hinvSTE6 -> Write();  hinvSTW6 -> Write();  hinvST6 -> Write();

  c7 -> Write();
  d_outfile -> Close();
}
void loop_a_list_of_tree(char* inName,TString outName){
  cout << "Starting the invMassLambda procedure" << endl;
  cout << "Input argument: " << inName << endl;
  Vinh *ana = new Vinh();
  ifstream ifile(inName);
  char filename[200];
  int nfiles=0;
  while(ifile.getline(filename,200)) {
    cout << nfiles << ": processing " << filename << endl;
    ana->loop_a_file(filename);
    nfiles++;
  }
  ana->ana_end(outName.Data());
  cout << "Histfile written!" << endl;
}

