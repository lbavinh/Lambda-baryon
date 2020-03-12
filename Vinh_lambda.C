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

static const double SoL=299792458;// speed of light
static const double m2pion=pow(139.57061e-3,2);
static const double m2kaon=pow(493.677e-3,2);
static const double mkaon=493.677e-3;
static const double m2proton=0.880354511;
static const double pTmax=3.8;
static const double pTmin=0.3;
static const double pi=3.141592654;



TH1F *hinvLambdaE = new TH1F("hinvLambdaE","East Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambdaW = new TH1F("hinvLambdaW","West Arm #Lambda-baryon",300,0.98,1.28);
TH1F *hinvLambda = new TH1F("hinvLambda","#Lambda-baryon",300,0.98,1.28);
TH1F *hpE = new TH1F("hpE","Proton East",200,-0.2,1.5);
TH1F *hpW = new TH1F("hpW","Proton West",200,-0.2,1.5);
TH1F *hm2E = new TH1F("hm2E","mass^{2} for Tof.East",200,-0.2,1.5);
TH1F *hm2W = new TH1F("hm2W","mass^{2} for Tof.West",200,-0.2,1.5);
TH1F *hpiW = new TH1F("hpiW","West Arm #pi^{-}",200,-0.2,1.5);
TH1F *hpiE = new TH1F("hpiE","East Arm #pi^{-}",200,-0.2,1.5);
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
  cout <<"one file processed"<<endl;
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
    if(fabs(bbcz)<30 && cent>=20 && cent<=60){//event selection, centrality selection
        if(ientry%100000==0) cout << ientry << endl;
        vector<trk> pW,pE,piW,piE;
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
            }// end of Tof West selection   
          }//track's momentum selection
        }// end of track loop

        if(piE.size()>0 && pE.size()>0) {//loop for east lambda
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
              if(pt>1.8 && 3*pTpi<pTp) {hinvLambdaE->Fill(Minv);hinvLambda->Fill(Minv);}
          }
        }//end of loop for east lambda
        if(piW.size()>0 && pW.size()>0) {//loop for west lambda
          for(unsigned int p = 0; p<pW.size(); p++)
	          for(unsigned int pi = 0; pi<piW.size(); pi++){
	            float phi0p = pW[p].phi0;
              float phi0pi = piW[pi].phi0;
              float the0p = pW[p].the0;
              float the0pi = piW[pi].the0;
              float pTp = pW[p].pt;
              float pTpi = piW[pi].pt;
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
              if(pt>1.8 && 3*pTpi<pTp) {hinvLambdaW->Fill(Minv);hinvLambda->Fill(Minv);}
          }
        }//end of loop for east lambda           
    }// end of event selection, end of centrality selection    
}
void Vinh::ana_end(TString outFile) {

  TCanvas*MyC1 = new TCanvas ("c1","My canvas",200,10,1600,1200);
  MyC1 -> Divide(3,1);
  MyC1 -> cd(1);
  hinvLambdaW -> GetXaxis() -> SetTitle ("m_{#pi^{-}p}, GeV/c^{2}");
  hinvLambdaW -> GetYaxis() -> SetTitle ("Entries");
  hinvLambdaW -> Draw();

  MyC1 -> cd(2);
  hinvLambdaE -> GetXaxis() -> SetTitle ("m_{#pi^{-}p}, GeV/c^{2}");
  hinvLambdaE -> GetYaxis() -> SetTitle ("Entries");
  hinvLambdaE -> Draw();

  MyC1 -> cd(3);
  hinvLambda -> GetXaxis() -> SetTitle ("m_{#pi^{-}p}, GeV/c^{2}");
  hinvLambda -> GetYaxis() -> SetTitle ("Entries");
  hinvLambda -> Draw();

  TCanvas*MyC2 = new TCanvas ("c2","Test canvas2",200,10,1600,1200);
  MyC2 -> Divide(3,2);
  MyC2 -> cd(1);
  hpW -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpW -> GetYaxis() -> SetTitle ("Entries");
  hpW -> Draw();  

  MyC2 -> cd(4);
  hpE -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpE -> GetYaxis() -> SetTitle ("Entries");
  hpE -> Draw();

  MyC2 -> cd(2);
  hpiW -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpiW -> GetYaxis() -> SetTitle ("Entries");
  hpiW -> Draw();  

  MyC2 -> cd(5);
  hpiE -> SetTitle("");
  hpiE -> SetTitleSize(0.03);
  hpiE -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hpiE -> GetYaxis() -> SetTitle ("Entries");
  hpiE -> Draw();

  MyC2->cd(3);
  hm2E -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hm2E -> GetYaxis() -> SetTitle ("Entries");
  hm2E -> Draw();  

  MyC2->cd(6);
  hm2W -> GetXaxis() -> SetTitle ("m^{2}, (GeV/c^{2})^{2}");
  hm2W -> GetYaxis() -> SetTitle ("Entries");
  hm2W -> Draw();

  TCanvas*c3 = new TCanvas ("c3","Some thing must be here",200,10,1600,1200);
  hinvLambda -> Draw();
  TLatex *latex = new TLatex(1,hinvLambda->GetMaximum(123123123123.),"cent: 20-60%, p_{T}>2.4 GeV/c^{2}, E_{p}>5E_{#pi^{-}}");
  latex -> SetTextFont(42);//https://root.cern.ch/doc/master/classTAttText.html#T5
  latex -> Draw();

  TFile *d_outfile = new TFile(outFile.Data(),"recreate");
  hpiE-> Write();
  hpiW-> Write();

  hpE -> Write();
  hpW -> Write();

  hm2E -> Write();
  hm2W -> Write();

  hinvLambdaE -> Write();
  hinvLambdaW -> Write();
  hinvLambda -> Write();
  MyC1->Write();
  MyC2->Write();
  c3->Write();
  d_outfile->Close();
}
void loop_a_list_of_tree(){
  Vinh *ana = new Vinh();
  ifstream ifile("/mnt/pool/2/lbavinh/runlist.list");
  char filename[200];
  int nfiles=0;
  while(ifile.getline(filename,200)) {
    cout << nfiles<<": processing "<<filename <<endl;
    ana->loop_a_file(filename);
    nfiles++;
  }
  cout<<"Done. "<<nfiles<<" files are processed"<<endl;
  ana->ana_end("/mnt/pool/2/lbavinh/lambda-2.4-inf-5-20-60.root");
  cout<<"Histfile written. Congratulations!"<<endl;
}