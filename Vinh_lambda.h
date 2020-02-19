//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 14 18:05:45 2020 by ROOT version 6.19/01
// from TTree mtree/Hadron EMC + TOF tree
// found on file: /home/lbavinh/Documents/NICA/Task/311032.root
//////////////////////////////////////////////////////////

#ifndef Vinh_h
#define Vinh_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Vinh {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Float_t         bbcz;
   Float_t         cent;
   Int_t           rh;
   Float_t         phir[47];   //[rh]
   Float_t         time[47];   //[rh]
   Float_t         qr0[47];   //[rh]
   Float_t         etar[47];   //[rh]
   Short_t         armr[47];   //[rh]
   Short_t         ring[47];   //[rh]
   Int_t           chid[47];   //[rh]
   Int_t           mh;
   Float_t         alpha[55];   //[mh]
   Short_t         dcarm[55];   //[mh]
   Float_t         p[55];   //[mh]
   Short_t         charge[55];   //[mh]
   Float_t         phi0[55];   //[mh]
   Float_t         the0[55];   //[mh]
   Float_t         phi[55];   //[mh]
   Float_t         ecore[55];   //[mh]
   Float_t         plemc[55];   //[mh]
   Float_t         ecent[55];   //[mh]
   Float_t         temc[55];   //[mh]
   Float_t         temcpi[55];   //[mh]
   Float_t         temcp[55];   //[mh]
   Float_t         temck[55];   //[mh]
   Short_t         sect[55];   //[mh]
   Float_t         isPiemc[55];   //[mh]
   Float_t         isPemc[55];   //[mh]
   Float_t         isKemc[55];   //[mh]
   Int_t           idtwr[55];   //[mh]
   Float_t         sigtof[55];   //[mh]
   Float_t         sigpc3[55];   //[mh]
   Float_t         sigemc[55];   //[mh]
   Float_t         res[55];   //[mh]
   Float_t         ttof[55];   //[mh]
   Int_t           slat[55];   //[mh]
   Float_t         pltof[55];   //[mh]
   Float_t         etof[55];   //[mh]
   Float_t         isPi[55];   //[mh]
   Float_t         isP[55];   //[mh]
   Float_t         isK[55];   //[mh]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_bbcz;   //!
   TBranch        *b_cent;   //!
   TBranch        *b_rh;   //!
   TBranch        *b_phir;   //!
   TBranch        *b_time;   //!
   TBranch        *b_qr0;   //!
   TBranch        *b_etar;   //!
   TBranch        *b_armr;   //!
   TBranch        *b_ring;   //!
   TBranch        *b_chid;   //!
   TBranch        *b_mh;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_dcarm;   //!
   TBranch        *b_p;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_phi0;   //!
   TBranch        *b_the0;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_ecore;   //!
   TBranch        *b_plemc;   //!
   TBranch        *b_ecent;   //!
   TBranch        *b_temc;   //!
   TBranch        *b_temcpi;   //!
   TBranch        *b_temcp;   //!
   TBranch        *b_temck;   //!
   TBranch        *b_sect;   //!
   TBranch        *b_isPiemc;   //!
   TBranch        *b_isPemc;   //!
   TBranch        *b_isKemc;   //!
   TBranch        *b_idtwr;   //!
   TBranch        *b_sigtof;   //!
   TBranch        *b_sigpc3;   //!
   TBranch        *b_sigemc;   //!
   TBranch        *b_res;   //!
   TBranch        *b_ttof;   //!
   TBranch        *b_slat;   //!
   TBranch        *b_pltof;   //!
   TBranch        *b_etof;   //!
   TBranch        *b_isPi;   //!
   TBranch        *b_isP;   //!
   TBranch        *b_isK;   //!

   Vinh(TTree *tree=0);
   virtual ~Vinh();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   // additional member functions
   void ana_event(int jevent, int ievent);
   void loop_a_file(char *file);
   void ana_end(char *outFile);
   //void BookHist(char *outFile);
   static float IsKaonW(float m2tof, float pt);
   static float IsKaonE(float m2tof, float pt);

   static float IsPionW(float m2tof, float pt);
   static float IsPionE(float m2tof, float pt);
   static float IsProtonW(float m2tof, float pt);
   static float IsProtonE(float m2tof, float pt);
};

#endif

#ifdef Vinh_cxx
Vinh::Vinh(TTree *tree) : fChain(0) 
{
 // if parameter tree is not specified (or zero), connect the file
 // used to generate this class and read the Tree.
   /*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/media/lbavinh/My Passport/ROOT_FILE/runlist.list");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/lbavinh/Documents/NICA/Task/asdfgh.root");
      }
      f->GetObject("mtree",tree);

   }
   */
   if (tree) Init(tree);
}

Vinh::~Vinh()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Vinh::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Vinh::LoadTree(Long64_t entry){
 // Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Vinh::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("bbcz", &bbcz, &b_bbcz);
   fChain->SetBranchAddress("cent", &cent, &b_cent);
   fChain->SetBranchAddress("rh", &rh, &b_rh);
   fChain->SetBranchAddress("phir", phir, &b_phir);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("qr0", qr0, &b_qr0);
   fChain->SetBranchAddress("etar", etar, &b_etar);
   fChain->SetBranchAddress("armr", armr, &b_armr);
   fChain->SetBranchAddress("ring", ring, &b_ring);
   fChain->SetBranchAddress("chid", chid, &b_chid);
   fChain->SetBranchAddress("mh", &mh, &b_mh);
   fChain->SetBranchAddress("alpha", alpha, &b_alpha);
   fChain->SetBranchAddress("dcarm", dcarm, &b_dcarm);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("phi0", phi0, &b_phi0);
   fChain->SetBranchAddress("the0", the0, &b_the0);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("ecore", ecore, &b_ecore);
   fChain->SetBranchAddress("plemc", plemc, &b_plemc);
   fChain->SetBranchAddress("ecent", ecent, &b_ecent);
   fChain->SetBranchAddress("temc", temc, &b_temc);
   fChain->SetBranchAddress("temcpi", temcpi, &b_temcpi);
   fChain->SetBranchAddress("temcp", temcp, &b_temcp);
   fChain->SetBranchAddress("temck", temck, &b_temck);
   fChain->SetBranchAddress("sect", sect, &b_sect);
   fChain->SetBranchAddress("isPiemc", isPiemc, &b_isPiemc);
   fChain->SetBranchAddress("isPemc", isPemc, &b_isPemc);
   fChain->SetBranchAddress("isKemc", isKemc, &b_isKemc);
   fChain->SetBranchAddress("idtwr", idtwr, &b_idtwr);
   fChain->SetBranchAddress("sigtof", sigtof, &b_sigtof);
   fChain->SetBranchAddress("sigpc3", sigpc3, &b_sigpc3);
   fChain->SetBranchAddress("sigemc", sigemc, &b_sigemc);
   fChain->SetBranchAddress("res", res, &b_res);
   fChain->SetBranchAddress("ttof", ttof, &b_ttof);
   fChain->SetBranchAddress("slat", slat, &b_slat);
   fChain->SetBranchAddress("pltof", pltof, &b_pltof);
   fChain->SetBranchAddress("etof", etof, &b_etof);
   fChain->SetBranchAddress("isPi", isPi, &b_isPi);
   fChain->SetBranchAddress("isP", isP, &b_isP);
   fChain->SetBranchAddress("isK", isK, &b_isK);
   Notify();
}

Bool_t Vinh::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Vinh::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Vinh::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Vinh_cxx
