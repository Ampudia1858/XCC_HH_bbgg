#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TEllipse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TF1.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TMatrixDSymEigen.h"
#include "EventShape/Class/interface/EventShape.h"
#include "vector"
#include <vector>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <cstdlib>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#endif

#include "epConstrainHH.h"
#include "LorentzVectorWithErrors.h"

using namespace std;

//------------------------------------------------------------------------------
//////////////TIM fit 
Double_t funPoly4(Double_t* x, Double_t* par) {

  Double_t xx=x[0];
  Double_t yy=x[1];
  
  Double_t result = par[0]*exp(par[1]+par[2]*xx+par[3]*yy+par[4]*pow(xx,2)+par[5]*pow(yy,2)+par[6]*xx*yy+par[7]*pow(xx,3)+par[8]*pow(yy,3)+par[9]*pow(xx,2)*yy+par[10]*xx*pow(yy,2)+par[11]*pow(xx,4)+par[12]*pow(yy,4)+par[13]*pow(xx,3)*yy+par[14]*xx*pow(yy,3)+par[15]*pow(xx*yy,2));
  cout << " funPoly4 " << " x= " << xx << " " << yy 
       << " result= " << result << endl;
  return result;
}

Double_t funPoly2(Double_t* x, Double_t* par) {
  
  Double_t result = par[0]+par[1]*x[0]+par[2]*x[1]+par[3]*pow(x[0],2)+par[4]*pow(x[1],2)+par[5]*x[0]*x[1];
  cout << " funPoly2 " << " x= " << x[0] << " " << x[1] 
       << " par= "  << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5]
       << " result= " << result << endl;
  return result;
}
////////////////TIM fit

void gammaGammaHHGenTreeggbb(const char *inputFile, const char *inputFileBack, const char *inputFileBacktt)
{
	gSystem->Load("libDelphes");
	
	////////creation of File for TMVA
  	TFile *outputTrees = new TFile("outputTreesggbbttESpreadNoFit.root", "recreate");
  	
  	// Create chain of root trees
  	TChain chain("Delphes");
  	chain.Add(inputFile);
  	
  	TChain chainBack("Delphes");
  	chainBack.Add(inputFileBack);
  	
  	TChain chainBacktt("Delphes");
  	chainBacktt.Add(inputFileBacktt);
  	
  	ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  	Long64_t numberOfEntries = treeReader->GetEntries();
  	
  	ExRootTreeReader *treeReaderBack = new ExRootTreeReader(&chainBack);
  	Long64_t numberOfEntriesBack = treeReaderBack->GetEntries();
  	
  	ExRootTreeReader *treeReaderBacktt = new ExRootTreeReader(&chainBacktt);
  	Long64_t numberOfEntriesBacktt = treeReaderBacktt->GetEntries();
  	
  	// Get pointers to branches used in this analysis
	TClonesArray *branchJet = treeReader->UseBranch("Jet");
	TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
	TClonesArray *branchParticle = treeReader->UseBranch("Particle");
	TClonesArray *branchEvent = treeReader->UseBranch("Event");
	TClonesArray *branchElectron = treeReader->UseBranch("Electron");
	TClonesArray *branchMuon = treeReader->UseBranch("Muon");
	
	///////////////TIM  fit
	  constexpr int nParPoly2=6;
	  Double_t paramPoly2[nParPoly2] = {1,0,0,0,0,0};

	  constexpr int nParPoly4=16;
	  Double_t paramPoly4[nParPoly4] = {1,0.5,-.01,-.01,0,0,0,0,0,0,0,0,0,0,0};


	  /* int nBins2dMass=380;
	  double min2dMass=0.;
	  double max2dMass=380.; */
	  
	  int nBins2dMass=180;
	  double min2dMass=70.;
	  double max2dMass=250.;
	  
	  /*int nBins2dMass=36;
	  double min2dMass=114.;
	  double max2dMass=150.; */
	  
	  /*  int nBins2dMass=80;
	  double min2dMass=110.;
	  double max2dMass=150.; */
	  /////////////////////TIM fit
	  
	  
	 /* TH2F *histInvMass2D = new TH2F("inv_mass2d", "Invariant mass of b-tagged vs. non-b-taged jet pairs", nBins2dMass,min2dMass,max2dMass,nBins2dMass,min2dMass,max2dMass);
  	  TH2F *histInvMass2DFit = new TH2F("inv_mass2dFit", "Fitted mass of b-tagged vs. non-b-taged jet pairs", nBins2dMass,min2dMass,max2dMass,nBins2dMass,min2dMass,max2dMass);
  	  TH2F *histInvMass2DTest = new TH2F("inv_mass2dTest", "Fitted mass of b-tagged vs. non-b-taged jet pairs TEST", nBins2dMass,min2dMass,max2dMass,nBins2dMass,min2dMass,max2dMass);
  	  TH1 *histInvMassB1 = new TH1F("inv_massb1", "Invariant mass of b-tagged jet 1 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB1qq = new TH1F("inv_massb1qq", "Invariant mass of b-tagged jet 1 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB1tt = new TH1F("inv_massb1tt", "Invariant mass of b-tagged jet 1 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB2 = new TH1F("inv_massb2", "Invariant mass of b-tagged jet 2 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB2qq = new TH1F("inv_massb2qq", "Invariant mass of b-tagged jet 2 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB2tt = new TH1F("inv_massb2tt", "Invariant mass of b-tagged jet 2 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB1 = new TH1F("inv_massnb1", "Invariant mass of non-b-tagged jet 1 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB1qq = new TH1F("inv_massnb1qq", "Invariant mass of non-b-tagged jet 1 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB1tt = new TH1F("inv_massnb1tt", "Invariant mass of non-b-tagged jet 1 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB2 = new TH1F("inv_massnb2", "Invariant mass of non-b-tagged jet 2 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB2qq = new TH1F("inv_massnb2qq", "Invariant mass of non-b-tagged jet 2 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB2tt = new TH1F("inv_massnb2tt", "Invariant mass of non-b-tagged jet 2 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	   TH1 *histInvMassB1g = new TH1F("inv_massb1g", "(GEN) Invariant mass of b-tagged jet 1 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB1qqg = new TH1F("inv_massb1qqg", "(GEN) Invariant mass of b-tagged jet 1 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB1ttg = new TH1F("inv_massb1ttg", "(GEN) Invariant mass of b-tagged jet 1 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB2g = new TH1F("inv_massb2g", "(GEN) Invariant mass of b-tagged jet 2 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB2qqg = new TH1F("inv_massb2qqg", "(GEN) Invariant mass of b-tagged jet 2 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassB2ttg = new TH1F("inv_massb2ttg", "(GEN) Invariant mass of b-tagged jet 2 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB1g = new TH1F("inv_massnb1g", "(GEN) Invariant mass of non-b-tagged jet 1 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB1qqg = new TH1F("inv_massnb1qqg", "(GEN) Invariant mass of non-b-tagged jet 1 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB1ttg = new TH1F("inv_massnb1ttg", "(GEN) Invariant mass of non-b-tagged jet 1 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB2g = new TH1F("inv_massnb2g", "(GEN) Invariant mass of non-b-tagged jet 2 (HH events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB2qqg = new TH1F("inv_massnb2qqg", "(GEN) Invariant mass of non-b-tagged jet 2 (qq events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histInvMassNB2ttg = new TH1F("inv_massnb2ttg", "(GEN) Invariant mass of non-b-tagged jet 2 (tt events where HH goes to 2 b-tagged and 2-non-b-tagged jets)", 100, 0.0, 30.0);
  	  TH1 *histB1DR04 = new TH1F("b1_dr04", "(HH) Particles with PID 4 or 5 at DR<0.4 from b-jet 1", 100, 0.0, 0.5);
  	  TH1 *histB1DR04qq = new TH1F("b1_dr04qq", "(qq) Particles with PID 4 or 5 at DR<0.4 from b-jet 1", 100, 0.0, 0.5);
  	  TH1 *histB1DR04tt = new TH1F("b1_dr04tt", "(tt) Particles with PID 4 or 5 at DR<0.4 from b-jet 1", 100, 0.0, 0.5);
  	  TH1 *histB2DR04 = new TH1F("b2_dr04", "(HH) Particles with PID 4 or 5 at DR<0.4 from b-jet 2", 100, 0.0, 0.5);
  	  TH1 *histB2DR04qq = new TH1F("b2_dr04qq", "(qq) Particles with PID 4 or 5 at DR<0.4 from b-jet 2", 100, 0.0, 0.5);
  	  TH1 *histB2DR04tt = new TH1F("b2_dr04tt", "(tt) Particles with PID 4 or 5 at DR<0.4 from b-jet 2", 100, 0.0, 0.5);
  	  TH1 *histNB1DR04 = new TH1F("nb1_dr04", "(HH) Particles with PID 4 or 5 at DR<0.4 from non-b-jet 1", 100, 0.0, 0.5);
  	  TH1 *histNB1DR04qq = new TH1F("nb1_dr04qq", "(qq) Particles with PID 4 or 5 at DR<0.4 from non-b-jet 1", 100, 0.0, 0.5);
  	  TH1 *histNB1DR04tt = new TH1F("nb1_dr04tt", "(tt) Particles with PID 4 or 5 at DR<0.4 from non-b-jet 1", 100, 0.0, 0.5);
  	  TH1 *histNB2DR04 = new TH1F("nb2_dr04", "(HH) Particles with PID 4 or 5 at DR<0.4 from non-b-jet 2", 100, 0.0, 0.5);
  	  TH1 *histNB2DR04qq = new TH1F("nb2_dr04qq", "(qq) Particles with PID 4 or 5 at DR<0.4 from non-b-jet 2", 100, 0.0, 0.5);
  	  TH1 *histNB2DR04tt = new TH1F("nb2_dr04tt", "(tt) Particles with PID 4 or 5 at DR<0.4 from non-b-jet 2", 100, 0.0, 0.5);
  	  
  	  
  	histInvMassB1->SetLineColor(kBlue);
	histInvMassB2->SetLineColor(kCyan);
	histInvMassNB1->SetLineColor(kGreen+2);
	histInvMassNB2->SetLineColor(kGreen);
	
	histInvMassB1qq->SetLineColor(kBlue);
	histInvMassB2qq->SetLineColor(kCyan);
	histInvMassNB1qq->SetLineColor(kGreen+2);
	histInvMassNB2qq->SetLineColor(kGreen);
	
	histInvMassB1tt->SetLineColor(kBlue);
	histInvMassB2tt->SetLineColor(kCyan);
	histInvMassNB1tt->SetLineColor(kGreen+2);
	histInvMassNB2tt->SetLineColor(kGreen);
	
	
	histInvMassB1g->SetLineColor(kBlue);
	histInvMassB2g->SetLineColor(kCyan);
	histInvMassNB1g->SetLineColor(kGreen+2);
	histInvMassNB2g->SetLineColor(kGreen);
	
	histInvMassB1qqg->SetLineColor(kBlue);
	histInvMassB2qqg->SetLineColor(kCyan);
	histInvMassNB1qqg->SetLineColor(kGreen+2);
	histInvMassNB2qqg->SetLineColor(kGreen);
	
	histInvMassB1ttg->SetLineColor(kBlue);
	histInvMassB2ttg->SetLineColor(kCyan);
	histInvMassNB1ttg->SetLineColor(kGreen+2);
	histInvMassNB2ttg->SetLineColor(kGreen);*/
	  
	  
	//ALL variables used
	int nJets=0;
	int nParticles=0, nElectrons=0, nMuons=0;
	int contG=0, contB=0;
	float contEventsHH=0, contEventsBBGGBack=0, contEventsBBGGBacktt=0;
	int contBJets=0, contNBJets=0, contPlus4Jets=0, contJets=0;
	float cont2B2NB=0, cont2B2NBBack=0, cont2B2NBBacktt=0;
	float ptSum=0;
	float contMassB=0, contMassNB=0, contMassBBack=0, contMassNBBack=0, contMassBBacktt=0, contMassNBBacktt=0;
	float contMassFilter=0, contMassFilterBack=0, contMassFilterBacktt=0;
	float eta=0, eEta=0, theta=0, jetEta=0, cosTheta=0;
	float contEtaCut=0, contCosThetaCut=0;
	float pt1=0, pt2=0, pt3=0, pt4=0;
	TLorentzVector jetB1, jetB2, jetPairB, jetNB1, jetNB2, jetPairNB, jetSum, jetPairBNB1, jetPairBNB2, jetPairBNB3, jetPairBNB4;
	TLorentzVector jetB1PF, jetB2PF, jetPairBPF, jetNB1PF, jetNB2PF, jetPairNBPF, jetSumPF, jetPairBNB1PF, jetPairBNB2PF, jetPairBNB3PF, jetPairBNB4PF;
	TLorentzVector jetB1Fit, jetB2Fit, jetPairBFit, jetNB1Fit, jetNB2Fit, jetPairNBFit, jetSumFit, jetPairBNB1Fit, jetPairBNB2Fit, jetPairBNB3Fit, jetPairBNB4Fit;  //////TIM fit
	int XXIIIIndex=0;
	float aplanarity=0, sphericity=0;
	float invMassB=0, invMassNB=0, minMass=0;  
	float invMassBFit=0, invMassNBFit=0, minMassFit=0;  
	float invMassBPF=0, invMassNBPF=0, minMassPF=0;  
	Double_t random;
	int conttt=0;
	float jetB1M=0, jetB2M=0, jetNB1M=0, jetNB2M=0;
	float jetB1MFit=0, jetB2MFit=0, jetNB1MFit=0, jetNB2MFit=0;
	float jetB1Mg=0, jetB2Mg=0, jetNB1Mg=0, jetNB2Mg=0;
	int indexb1, indexb2, indexnb1, indexnb2, indexgenb1=0, indexgenb2=0, indexgennb1=0, indexgennb2=0;
	int contDiffJets=0, contDiffJetsqq=0, contDiffJetstt=0, cont5Jets=0, cont5Jetsqq=0, cont5Jetstt=0, contLargeBMass=0, contLargeBMassqq=0, contLargeBMasstt=0;
	double fiveJets[10000];
	double fiveJetsqq[10000];
	double fiveJetstt[10000];
	double largeBMass[100000];
	double largeBMassqq[100000];
	double largeBMasstt[100000];
	int elecCont=0, elecContTotal=0, muonCont=0, muonContTotal=0, elecContqq=0, elecContTotalqq=0, muonContqq=0, muonContTotalqq=0, elecConttt=0, elecContTotaltt=0, muonConttt=0, muonContTotaltt=0, elecContBranch=0, elecContBranchqq=0, elecContBranchtt=0, muonContBranch=0, muonContBranchqq=0, muonContBranchtt=0, elecPt, elecPtqq, elecPttt, muonPt, muonPtqq, muonPttt; ///leptonic veto
	int contG1qq=0, contG123qq=0, contG1tt=0, contG123tt=0, contB1qq=0, contB123qq=0, contB1tt=0, contB123tt=0;
	int contG2qq=0, contG223qq=0, contG2tt=0, contG223tt=0, contB2qq=0, contB223qq=0, contB2tt=0, contB223tt=0;
	int contNG1qq=0, contNG123qq=0, contNG1tt=0, contNG123tt=0, contNB1qq=0, contNB123qq=0, contNB1tt=0, contNB123tt=0;
	int contNG2qq=0, contNG223qq=0, contNG2tt=0, contNG223tt=0, contNB2qq=0, contNB223qq=0, contNB2tt=0, contNB223tt=0;
	int contG1HH=0, contG123HH=0, contG2HH=0, contG223HH=0, contB1HH=0, contB123HH=0, contB2HH=0, contB223HH=0, contNG1HH=0, contNG123HH=0, contNG2HH=0, contNG223HH=0, contNB1HH=0, contNB123HH=0, contNB2HH=0, contNB223HH=0;	
	
	//double invMassBNonFit=0, invMassNBNonFit=0, invMassBFit=0, invMassNBFit=0;
	
	/////////////TIM fit 
	//  double backWt=0.1202325581;
  double backWt=0.04;
  //  double sigWt=0.00287;
  double sigWt=0.00016;
  cout << " sigWt= " << sigWt << " backWt= " << backWt << endl;


  double Ecm=380.;
  double sigmaEcm=1.e-2;  // constraint tolerance for Ecm; to be made larger later to reflect spread in gam gam luminosity near peak
  double sigmaP=1.e-2;  //  constraint tolerance for total Px, Py, Pz
  bool enableExtraTries=false;  //   if true then fit is performed several times with different intial betaX, betaY, betaZ valeus
  double nSigVar=3.;   // if enableExtraTries=true, controls spread in initial beta values w.r.t. values given by jetB1,jetB2,jetNB1,jetNB2

  /*  bool testFittedMass=true;
  double massMaxB=130.15;
  double massMinB=119.85;
  double massMaxNB=130.15;
  double massMinNB=119.85; */   // eff=0.50
  
  bool testFittedMass=true;
  double massMaxB1=134;
  double massMinB1=120;
  double massMaxB2=134;
  double massMinB2=120;  // eff=
  
  /*  bool testFittedMass=true;
  double massMaxB=127.9;
  double massMinB=122.1;
  double massMaxNB=127.9;
  double massMinNB=122.1; */   // eff=0.33
  
  /*  bool testFittedMass=false;
  double massMaxB=130;
  double massMinB=105.5;
  double massMaxNB=130;
  double massMinNB=115.5; */    // eff=0.50

  /*  bool testFittedMass=false;
  double massMaxB=130;
  double massMinB=113.4;
  double massMaxNB=130;
  double massMinNB=117.3; */   // eff=0.33

  cout << " testFittedMass= " <<  testFittedMass << " massMaxB1= " << massMaxB1 << " massMinB1= " << massMinB1 << " massMaxB2= " << massMaxB2 << " massMinB2= " << massMinB2 << endl;
		      
  //  double scaleMaxEntries=0.1;
  //  cout << " scaleMaxEntries= " << scaleMaxEntries << endl;

  //  maxEntries *= scaleMaxEntries;
  ///////////////////////TIM fit
	
  	TTree TreeS("TreeS","a simple Tree with simple variables");
	TreeS.Branch("aplanarity",&aplanarity,"aplanarity/F");
   	TreeS.Branch("invMassB",&invMassB,"invMassB/F");
   	TreeS.Branch("invMassNB",&invMassNB,"invMassNB/F");
   	TreeS.Branch("minMass",&minMass,"minMass/F");  
	TreeS.Branch("sphericity",&sphericity,"sphericity/F");
	TreeS.Branch("cosTheta",&cosTheta,"cosTheta/F");
	TreeS.Branch("nJets",&nJets,"nJets/F");
	TreeS.Branch("ptSum",&ptSum,"ptSum/F");
	TreeS.Branch("pt1",&pt1,"pt1/F");
	TreeS.Branch("pt2",&pt2,"pt2/F");
	TreeS.Branch("pt3",&pt3,"pt3/F");
	TreeS.Branch("pt4",&pt4,"pt4/F");
	TreeS.Branch("jetB1M",&jetB1M,"jetB1M/F");
	TreeS.Branch("jetB2M",&jetB2M,"jetB2M/F");
	TreeS.Branch("jetNB1M",&jetNB1M,"jetNB1M/F");
	TreeS.Branch("jetNB2M",&jetNB2M,"jetNB2M/F");
	  
	TTree TreeB("TreeB","a Bimple Tree with Bimple variables");
	TreeB.Branch("aplanarity",&aplanarity,"aplanarity/F");
   	TreeB.Branch("invMassB",&invMassB,"invMassB/F");
   	TreeB.Branch("invMassNB",&invMassNB,"invMassNB/F");
   	TreeB.Branch("minMass",&minMass,"minMass/F");
   	TreeB.Branch("sphericity",&sphericity,"sphericity/F"); 
   	TreeB.Branch("cosTheta",&cosTheta,"cosTheta/F");
	TreeB.Branch("nJets",&nJets,"nJets/F");
	TreeB.Branch("ptSum",&ptSum,"ptSum/F");
	TreeB.Branch("pt1",&pt1,"pt1/F");
	TreeB.Branch("pt2",&pt2,"pt2/F");
	TreeB.Branch("pt3",&pt3,"pt3/F");
	TreeB.Branch("pt4",&pt4,"pt4/F");
	TreeB.Branch("jetB1M",&jetB1M,"jetB1M/F");
	TreeB.Branch("jetB2M",&jetB2M,"jetB2M/F");
	TreeB.Branch("jetNB1M",&jetNB1M,"jetNB1M/F");
	TreeB.Branch("jetNB2M",&jetNB2M,"jetNB2M/F"); 
   	
   	   cout<<endl<<endl<<"For signal: "<<endl<<endl;
   	   
   	// fill the tree
   for(Int_t entry=0; entry<numberOfEntries; entry++) 
   {   
      treeReader->ReadEntry(entry);
          
      if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
      else nJets=0;
      
      contBJets=0;
      contNBJets=0;
      contG=0;
      contB=0;
      contEtaCut=0;
      contCosThetaCut=0;
      ptSum=0;
      
      if(nJets==4)
      {
      	  for(int i=0;i<nJets;i++)
	  {
		Jet *jet = (Jet*) branchJet->At(i);
		if(jet->BTag==1)
		{
			contBJets++;
			if(contBJets==1)
			{
				jetB1=jet->P4();
				jetB1M=jet->Mass;
				indexb1=i;
			}	
			else if(contBJets==2)
			{
				jetB2=jet->P4();
				jetB2M=jet->Mass;
				indexb2=i;
			}
		}
		else
		{
			contNBJets++;
			if(contNBJets==1)
			{
				jetNB1=jet->P4();
				jetNB1M=jet->Mass;
				indexnb1=i;
			}
			else if(contNBJets==2)
			{
				jetNB2=jet->P4();
				jetNB2M=jet->Mass;
				indexnb2=i;
			}
		}
		if(contBJets>2 || contNBJets>2) 
		{
			contPlus4Jets++;
			break;
		}
	  }
	  if(contBJets==2 && contNBJets==2)
	  {
	  	if(branchParticle->GetEntries() > 0) nParticles =  branchParticle->GetEntries();
    		else nParticles=0;
    		bool elecFlag=false;
    		int elecContInd=0;
    		for(int j=0;j<nParticles;j++)
    		{
  			GenParticle *particle = (GenParticle*) branchParticle->At(j); 
  			if(particle->Status==23) XXIIIIndex=j;
		  	int d1Index=particle->D1;
		  	int d2Index=particle->D2;
		  	if (d1Index >= 0 && d2Index >= 0)
		  	{
			  int d1=abs(static_cast<GenParticle*>(branchParticle->At(d1Index))->PID);
			  int d2=abs(static_cast<GenParticle*>(branchParticle->At(d2Index))->PID);
			  if ((particle->PID == 25 || particle->PID == 36) && (d1==21 && d2==21)) contG++;
			  if ((particle->PID == 25 || particle->PID == 36) && (d1==5 && d2==5)) contB++; 
  			}
     		}
     		//cout<<entry<<" has electrons: "<<elecContInd<<endl<<endl;
     		
     		if(contG>=1 && contB>=1)
     		{
			///////lepton veto
			if(branchElectron->GetEntries() > 0) nElectrons =  branchElectron->GetEntries();
	    		else nElectrons=0;
	     		Electron *electron = (Electron*) branchElectron->At(0);
	     		if (nElectrons!=0)
	     		{ 
	     			elecContBranch+=branchElectron->GetEntries();
	     			elecPt=electron->PT;
	     		}
	     		
	     		if(branchMuon->GetEntries() > 0) nMuons =  branchMuon->GetEntries();
	    		else nMuons=0;
	     		Muon *muon = (Muon*) branchMuon->At(0);
	     		if (nMuons!=0)
	     		{
	     			muonContBranch+=branchMuon->GetEntries();
	     			muonPt=muon->PT;
	     		}
			//////lepton veto
			
			contEventsHH++;
			//cout<<i<<", ";
			for(int i=0;i<nJets;i++)
			{
				Jet *jet = (Jet*) branchJet->At(i);
				
				//////////jet pt
				ptSum+=jet->PT;
				if(i==0)
				{ 
					pt1 = jet->PT;
				}
				else if(i==1)
				{
					pt2 = jet->PT;
				}
				else if(i==2)
				{ 
					pt3 = jet->PT;
				}
				else if(i==3)
				{
					pt4 = jet->PT;
				}
				/////////jet pt
				
				////////angles
				eta = jet->Eta;
				eEta = pow(TMath::E(), -eta);
				theta = 2*TMath::ATan(eEta);
				cosTheta = TMath::Cos(theta);
        			jetEta=abs(jet->Eta);
				if(jetEta<1.5) contEtaCut++;
        			if(cosTheta>-0.91 && cosTheta<0.91) contCosThetaCut++;
        			////////angles
			}
			if(contBJets==2 && contNBJets==2 && contCosThetaCut==4 && pt1>80 && pt2>60 && pt3>15 && pt4>15 && ptSum>240 && elecPt<6 && muonPt<8)
			{
				cont2B2NB++;
				
				/////////event shape
				Jet *jet1 = (Jet*) branchJet->At(0);
				Jet *jet2 = (Jet*) branchJet->At(1);
				Jet *jet3 = (Jet*) branchJet->At(2);
				Jet *jet4 = (Jet*) branchJet->At(3);
				TLorentzVector j1 = jet1->P4(), j2=jet2->P4(), j3=jet3->P4(), j4 = jet4->P4();
				vector<TVector3> vectorCollection = {
				  TVector3(j1.Px(), j1.Py(), j1.Pz()),
				  TVector3(j2.Px(), j2.Py(), j2.Pz()),
				  TVector3(j3.Px(), j3.Py(), j3.Pz()),
				  TVector3(j4.Px(), j4.Py(), j4.Pz())
				};
				TMatrixDSym momentumTensor(3);
				for(std::vector<TVector3>::const_iterator p=vectorCollection.begin(); p!=vectorCollection.end(); p++){
				    for(int k=0; k<3; k++){
				      for(int m=0; m<=k; m++){
					momentumTensor[k][m] += (*p)[k]*(*p)[m];
				      }
				    }
				}
				momentumTensor *= 1/(momentumTensor[0][0] + momentumTensor[1][1] + momentumTensor[2][2]);
				TMatrixDSymEigen eigenSystem(momentumTensor);
				TVectorD eigenvalues = eigenSystem.GetEigenValues();
				float eigenvalue1_ = eigenvalues[0];
				float eigenvalue2_ = eigenvalues[1];
				float eigenvalue3_ = eigenvalues[2];
				sphericity = 1.5*(eigenvalue2_ + eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
				aplanarity = 1.5*(eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
				/////////event shape
				
				////////////////////TIM fit
				/////// Beam energy-momentum constrained fit for jetB1,jetB2,jetNB1,jetNB2

		      vector<LorentzVectorWithErrors>* inputLVWE=new vector<LorentzVectorWithErrors>(4);
		      // for now just assign uniform 100% errors
		      inputLVWE->at(0)=LorentzVectorWithErrors(jetB1,jetB1.Energy(),1.,1.,1.);
		      inputLVWE->at(1)=LorentzVectorWithErrors(jetB2,jetB2.Energy(),1.,1.,1.);
		      inputLVWE->at(2)=LorentzVectorWithErrors(jetNB1,jetNB1.Energy(),1.,1.,1.);
		      inputLVWE->at(3)=LorentzVectorWithErrors(jetNB2,jetNB2.Energy(),1.,1.,1.);

		      epConstrainHH eCHH(inputLVWE,Ecm,sigmaEcm,sigmaP,enableExtraTries,nSigVar);
		      vector<LorentzVectorWithErrors>* outputLVWE=eCHH.NumericalMinimization();


		      if(cont2B2NB < 20) { 
			for(int i=0; i<4; i++) {
			  TLorentzVector tLinput=inputLVWE->at(i).getLV();
			  TLorentzVector tLoutput=outputLVWE->at(i).getLV();
			  //cout << " i= " << i << "  input " << " E,px,py,pz= " << tLinput.E() << " " << tLinput.Px()  << " " << tLinput.Py()  << " " << tLinput.Pz() << endl;
			 // cout << " i= " << i << " output " << " E,px,py,pz= " << tLoutput.E() << " " << tLoutput.Px()  << " " << tLoutput.Py()  << " " << tLoutput.Pz() << endl;
			}
		      }

		      jetB1Fit=outputLVWE->at(0).getLV();
		      jetB2Fit=outputLVWE->at(1).getLV();
		      jetNB1Fit=outputLVWE->at(2).getLV();
		      jetNB2Fit=outputLVWE->at(3).getLV();

		      delete inputLVWE;
		      delete outputLVWE;

			////////////////////TIM fit
			
			
			////////W/O fit
				/*jetB1M=static_cast<float>(jetB1.M());
				jetB2M=static_cast<float>(jetB2.M());
				jetNB1M=static_cast<float>(jetNB1.M());
				jetNB2M=static_cast<float>(jetNB2.M());*/
				int nJetsg;
				if(branchGenJet->GetEntries() > 0)  nJetsg =  branchGenJet->GetEntries();
      				else nJetsg=0;
      				if(nJetsg!=4) contDiffJets++;
      				if(nJetsg==5) fiveJets[cont5Jets]=entry;
      				if(nJetsg==5) cont5Jets++;
      				if(jetNB1M>25) largeBMass[contLargeBMass]=entry;
				if(jetNB1M>25) contLargeBMass++;
				//cout<<"nJetsg: "<<nJetsg<<endl<<"indexb1: "<<indexb1<<endl<<"indexb2: "<<indexb2<<endl<<"indexnb1: "<<indexnb1<<endl<<"indexnb2: "<<indexnb2<<endl<<endl;
				if (nJetsg==4)
				{
					Jet *jetb1 = (Jet*) branchJet->At(indexb1);
					Jet *jetb2 = (Jet*) branchJet->At(indexb2);
					Jet *jetnb1 = (Jet*) branchJet->At(indexnb1);
					Jet *jetnb2 = (Jet*) branchJet->At(indexnb2);
					
					////////cont parts
					int jetFlav=jetb1->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==5) contB1HH++;
					else contNB1HH++;
					if(jetFlav==5 && (jetb1->Mass)>23) contB123HH++;
					if(jetFlav!=5 && (jetb1->Mass)>23) contNB123HH++;
					
					jetFlav=jetb2->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==5) contB2HH++;
					else contNB2HH++;
					if(jetFlav==5 && (jetb2->Mass)>23) contB223HH++;
					if(jetFlav!=5 && (jetb2->Mass)>23) contNB223HH++;
					
					jetFlav=jetnb1->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==21) contG1HH++;
					else contNG1HH++;
					if(jetFlav==21 && (jetnb1->Mass)>23) contG123HH++;
					if(jetFlav!=21 && (jetnb1->Mass)>23) contNG123HH++;
					
					jetFlav=jetnb2->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==21) contG2HH++;
					else contNG2HH++;
					if(jetFlav==21 && (jetnb2->Mass)>23) contG223HH++;
					if(jetFlav!=21 && (jetnb2->Mass)>23) contNG223HH++;
					///////cont parts
					
					double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2, deltaRMin;
					
					/////gen jet masses
					for(int j=0;j<nJetsg;j++)
					{
						Jet *genjet = (Jet*) branchGenJet->At(j);
						double genjetEta = genjet->Eta;
     						double genjetPhi = genjet->Phi;
     			
			     			double jetB1Eta = jetb1->Eta;
			     			double jetB1Phi = jetb1->Phi;
			     			deltaRB1 = sqrt(pow(genjetEta-jetB1Eta, 2) + pow(genjetPhi-jetB1Phi, 2));
			     			
			     			double jetB2Eta = jetb2->Eta;
			     			double jetB2Phi = jetb2->Phi;
			     			deltaRB2 = sqrt(pow(genjetEta-jetB2Eta, 2) + pow(genjetPhi-jetB2Phi, 2));
			     			
			     			double jetNB1Eta = jetnb1->Eta;
			     			double jetNB1Phi = jetnb1->Phi;
			     			deltaRNB1 = sqrt(pow(genjetEta-jetNB1Eta, 2) + pow(genjetPhi-jetNB1Phi, 2));
			     			
			     			double jetNB2Eta = jetnb2->Eta;
			     			double jetNB2Phi = jetnb2->Phi;
			     			deltaRNB2 = sqrt(pow(genjetEta-jetNB2Eta, 2) + pow(genjetPhi-jetNB2Phi, 2));
			     			
			     			if(deltaRB1<deltaRB2 && deltaRB1<deltaRNB1 && deltaRB1<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassB1g->Fill(genjet->Mass, 0.001225);
			     				indexgenb1=j;
			     			}
			     			else if(deltaRB2<deltaRB1 && deltaRB2<deltaRNB1 && deltaRB2<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassB2g->Fill(genjet->Mass, 0.001225);
			     				indexgenb2=j;
			     			}
			     			else if(deltaRNB1<deltaRB2 && deltaRNB1<deltaRB1 && deltaRNB1<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassNB1g->Fill(genjet->Mass, 0.001225);
			     				indexgennb1=j;
			     			}
			     			else if(deltaRNB2<deltaRB2 && deltaRNB2<deltaRB1 && deltaRNB2<deltaRNB1)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassNB2g->Fill(genjet->Mass, 0.001225);
			     				indexgennb2=j;
			     			}
					}    ////gen jet masses

					/*cout<<"jetB1M: "<<jetB1M<<"      b-tag: "<<jetb1->BTag<<"      Flavor: "<<jetb1->Flavor<<endl<<"jetB1Mg: "<<jetB1Mg<<endl;
					cout<<"jetB2M: "<<jetB2M<<"      b-tag: "<<jetb2->BTag<<"      Flavor: "<<jetb2->Flavor<<endl<<"jetB2Mg: "<<jetB2Mg<<endl;
					cout<<"jetNB1M: "<<jetNB1M<<"      b-tag: "<<jetnb1->BTag<<"      Flavor: "<<jetnb1->Flavor<<endl<<"jetNB1Mg: "<<jetNB1Mg<<endl;
					cout<<"jetNB2M: "<<jetNB2M<<"      b-tag: "<<jetnb2->BTag<<"      Flavor: "<<jetnb2->Flavor<<endl<<"jetNB2Mg: "<<jetNB2Mg<<endl<<endl;*/
					
					/////RecoJets Masses
					/*if ((jetb1->PT)>60 && (jetb1->PT)<70) histInvMassB1->Fill(jetB1M, 0.001225);
					if ((jetb2->PT)>60 && (jetb2->PT)<70) histInvMassB2->Fill(jetB2M, 0.001225);
					if ((jetnb1->PT)>60 && (jetnb1->PT)<70) histInvMassNB1->Fill(jetNB1M, 0.001225);
					if ((jetnb2->PT)>60 && (jetnb2->PT)<70) histInvMassNB2->Fill(jetNB2M, 0.001225);*/
					
					Jet *genjetb1 = (Jet*) branchGenJet->At(indexgenb1);
					Jet *genjetb2 = (Jet*) branchGenJet->At(indexgenb2);
					Jet *genjetnb1 = (Jet*) branchGenJet->At(indexgennb1);
					Jet *genjetnb2 = (Jet*) branchGenJet->At(indexgennb2);
					
					int nParticlesDR=0;
					if(branchParticle->GetEntries() > 0) nParticlesDR =  branchParticle->GetEntries();
			    		else nParticlesDR=0;
										
					for(int j=0;j<nParticlesDR;j++)
			    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(j); 
			  			int pid=abs(particle->PID);
			  			 if (pid==5 || pid==4)
		  				{
							double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2;
							
							double partEta = particle->Eta;
					     		double partPhi = particle->Phi;
					     		
					     		double genjetB1Eta = genjetb1->Eta;
				     			double genjetB1Phi = genjetb1->Phi;
				     			deltaRB1 = sqrt(pow(partEta-genjetB1Eta, 2) + pow(partPhi-genjetB1Phi, 2));
				     			
				     			double genjetB2Eta = genjetb2->Eta;
				     			double genjetB2Phi = genjetb2->Phi;
				     			deltaRB2 = sqrt(pow(partEta-genjetB2Eta, 2) + pow(partPhi-genjetB2Phi, 2));
				     			
				     			double genjetNB1Eta = genjetnb1->Eta;
				     			double genjetNB1Phi = genjetnb1->Phi;
				     			deltaRNB1 = sqrt(pow(partEta-genjetNB1Eta, 2) + pow(partPhi-genjetNB1Phi, 2));
				     			
				     			double genjetNB2Eta = genjetnb2->Eta;
				     			double genjetNB2Phi = genjetnb2->Phi;
				     			deltaRNB2 = sqrt(pow(partEta-genjetNB2Eta, 2) + pow(partPhi-genjetNB2Phi, 2));
				     			
				     			/*if (deltaRB1<0.4) histB1DR04->Fill(deltaRB1, 0.001225);
				     			if (deltaRB2<0.4) histB2DR04->Fill(deltaRB2, 0.001225);
				     			if (deltaRNB1<0.4) histNB1DR04->Fill(deltaRNB1, 0.001225);
				     			if (deltaRNB2<0.4) histNB2DR04->Fill(deltaRNB2, 0.001225);*/
  		 				 }
					}
					
				}
				jetPairB=jetB1+jetB2;
				jetPairNB=jetNB1+jetNB2;
				jetSum=jetB1+jetB2+jetNB1+jetNB2;
				jetPairBNB1=jetB1+jetNB1;
				jetPairBNB2=jetB2+jetNB2;
				jetPairBNB3=jetB1+jetNB2;
				jetPairBNB4=jetB2+jetNB1;
				minMass = TMath::Min(static_cast<float>(jetPairBNB1.M()), static_cast<float>(jetPairBNB2.M()));
				minMass = TMath::Min(minMass, static_cast<float>(jetPairBNB3.M()));
				minMass = TMath::Min(minMass, static_cast<float>(jetPairBNB4.M()));
				invMassB = jetPairB.M();
				invMassNB = jetPairNB.M();
				////////W/O fit
				/////With fit
				jetB1MFit=static_cast<float>(jetB1Fit.M());
				jetB2MFit=static_cast<float>(jetB2Fit.M());
				jetNB1MFit=static_cast<float>(jetNB1Fit.M());
				jetNB2MFit=static_cast<float>(jetNB2Fit.M());
				jetPairBFit=jetB1Fit+jetB2Fit;
				jetPairNBFit=jetNB1Fit+jetNB2Fit;
				jetSumFit=jetB1Fit+jetB2Fit+jetNB1Fit+jetNB2Fit;
				jetPairBNB1Fit=jetB1Fit+jetNB1Fit;
				jetPairBNB2Fit=jetB2Fit+jetNB2Fit;
				jetPairBNB3Fit=jetB1Fit+jetNB2Fit;
				jetPairBNB4Fit=jetB2Fit+jetNB1Fit;
				minMassFit = TMath::Min(static_cast<float>(jetPairBNB1Fit.M()), static_cast<float>(jetPairBNB2Fit.M()));
				minMassFit = TMath::Min(minMassFit, static_cast<float>(jetPairBNB3Fit.M()));
				minMassFit = TMath::Min(minMassFit, static_cast<float>(jetPairBNB4Fit.M()));
				invMassBFit = jetPairBFit.M();
				invMassNBFit = jetPairNBFit.M();
				///////////With fit
				float invMassSum = jetSum.M();
				float invMassSumFit = jetSumFit.M(); //////////////TIM fit 
				////////////invariant mass
				
				
				/*if(invMassBFit>30 && invMassNBFit>30 && invMassBFit<=380 && invMassNBFit<=380)
				{
					invMassB=invMassBFit;
					invMassNB=invMassNBFit;
					if (minMassFit>0) minMass=minMassFit;
					jetB1M=jetB1MFit;
					jetB2M=jetB2MFit;
					jetNB1M=jetNB1MFit;
					jetNB2M=jetNB2MFit;
				}*/
				
				if(jetB1M<0) jetB1M=0;
				if(jetB2M<0) jetB2M=0;
				if(jetNB1M<0) jetNB1M=0;
				if(jetNB2M<0) jetNB2M=0;
				
				/*histInvMass2D->Fill(invMassBNonFit, invMassNBNonFit, sigWt);  /////TIM fit
		      		histInvMass2DFit->Fill(invMassBFit, invMassNBFit, sigWt);   /////TIM Fit
		      		histInvMass2DTest->Fill(invMassB, invMassNB, sigWt);  /////TIM fit*/
				////////////invariant mass
				
				
				if(aplanarity<0) cout<<"Event: "<<entry+1<<endl<<"aplanarity: "<<aplanarity<<endl;
				if(invMassB<0) cout<<"Event: "<<entry+1<<endl<<"invMassB: "<<invMassB<<endl;
				if(invMassNB<0) cout<<"Event: "<<entry+1<<endl<<"invmassNB: "<<invMassNB<<endl;
				if(minMass<0) cout<<"Event: "<<entry+1<<endl<<"minMass: "<<minMass<<endl;
				if(sphericity<0) cout<<"Event: "<<entry+1<<endl<<"sohericity: "<<sphericity<<endl;
				if(ptSum<0) cout<<"Event: "<<entry+1<<endl<<"ptSum: "<<ptSum<<endl;
				if(pt1<0) cout<<"Event: "<<entry+1<<endl<<"pt1: "<<pt1<<endl;
				if(pt2<0) cout<<"Event: "<<entry+1<<endl<<"pt2: "<<pt2<<endl;
				if(pt3<0) cout<<"Event: "<<entry+1<<endl<<"pt3: "<<pt3<<endl;
				if(pt4<0) cout<<"Event: "<<entry+1<<endl<<"pt4: "<<pt4<<endl;
				if(jetB1M<0) cout<<"Event: "<<entry+1<<endl<<"jetB1M: "<<jetB1M<<endl;
				if(jetB2M<0) cout<<"Event: "<<entry+1<<endl<<"jetB2M: "<<jetB2M<<endl;
				if(jetNB1M<0) cout<<"Event: "<<entry+1<<endl<<"jetNB1M: "<<jetNB1M<<endl;
				if(jetNB2M<0) cout<<"Event: "<<entry+1<<endl<<"jetNB2M: "<<jetNB2M<<endl<<endl;
				
				//if (entry!=150983 && entry!=450261 && entry!=956345 && entry!=963429) TreeS.Fill();
				
				if(invMassB>=400) cout<<entry<<": "<<invMassB<<endl;
				if(invMassNB>=400) cout<<entry<<": "<<invMassNB<<endl<<endl;
				TreeS.Fill();
		       	}
	  	}
	  }
	  
      }
      
      
   }   
   //////////////////Signal
   
   
   
   
   /////////////////Background
   // Get pointers to branches used in this analysis
   branchJet = treeReaderBack->UseBranch("Jet");
   branchGenJet = treeReaderBack->UseBranch("GenJet");
   branchParticle = treeReaderBack->UseBranch("Particle");
   branchEvent = treeReaderBack->UseBranch("Event");
   branchElectron = treeReaderBack->UseBranch("Electron");
   branchMuon = treeReaderBack->UseBranch("Muon");
   
   cout<<endl<<endl<<"For back qq: "<<endl<<endl;
   
   for(Int_t entry=0; entry<numberOfEntriesBack; entry++) 
   {   
      treeReaderBack->ReadEntry(entry);
          
      if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
      else nJets=0;
      
      contBJets=0;
      contNBJets=0;
      contG=0;
      contB=0;
      contEtaCut=0;
      contCosThetaCut=0;
      ptSum=0;
      
      if(nJets==4)
      {
      	  for(int i=0;i<nJets;i++)
	  {
		Jet *jet = (Jet*) branchJet->At(i);
		
		if(jet->BTag==1)
		{
			contBJets++;
			if(contBJets==1)
			{
				jetB1=jet->P4();
				jetB1M=jet->Mass;
				indexb1=i;
			}	
			else if(contBJets==2)
			{
				jetB2=jet->P4();
				jetB2M=jet->Mass;
				indexb2=i;
			}
		}
		else
		{
			contNBJets++;
			if(contNBJets==1)
			{
				jetNB1=jet->P4();
				jetNB1M=jet->Mass;
				indexnb1=i;
			}
			else if(contNBJets==2)
			{
				jetNB2=jet->P4();
				jetNB2M=jet->Mass;
				indexnb2=i;
			}
		}
		if(contBJets>2 || contNBJets>2) 
		{
			contPlus4Jets++;
			break;
		}
	  }
	  if(contBJets==2 && contNBJets==2)
	  {
			///////lepton veto
			if(branchElectron->GetEntries() > 0) nElectrons =  branchElectron->GetEntries();
	    		else nElectrons=0;
	     		Electron *electron = (Electron*) branchElectron->At(0);
	     		if (nElectrons!=0)
	     		{ 
	     			elecContBranchqq+=branchElectron->GetEntries();
	     			elecPt=electron->PT;
	     		}
	     		
	     		if(branchMuon->GetEntries() > 0) nMuons =  branchMuon->GetEntries();
	    		else nMuons=0;
	     		Muon *muon = (Muon*) branchMuon->At(0);
	     		if (nMuons!=0)
	     		{
	     			muonContBranchqq+=branchMuon->GetEntries();
	     			muonPt=muon->PT;
	     		}
			//////lepton veto
			
			contEventsBBGGBack++;
			//cout<<i<<", ";
			for(int i=0;i<nJets;i++)
			{
				Jet *jet = (Jet*) branchJet->At(i);
				
				//////////jet pt
				ptSum+=jet->PT;
				if(i==0)
				{ 
					pt1 = jet->PT;
				}
				else if(i==1)
				{
					pt2 = jet->PT;
				}
				else if(i==2)
				{ 
					pt3 = jet->PT;
				}
				else if(i==3)
				{
					pt4 = jet->PT;
				}
				/////////jet pt
				
				////////angles
				eta = jet->Eta;
				eEta = pow(TMath::E(), -eta);
				theta = 2*TMath::ATan(eEta);
				cosTheta = TMath::Cos(theta);
        			jetEta=abs(jet->Eta);
				if(jetEta<1.5) contEtaCut++;
        			if(cosTheta>-0.91 && cosTheta<0.91) contCosThetaCut++;
        			////////angles
			}
			if(contBJets==2 && contNBJets==2 && contCosThetaCut==4 && pt1>80 && pt2>60 && pt3>15 && pt4>15 && ptSum>240 && elecPt<6 && muonPt<8)
			{
				cont2B2NBBack++;
				
				/////////event shape
				Jet *jet1 = (Jet*) branchJet->At(0);
				Jet *jet2 = (Jet*) branchJet->At(1);
				Jet *jet3 = (Jet*) branchJet->At(2);
				Jet *jet4 = (Jet*) branchJet->At(3);
				TLorentzVector j1 = jet1->P4(), j2=jet2->P4(), j3=jet3->P4(), j4 = jet4->P4();
				vector<TVector3> vectorCollection = {
				  TVector3(j1.Px(), j1.Py(), j1.Pz()),
				  TVector3(j2.Px(), j2.Py(), j2.Pz()),
				  TVector3(j3.Px(), j3.Py(), j3.Pz()),
				  TVector3(j4.Px(), j4.Py(), j4.Pz())
				};
				TMatrixDSym momentumTensor(3);
				for(std::vector<TVector3>::const_iterator p=vectorCollection.begin(); p!=vectorCollection.end(); p++){
				    for(int k=0; k<3; k++){
				      for(int m=0; m<=k; m++){
					momentumTensor[k][m] += (*p)[k]*(*p)[m];
				      }
				    }
				}
				momentumTensor *= 1/(momentumTensor[0][0] + momentumTensor[1][1] + momentumTensor[2][2]);
				TMatrixDSymEigen eigenSystem(momentumTensor);
				TVectorD eigenvalues = eigenSystem.GetEigenValues();
				float eigenvalue1_ = eigenvalues[0];
				float eigenvalue2_ = eigenvalues[1];
				float eigenvalue3_ = eigenvalues[2];
				sphericity = 1.5*(eigenvalue2_ + eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
				aplanarity = 1.5*(eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
				/////////event shape
				
				
				////////////////////TIM fit
				/////// Beam energy-momentum constrained fit for jetB1,jetB2,jetNB1,jetNB2

		      vector<LorentzVectorWithErrors>* inputLVWE=new vector<LorentzVectorWithErrors>(4);
		      // for now just assign uniform 100% errors
		      inputLVWE->at(0)=LorentzVectorWithErrors(jetB1,jetB1.Energy(),1.,1.,1.);
		      inputLVWE->at(1)=LorentzVectorWithErrors(jetB2,jetB2.Energy(),1.,1.,1.);
		      inputLVWE->at(2)=LorentzVectorWithErrors(jetNB1,jetNB1.Energy(),1.,1.,1.);
		      inputLVWE->at(3)=LorentzVectorWithErrors(jetNB2,jetNB2.Energy(),1.,1.,1.);

		      epConstrainHH eCHH(inputLVWE,Ecm,sigmaEcm,sigmaP,enableExtraTries,nSigVar);
		      vector<LorentzVectorWithErrors>* outputLVWE=eCHH.NumericalMinimization();


		      if(cont2B2NB < 20) { 
			for(int i=0; i<4; i++) {
			  TLorentzVector tLinput=inputLVWE->at(i).getLV();
			  TLorentzVector tLoutput=outputLVWE->at(i).getLV();
			 // cout << " i= " << i << "  input " << " E,px,py,pz= " << tLinput.E() << " " << tLinput.Px()  << " " << tLinput.Py()  << " " << tLinput.Pz() << endl;
			  //cout << " i= " << i << " output " << " E,px,py,pz= " << tLoutput.E() << " " << tLoutput.Px()  << " " << tLoutput.Py()  << " " << tLoutput.Pz() << endl;
			}
		      }

		      jetB1Fit=outputLVWE->at(0).getLV();
		      jetB2Fit=outputLVWE->at(1).getLV();
		      jetNB1Fit=outputLVWE->at(2).getLV();
		      jetNB2Fit=outputLVWE->at(3).getLV();

		      delete inputLVWE;
		      delete outputLVWE;

			////////////////////TIM fit
			
			////////W/O fit
				/*jetB1M=static_cast<float>(jetB1.M());
				jetB2M=static_cast<float>(jetB2.M());
				jetNB1M=static_cast<float>(jetNB1.M());
				jetNB2M=static_cast<float>(jetNB2.M());*/
				int nJetsg;
				if(branchGenJet->GetEntries() > 0)  nJetsg =  branchGenJet->GetEntries();
      				else nJetsg=0;
      				if(nJetsg!=4) contDiffJetsqq++;
      				if(nJetsg==5) fiveJetsqq[cont5Jetsqq]=entry;
      				if(nJetsg==5) cont5Jetsqq++;
      				//if(nJetsg==5) cout<<"entry: "<<entry<<"    njetsg: "<<nJetsg<<"     nJets: "<<nJets<<endl;
      				if(jetNB1M>25) largeBMassqq[contLargeBMassqq]=entry;
				if(jetNB1M>25) contLargeBMassqq++;
				//cout<<"nJetsg: "<<nJetsg<<endl<<"indexb1: "<<indexb1<<endl<<"indexb2: "<<indexb2<<endl<<"indexnb1: "<<indexnb1<<endl<<"indexnb2: "<<indexnb2<<endl<<endl;
				if (nJetsg==4)
				{
					Jet *jetb1 = (Jet*) branchJet->At(indexb1);
					Jet *jetb2 = (Jet*) branchJet->At(indexb2);
					Jet *jetnb1 = (Jet*) branchJet->At(indexnb1);
					Jet *jetnb2 = (Jet*) branchJet->At(indexnb2);
					
					////////cont parts
					int jetFlav=jetb1->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==5) contB1qq++;
					else contNB1qq++;
					if(jetFlav==5 && (jetb1->Mass)>23) contB123qq++;
					if(jetFlav!=5 && (jetb1->Mass)>23) contNB123qq++;
					
					jetFlav=jetb2->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==5) contB2qq++;
					else contNB2qq++;
					if(jetFlav==5 && (jetb2->Mass)>23) contB223qq++;
					if(jetFlav!=5 && (jetb2->Mass)>23) contNB223qq++;
					
					jetFlav=jetnb1->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==21) contG1qq++;
					else contNG1qq++;
					if(jetFlav==21 && (jetnb1->Mass)>23) contG123qq++;
					if(jetFlav!=21 && (jetnb1->Mass)>23) contNG123qq++;
					
					jetFlav=jetnb2->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==21) contG2qq++;
					else contNG2qq++;
					if(jetFlav==21 && (jetnb2->Mass)>23) contG223qq++;
					if(jetFlav!=21 && (jetnb2->Mass)>23) contNG223qq++;
					///////cont parts
					
					double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2, deltaRMin;
					
					////gen jet masses
					for(int j=0;j<nJetsg;j++)
					{
						Jet *genjet = (Jet*) branchGenJet->At(j);
						double genjetEta = genjet->Eta;
     						double genjetPhi = genjet->Phi;
     			
			     			double jetB1Eta = jetb1->Eta;
			     			double jetB1Phi = jetb1->Phi;
			     			deltaRB1 = sqrt(pow(genjetEta-jetB1Eta, 2) + pow(genjetPhi-jetB1Phi, 2));
			     			
			     			double jetB2Eta = jetb2->Eta;
			     			double jetB2Phi = jetb2->Phi;
			     			deltaRB2 = sqrt(pow(genjetEta-jetB2Eta, 2) + pow(genjetPhi-jetB2Phi, 2));
			     			
			     			double jetNB1Eta = jetnb1->Eta;
			     			double jetNB1Phi = jetnb1->Phi;
			     			deltaRNB1 = sqrt(pow(genjetEta-jetNB1Eta, 2) + pow(genjetPhi-jetNB1Phi, 2));
			     			
			     			double jetNB2Eta = jetnb2->Eta;
			     			double jetNB2Phi = jetnb2->Phi;
			     			deltaRNB2 = sqrt(pow(genjetEta-jetNB2Eta, 2) + pow(genjetPhi-jetNB2Phi, 2));
			     			
			     			if(deltaRB1<deltaRB2 && deltaRB1<deltaRNB1 && deltaRB1<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassB1qqg->Fill(genjet->Mass, 0.3062);
			     				indexgenb1=j;
			     			}
			     			else if(deltaRB2<deltaRB1 && deltaRB2<deltaRNB1 && deltaRB2<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassB2qqg->Fill(genjet->Mass, 0.3062);
			     				indexgenb2=j;
			     			}
			     			else if(deltaRNB1<deltaRB2 && deltaRNB1<deltaRB1 && deltaRNB1<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassNB1qqg->Fill(genjet->Mass, 0.3062);
			     				indexgennb1=j;
			     			}
			     			else if(deltaRNB2<deltaRB2 && deltaRNB2<deltaRB1 && deltaRNB2<deltaRNB1)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassNB2qqg->Fill(genjet->Mass, 0.3062);
			     				indexgennb2=j;
			     			}
					} /////gen jet masses

					/*cout<<"jetB1M: "<<jetB1M<<"      b-tag: "<<jetb1->BTag<<"      Flavor: "<<jetb1->Flavor<<endl<<"jetB1Mg: "<<jetB1Mg<<endl;
					cout<<"jetB2M: "<<jetB2M<<"      b-tag: "<<jetb2->BTag<<"      Flavor: "<<jetb2->Flavor<<endl<<"jetB2Mg: "<<jetB2Mg<<endl;
					cout<<"jetNB1M: "<<jetNB1M<<"      b-tag: "<<jetnb1->BTag<<"      Flavor: "<<jetnb1->Flavor<<endl<<"jetNB1Mg: "<<jetNB1Mg<<endl;
					cout<<"jetNB2M: "<<jetNB2M<<"      b-tag: "<<jetnb2->BTag<<"      Flavor: "<<jetnb2->Flavor<<endl<<"jetNB2Mg: "<<jetNB2Mg<<endl<<endl;*/
					
					/////RecoJets Masses
					/*if ((jetb1->PT)>60 && (jetb1->PT)<70) histInvMassB1qq->Fill(jetB1M, 0.3062);
					if ((jetb2->PT)>60 && (jetb2->PT)<70) histInvMassB2qq->Fill(jetB2M, 0.3062);
					if ((jetnb1->PT)>60 && (jetnb1->PT)<70) histInvMassNB1qq->Fill(jetNB1M, 0.3062);
					if ((jetnb2->PT)>60 && (jetnb2->PT)<70) histInvMassNB2qq->Fill(jetNB2M, 0.3062);*/
					
					Jet *genjetb1 = (Jet*) branchGenJet->At(indexgenb1);
					Jet *genjetb2 = (Jet*) branchGenJet->At(indexgenb2);
					Jet *genjetnb1 = (Jet*) branchGenJet->At(indexgennb1);
					Jet *genjetnb2 = (Jet*) branchGenJet->At(indexgennb2);
					
					int nParticlesDR=0;
					if(branchParticle->GetEntries() > 0) nParticlesDR =  branchParticle->GetEntries();
			    		else nParticlesDR=0;
										
					for(int j=0;j<nParticlesDR;j++)
			    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(j); 
			  			int pid=abs(particle->PID);
			  			 if (pid==5 || pid==4)
		  				{
							double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2;
							
							double partEta = particle->Eta;
					     		double partPhi = particle->Phi;
					     		
					     		double genjetB1Eta = genjetb1->Eta;
				     			double genjetB1Phi = genjetb1->Phi;
				     			deltaRB1 = sqrt(pow(partEta-genjetB1Eta, 2) + pow(partPhi-genjetB1Phi, 2));
				     			
				     			double genjetB2Eta = genjetb2->Eta;
				     			double genjetB2Phi = genjetb2->Phi;
				     			deltaRB2 = sqrt(pow(partEta-genjetB2Eta, 2) + pow(partPhi-genjetB2Phi, 2));
				     			
				     			double genjetNB1Eta = genjetnb1->Eta;
				     			double genjetNB1Phi = genjetnb1->Phi;
				     			deltaRNB1 = sqrt(pow(partEta-genjetNB1Eta, 2) + pow(partPhi-genjetNB1Phi, 2));
				     			
				     			double genjetNB2Eta = genjetnb2->Eta;
				     			double genjetNB2Phi = genjetnb2->Phi;
				     			deltaRNB2 = sqrt(pow(partEta-genjetNB2Eta, 2) + pow(partPhi-genjetNB2Phi, 2));
				     			
				     			/*if (deltaRB1<0.4) histB1DR04qq->Fill(deltaRB1, 0.3062);
				     			if (deltaRB2<0.4) histB2DR04qq->Fill(deltaRB2, 0.3062);
				     			if (deltaRNB1<0.4) histNB1DR04qq->Fill(deltaRNB1, 0.3062);
				     			if (deltaRNB2<0.4) histNB2DR04qq->Fill(deltaRNB2, 0.3062);*/
  		 				 }
					}
				}
				jetPairB=jetB1+jetB2;
				jetPairNB=jetNB1+jetNB2;
				jetSum=jetB1+jetB2+jetNB1+jetNB2;
				jetPairBNB1=jetB1+jetNB1;
				jetPairBNB2=jetB2+jetNB2;
				jetPairBNB3=jetB1+jetNB2;
				jetPairBNB4=jetB2+jetNB1;
				minMass = TMath::Min(static_cast<float>(jetPairBNB1.M()), static_cast<float>(jetPairBNB2.M()));
				minMass = TMath::Min(minMass, static_cast<float>(jetPairBNB3.M()));
				minMass = TMath::Min(minMass, static_cast<float>(jetPairBNB4.M()));
				invMassB = jetPairB.M();
				invMassNB = jetPairNB.M();
				////////W/O fit
				/////With fit
				jetB1MFit=static_cast<float>(jetB1Fit.M());
				jetB2MFit=static_cast<float>(jetB2Fit.M());
				jetNB1MFit=static_cast<float>(jetNB1Fit.M());
				jetNB2MFit=static_cast<float>(jetNB2Fit.M());
				jetPairBFit=jetB1Fit+jetB2Fit;
				jetPairNBFit=jetNB1Fit+jetNB2Fit;
				jetSumFit=jetB1Fit+jetB2Fit+jetNB1Fit+jetNB2Fit;
				jetPairBNB1Fit=jetB1Fit+jetNB1Fit;
				jetPairBNB2Fit=jetB2Fit+jetNB2Fit;
				jetPairBNB3Fit=jetB1Fit+jetNB2Fit;
				jetPairBNB4Fit=jetB2Fit+jetNB1Fit;
				minMassFit = TMath::Min(static_cast<float>(jetPairBNB1Fit.M()), static_cast<float>(jetPairBNB2Fit.M()));
				minMassFit = TMath::Min(minMassFit, static_cast<float>(jetPairBNB3Fit.M()));
				minMassFit = TMath::Min(minMassFit, static_cast<float>(jetPairBNB4Fit.M()));
				invMassBFit = jetPairBFit.M();
				invMassNBFit = jetPairNBFit.M();
				///////////With fit
				float invMassSum = jetSum.M();
				float invMassSumFit = jetSumFit.M(); //////////////TIM fit 
				////////////invariant mass
				
				
				/*if(invMassBFit>30 && invMassNBFit>30 && invMassBFit<=380 && invMassNBFit<=380)
				{
					invMassB=invMassBFit;
					invMassNB=invMassNBFit;
					if (minMassFit>0) minMass=minMassFit;
					jetB1M=jetB1MFit;
					jetB2M=jetB2MFit;
					jetNB1M=jetNB1MFit;
					jetNB2M=jetNB2MFit;
				}*/
				
				if(jetB1M<0) jetB1M=0;
				if(jetB2M<0) jetB2M=0;
				if(jetNB1M<0) jetNB1M=0;
				if(jetNB2M<0) jetNB2M=0;
				
				
				if(aplanarity<0) cout<<"Event: "<<entry+1<<endl<<"aplanarity: "<<aplanarity<<endl;
				if(invMassB<0) cout<<"Event: "<<entry+1<<endl<<"invMassB: "<<invMassB<<endl;
				if(invMassNB<0) cout<<"Event: "<<entry+1<<endl<<"invmassNB: "<<invMassNB<<endl;
				if(minMass<0) cout<<"Event: "<<entry+1<<endl<<"minMass: "<<minMass<<endl;
				if(sphericity<0) cout<<"Event: "<<entry+1<<endl<<"sohericity: "<<sphericity<<endl;
				if(ptSum<0) cout<<"Event: "<<entry+1<<endl<<"ptSum: "<<ptSum<<endl;
				if(pt1<0) cout<<"Event: "<<entry+1<<endl<<"pt1: "<<pt1<<endl;
				if(pt2<0) cout<<"Event: "<<entry+1<<endl<<"pt2: "<<pt2<<endl;
				if(pt3<0) cout<<"Event: "<<entry+1<<endl<<"pt3: "<<pt3<<endl;
				if(pt4<0) cout<<"Event: "<<entry+1<<endl<<"pt4: "<<pt4<<endl;
				if(jetB1M<0) cout<<"Event: "<<entry+1<<endl<<"jetB1M: "<<jetB1M<<endl;
				if(jetB2M<0) cout<<"Event: "<<entry+1<<endl<<"jetB2M: "<<jetB2M<<endl;
				if(jetNB1M<0) cout<<"Event: "<<entry+1<<endl<<"jetNB1M: "<<jetNB1M<<endl;
				if(jetNB2M<0) cout<<"Event: "<<entry+1<<endl<<"jetNB2M: "<<jetNB2M<<endl<<endl;
				
				//if (entry!=150983 && entry!=450261 && entry!=956345 && entry!=963429) TreeB.Fill();
				if(invMassB>=400) cout<<entry<<": "<<invMassB<<endl;
				if(invMassNB>=400) cout<<entry<<": "<<invMassNB<<endl<<endl;
				//TreeB.Fill();
		       	}
	  	}
	  }
	  
      }
      
      
      
      
      
      
      
      
      
      
      
      
      
      /////////////////Background tt
   // Get pointers to branches used in this analysis
   branchJet = treeReaderBacktt->UseBranch("Jet");
   branchGenJet = treeReaderBacktt->UseBranch("GenJet");
   branchParticle = treeReaderBacktt->UseBranch("Particle");
   branchEvent = treeReaderBacktt->UseBranch("Event");
   branchElectron = treeReaderBacktt->UseBranch("Electron");
   branchMuon = treeReaderBacktt->UseBranch("Muon");
   
   cout<<endl<<endl<<"For back tt: "<<endl<<endl;
   
   for(Int_t entry=0; entry<numberOfEntriesBacktt; entry++) 
   {   
      treeReaderBacktt->ReadEntry(entry);
          
      if(branchJet->GetEntries() > 0)  nJets =  branchJet->GetEntries();
      else nJets=0;
      
      contBJets=0;
      contNBJets=0;
      contG=0;
      contB=0;
      contEtaCut=0;
      contCosThetaCut=0;
      ptSum=0;
      
      if(nJets==4)
      {
      	  for(int i=0;i<nJets;i++)
	  {
		Jet *jet = (Jet*) branchJet->At(i);
		if(jet->BTag==1)
		{
			contBJets++;
			if(contBJets==1)
			{
				jetB1=jet->P4();
				jetB1M=jet->Mass;
				indexb1=i;
			}	
			else if(contBJets==2)
			{
				jetB2=jet->P4();
				jetB2M=jet->Mass;
				indexb2=i;
			}
		}
		else
		{
			contNBJets++;
			if(contNBJets==1)
			{
				jetNB1=jet->P4();
				jetNB1M=jet->Mass;
				indexnb1=i;
			}
			else if(contNBJets==2)
			{
				jetNB2=jet->P4();
				jetNB2M=jet->Mass;
				indexnb2=i;
			}
		}
		if(contBJets>2 || contNBJets>2) 
		{
			contPlus4Jets++;
			break;
		}
	  }
	  if(contBJets==2 && contNBJets==2)
	  {
			///////lepton veto
			if(branchElectron->GetEntries() > 0) nElectrons =  branchElectron->GetEntries();
	    		else nElectrons=0;
	     		Electron *electron = (Electron*) branchElectron->At(0);
	     		if (nElectrons!=0)
	     		{ 
	     			elecContBranchtt+=branchElectron->GetEntries();
	     			elecPt=electron->PT;
	     		}
	     		
	     		if(branchMuon->GetEntries() > 0) nMuons =  branchMuon->GetEntries();
	    		else nMuons=0;
	     		Muon *muon = (Muon*) branchMuon->At(0);
	     		if (nMuons!=0)
	     		{
	     			muonContBranchtt+=branchMuon->GetEntries();
	     			muonPt=muon->PT;
	     		}
			//////lepton veto
			
			contEventsBBGGBacktt++;
			//cout<<i<<", ";
			for(int i=0;i<nJets;i++)
			{
				Jet *jet = (Jet*) branchJet->At(i);
				
				//////////jet pt
				ptSum+=jet->PT;
				if(i==0)
				{ 
					pt1 = jet->PT;
				}
				else if(i==1)
				{
					pt2 = jet->PT;
				}
				else if(i==2)
				{ 
					pt3 = jet->PT;
				}
				else if(i==3)
				{
					pt4 = jet->PT;
				}
				/////////jet pt
				
				////////angles
				eta = jet->Eta;
				eEta = pow(TMath::E(), -eta);
				theta = 2*TMath::ATan(eEta);
				cosTheta = TMath::Cos(theta);
        			jetEta=abs(jet->Eta);
				if(jetEta<1.5) contEtaCut++;
        			if(cosTheta>-0.91 && cosTheta<0.91) contCosThetaCut++;
        			////////angles
			}
			if(contBJets==2 && contNBJets==2 && contCosThetaCut==4 && pt1>80 && pt2>60 && pt3>15 && pt4>15 && ptSum>240 && elecPt<6 && muonPt<8)
			{
				cont2B2NBBacktt++;
				
				/////////event shape
				Jet *jet1 = (Jet*) branchJet->At(0);
				Jet *jet2 = (Jet*) branchJet->At(1);
				Jet *jet3 = (Jet*) branchJet->At(2);
				Jet *jet4 = (Jet*) branchJet->At(3);
				TLorentzVector j1 = jet1->P4(), j2=jet2->P4(), j3=jet3->P4(), j4 = jet4->P4();
				vector<TVector3> vectorCollection = {
				  TVector3(j1.Px(), j1.Py(), j1.Pz()),
				  TVector3(j2.Px(), j2.Py(), j2.Pz()),
				  TVector3(j3.Px(), j3.Py(), j3.Pz()),
				  TVector3(j4.Px(), j4.Py(), j4.Pz())
				};
				TMatrixDSym momentumTensor(3);
				for(std::vector<TVector3>::const_iterator p=vectorCollection.begin(); p!=vectorCollection.end(); p++){
				    for(int k=0; k<3; k++){
				      for(int m=0; m<=k; m++){
					momentumTensor[k][m] += (*p)[k]*(*p)[m];
				      }
				    }
				}
				momentumTensor *= 1/(momentumTensor[0][0] + momentumTensor[1][1] + momentumTensor[2][2]);
				TMatrixDSymEigen eigenSystem(momentumTensor);
				TVectorD eigenvalues = eigenSystem.GetEigenValues();
				float eigenvalue1_ = eigenvalues[0];
				float eigenvalue2_ = eigenvalues[1];
				float eigenvalue3_ = eigenvalues[2];
				sphericity = 1.5*(eigenvalue2_ + eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
				aplanarity = 1.5*(eigenvalue3_)/(eigenvalue1_ + eigenvalue2_ + eigenvalue3_);
				/////////event shape
				
				
					////////////////////TIM fit
				/////// Beam energy-momentum constrained fit for jetB1,jetB2,jetNB1,jetNB2

		      vector<LorentzVectorWithErrors>* inputLVWE=new vector<LorentzVectorWithErrors>(4);
		      // for now just assign uniform 100% errors
		      inputLVWE->at(0)=LorentzVectorWithErrors(jetB1,jetB1.Energy(),1.,1.,1.);
		      inputLVWE->at(1)=LorentzVectorWithErrors(jetB2,jetB2.Energy(),1.,1.,1.);
		      inputLVWE->at(2)=LorentzVectorWithErrors(jetNB1,jetNB1.Energy(),1.,1.,1.);
		      inputLVWE->at(3)=LorentzVectorWithErrors(jetNB2,jetNB2.Energy(),1.,1.,1.);

		      epConstrainHH eCHH(inputLVWE,Ecm,sigmaEcm,sigmaP,enableExtraTries,nSigVar);
		      vector<LorentzVectorWithErrors>* outputLVWE=eCHH.NumericalMinimization();


		      if(cont2B2NB < 20) { 
			for(int i=0; i<4; i++) {
			  TLorentzVector tLinput=inputLVWE->at(i).getLV();
			  TLorentzVector tLoutput=outputLVWE->at(i).getLV();
			  //cout << " i= " << i << "  input " << " E,px,py,pz= " << tLinput.E() << " " << tLinput.Px()  << " " << tLinput.Py()  << " " << tLinput.Pz() << endl;
			  //cout << " i= " << i << " output " << " E,px,py,pz= " << tLoutput.E() << " " << tLoutput.Px()  << " " << tLoutput.Py()  << " " << tLoutput.Pz() << endl;
			}
		      }

		      jetB1Fit=outputLVWE->at(0).getLV();
		      jetB2Fit=outputLVWE->at(1).getLV();
		      jetNB1Fit=outputLVWE->at(2).getLV();
		      jetNB2Fit=outputLVWE->at(3).getLV();

		      delete inputLVWE;
		      delete outputLVWE;

			////////////////////TIM fit
			
			////////W/O fit
				/*jetB1M=static_cast<float>(jetB1.M());
				jetB2M=static_cast<float>(jetB2.M());
				jetNB1M=static_cast<float>(jetNB1.M());
				jetNB2M=static_cast<float>(jetNB2.M());*/
				int nJetsg;
				if(branchGenJet->GetEntries() > 0)  nJetsg =  branchGenJet->GetEntries();
      				else nJetsg=0;
      				if(nJetsg!=4) contDiffJetstt++;
      				if(nJetsg==5) fiveJetstt[cont5Jetstt]=entry;
      				if(nJetsg==5) cont5Jetstt++;
      				if(jetNB1M>25) largeBMasstt[contLargeBMasstt]=entry;
				if(jetNB1M>25) contLargeBMasstt++;
				//cout<<"nJetsg: "<<nJetsg<<endl<<"indexb1: "<<indexb1<<endl<<"indexb2: "<<indexb2<<endl<<"indexnb1: "<<indexnb1<<endl<<"indexnb2: "<<indexnb2<<endl<<endl;
				if (nJetsg==4)
				{
					Jet *jetb1 = (Jet*) branchJet->At(indexb1);
					Jet *jetb2 = (Jet*) branchJet->At(indexb2);
					Jet *jetnb1 = (Jet*) branchJet->At(indexnb1);
					Jet *jetnb2 = (Jet*) branchJet->At(indexnb2);
					
					////////cont parts
					int jetFlav=jetb1->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==5) contB1tt++;
					else contNB1tt++;
					if(jetFlav==5 && (jetb1->Mass)>23) contB123tt++;
					if(jetFlav!=5 && (jetb1->Mass)>23) contNB123tt++;
					
					jetFlav=jetb2->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==5) contB2tt++;
					else contNB2tt++;
					if(jetFlav==5 && (jetb2->Mass)>23) contB223tt++;
					if(jetFlav!=5 && (jetb2->Mass)>23) contNB223tt++;
					
					jetFlav=jetnb1->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==21) contG1tt++;
					else contNG1tt++;
					if(jetFlav==21 && (jetnb1->Mass)>23) contG123tt++;
					if(jetFlav!=21 && (jetnb1->Mass)>23) contNG123tt++;
					
					jetFlav=jetnb2->Flavor;
					jetFlav=abs(jetFlav);
					if(jetFlav==21) contG2tt++;
					else contNG2tt++;
					if(jetFlav==21 && (jetnb2->Mass)>23) contG223tt++;
					if(jetFlav!=21 && (jetnb2->Mass)>23) contNG223tt++;
					///////cont parts
					
					double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2, deltaRMin;
					
					/////gen jet masses
					for(int j=0;j<nJetsg;j++)
					{
						Jet *genjet = (Jet*) branchGenJet->At(j);
						double genjetEta = genjet->Eta;
     						double genjetPhi = genjet->Phi;
     			
			     			double jetB1Eta = jetb1->Eta;
			     			double jetB1Phi = jetb1->Phi;
			     			deltaRB1 = sqrt(pow(genjetEta-jetB1Eta, 2) + pow(genjetPhi-jetB1Phi, 2));
			     			
			     			double jetB2Eta = jetb2->Eta;
			     			double jetB2Phi = jetb2->Phi;
			     			deltaRB2 = sqrt(pow(genjetEta-jetB2Eta, 2) + pow(genjetPhi-jetB2Phi, 2));
			     			
			     			double jetNB1Eta = jetnb1->Eta;
			     			double jetNB1Phi = jetnb1->Phi;
			     			deltaRNB1 = sqrt(pow(genjetEta-jetNB1Eta, 2) + pow(genjetPhi-jetNB1Phi, 2));
			     			
			     			double jetNB2Eta = jetnb2->Eta;
			     			double jetNB2Phi = jetnb2->Phi;
			     			deltaRNB2 = sqrt(pow(genjetEta-jetNB2Eta, 2) + pow(genjetPhi-jetNB2Phi, 2));
			     			
			     			if(deltaRB1<deltaRB2 && deltaRB1<deltaRNB1 && deltaRB1<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassB1ttg->Fill(genjet->Mass, 3.85109375);
			     				indexgenb1=j;
			     			}
			     			else if(deltaRB2<deltaRB1 && deltaRB2<deltaRNB1 && deltaRB2<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassB2ttg->Fill(genjet->Mass, 3.85109375);
			     				indexgenb2=j;
			     			}
			     			else if(deltaRNB1<deltaRB2 && deltaRNB1<deltaRB1 && deltaRNB1<deltaRNB2)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassNB1ttg->Fill(genjet->Mass, 3.85109375);
			     				indexgennb1=j;
			     			}
			     			else if(deltaRNB2<deltaRB2 && deltaRNB2<deltaRB1 && deltaRNB2<deltaRNB1)
			     			{
			     				//if ((genjet->PT)>60 && (genjet->PT)<70) histInvMassNB2ttg->Fill(genjet->Mass, 3.85109375);
			     				indexgennb2=j;
			     			}
					} /////gen jet masses

					/*cout<<"jetB1M: "<<jetB1M<<"      b-tag: "<<jetb1->BTag<<"      Flavor: "<<jetb1->Flavor<<endl<<"jetB1Mg: "<<jetB1Mg<<endl;
					cout<<"jetB2M: "<<jetB2M<<"      b-tag: "<<jetb2->BTag<<"      Flavor: "<<jetb2->Flavor<<endl<<"jetB2Mg: "<<jetB2Mg<<endl;
					cout<<"jetNB1M: "<<jetNB1M<<"      b-tag: "<<jetnb1->BTag<<"      Flavor: "<<jetnb1->Flavor<<endl<<"jetNB1Mg: "<<jetNB1Mg<<endl;
					cout<<"jetNB2M: "<<jetNB2M<<"      b-tag: "<<jetnb2->BTag<<"      Flavor: "<<jetnb2->Flavor<<endl<<"jetNB2Mg: "<<jetNB2Mg<<endl<<endl;*/
					
					/////RecoJets Masses
					/*if ((jetb1->PT)>60 && (jetb1->PT)<70) histInvMassB1tt->Fill(jetB1M, 3.85109375);
					if ((jetb2->PT)>60 && (jetb2->PT)<70) histInvMassB2tt->Fill(jetB2M, 3.85109375);
					if ((jetnb1->PT)>60 && (jetnb1->PT)<70) histInvMassNB1tt->Fill(jetNB1M, 3.85109375);
					if ((jetnb2->PT)>60 && (jetnb2->PT)<70) histInvMassNB2tt->Fill(jetNB2M, 3.85109375);*/
					
					Jet *genjetb1 = (Jet*) branchGenJet->At(indexgenb1);
					Jet *genjetb2 = (Jet*) branchGenJet->At(indexgenb2);
					Jet *genjetnb1 = (Jet*) branchGenJet->At(indexgennb1);
					Jet *genjetnb2 = (Jet*) branchGenJet->At(indexgennb2);
					
					int nParticlesDR=0;
					if(branchParticle->GetEntries() > 0) nParticlesDR =  branchParticle->GetEntries();
			    		else nParticlesDR=0;
										
					for(int j=0;j<nParticlesDR;j++)
			    		{
			  			GenParticle *particle = (GenParticle*) branchParticle->At(j); 
			  			int pid=abs(particle->PID);
			  			 if (pid==5 || pid==4)
		  				{
							double deltaRB1, deltaRB2, deltaRNB1, deltaRNB2;
							
							double partEta = particle->Eta;
					     		double partPhi = particle->Phi;
					     		
					     		double genjetB1Eta = genjetb1->Eta;
				     			double genjetB1Phi = genjetb1->Phi;
				     			deltaRB1 = sqrt(pow(partEta-genjetB1Eta, 2) + pow(partPhi-genjetB1Phi, 2));
				     			
				     			double genjetB2Eta = genjetb2->Eta;
				     			double genjetB2Phi = genjetb2->Phi;
				     			deltaRB2 = sqrt(pow(partEta-genjetB2Eta, 2) + pow(partPhi-genjetB2Phi, 2));
				     			
				     			double genjetNB1Eta = genjetnb1->Eta;
				     			double genjetNB1Phi = genjetnb1->Phi;
				     			deltaRNB1 = sqrt(pow(partEta-genjetNB1Eta, 2) + pow(partPhi-genjetNB1Phi, 2));
				     			
				     			double genjetNB2Eta = genjetnb2->Eta;
				     			double genjetNB2Phi = genjetnb2->Phi;
				     			deltaRNB2 = sqrt(pow(partEta-genjetNB2Eta, 2) + pow(partPhi-genjetNB2Phi, 2));
				     			
				     			/*if (deltaRB1<0.4) histB1DR04tt->Fill(deltaRB1, 3.85109375);
				     			if (deltaRB2<0.4) histB2DR04tt->Fill(deltaRB2, 3.85109375);
				     			if (deltaRNB1<0.4) histNB1DR04tt->Fill(deltaRNB1, 3.85109375);
				     			if (deltaRNB2<0.4) histNB2DR04tt->Fill(deltaRNB2, 3.85109375);*/
  		 				 }
					}
				}
				
				jetPairB=jetB1+jetB2;
				jetPairNB=jetNB1+jetNB2;
				jetSum=jetB1+jetB2+jetNB1+jetNB2;
				jetPairBNB1=jetB1+jetNB1;
				jetPairBNB2=jetB2+jetNB2;
				jetPairBNB3=jetB1+jetNB2;
				jetPairBNB4=jetB2+jetNB1;
				minMass = TMath::Min(static_cast<float>(jetPairBNB1.M()), static_cast<float>(jetPairBNB2.M()));
				minMass = TMath::Min(minMass, static_cast<float>(jetPairBNB3.M()));
				minMass = TMath::Min(minMass, static_cast<float>(jetPairBNB4.M()));
				invMassB = jetPairB.M();
				invMassNB = jetPairNB.M();
				////////W/O fit
				/////With fit
				jetB1MFit=static_cast<float>(jetB1Fit.M());
				jetB2MFit=static_cast<float>(jetB2Fit.M());
				jetNB1MFit=static_cast<float>(jetNB1Fit.M());
				jetNB2MFit=static_cast<float>(jetNB2Fit.M());
				jetPairBFit=jetB1Fit+jetB2Fit;
				jetPairNBFit=jetNB1Fit+jetNB2Fit;
				jetSumFit=jetB1Fit+jetB2Fit+jetNB1Fit+jetNB2Fit;
				jetPairBNB1Fit=jetB1Fit+jetNB1Fit;
				jetPairBNB2Fit=jetB2Fit+jetNB2Fit;
				jetPairBNB3Fit=jetB1Fit+jetNB2Fit;
				jetPairBNB4Fit=jetB2Fit+jetNB1Fit;
				minMassFit = TMath::Min(static_cast<float>(jetPairBNB1Fit.M()), static_cast<float>(jetPairBNB2Fit.M()));
				minMassFit = TMath::Min(minMassFit, static_cast<float>(jetPairBNB3Fit.M()));
				minMassFit = TMath::Min(minMassFit, static_cast<float>(jetPairBNB4Fit.M()));
				invMassBFit = jetPairBFit.M();
				invMassNBFit = jetPairNBFit.M();
				///////////With fit
				float invMassSum = jetSum.M();
				float invMassSumFit = jetSumFit.M(); //////////////TIM fit 
				////////////invariant mass
				
				
				/*if(invMassBFit>30 && invMassNBFit>30 && invMassBFit<=380 && invMassNBFit<=380)
				{
					invMassB=invMassBFit;
					invMassNB=invMassNBFit;
					if (minMassFit>0) minMass=minMassFit;
					jetB1M=jetB1MFit;
					jetB2M=jetB2MFit;
					jetNB1M=jetNB1MFit;
					jetNB2M=jetNB2MFit;
				}*/
				
				
				if(jetB1M<0) jetB1M=0;
				if(jetB2M<0) jetB2M=0;
				if(jetNB1M<0) jetNB1M=0;
				if(jetNB2M<0) jetNB2M=0;
				////////////invariant mass
				
				
				if(aplanarity<0) cout<<"Event: "<<entry+1<<endl<<"aplanarity: "<<aplanarity<<endl;
				if(invMassB<0) cout<<"Event: "<<entry+1<<endl<<"invMassB: "<<invMassB<<endl;
				if(invMassNB<0) cout<<"Event: "<<entry+1<<endl<<"invmassNB: "<<invMassNB<<endl;
				if(minMass<0) cout<<"Event: "<<entry+1<<endl<<"minMass: "<<minMass<<endl;
				if(sphericity<0) cout<<"Event: "<<entry+1<<endl<<"sohericity: "<<sphericity<<endl;
				if(ptSum<0) cout<<"Event: "<<entry+1<<endl<<"ptSum: "<<ptSum<<endl;
				if(pt1<0) cout<<"Event: "<<entry+1<<endl<<"pt1: "<<pt1<<endl;
				if(pt2<0) cout<<"Event: "<<entry+1<<endl<<"pt2: "<<pt2<<endl;
				if(pt3<0) cout<<"Event: "<<entry+1<<endl<<"pt3: "<<pt3<<endl;
				if(pt4<0) cout<<"Event: "<<entry+1<<endl<<"pt4: "<<pt4<<endl;
				if(jetB1M<0) cout<<"Event: "<<entry+1<<endl<<"jetB1M: "<<jetB1M<<endl;
				if(jetB2M<0) cout<<"Event: "<<entry+1<<endl<<"jetB2M: "<<jetB2M<<endl;
				if(jetNB1M<0) cout<<"Event: "<<entry+1<<endl<<"jetNB1M: "<<jetNB1M<<endl;
				if(jetNB2M<0) cout<<"Event: "<<entry+1<<endl<<"jetNB2M: "<<jetNB2M<<endl<<endl;
				
				/*if (entry==150983 || entry==450260 || entry==956345 || entry==963429)
				{
					cout<<"Entry: "<<entry<<endl<<"aplanarity: "<<aplanarity<<endl;
					cout<<"invMassB: "<<invMassB<<endl;
					cout<<"invmassNB: "<<invMassNB<<endl;
					cout<<"minMass: "<<minMass<<endl;
					cout<<"sphericity: "<<sphericity<<endl;
					cout<<"ptSum: "<<ptSum<<endl;
					cout<<"pt1: "<<pt1<<endl;
					cout<<"pt2: "<<pt2<<endl;
					cout<<"pt3: "<<pt3<<endl;
					cout<<"pt4: "<<pt4<<endl<<endl;
				}*/
				//if (entry!=150983 && entry!=450260 && entry!=956345 && entry!=963429) TreeB.Fill();
				//if (entry!=150983 && entry!=956345 && entry!=963429) TreeB.Fill();
				//if (entry!=963429) TreeB.Fill();
				
				TreeB.Fill();
				
				/*if(conttt<2788)
				{
					TreeB.Fill();
					conttt++;
				}*/
		       	}
	  	}
	  }
	  
      }
      
      
     
   
   
   
   /* TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    TCanvas *c3 = new TCanvas();
    TCanvas *c4 = new TCanvas();
    TCanvas *c5 = new TCanvas();
    TCanvas *c6 = new TCanvas();
    TCanvas *c7 = new TCanvas();
    TCanvas *c8 = new TCanvas();
    TCanvas *c9 = new TCanvas();
    TCanvas *c10 = new TCanvas();
    TCanvas *c11 = new TCanvas();
    TCanvas *c12 = new TCanvas();
    TCanvas *c13 = new TCanvas();
    TCanvas *c14 = new TCanvas();
    TCanvas *c15 = new TCanvas();
    TCanvas *c16 = new TCanvas();
    TCanvas *c17 = new TCanvas();
    TCanvas *c18 = new TCanvas();
    TCanvas *c19 = new TCanvas();
    TCanvas *c20 = new TCanvas();
    TCanvas *c21 = new TCanvas();
    TCanvas *c22 = new TCanvas();
    TCanvas *c23 = new TCanvas();
    TCanvas *c24 = new TCanvas();
    TCanvas *c25 = new TCanvas();
    TCanvas *c26 = new TCanvas();
    TCanvas *c27 = new TCanvas();
    TCanvas *c28 = new TCanvas();
    TCanvas *c29 = new TCanvas();
    TCanvas *c30 = new TCanvas();
    TCanvas *c31 = new TCanvas();
    TCanvas *c32 = new TCanvas();
    TCanvas *c33 = new TCanvas();
    TCanvas *c34 = new TCanvas();
    TCanvas *c35 = new TCanvas();
    TCanvas *c36 = new TCanvas();
    
    c1->cd();
    histInvMassB1->Draw("HIST");
    
    c2->cd();
    histInvMassB1qq->Draw("HIST");
    
    c3->cd();
    histInvMassB1tt->Draw("HIST");
    
    c4->cd();
    histInvMassB2->Draw("HIST");
    
    c5->cd();
    histInvMassB2qq->Draw("HIST");
    
    c6->cd();
    histInvMassB2tt->Draw("HIST");
    
    c7->cd();
    histInvMassNB1->Draw("HIST");
    
    c8->cd();
    histInvMassNB1qq->Draw("HIST");
    
    c9->cd();
    histInvMassNB1tt->Draw("HIST");
   
    c10->cd();
    histInvMassNB2->Draw("HIST");
    
    c11->cd();
    histInvMassNB2qq->Draw("HIST");
    
    c12->cd();
    histInvMassNB2tt->Draw("HIST");
    
    c13->cd();
    histInvMassB1g->Draw("HIST");
    
    c14->cd();
    histInvMassB1qqg->Draw("HIST");
    
    c15->cd();
    histInvMassB1ttg->Draw("HIST");
    
    c16->cd();
    histInvMassB2g->Draw("HIST");
    
    c17->cd();
    histInvMassB2qqg->Draw("HIST");
    
    c18->cd();
    histInvMassB2ttg->Draw("HIST");
    
    c19->cd();
    histInvMassNB1g->Draw("HIST");
    
    c20->cd();
    histInvMassNB1qqg->Draw("HIST");
    
    c21->cd();
    histInvMassNB1ttg->Draw("HIST");
   
    c22->cd();
    histInvMassNB2g->Draw("HIST");
    
    c23->cd();
    histInvMassNB2qqg->Draw("HIST");
    
    c24->cd();
    histInvMassNB2ttg->Draw("HIST");
    
    c25->cd();
    histB1DR04->Draw("HIST");
    
    c26->cd();
    histB1DR04qq->Draw("HIST");
    
    c27->cd();
    histB1DR04tt->Draw("HIST");
    
    c28->cd();
    histB2DR04->Draw("HIST");
    
    c29->cd();
    histB2DR04qq->Draw("HIST");
    
    c30->cd();
    histB2DR04tt->Draw("HIST");
    
    c31->cd();
    histNB1DR04->Draw("HIST");
    
    c32->cd();
    histNB1DR04qq->Draw("HIST");
    
    c33->cd();
    histNB1DR04tt->Draw("HIST");
    
    c34->cd();
    histNB2DR04->Draw("HIST");
    
    c35->cd();
    histNB2DR04qq->Draw("HIST");
    
    c36->cd();
    histNB2DR04tt->Draw("HIST");*/
    
    
   
   
   
   
   
   
   
   
   
   
   cout<<endl<<"HH to bbgg events: "<<contEventsHH<<endl;
   cout<<"HH events that pass the preliminary filter: "<<cont2B2NB<<endl;
   cout<<endl<<"Back to bbgg events: "<<contEventsBBGGBack<<endl;
   cout<<"Back events that pass the preliminary filter: "<<cont2B2NBBack<<endl;
   cout<<endl<<"Back tt to bbgg events: "<<contEventsBBGGBacktt<<endl;
   cout<<"Back tt events that pass the preliminary filter: "<<cont2B2NBBacktt<<endl<<endl;
   
   cout<<"Events where there are 4 jets but not 4 genJets (HH): "<<contDiffJets<<endl<<"Events where there are 4 jets but not 4 genJets (qqbar): "<<contDiffJetsqq<<endl<<"Events where there are 4 jets but not 4 genJets (ttbar): "<<contDiffJetstt<<endl<<endl;
    cout<<"Events where there are 4 jets but 5 genJets (HH): "<<cont5Jets<<endl<<"Events where there are 4 jets but 5 genJets (qqbar): "<<cont5Jetsqq<<endl<<"Events where there are 4 jets but 5 genJets (ttbar): "<<cont5Jetstt<<endl;
   
   
   cout<<endl<<"Signal events with 5 jets: ";
   for(int i=0;i<cont5Jets;i++)
   {
   	cout<<fiveJets[i]<<", ";
   }
   cout<<endl<<endl<<"qqbar events with 5 jets: ";
   for(int i=0;i<cont5Jetsqq;i++)
   {
   	cout<<fiveJetsqq[i]<<", ";
   }
   cout<<endl<<endl<<"ttbar events with 5 jets: ";
   for(int i=0;i<cont5Jetstt;i++)
   {
   	cout<<fiveJetstt[i]<<", ";
   }
   cout<<endl<<endl<<"Signal events with jetNB1M>25: ";
   for(int i=0;i<contLargeBMass;i++)
   {
   	cout<<largeBMass[i]<<", ";
   }
   cout<<endl<<endl<<"qqbar events with jetNB1M>25: ";
   for(int i=0;i<contLargeBMassqq;i++)
   {
   	cout<<largeBMassqq[i]<<", ";
   }
   cout<<endl<<endl<<"ttbar events with jetNB1M>25: ";
   for(int i=0;i<contLargeBMasstt;i++)
   {
   	cout<<largeBMasstt[i]<<", ";
   }
   
   cout<<endl;
   
  cout<<"Signal electrons: "<<elecContBranch<<endl; 
  cout<<"qq electrons: "<<elecContBranchqq<<endl;
  cout<<"tt electrons: "<<elecContBranchtt<<endl<<endl; 
  
  cout<<"Signal muons: "<<muonContBranch<<endl; 
  cout<<"qq muons: "<<muonContBranchqq<<endl;
  cout<<"tt muons: "<<muonContBranchtt<<endl<<endl; 
  
  cout<<"contB1HH: "<<contB1HH<<endl;
  cout<<"contNB1HH: "<<contNB1HH<<endl;
  cout<<"contB123HH: "<<contB123HH<<endl;
  cout<<"contNB123HH: "<<contNB123HH<<endl<<endl;
  
  cout<<"contB2HH: "<<contB2HH<<endl;
  cout<<"contNB2HH: "<<contNB2HH<<endl;
  cout<<"contB223HH: "<<contB223HH<<endl;
  cout<<"contNB223HH: "<<contNB223HH<<endl<<endl;
  
  cout<<"contG1HH: "<<contG1HH<<endl;
  cout<<"contNG1HH: "<<contNG1HH<<endl;
  cout<<"contG123HH: "<<contG123HH<<endl;
  cout<<"contNG123HH: "<<contNG123HH<<endl<<endl;
  
  cout<<"contG2HH: "<<contG2HH<<endl;
  cout<<"contNG2HH: "<<contNG2HH<<endl;
  cout<<"contG223HH: "<<contG223HH<<endl;
  cout<<"contNG223HH: "<<contNG223HH<<endl<<endl<<endl;
  
  cout<<"contB1qq: "<<contB1qq<<endl;
  cout<<"contNB1qq: "<<contNB1qq<<endl;
  cout<<"contB123qq: "<<contB123qq<<endl;
  cout<<"contNB123qq: "<<contNB123qq<<endl<<endl;
  
  cout<<"contB2qq: "<<contB2qq<<endl;
  cout<<"contNB2qq: "<<contNB2qq<<endl;
  cout<<"contB223qq: "<<contB223qq<<endl;
  cout<<"contNB223qq: "<<contNB223qq<<endl<<endl;
  
  cout<<"contG1qq: "<<contG1qq<<endl;
  cout<<"contNG1qq: "<<contNG1qq<<endl;
  cout<<"contG123qq: "<<contG123qq<<endl;
  cout<<"contNG123qq: "<<contNG123qq<<endl<<endl;
  
  cout<<"contG2qq: "<<contG2qq<<endl;
  cout<<"contNG2qq: "<<contNG2qq<<endl;
  cout<<"contG223qq: "<<contG223qq<<endl;
  cout<<"contNG223qq: "<<contNG223qq<<endl<<endl<<endl;
  
  cout<<"contB1tt: "<<contB1tt<<endl;
  cout<<"contNB1tt: "<<contNB1tt<<endl;
  cout<<"contB123tt: "<<contB123tt<<endl;
  cout<<"contNB123tt: "<<contNB123tt<<endl<<endl;
  
  cout<<"contB2tt: "<<contB2tt<<endl;
  cout<<"contNB2tt: "<<contNB2tt<<endl;
  cout<<"contB223tt: "<<contB223tt<<endl;
  cout<<"contNB223tt: "<<contNB223tt<<endl<<endl;
  
  cout<<"contG1tt: "<<contG1tt<<endl;
  cout<<"contNG1tt: "<<contNG1tt<<endl;
  cout<<"contG123tt: "<<contG123tt<<endl;
  cout<<"contNG123tt: "<<contNG123tt<<endl<<endl;
  
  cout<<"contG2tt: "<<contG2tt<<endl;
  cout<<"contNG2tt: "<<contNG2tt<<endl;
  cout<<"contG223tt: "<<contG223tt<<endl;
  cout<<"contNG223tt: "<<contNG223tt<<endl<<endl<<endl;
    
    
    
    
   
   // save the Tree heade; the file will be automatically closed
   // when going out of the function scope
   TreeS.Write();
   TreeB.Write();
   
	  
	  
	  
}







