#include "epConstrainHH.h"
#include "Fit/ParameterSettings.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
//------------------------------------------------------------------------------------------------
epConstrainHH::epConstrainHH(vector<LorentzVectorWithErrors>* mLV, double ecm, double sigmaEcm, double sigmaP, bool enableExtraTries, double nSigVar)
  : m_mLV(mLV),m_num_vars(16),m_xx0(new double[m_num_vars]),m_xx(new double[m_num_vars]),m_minVal(new double[m_num_vars]),m_maxVal(new double[m_num_vars]),
    m_step(new double[m_num_vars]),m_name(new TString[m_num_vars]),m_stepFraction(1.e-4),
    m_yy0(new double[m_num_vars]),m_sig0(new double[m_num_vars]),m_xe(new double[m_num_vars]),m_cov_matrix(0),m_number_freeparam(0),
    m_ecm(ecm),m_sigmaEcm(sigmaEcm),m_sigmaP(sigmaP),m_enableExtraTries(enableExtraTries),m_nSigVar(nSigVar),m_tolerance_factor(1.),
    m_bM(4,4),m_energyM(4,1),m_sqrtsM(4,1),m_vecStringInd(4),m_chi2(0.),m_status(-1)
{
  for(int i=0; i<4; i++) {
    m_vecStringInd[i]= new map<TString,int>();
    m_vecStringInd[i]->insert(std::pair<TString,int>("Energy",4*i));
    m_vecStringInd[i]->insert(std::pair<TString,int>("BetaX",4*i+1));
    m_vecStringInd[i]->insert(std::pair<TString,int>("BetaY",4*i+2));
    m_vecStringInd[i]->insert(std::pair<TString,int>("BetaZ",4*i+3));
    double energy=m_mLV->at(i).getLV().E();
    m_bM(0,i)=1;
    m_bM(1,i)=m_mLV->at(i).getLV().Px()/energy;
    m_bM(2,i)=m_mLV->at(i).getLV().Py()/energy;
    m_bM(3,i)=m_mLV->at(i).getLV().Pz()/energy;
  }
  m_sqrtsM(0,0)=m_ecm;
  m_sqrtsM(1,0)=0;
  m_sqrtsM(2,0)=0;
  m_sqrtsM(3,0)=0;
  m_energyM=TMatrixD(TMatrixD(TMatrixD::kInverted,m_bM),TMatrixD::kMult,m_sqrtsM);
  TMatrixD testConstraint(m_bM,TMatrixD::kMult,m_energyM);
  for(int i=0; i<4; i++) {
    m_xx0[m_vecStringInd[i]->at("Energy")]=m_energyM(i,0);
    m_xx0[m_vecStringInd[i]->at("BetaX")]=m_bM(1,i);
    m_xx0[m_vecStringInd[i]->at("BetaY")]=m_bM(2,i);
    m_xx0[m_vecStringInd[i]->at("BetaZ")]=m_bM(3,i);
    
    m_yy0[m_vecStringInd[i]->at("Energy")]=m_mLV->at(i).getLV().E();
    m_yy0[m_vecStringInd[i]->at("BetaX")]=m_bM(1,i);
    m_yy0[m_vecStringInd[i]->at("BetaY")]=m_bM(2,i);
    m_yy0[m_vecStringInd[i]->at("BetaZ")]=m_bM(3,i);
    
    m_sig0[m_vecStringInd[i]->at("Energy")]=m_mLV->at(i).getSigEnergy();
    m_sig0[m_vecStringInd[i]->at("BetaX")]=m_mLV->at(i).getSigBetaX();
    m_sig0[m_vecStringInd[i]->at("BetaY")]=m_mLV->at(i).getSigBetaY();
    m_sig0[m_vecStringInd[i]->at("BetaZ")]=m_mLV->at(i).getSigBetaZ();
    
    m_minVal[m_vecStringInd[i]->at("Energy")]=1.e-4;
    m_minVal[m_vecStringInd[i]->at("BetaX")]=-1.;
    m_minVal[m_vecStringInd[i]->at("BetaY")]=-1.;
    m_minVal[m_vecStringInd[i]->at("BetaZ")]=-1.;
    
    m_maxVal[m_vecStringInd[i]->at("Energy")]=0.5*m_ecm;
    m_maxVal[m_vecStringInd[i]->at("BetaX")]=1.;
    m_maxVal[m_vecStringInd[i]->at("BetaY")]=1.;
    m_maxVal[m_vecStringInd[i]->at("BetaZ")]=1.;
    
    m_step[m_vecStringInd[i]->at("Energy")]=m_stepFraction*m_sig0[m_vecStringInd[i]->at("Energy")];
    m_step[m_vecStringInd[i]->at("BetaX")]=m_stepFraction*m_sig0[m_vecStringInd[i]->at("BetaX")];
    m_step[m_vecStringInd[i]->at("BetaY")]=m_stepFraction*m_sig0[m_vecStringInd[i]->at("BetaY")];
    m_step[m_vecStringInd[i]->at("BetaZ")]=m_stepFraction*m_sig0[m_vecStringInd[i]->at("BetaZ")];
    
    m_name[m_vecStringInd[i]->at("Energy")]=TString("Energy")+TString::Format("%d",i+1);
    m_name[m_vecStringInd[i]->at("BetaX")]=TString("BetaX")+TString::Format("%d",i+1);
    m_name[m_vecStringInd[i]->at("BetaY")]=TString("BetaY")+TString::Format("%d",i+1);
    m_name[m_vecStringInd[i]->at("BetaZ")]=TString("BetaZ")+TString::Format("%d",i+1);
    
  }
}

//------------------------------------------------------------------------------------------------
epConstrainHH::~epConstrainHH()
{
  //  cout << " ~epConstrainHH point 000 " << endl;
}


//------------------------------------------------------------------------------------------------
vector<LorentzVectorWithErrors>* epConstrainHH::NumericalMinimization()
{
  
  
  vector<LorentzVectorWithErrors>* result=new vector<LorentzVectorWithErrors>(4);
  
  for(int i=0; i<4; i++) {
    Double_t ee=m_xx0[m_vecStringInd[i]->at("Energy")];
    Double_t px=ee*m_xx0[m_vecStringInd[i]->at("BetaX")];
    Double_t py=ee*m_xx0[m_vecStringInd[i]->at("BetaY")];
    Double_t pz=ee*m_xx0[m_vecStringInd[i]->at("BetaZ")];
    TLorentzVector tL(px,py,pz,ee);
    Double_t sigEnergy=m_sig0[m_vecStringInd[i]->at("Energy")];
    Double_t sigBetaX=m_sig0[m_vecStringInd[i]->at("BetaX")];
    Double_t sigBetaY=m_sig0[m_vecStringInd[i]->at("BetaY")];
    Double_t sigBetaZ=m_sig0[m_vecStringInd[i]->at("BetaZ")];
    result->at(i)=LorentzVectorWithErrors(tL,sigEnergy,sigBetaX,sigBetaY,sigBetaZ);
  } 

  return result;
  
}


//------------------------------------------------------------------------------------------------
double epConstrainHH::ChiSquared( const double* xx ) 
{
  double result(0.); 

  for(int i=0; i<m_num_vars; i++) {
    result += pow((xx[i]-m_yy0[i])/m_sig0[i],2);
  }
  
  for(int i=0; i<4; i++) {
    m_energyM(i,0)=xx[m_vecStringInd[i]->at("Energy")];
    m_bM(0,i)=1;
    m_bM(1,i)=xx[m_vecStringInd[i]->at("BetaX")];
    m_bM(2,i)=xx[m_vecStringInd[i]->at("BetaY")];
    m_bM(3,i)=xx[m_vecStringInd[i]->at("BetaZ")];
  }
  m_sqrtsM=TMatrixD(m_bM,TMatrixD::kMult,m_energyM);
  result += pow((m_sqrtsM(0,0)-m_ecm)/m_sigmaEcm,2);
  result += pow(m_sqrtsM(1,0)/m_sigmaP,2);
  result += pow(m_sqrtsM(2,0)/m_sigmaP,2);
  result += pow(m_sqrtsM(3,0)/m_sigmaP,2);

  m_chi2 = result;

  return result;
}
void epConstrainHH::SetInitialFitParameters( ROOT::Math::Minimizer* min, double* variableValuesBest, double factorSig12Bx, double factorSig34Bx, double factorSig12By, double factorSig34By,
					     double factorSig12Bz, double factorSig34Bz)
{



  double stepFraction=1.e-4;
  

  if(variableValuesBest) {
    m_bM(0,0)=1;
    m_bM(1,0)=variableValuesBest[m_vecStringInd[0]->at("BetaX")];
    m_bM(2,0)=variableValuesBest[m_vecStringInd[0]->at("BetaY")];
    m_bM(3,0)=variableValuesBest[m_vecStringInd[0]->at("BetaZ")];
    m_bM(0,1)=1;
    m_bM(1,1)=variableValuesBest[m_vecStringInd[1]->at("BetaX")];
    m_bM(2,1)=variableValuesBest[m_vecStringInd[1]->at("BetaY")];
    m_bM(3,1)=variableValuesBest[m_vecStringInd[1]->at("BetaZ")];
  
    m_bM(0,2)=1;
    m_bM(1,2)=variableValuesBest[m_vecStringInd[2]->at("BetaX")];
    m_bM(2,2)=variableValuesBest[m_vecStringInd[2]->at("BetaY")];
    m_bM(3,2)=variableValuesBest[m_vecStringInd[2]->at("BetaZ")];
    m_bM(0,3)=1;
    m_bM(1,3)=variableValuesBest[m_vecStringInd[3]->at("BetaX")];
    m_bM(2,3)=variableValuesBest[m_vecStringInd[3]->at("BetaY")];
    m_bM(3,3)=variableValuesBest[m_vecStringInd[3]->at("BetaZ")];
  }
  else {
    m_bM(0,0)=1;
    m_bM(1,0)=m_xx0[m_vecStringInd[0]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[0]->at("BetaX")]*(factorSig12Bx-1.);
    m_bM(2,0)=m_xx0[m_vecStringInd[0]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[0]->at("BetaY")]*(factorSig12By-1.);
    m_bM(3,0)=m_xx0[m_vecStringInd[0]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[0]->at("BetaZ")]*(factorSig12Bz-1.);
    m_bM(0,1)=1;
    m_bM(1,1)=m_xx0[m_vecStringInd[1]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[1]->at("BetaX")]*(factorSig12Bx-1.);
    m_bM(2,1)=m_xx0[m_vecStringInd[1]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[1]->at("BetaY")]*(factorSig12By-1.);
    m_bM(3,1)=m_xx0[m_vecStringInd[1]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[1]->at("BetaZ")]*(factorSig12Bz-1.);
  
    m_bM(0,2)=1;
    m_bM(1,2)=m_xx0[m_vecStringInd[2]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[2]->at("BetaX")]*(factorSig34Bx-1.);
    m_bM(2,2)=m_xx0[m_vecStringInd[2]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[2]->at("BetaY")]*(factorSig34By-1.);
    m_bM(3,2)=m_xx0[m_vecStringInd[2]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[2]->at("BetaZ")]*(factorSig34Bz-1.);
    m_bM(0,3)=1;
    m_bM(1,3)=m_xx0[m_vecStringInd[3]->at("BetaX")]+m_nSigVar*m_sig0[m_vecStringInd[3]->at("BetaX")]*(factorSig34Bx-1.);
    m_bM(2,3)=m_xx0[m_vecStringInd[3]->at("BetaY")]+m_nSigVar*m_sig0[m_vecStringInd[3]->at("BetaY")]*(factorSig34By-1.);
    m_bM(3,3)=m_xx0[m_vecStringInd[3]->at("BetaZ")]+m_nSigVar*m_sig0[m_vecStringInd[3]->at("BetaZ")]*(factorSig34Bz-1.);
  }
  
    

  m_energyM=TMatrixD(TMatrixD(TMatrixD::kInverted,m_bM),TMatrixD::kMult,m_sqrtsM);
  for(int i=0; i<4; i++) {
    m_xx[m_vecStringInd[i]->at("Energy")]=m_energyM(i,0);
    m_xx[m_vecStringInd[i]->at("BetaX")]=m_bM(1,i);
    m_xx[m_vecStringInd[i]->at("BetaY")]=m_bM(2,i);
    m_xx[m_vecStringInd[i]->at("BetaZ")]=m_bM(3,i);
  }

  for(int i=0; i<m_num_vars; i++) {
    min -> SetLimitedVariable( i, m_name[i].Data(), m_xx[i], m_step[i], m_minVal[i], m_maxVal[i] );
  }
      

}




