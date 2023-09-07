#ifndef EPCONSTRAINHH
#define EPCONSTRAINHH

#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include "TMatrixD.h"
#include "TString.h"

#include "Math/Minimizer.h"
#include "LorentzVectorWithErrors.h"



using namespace std;

class epConstrainHH {

public:
  //
  // Constructor and destructor
  //
  epConstrainHH(vector<LorentzVectorWithErrors>* mLV, double ecm, double sigmaEcm, double sigmaP, bool enableExtraTries, double nSigVar);
  ~epConstrainHH(); 

  //
  //
  // Minimization and chi squared
  //
  vector<LorentzVectorWithErrors>* NumericalMinimization();
    

private:   

  double ChiSquared(const double* xx);
  void SetInitialFitParameters( ROOT::Math::Minimizer* min, double* variableValuesBest, double factorSig12Bx=0, double factorSig34Bx=0, double factorSig12By=0, double factorSig34By=0,
					     double factorSig12Bz=0, double factorSig34Bz=0);
  vector<LorentzVectorWithErrors>* m_mLV;
  int m_num_vars;
  double* m_xx0;
  double* m_xx;
  double* m_minVal;
  double* m_maxVal;
  double* m_step;
  TString* m_name;
  double m_stepFraction;
  double* m_yy0;
  double* m_sig0;
  const double* m_xe;
  double* m_cov_matrix;
  int m_number_freeparam;
  double m_ecm;
  double m_sigmaEcm;
  double m_sigmaP;
  bool m_enableExtraTries;
  double m_nSigVar;
  double m_tolerance_factor;
  TMatrixD m_bM;
  TMatrixD m_energyM;
  TMatrixD m_sqrtsM;
  vector<map<TString,int>*> m_vecStringInd;
  double m_chi2;
  int m_status;
  
};

#endif
