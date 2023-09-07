#ifndef LORENTZVECTORWITHERRORS
#define LORENTZVECTORWITHERRORS

#include "TLorentzVector.h"
using namespace std;


class LorentzVectorWithErrors {
public:

  LorentzVectorWithErrors() : m_LV(TLorentzVector()),
			      m_sigEnergy(0.),
			      m_sigBetaX(0.), 
			      m_sigBetaY(0.), 
			      m_sigBetaZ(0.) {}

  LorentzVectorWithErrors(const TLorentzVector& LV, Double_t sigEnergy,
			  Double_t sigBetaX, Double_t sigBetaY, Double_t sigBetaZ) : m_LV(LV),
										     m_sigEnergy(sigEnergy),
										     m_sigBetaX(sigBetaX), 
										     m_sigBetaY(sigBetaY), 
										     m_sigBetaZ(sigBetaZ) {}

  TLorentzVector getLV() {return m_LV;}
  Double_t getSigEnergy() {return m_sigEnergy;}
  Double_t getSigBetaX() {return m_sigBetaX;}
  Double_t getSigBetaY() {return m_sigBetaY;}
  Double_t getSigBetaZ() {return m_sigBetaZ;}


protected:

  
  TLorentzVector m_LV;
  Double_t m_sigEnergy;
  Double_t m_sigBetaX;
  Double_t m_sigBetaY;
  Double_t m_sigBetaZ;


};

#endif
