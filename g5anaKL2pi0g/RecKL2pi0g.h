#ifndef RECKLPI0GG_H
#define RECKLPI0GG_H

#include <vector>
#include <list>

#include "CLHEP/Vector/ThreeVector.h"

#include "gamma/Gamma.h"
#include "pi0/Pi0.h"

#include "g5anaKL2pi0g/KL2pi0g.h"

class RecKL2pi0g
{
 public:
  RecKL2pi0g();
  ~RecKL2pi0g();

  std::vector<KL2pi0g> recK2pi0g( const std::list<Gamma> &glist, int userFlag = 20180601 ) const;
  std::vector<KL2pi0g> recK2pi0g( const std::list<Pi0>& pi0list, const std::list<Gamma> &glist, int userFlag = 20180601, double pi0sig2cut=40000. ) const;
  Gamma GetResGamma( const Pi0& p1, const Pi0& p2, const std::list<Gamma> &glist ) const;

  void SetKlongTime( KL2pi0g &kl ) const;
  void SetDebugLevel( Int_t debuglevel ) { m_debugLevel = debuglevel; } 
  
 

 private:
  Int_t    m_debugLevel;


};


#endif
