#include <iostream>
#include <iomanip>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "rec2g/Rec2g.h"
#include "pi0/Pi0.h"
#include "MTAnalysisLibrary/MTBasicParameters.h"
#include "MTAnalysisLibrary/MTFunction.h"

#include "g5anaKL2pi0g/RecKL2pi0g.h"

RecKL2pi0g::RecKL2pi0g()
   : m_debugLevel(0)
{
   ;
}

RecKL2pi0g::~RecKL2pi0g()
{
   ;
}


Gamma RecKL2pi0g::GetResGamma( const Pi0& p1, const Pi0& p2, const std::list<Gamma> &glist ) const
{
  Gamma g5;
  for(std::list<Gamma>::const_iterator ig=glist.begin(); ig!=glist.end(); ig++){
    if(ig->id()==p1.g1().id() || ig->id()==p1.g2().id()) continue;
    if(ig->id()==p2.g1().id() || ig->id()==p2.g2().id()) continue;
    g5 = *ig;
    break;
  }
    return g5;
}

std::vector<KL2pi0g>
RecKL2pi0g::recK2pi0g( const std::list<Gamma> &glist, int userFlag ) const
{

  Rec2g rec2g; // pi0 reconstruction
  const std::list<Pi0> pi0list = rec2g.recPi0withConstM( glist );  

  return recK2pi0g(pi0list, glist, userFlag); 
}

std::vector<KL2pi0g>
RecKL2pi0g::recK2pi0g( const std::list<Pi0>& pi0list, const std::list<Gamma>& glist, int userFlag, double pi0sig2cut ) const
{

  std::vector<KL2pi0g> klvec;
  for( std::list<Pi0>::const_iterator p1=pi0list.begin(); p1!=pi0list.end(); p1++ ) {
    for( std::list<Pi0>::const_iterator p2=p1; p2!=pi0list.end(); p2++ ) {
      if ( p2->g1().id() == p1->g1().id() || p2->g1().id() == p1->g2().id() ||
           p2->g2().id() == p1->g1().id() || p2->g2().id() == p1->g2().id()  ){
        continue; // same gamma.
      }
      else {
    //comment out by Qisen 2020/05/04
        if( //p1->status() == 1 &&
        //p2->status() == 1 &&
          p1->recZsig2() <= pi0sig2cut &&
          p2->recZsig2() <= pi0sig2cut ) {

            KL2pi0g kl2pi0g;
            Gamma g5 = GetResGamma(*p1,*p2,glist);
            kl2pi0g.SetPi0( (*p1),(*p2) );
            kl2pi0g.SetGamma(g5); 
            //std::cout << "level: " << m_debugLevel << std::endl;
            if( p1->id()==0 && p2->id()==5 && m_debugLevel==1){
              std::cout << "debug"                    << std::endl;
              kl2pi0g.SetDebugLevel(m_debugLevel);
/*
              std::cout << kl2pi0g.pi0()[0].g1().id() << std::endl;
              std::cout << kl2pi0g.pi0()[0].g2().id() << std::endl;
              std::cout << kl2pi0g.pi0()[1].g1().id() << std::endl;
              std::cout << kl2pi0g.pi0()[1].g2().id() << std::endl;
              std::cout << kl2pi0g.g5().id() << std::endl;
*/
/*
              std::cout << kl2pi0g.pi0()[0].g1() << std::endl;
              std::cout << kl2pi0g.pi0()[0].g2() << std::endl;
              std::cout << kl2pi0g.pi0()[1].g1() << std::endl;
              std::cout << kl2pi0g.pi0()[1].g2() << std::endl;
              std::cout << kl2pi0g.g5() << std::endl;
*/
            }
            Int_t status = kl2pi0g.UpdateVars(); 
            if(status!=1) continue;
            kl2pi0g.SetSortFlag(KL2pi0g::SORT_BY_CHISQ);
            kl2pi0g.SetUserFlag(userFlag);
            SetKlongTime(kl2pi0g);        

            klvec.push_back( kl2pi0g );
        }
      }
    }
  }

  /// sorting
  std::sort( klvec.begin(), klvec.end() );
  Int_t kl_id = 0;
  for(std::vector<KL2pi0g>::iterator it=klvec.begin(); it!=klvec.end(); ++it )
    it->SetId( kl_id++ );

  return klvec; 
}

void RecKL2pi0g::SetKlongTime( KL2pi0g &kl ) const
{
  std::vector<Gamma> gvec;

  for(std::vector<Pi0>::const_iterator ip=kl.pi0().begin();ip!=kl.pi0().end();ip++){
    gvec.push_back(ip->g1());
    gvec.push_back(ip->g2());
  } 

  gvec.push_back(kl.g5() );
  if(gvec.size()!=5){
    std::cout << "error: kl size mismtach" << std::endl;
    return; 
  }

  std::vector<Double_t> tvec, tvec_csi, wvec;

  for(std::vector<Gamma>::const_iterator ig = gvec.begin(); ig!=gvec.end(); ++ig )
  {
    Double_t t, t_csi, tsigma;
    Double_t distance =   TMath::Power(ig->pos().x()-kl.v().x(),2)
                        + TMath::Power(ig->pos().y()-kl.v().y(),2)
                        + TMath::Power(ig->pos().z()-kl.v().z(),2);
    distance = TMath::Sqrt( distance );
    Double_t tof = distance / (TMath::C()/1.E6);

    t     = ig->t() - tof;
    t_csi = ig->t();

    tsigma = MTFUNC::TResolutionOfCluster( ig->e() );

    tvec.push_back(t);
    tvec_csi.push_back(t_csi);

    wvec.push_back(1./tsigma/tsigma);
  }

  const Double_t kl_t   = TMath::Mean(tvec.begin(), tvec.end(), wvec.begin() );
  const Double_t clus_t = TMath::Mean(tvec_csi.begin(), tvec_csi.end(), wvec.begin() );

  kl.SetTime( kl_t );
  kl.SetCsiTime( clus_t );
}



