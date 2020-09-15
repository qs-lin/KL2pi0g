#include <iostream>
#include <list>

#include "TObject.h"
#include "TMath.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

#include "g5anaKL2pi0g/KL2pi0g.h"
#include "g5anaKL2pi0g/E14GNAnaDataContainer.h"
#include "g5anaKL2pi0g/Veto.h"

Int_t KinematicCut( E14GNAnaDataContainer &data, std::vector<KL2pi0g> const &klvec);
Int_t VetoCut(std::vector<Veto*> const &vetovec);

void user_cut( E14GNAnaDataContainer &data, std::vector<KL2pi0g> const &klvec,
                                    std::vector<Veto*> const &vetovec)
{
   
  data.CutCondition  = KinematicCut(data, klvec);
  data.VetoCondition = VetoCut(vetovec);
}


Int_t KinematicCut( E14GNAnaDataContainer &data, std::vector<KL2pi0g> const &klvec)
{
  Int_t cutval = 0;

  /// cut threshold ///
  const Double_t k_gEminThreshold         = 50; 
  const Double_t k_minXYThreshold         = 150.;
  const Double_t k_maxRThreshold          = 850.;
  const Double_t k_gDistMinThreshold      = 150.; 
  const Double_t k_EtThreshold            = 600.;
  //const Double_t k_deltaTimeThreshold     = 3.  ; 

  enum{ EGAMMA_MIN=0,    FIDUCIAL_XY,    FIDUCIAL_R,   GAMMA_D_MIN,
        CSI_ET_MIN,      DELTA_VTX_T };

  ///// cut anlaysis /////
  std::list<Gamma> glist;
  glist.push_back( klvec[0].pi0()[0].g1() );
  glist.push_back( klvec[0].pi0()[0].g2() );
  glist.push_back( klvec[0].pi0()[1].g1() );
  glist.push_back( klvec[0].pi0()[1].g2() );
  glist.push_back( klvec[0].g5() );

  Double_t gEmin(HUGE), minXY(HUGE), gDmin(HUGE), esum(0.);
  Double_t maxR = -1.;
  std::vector<Double_t> tvec;

  for(std::list<Gamma>::const_iterator ig=glist.begin(); ig!=glist.end(); ++ig )
  {
    if( ig->e()<gEmin ) gEmin = ig->e();
    esum += ig->e();

    Double_t csi_xy =   (TMath::Abs(ig->x()) > TMath::Abs(ig->y()) ) 
                        ? TMath::Abs(ig->x()) : TMath::Abs(ig->y());

    if( csi_xy < minXY ) minXY = csi_xy;
    if( ig->pos().perp() > maxR ) maxR = ig->pos().perp();

    Double_t dist_xy = TMath::Hypot(ig->x() - klvec[0].v().x(), ig->y() - klvec[0].v().y() );
    Double_t dist    = TMath::Hypot(ig->z() - klvec[0].v().z(), dist_xy );

    Double_t tof     = dist / (TMath::C()/1.e6);
    tvec.push_back( ig->t() - tof );
  }

  std::sort( tvec.begin(), tvec.end() );
  Double_t delT = *(tvec.rbegin() ) - *(tvec.begin() );

  for(std::list<Gamma>::const_iterator g1=glist.begin(); g1!=glist.end(); ++g1 )
  {
    for(std::list<Gamma>::const_iterator g2=glist.begin(); g2!=glist.end(); ++g2 )
    {
      if( g1->id()>=g2->id() ) continue;
      Double_t gD = TMath::Hypot( g1->x()-g2->x(), g1->y()-g2->y() );
      if( gD < gDmin ) gDmin = gD;
    }
  }

/*
  Double_t minK2pi0Chisq(HUGE);
  for(std::vector<KL2pi0g>::const_iterator kl=klvec.begin(); kl!=klvec.end(); ++kl )
    minK2pi0Chisq = ( kl->k2pi0Chisq() < minK2pi0Chisq ) ? kl->k2pi0Chisq() : minK2pi0Chisq;
*/

  //Double_t CoeR   = TMath::Hypot(klvec[0].coe().x(), klvec[0].coe().y());
  //Double_t delM12 = TMath::Abs(klvec[0].pi0().m() - MTBP::Pi0_MASS);

  /// save vars for analysis ///
  data.MinGammaE          = gEmin;
  data.MinFiducialXY      = minXY;
  data.MaxFiducialR       = maxR;
  data.MinClusterDistance = gDmin;
  data.MaxDeltaVertexTime = delT;

  /// prepara cut condition ///
  if( gEmin         < k_gEminThreshold         ) cutval |= ( 1 << EGAMMA_MIN      );
  if( minXY         < k_minXYThreshold         ) cutval |= ( 1 << FIDUCIAL_XY     );
  if( maxR          > k_maxRThreshold          ) cutval |= ( 1 << FIDUCIAL_R      );
  if( gDmin         < k_gDistMinThreshold      ) cutval |= ( 1 << GAMMA_D_MIN     ); 
  if( esum          < k_EtThreshold            ) cutval |= ( 1 << CSI_ET_MIN      );
  //if( delM12        > k_m12MaxDiffThreshold    ) cutval |= ( 1 << M12_PI_MASS     );
  //if( CoeR          > k_CoeRMaxThreshold       ) cutval |= ( 1 << COE_CUT         );
  //if( delT          > k_deltaTimeThreshold     ) cutval |= ( 1 << DELTA_VTX_T     );   
  //if( minK2pi0Chisq < k_minK2pi0ChisqThreshold ) cutval |= ( 1 << MIN_KPIPI_CHISQ ); 

  return cutval;  
}

///
Int_t VetoCut(std::vector<Veto*> const &vetovec)
{
  Int_t vetoval = 0;

  for( std::vector<Veto*>::const_iterator iv=vetovec.begin(); 
       iv!=vetovec.end(); ++iv )
  {
    if( (*iv)->IsVeto() ) vetoval += ( 1 << (*iv)->GetId() );
  }

  return vetoval;
}

