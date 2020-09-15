#include "g5anaKL2pi0g/VetoCSI.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"
#include "csimap/CsiMap.h"

VetoCSI::VetoCSI()
   : Veto125MHz("CSI",20180601),
     m_esupp(20180601)
{
    Init();
}

VetoCSI::VetoCSI( const Int_t userflag )
   : Veto125MHz("CSI",userflag),
     m_esupp(userflag)
{
    Init();
}

VetoCSI::~VetoCSI()
{
   ;
}

void VetoCSI::Init()
{
   const Double_t t0 = 0.;
   SetT0(t0);

   SetLooseVetoEneThreshold(3.);
   SetVetoEneThreshold(3.);
   SetTightVetoEneThreshold(3.);

   SetCandidateWindow( t0 - 20., t0 + 20.);
   SetVetoWindow(t0-3., t0+3.);
}

Double_t VetoCSI::Get3SigmaSuppEneThreshold( const Int_t modId ) const
{
   return m_esupp.GetEneSuppThreshold(modId);
}

void VetoCSI::Reset()
{
   CommonReset();
   for( Int_t i=0; i<s_arrSize; ++i )
      m_candDist[i] = 0.;
}

bool VetoCSI::SetBranchAddresses( TTree *tr )
{
   if( tr==NULL ) return false;
   tr->SetBranchAddress("CSINumber", &m_nmod  );
   tr->SetBranchAddress("CSIModID" ,  m_modId );
   tr->SetBranchAddress("CSIEne"   ,  m_ene   );
   tr->SetBranchAddress("CSITime"  ,  m_time  );
   return true;
}

bool VetoCSI::AddBranches( TTree *tr )
{
   if( !AddCommonBranches(tr) ) return false;
   tr->Branch("CSICandDist",m_candDist,"CSICandDist[CSICandNumber]/D");
   return true;
}

bool VetoCSI::EvalCsiVeto( const Double_t dist,
                           const Double_t ene,
                           const Double_t ene_thr ) const
{
   const Double_t dist_UL = 400 * (10 - ene_thr) / 7. + 200;

   if( dist < 200. ){
      if( ene > 10. && ene > ene_thr ) return true;
   }else if( dist > 200. && dist < dist_UL ){
      if( ene > 10. - (dist - 200. ) / 400. * 7 ) return true;
   }else{
      if( ene > ene_thr ) return true;
   }

   return false;
}

void VetoCSI::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const std::map<Int_t, Int_t> idmap = GetCsiVetoMap(kl);
   const Int_t ng = 5;
   Double_t cx[ng], cy[ng], ct[ng];
   GetClusterInfo(kl, cx, cy, ct); 

   for(std::map<Int_t, Int_t>::const_iterator it = idmap.begin(); it!= idmap.end(); ++it )
   {
      Double_t dist    = 1000.;
      Double_t delta_t = -1.;

      for( Int_t ig=0; ig<ng; ++ig )
      {
         const Double_t this_dist = CalcDistFromCsi(cx[ig], cy[ig], it->first );
         if( this_dist < dist )
         {
            dist = this_dist;
            delta_t = TMath::Abs( ct[ig] - m_time[it->second] );
         }
      }
   
      /// output vars ///
      m_candId[m_ncand]   = it->first;
      m_candEne[m_ncand]  = m_ene[it->second];
      m_candTime[m_ncand] = delta_t;
      m_candDist[m_ncand] = dist;

      m_ncand++;
   }
}

std::map<Int_t, Int_t> VetoCSI::GetCsiVetoMap( const KL2pi0g &kl ) const
{
   std::map<Int_t, Int_t> idmap;
   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      const Int_t id = m_modId[imod];
      if( m_ene[imod] < Get3SigmaSuppEneThreshold(id) ) continue;
      idmap.insert( std::make_pair(id, imod) );
   }

   ExcludeCsiInCluster(idmap, kl.pi0()[0].g1().cluster().clusterIdVec() );
   ExcludeCsiInCluster(idmap, kl.pi0()[0].g2().cluster().clusterIdVec() );
   ExcludeCsiInCluster(idmap, kl.pi0()[1].g1().cluster().clusterIdVec() );
   ExcludeCsiInCluster(idmap, kl.pi0()[1].g2().cluster().clusterIdVec() );
   ExcludeCsiInCluster(idmap, kl.g5().cluster().clusterIdVec()       );

   return idmap;
}

void VetoCSI::ExcludeCsiInCluster( std::map<Int_t, Int_t> &idmap,
                                   std::vector<Int_t> const& cls_vec ) const
{
   for(std::vector<Int_t>::const_iterator it = cls_vec.begin(); it!= cls_vec.end(); ++it )
      idmap.erase(*it);
}

Double_t VetoCSI::CalcDistFromCsi( const Double_t x, const Double_t y, const Int_t id ) const
{
   const Double_t csi_x = CsiMap::getCsiMap()->getX(id);
   const Double_t csi_y = CsiMap::getCsiMap()->getY(id);
   return TMath::Hypot(csi_x - x, csi_y - y );
}

void VetoCSI::GetClusterInfo( const KL2pi0g &kl, Double_t *x, Double_t *y, Double_t *t ) const
{
   x[0] = kl.pi0()[0].g1().x();
   y[0] = kl.pi0()[0].g1().y();
   t[0] = kl.pi0()[0].g1().t();

   x[1] = kl.pi0()[0].g2().x();
   y[1] = kl.pi0()[0].g2().y();
   t[1] = kl.pi0()[0].g2().t();

   x[2] = kl.pi0()[1].g1().x();
   y[2] = kl.pi0()[1].g1().y();
   t[2] = kl.pi0()[1].g1().t();

   x[3] = kl.pi0()[1].g2().x();
   y[3] = kl.pi0()[1].g2().y();
   t[3] = kl.pi0()[1].g2().t();


   x[4] = kl.g5().x();
   y[4] = kl.g5().y();
   t[4] = kl.g5().t();

}

bool VetoCSI::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      const Int_t id = m_candId[imod];
      const Double_t ene_thr = ( Get3SigmaSuppEneThreshold(id) > 3. ) ? 
                                Get3SigmaSuppEneThreshold(id) : 3.;
      if( !IsInsideVetoWindow(m_candTime[imod]) ) continue;
      if( EvalCsiVeto(m_candDist[imod], m_candEne[imod], ene_thr) ) 
         return true;
   }
   return false;
}

bool VetoCSI::IsLooseVeto()
{
   return false;
}

bool VetoCSI::IsTightVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      const Int_t id = m_candId[imod];
      if( !IsInsideVetoWindow(m_candTime[imod]) ) continue;
      if( EvalCsiVeto(m_candDist[imod], m_candEne[imod], Get3SigmaSuppEneThreshold(id)) )
         return true;
   }
   return false;
}

