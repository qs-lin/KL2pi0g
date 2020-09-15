#include "g5anaKL2pi0g/VetoBHPV.h"

#include <map>

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoBHPV::VetoBHPV()
   : Veto500MHz("BHPV",20180601)
{
    Init();
}

VetoBHPV::VetoBHPV( const Int_t userflag )
   : Veto500MHz("BHPV",userflag)
{
    Init();
}

VetoBHPV::~VetoBHPV()
{
   ;
}

void VetoBHPV::Init()
{
   const Double_t t0 = -76.2;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
   SetLooseVetoEneThreshold(2.5);
   SetVetoEneThreshold(2.5);
   SetTightVetoEneThreshold(2.5);

   const Double_t cand_w = 30.;
   const Double_t veto_w = 7.5;
   SetCandidateWindow( t0 - cand_w, t0 + cand_w );
   SetVetoWindow( t0 - veto_w, t0 + veto_w );
   SetNcoinThreshold(3);
   SetNcoinTightThreshold(2);

   SetTofVec();
}

void VetoBHPV::Reset()
{
   CommonReset();
   for( Int_t i=0; i<s_arrSize; ++i )
   {
      m_candNcoin[i] = 0; 
      m_candTimeSpread[i] = 0.;
   }

   m_proper2CoinTime = m_proper3CoinTime = -9999.;
}

bool VetoBHPV::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoBHPV::AddBranches( TTree *tr )
{
   if( !AddCommonBranches(tr) ) return false;
   tr->Branch("BHPVCandNcoin",m_candNcoin,"BHPVCandNcoin[BHPVCandNumber]/I");
   tr->Branch("BHPVProper2CoinTime",&m_proper2CoinTime,"BHPVProper2CoinTime/D");
   tr->Branch("BHPVProper3CoinTime",&m_proper3CoinTime,"BHPVProper3CoinTime/D");
   return true;
}

std::map<Int_t, std::list<VetoBHPV::ModHit> >
VetoBHPV::GetModHitMap() const
{
   VetoBHPV::ChMap_t nmap, smap;

   for( Int_t ich=0; ich<m_nch; ++ich )
   {
      const Int_t modId = m_chId[ich] / 2 ;
      if( modId >= MTBP::nModBHPV ) continue; 

      const bool isNorth = ( (ich%2)==0 ) ? true : false; 
   
      std::list<Int_t> hitlist;
      for( Int_t ihit=0; ihit<m_nhit[ich]; ++ihit )
      { 
         /// to go back to the original algorithm, use following statement instead ///
         // if( m_ene[ich][ihit] < GetVetoThreshold() ) continue;
         if( m_ene[ich][ihit] < GetSuppEneThreshold() ) continue;
         hitlist.push_back(ihit);
      }   

      if( hitlist.empty() ) continue;
      HitId_t hitId(ich, hitlist);

      if( isNorth ){
         nmap.insert( std::make_pair(modId, hitId) );
      }else{
         smap.insert( std::make_pair(modId, hitId) );
      }
   }

   return GetModHitMap(nmap, smap);
}

std::map<Int_t, std::list<VetoBHPV::ModHit> > 
VetoBHPV::GetModHitMap( VetoBHPV::ChMap_t &nmap, VetoBHPV::ChMap_t &smap) const
{
   std::map<Int_t, std::list<VetoBHPV::ModHit> > modmap;

   /// north module loop ///
   for(VetoBHPV::ChMap_t::iterator ch_n = nmap.begin(); ch_n != nmap.end(); ++ch_n )
   {
      std::list<VetoBHPV::ModHit> hlist_mod;

      const Int_t modId  = ch_n->first;
      HitId_t &hitId_n = (ch_n->second);
      const Int_t idx_n = hitId_n.first;
      std::list<Int_t> &hlist_n = hitId_n.second;
      
      ChMap_t::iterator  ch_s = smap.find(modId);
      bool isSouthExist = (ch_s != smap.end() );

      for(std::list<Int_t>::iterator hit_n  = hlist_n.begin();
                                     hit_n != hlist_n.end(); ++hit_n )
      {
         ModHit hit;
         hit.id = modId;

         const Double_t time_n = m_time[idx_n][*hit_n];
         bool isCoinHit = false;

         if( isSouthExist ){
            HitId_t &hitId_s = (ch_s->second);
            const Int_t idx_s = hitId_s.first;
            std::list<Int_t> &hlist_s = hitId_s.second;

            for(std::list<Int_t>::iterator hit_s  = hlist_s.begin();
                                           hit_s != hlist_s.end(); ++hit_s )
            {
               const Double_t time_s = m_time[idx_s][*hit_s];
               if( TMath::Abs(time_n - time_s) < MTBP::BHPVCoincidenceTimeWindow ){
                  isCoinHit = true;
                  hit.e = m_ene[idx_n][*hit_n] + m_ene[idx_s][*hit_s];
                  hit.t = (time_n + time_s) / 2.;
                  hlist_s.erase(hit_s);
                  break;
               }
            }

            /// if all hits in south channel are analyzed, treat it as vacancy later on.
            if( hlist_s.empty() ) isSouthExist = false;
         }

         if( !isCoinHit ){
            hit.e = m_ene[idx_n][*hit_n];
            hit.t = time_n;
         }

         if( hit.e < GetVetoEneThreshold() ) continue;
         hlist_mod.push_back(hit);
      } // north hit loop

      if( isSouthExist ){
         HitId_t &hitId_s = (ch_s->second);
         const Int_t idx_s = hitId_s.first;
         std::list<Int_t> &hlist_s = hitId_s.second;

         for(std::list<Int_t>::iterator hit_s  = hlist_s.begin();
                                        hit_s != hlist_s.end(); ++hit_s )
         {
            VetoBHPV::ModHit hit;
            hit.e = m_ene[idx_s][*hit_s];
            if( hit.e < GetVetoEneThreshold() ) continue; 
            hit.t = m_time[idx_s][*hit_s];
            hlist_mod.push_back(hit);
         }
      }
 
      smap.erase(modId);
 
      hlist_mod.sort(); 
      modmap.insert( std::make_pair(modId, hlist_mod) ); 
      
   } /// module loop

   /// remaining south module loop ///
   for(ChMap_t::iterator ch_s = smap.begin(); ch_s != smap.end(); ++ch_s ) 
   {
     std::list<VetoBHPV::ModHit> hlist_mod;

      const Int_t modId  = ch_s->first;
      HitId_t &hitId_s = (ch_s->second);
      const Int_t idx_s = hitId_s.first;
      std::list<Int_t> &hlist_s = hitId_s.second; 

      for(std::list<Int_t>::iterator hit_s  = hlist_s.begin();
                                     hit_s != hlist_s.end(); ++hit_s )
      {
         if( m_ene[idx_s][*hit_s] < GetVetoEneThreshold() ) continue;
         VetoBHPV::ModHit hit;
         hit.id = modId;
         hit.e = m_ene[idx_s][*hit_s];
         hit.t = m_time[idx_s][*hit_s];
         hlist_mod.push_back(hit);
      }

      if( hlist_mod.empty() ) continue;

      modmap.insert( std::make_pair(modId, hlist_mod) );
   } // END south loop

   return modmap; 
}

void VetoBHPV::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t detz = MTBP::BHPVZPosition + 20.;
   const Double_t common_tof = (detz - kl.v().z() ) / (TMath::C()/1.e6 );

   const UInt_t nMinCoinMod = 2;
   typedef std::map<Int_t, std::list<VetoBHPV::ModHit> > ModMap_t;
   ModMap_t modmap = GetModHitMap();
   if( modmap.size() < nMinCoinMod ) return;

   for(ModMap_t::iterator m1 = modmap.begin(); m1!=modmap.end(); ++m1 )
   {
      const Int_t id1 = m1->first;

      /// concidence timing analysis
      std::list<VetoBHPV::ModHit> &hlist1 = m1->second;
      for(std::list<VetoBHPV::ModHit>::iterator h1 = hlist1.begin(); h1!=hlist1.end(); ++h1 )
      {
         std::vector<Double_t> tvec;
         Double_t esum = 0.;         
    
         const Double_t t1 = h1->t - m_tofvec[id1] - common_tof - kl.t(); 
         if( TMath::Abs( t1 - GetT0() ) > MTBP::BHPVCoincidenceStoreTimeRange ) continue;
         tvec.push_back( t1 );
         esum += h1->e;

         for( Int_t id2 = id1 + tvec.size(); id2 < MTBP::nModBHPV; ++id2 )
         {
            ModMap_t::iterator m2 = modmap.find(id2);
            if( m2==modmap.end() ) break;

            const Double_t t_mean = TMath::Mean(tvec.size(), &tvec[0]);
            std::list<VetoBHPV::ModHit> &hlist2 = m2->second;

            bool isCoinHit = false;
            for(std::list<VetoBHPV::ModHit>::iterator h2 = hlist2.begin(); h2!=hlist2.end(); ++h2 )
            {
               const Double_t t2 = h2->t - m_tofvec[id2] - common_tof - kl.t(); 
               if( TMath::Abs( t2 - t_mean ) < MTBP::BHPVCoincidenceTimeWindow )
               {
                  tvec.push_back(t2);
                  esum += h2->e;
                  isCoinHit = true;
                  hlist2.erase(h2);
                  break;
               } 
            }

            if( !isCoinHit ) break;
         }

         if( tvec.size() < nMinCoinMod ) continue;
         
         /// output vars ///
         m_candId[m_ncand]        = id1;
         m_candNcoin[m_ncand]     = tvec.size();
         m_candEne[m_ncand]       = esum;
         m_candTime[m_ncand]       = TMath::Mean(tvec.size(), &tvec[0] );
         m_candTimeSpread[m_ncand] = TMath::RMS (tvec.size(), &tvec[0] );

         const Double_t tdiff = m_candTime[m_ncand] - GetT0();
         if( m_candNcoin[m_ncand]>=2 && TMath::Abs(tdiff)<TMath::Abs(m_proper2CoinTime)  )
            m_proper2CoinTime = tdiff;
         if( m_candNcoin[m_ncand]>=3 && TMath::Abs(tdiff)<TMath::Abs(m_proper3CoinTime)  )
            m_proper3CoinTime = tdiff;

         m_ncand++;
         if( m_ncand>=s_arrSize ){
            std::cout<<" [Error] <VetoBHPV::UpdateVars> Too many candidates : n =  " 
                     << m_ncand << std::endl;
            return;
         }
      } // first module hit loop

   } // first module loop

}

void VetoBHPV::SetTofVec()
{
   Double_t dist = 0;
   m_tofvec.clear();
   for( Int_t imod=0 ; imod<MTBP::nModBHPV ; ++imod )
   {
      m_tofvec.push_back( dist / (TMath::C()*1e-6) );
      dist += (MTBP::BHPVModZLength + MTBP::BHPVModGapZ[imod]);
   }
}

bool VetoBHPV::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candNcoin[imod] >= GetNcoinThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
         return true;
   }

   return false;
}

bool VetoBHPV::IsLooseVeto()
{
   return false;
}

bool VetoBHPV::IsTightVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candNcoin[imod] >= GetNcoinTightThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
         return true;
   }

   return false;
}
