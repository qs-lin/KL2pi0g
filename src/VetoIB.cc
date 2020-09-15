#include "g5anaKL2pi0g/VetoIB.h"

#include <map>

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoIB::VetoIB()
   : Veto500MHz("IB",20180601)
{
    Init();
}

VetoIB::VetoIB( const Int_t userflag )
   : Veto500MHz("IB",userflag)
{
    Init();
}

VetoIB::~VetoIB()
{
   ;
}

void VetoIB::Init()
{
   const Double_t t0 = -40.5;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
   //SetLooseVetoEneThreshold(5.);
   SetVetoEneThreshold(1.);
   //SetTightVetoEneThreshold(1.);

   m_isDeadChannel = false;
   if( m_userflag<20181231 ) m_isDeadChannel = true;

   SetCandidateWindow( t0-25., t0+70. );
   SetVetoWindow( t0-15., t0+35.);

   /// constant ///
   m_zcenter = MTBP::IBZPosition + MTBP::IBLength/2.;
}

void VetoIB::Reset()
{
   CommonReset();
   for( Int_t i=0; i<s_arrSize; ++i )
      m_candHitZ[i] = 0.;

   /// single ///
   m_nsingle = 0;
   for( Int_t i=0; i<s_hitSize; ++i )
   {
      m_singleId[i] = -1;
      m_singleEne[i] = 0.;
      m_singleTime[i] = -9999.;
   }
   m_singleProperTime = -9999.;
   m_singleProperArrId = -1;
}

bool VetoIB::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoIB::AddBranches( TTree *tr )
{
   if( !AddCommonBranches(tr) ) return false;
   tr->Branch("IBCandHitZ",m_candHitZ,"IBCandHitZ[IBCandNumber]/D");
   //tr->Branch("IBCandBackTime",m_candBackTime,"IBCandBackTime[IBCandNumber]/D");

   tr->Branch("IBSingleNumber",&m_nsingle,"IBSingleNumber/I");
   tr->Branch("IBSingleId",m_singleId,"IBSingleId[IBSingleNumber]/I");
   tr->Branch("IBSingleEne",m_singleEne,"IBSingleEne[IBSingleNumber]/D");
   tr->Branch("IBSingleTime",m_singleTime,"IBSingleTime[IBSingleNumber]/D");
   tr->Branch("IBSingleProperTime",&m_singleProperTime,"IBSingleProperTime/D");
   tr->Branch("IBSingleProperArrId",&m_singleProperArrId,"IBSingleProperArrId/I");
   return true;
}

void VetoIB::UpdateVars( const KL2pi0g &kl )
{
   Reset();

   std::vector<VetoIB::ModHit> hitvec = GetModHitVec();
   for(std::vector<VetoIB::ModHit>::const_iterator it = hitvec.begin(); it != hitvec.end(); ++it )
   {
      if( !it->single ){
         if( m_isDeadChannel && it->id==23 && it->e < GetSuppEneThreshold()/2. ) continue;
         else if( it->e < GetSuppEneThreshold() ) continue;

         const Double_t vtime  = GetVetoTime(it->z, it->t, kl.csi_t() );
         if( !IsInsideCandidateWindow(vtime) ) continue;

         /// output vars ///
         m_candId[m_ncand]   = it->id;
         m_candEne[m_ncand]  = it->e;
         m_candTime[m_ncand] = vtime;
         m_candHitZ[m_ncand] = it->z + m_zcenter;
         //m_candBackTime[m_ncand] = tback;

         m_ncand++;
      }else{
         if( it->e < GetSuppEneThreshold()/2. ) continue;

         const Double_t vtime = GetVetoTime(it->z, it->t, kl.csi_t() );
         if( !IsInsideCandidateWindow(vtime) ) continue;

         /// output vars ///
         m_singleId[m_nsingle] = it->id;
         m_singleEne[m_nsingle] = it->e;
         m_singleTime[m_nsingle] = vtime;
         
         m_nsingle++;
      }
   } // END mod loop

   EvalProperTime();
   EvalSingleProperTime();
   EvalMaxEneHit();

}

/*
void VetoIB::UpdateVars( const KL2pi0g &kl )
{
   Reset();

   std::vector<VetoIB::ModHit> hitvec = GetModHitVec(); 
   for(std::vector<VetoIB::ModHit>::const_iterator it = hitvec.begin(); it != hitvec.end(); ++it )
   {
      if( m_isDeadChannel && it->id==23 && it->e < 0.5 ) continue; 
      else if( it->e < 1. ) continue;

      const Double_t radius = MTBP::IBInnerInnerR;
      Double_t       hitz   = it->z;
      if( hitz < -MTBP::IBLength/2. ) hitz = -MTBP::IBLength/2.;
      if( hitz >  MTBP::IBLength/2. ) hitz =  MTBP::IBLength/2.; 
      hitz += zcenter;

      const Double_t delz = hitz - kl.v().z();

      const Double_t dist   = TMath::Hypot(radius,delz);
      const Double_t tof    = dist / (TMath::C()/1.E6);
      const Double_t vtime  = it->t - tof - kl.t();

      if( !IsInsideCandidateWindow(vtime) ) continue;

      /// output vars ///
      m_candId[m_ncand]   = it->id;
      m_candEne[m_ncand]  = it->e;
      m_candTime[m_ncand] = vtime;
      m_candHitZ[m_ncand] = hitz;

      m_ncand++;   
   } // END mod loop
}
*/
std::vector<VetoIB::ModHit> 
VetoIB::GetModHitVec() const
{
   const Double_t MaxDeltaT = 25.;
   std::vector<ModHit> hitvec;

   std::map<Int_t, Int_t> up_idmap;
   std::map<Int_t, std::vector<bool> > up_tagmap; 
   for( Int_t ich=0; ich<m_nch; ++ich )
   {
      if( m_nhit[ich]<=0 ) continue;
      const Int_t id = m_chId[ich];
      const bool isUp = ( id < 100 ) ? true : false;

      if( isUp ){
         if( m_isDeadChannel && id==23 ) continue; // run70s: treat ch23 as dead one
         up_idmap.insert( std::make_pair(id, ich) );
         std::vector<bool> tagvec(m_nhit[ich],false);
         up_tagmap.insert( std::make_pair(id, tagvec) );
      }else{
         const Int_t partner_id = id - 100;

         /// run70s: mod 23 special treatment ///
         if( m_isDeadChannel && partner_id==23 )
         {
            for( Int_t ihit=0; ihit<m_nhit[ich]; ++ihit )
            {
               if( m_ene[ich][ihit]<0.5 ) continue;
               VetoIB::ModHit hit;
               hit.id = partner_id;
               hit.t  = m_time[ich][ihit];
               hit.e  = m_ene[ich][ihit];
               hit.z  = 500.; // roughtly true z peak value
               hit.single = false;
               hitvec.push_back(hit); 
            }
         }

         std::map<Int_t, Int_t>::const_iterator it = up_idmap.find(partner_id);

         for( Int_t ihit_d=0; ihit_d<m_nhit[ich]; ++ihit_d )
         {
            bool isFindPartner = false;
            /// if there are hits in upstrem ///
            if( it!=up_idmap.end() ){
               const Int_t up_index = it->second;
               for( Int_t ihit_u=0; ihit_u<m_nhit[up_index]; ++ihit_u )
               {
                  if( TMath::Abs( m_time[ich][ihit_d] - m_time[up_index][ihit_u] ) < MaxDeltaT )
                  {
                     isFindPartner = true;
                     VetoIB::ModHit hit;
                     hit.id = partner_id;
                     hit.t  = GetHitTime( m_time[up_index][ihit_u], m_time[ich][ihit_d] );
                     hit.z  = GetHitZ   ( m_time[up_index][ihit_u], m_time[ich][ihit_d] );
                     hit.e  = GetHitEne ( m_ene [up_index][ihit_u], m_ene [ich][ihit_d],
                                          m_time[up_index][ihit_u], m_time[ich][ihit_d] );
                     hit.single = false;
                     hitvec.push_back(hit);
                     std::map<Int_t, std::vector<bool> >::iterator tagit 
                                                                = up_tagmap.find(partner_id);
                     if( tagit!=up_tagmap.end() )
                        (tagit->second)[ihit_u] = true;
                  }
               } // u loop

            } // END up_idmap not found IF
          
            if( !isFindPartner ){
               VetoIB::ModHit hit;
               hit.id = id;
               hit.t  = m_time[ich][ihit_d];
               hit.z  = -MTBP::IBLength/2.;
               hit.e  = m_ene[ich][ihit_d];
               hit.single = true;
               hitvec.push_back(hit);
            }  
         } // END d hit loop

      } // END else (channel loop)

   } // channel loop

   /// remaining up channel
   for(std::map<Int_t, Int_t>::const_iterator it = up_idmap.begin(); it!=up_idmap.end(); ++it )
   {
      const Int_t arrId = it->first;
      const Int_t modId = it->second;
      std::map<Int_t, std::vector<bool> >::const_iterator tagit = up_tagmap.find(modId);
      if( tagit!=up_tagmap.end() ){
         for(UInt_t ihit=0; ihit<(tagit->second).size(); ++ihit )
         {
            if( !(tagit->second)[ihit] ){
               VetoIB::ModHit hit;
               hit.id = modId;
               hit.t  = m_time[arrId][ihit];
               hit.z  = MTBP::IBLength/2.;
               hit.e  = m_ene[arrId][ihit];
               hit.single = true;
               hitvec.push_back(hit);           
            }
         }
      }
   } // END remaining up loop

   return hitvec;
}

Double_t VetoIB::GetHitTime( const Double_t utime, const Double_t dtime ) const
{
   return (utime + dtime) / 2.;
}

Double_t VetoIB::GetHitZ( const Double_t utime, const Double_t dtime ) const
{
   return MTBP::IBPropVelo * (utime - dtime) / 2.;
}

Double_t VetoIB::GetHitEne ( const Double_t uene , const Double_t dene  ,
                             const Double_t utime, const Double_t dtime ) const
{
   const Double_t hitz = GetHitZ(utime, dtime);
   if( hitz < - MTBP::IBLength / 2. ) return dene; // dowmstrem overflow
   if( hitz >   MTBP::IBLength / 2. ) return uene;
   const Double_t e =  uene / TMath::Exp( -hitz / ( MTBP::IBLAMBDA_U + MTBP::IBALPHA_U * hitz) )
                     + dene / TMath::Exp(  hitz / ( MTBP::IBLAMBDA_D - MTBP::IBALPHA_D * hitz) );
   return e;
}

Double_t VetoIB::GetVetoTime( const Double_t hitz, 
                              const Double_t hit_time,
                              const Double_t csi_time ) const
{
   const Double_t csiz    = MTBP::CSIZPosition + 20.;
   const Double_t radius = MTBP::IBInnerInnerR;
   
   Double_t real_hitz = hitz;

   if( hitz < -MTBP::IBLength/2. ) real_hitz = -MTBP::IBLength/2.;
   if( hitz >  MTBP::IBLength/2. ) real_hitz =  MTBP::IBLength/2.;
   real_hitz += m_zcenter;

   const Double_t delz   = csiz - real_hitz;
   const Double_t dist   = TMath::Hypot(radius,delz);
   const Double_t tback  = dist / (TMath::C()/1.E6);
   const Double_t vtime  = hit_time - csi_time + tback;

   return vtime;
}

void VetoIB::EvalSingleProperTime()
{
   Double_t min_tdiff = 1000.;

   m_singleProperTime = -9999.;
   m_singleProperArrId = -1;

   for( Int_t imod=0; imod<m_nsingle; ++imod )
   {
      if( m_singleEne[imod] < GetVetoEneThreshold()/2. ) continue;
      Double_t tdiff = TMath::Abs(m_singleTime[imod] - GetT0() );
      if( tdiff < min_tdiff ){
         m_singleProperArrId = imod;
         m_singleProperTime = m_singleTime[imod] - GetT0();
         min_tdiff = tdiff;
      }
   }
}

bool VetoIB::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candEne[imod]>GetVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) ){
         return true;
      }
   }

   return false;
}

bool VetoIB::IsLooseVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if(   m_candEne[imod]>GetLooseVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
      {
         return true;
      }
   }
   return false;
}

bool VetoIB::IsTightVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if(   m_candEne[imod]>GetTightVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
      {
         return true;
      }
   }
   return false;
}
