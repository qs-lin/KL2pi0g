#include "g5anaKL2pi0g/VetoBHGC.h"

#include <map>

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoBHGC::VetoBHGC()
   : Veto500MHz("BHGC",20180601)
{
    Init();
}

VetoBHGC::VetoBHGC( const Int_t userflag )
   : Veto500MHz("BHGC",userflag)
{
    Init();
}

VetoBHGC::~VetoBHGC()
{
   ;
}

void VetoBHGC::Init()
{
   const Double_t t0 = -77.4;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
//   SetLooseVetoEneThreshold(2.5);
   SetVetoEneThreshold(2.5);
//   SetTightVetoEneThreshold(2.5);

   const Double_t cand_w = 30.;
   const Double_t veto_w = 7.5;
   SetCandidateWindow( t0 - cand_w, t0 + cand_w );
   SetVetoWindow( t0 - veto_w, t0 + veto_w );

   SetCoinWindowWidth(5.);
}

void VetoBHGC::Reset()
{
   CommonReset();
}

bool VetoBHGC::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoBHGC::AddBranches( TTree *tr )
{
   if( !AddCommonBranches(tr) ) return false;
   return true;
}

std::map<Int_t, std::pair<Int_t, Int_t> > 
VetoBHGC::GetModIdMap() const
{
   std::map<Int_t, std::pair<Int_t, Int_t> > modIdMap;
   std::map<Int_t, Int_t> mapA, mapB;

   for( Int_t ich=0; ich<m_nch; ++ich )
   {
      const Int_t chId = m_chId[ich];
      bool isChannelA = (chId%2==0) ? true : false;
      const Int_t modId = chId / 2;

      if( isChannelA ) mapA.insert( std::make_pair(modId, ich) );
      else             mapB.insert( std::make_pair(modId, ich) );   
   }

   for(std::map<Int_t, Int_t>::const_iterator a = mapA.begin(); a != mapA.end(); ++a )
   {
      const Int_t modId = a->first;
      std::map<Int_t, Int_t>::const_iterator b = mapB.find(modId);
      if( b != mapB.end() ){
         std::pair<Int_t, Int_t> modpair(a->second, b->second);
         modIdMap.insert( std::make_pair(modId, modpair) );
      }else{
         continue;
      }
   }

   return modIdMap;
}

void VetoBHGC::UpdateVars( const KL2pi0g &kl )
{
   Reset();

   const Double_t detz = MTBP::BHGCZPosition + 20.;
   std::map<Int_t, std::pair<Int_t, Int_t> > modmap = GetModIdMap();
   
   for(std::map<Int_t, std::pair<Int_t, Int_t> >::const_iterator im = modmap.begin(); 
                                                                 im != modmap.end(); ++im )
   {
      const Int_t modId =  im->first;
      const Int_t id_A  = (im->second).first;
      const Int_t id_B  = (im->second).second;

      Double_t dist = detz - kl.v().z();
      dist += (modId<2) ? 0. : MTBP::BHGCModuleInterval;

      const Double_t tof = dist / (TMath::C()*1.e-6);

      for( Int_t ihit_A=0; ihit_A<m_nhit[id_A]; ++ihit_A )
      {
         for( Int_t ihit_B=0; ihit_B<m_nhit[id_B]; ++ihit_B )
         {
            if(    m_ene[id_A][ihit_A] < GetSuppEneThreshold() 
                || m_ene[id_B][ihit_B] < GetSuppEneThreshold() ) continue;

            const Double_t esum  = (m_ene[id_A][ihit_A] + m_ene[id_B][ihit_B]);
            if( esum < GetVetoEneThreshold() ) continue;

            if( !IsCoinHitPair(m_time[id_A][ihit_A], m_time[id_B][ihit_B]) ) continue;

            const Double_t htime = (m_time[id_A][ihit_A] + m_time[id_B][ihit_B]) / 2.;
            const Double_t vtime = htime - tof - kl.t();

            if( !IsInsideCandidateWindow(vtime) ) continue;

            /// save output vars ///
            m_candId[m_ncand]        = modId;
            m_candEne[m_ncand]       = esum;
            m_candTime[m_ncand]      = vtime;

            m_ncand++;
         }
      }
   }

   EvalProperTime();
   EvalMaxEneHit();

}

bool VetoBHGC::IsCoinHitPair( const Double_t t1, const Double_t t2 ) const
{
   return (TMath::Abs(t1-t2) < GetCoinWindowWidth() );
}

bool VetoBHGC::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if(   m_candEne[imod]>GetVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
      {
         return true;
      }
   }
   return false;
}

bool VetoBHGC::IsLooseVeto()
{
   return false;
}

bool VetoBHGC::IsTightVeto()
{
   return IsVeto();
}
