#include "g5anaKL2pi0g/VetoIBCV.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoIBCV::VetoIBCV()
   : Veto125MHz("IBCV",20180601)
{
    Init();
}

VetoIBCV::VetoIBCV( const Int_t userflag )
   : Veto125MHz("IBCV",userflag)
{
    Init();
}

VetoIBCV::~VetoIBCV()
{
   ;
}

void VetoIBCV::Init()
{
   const Double_t t0 = -10.8;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
//   SetLooseVetoEneThreshold(1.);
   SetVetoEneThreshold(1.);
//   SetTightVetoEneThreshold(9999.);

   m_isDeadChannel = false;
   if( m_userflag < 20181231 ) m_isDeadChannel = true;

   const Double_t cand_w = 40.;
   const Double_t veto_w = 10.;
   SetCandidateWindow( t0 - cand_w, t0 + cand_w);
   SetVetoWindow( t0 - veto_w, t0 + veto_w);
}

void VetoIBCV::Reset()
{
   CommonReset();
   for( Int_t i=0; i<s_arrSize; ++i )
      m_candHitZ[i] = 0.;
}

bool VetoIBCV::SetBranchAddresses( TTree *tr )
{
   if( !SetCommonBranches(tr) ) return false;
   tr->SetBranchAddress(Form("%sModuleHitZ",m_detname.c_str()), m_hitz);
   return true;
}

bool VetoIBCV::AddBranches( TTree *tr )
{
   if( !AddCommonBranches(tr) ) return false;
   tr->Branch("IBCVCandHitZ",m_candHitZ,"IBCVCandHitZ[IBCVCandNumber]/D");
   return true;
}


void VetoIBCV::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t zcenter = MTBP::IBCVZPosition + MTBP::IBCVLength/2.;
   const Double_t csiz    = MTBP::CSIZPosition + 20.;

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_isDeadChannel && m_modId[imod]==27 && m_ene[imod]<GetSuppEneThreshold()/2. ) continue;
      else if( m_ene[imod]<GetSuppEneThreshold() ) continue;

      const Double_t radius = MTBP::IBCVInnerR;

      /// overflow treatment
      Double_t hitz       = m_hitz[imod];
      if( hitz < -MTBP::IBCVLength/2. ) hitz = -MTBP::IBCVLength/2.;
      if( hitz >  MTBP::IBCVLength/2. ) hitz =  MTBP::IBCVLength/2.;

      hitz += zcenter;
      const Double_t delz   = csiz - hitz; 

      const Double_t dist   = TMath::Hypot(radius,delz); 
      const Double_t tback  = dist / (TMath::C()/1.E6);
      const Double_t vtime  = m_time[imod] - kl.csi_t() + tback;

      if( !IsInsideCandidateWindow(vtime) ) continue;

      /// output vars ///
      m_candId[m_ncand]   = m_modId[imod];
      m_candEne[m_ncand]  = m_ene[imod];
      m_candTime[m_ncand] = vtime;
      m_candHitZ[m_ncand] = hitz;

      m_ncand++;
   }

   EvalProperTime();
   EvalMaxEneHit();
}

bool VetoIBCV::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candEne[imod]>GetVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) ){
         return true;
      }
   }
   return false;
/*
   /// temporarily settings ///
   const Double_t k_t1 = GetT0() - 15.;
   const Double_t k_t2 = GetT0() + 15.;
   const Double_t k_ethreshold = 1.;

   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_isDeadChannel && m_candId[imod]==27 ){
         if(    m_candEne[imod] > k_ethreshold 
             && m_candTime[imod] > k_t1 - 20.
             && m_candTime[imod] < k_t2 + 20. )
         {
            return true;
         }
      }

      if(   m_candEne[imod] > k_ethreshold 
         && m_candTime[imod] > k_t1 
         && m_candTime[imod] < k_t2          )
      {
         return true;
      }
   }
   return false;
*/
}

bool VetoIBCV::IsLooseVeto()
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

bool VetoIBCV::IsTightVeto()
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
