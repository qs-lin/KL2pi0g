#include "g5anaKL2pi0g/VetoOEV.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoOEV::VetoOEV()
   : Veto125MHz("OEV",20180601)
{
    Init();
}

VetoOEV::VetoOEV( const Int_t userflag )
   : Veto125MHz("OEV",userflag)
{
    Init();
}

VetoOEV::~VetoOEV()
{
   ;
}

void VetoOEV::Init()
{
   const Double_t t0 = 19.5;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
//   SetLooseVetoEneThreshold(5.);
   SetVetoEneThreshold(1.);
//   SetTightVetoEneThreshold(1.);

   const Double_t cand_w = 30.;
   const Double_t veto_w = 10.;
   SetCandidateWindow( t0 - cand_w, t0 + cand_w);
   SetVetoWindow( t0 - veto_w, t0 + veto_w);
}

void VetoOEV::Reset()
{
   CommonReset();
}

bool VetoOEV::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoOEV::AddBranches( TTree *tr )
{
   return AddCommonBranches(tr);
}


void VetoOEV::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t z_surface = MTBP::OEVZPosition + 20.;   

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<1.) continue;

      const Double_t delz   = TMath::Abs(kl.v().z() - z_surface );
      const Double_t radius = 950.; // (rough) OEV radius

      const Double_t dist   = TMath::Hypot(delz, radius);
      const Double_t tof    = dist / (TMath::C()/1.E6);
      const Double_t vtime  = m_time[imod] - tof - kl.t();

      if( !IsInsideCandidateWindow(vtime) ) continue;

      /// output vars ///
      m_candId[m_ncand]   = m_modId[imod];
      m_candEne[m_ncand]  = m_ene[imod];
      m_candTime[m_ncand] = vtime;

      m_ncand++;
   }

   EvalProperTime();
   EvalMaxEneHit();
}

bool VetoOEV::IsVeto()
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
   const Double_t k_ethreshold = 2.;

   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
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

bool VetoOEV::IsLooseVeto()
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

bool VetoOEV::IsTightVeto()
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
