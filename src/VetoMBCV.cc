#include "g5anaKL2pi0g/VetoMBCV.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoMBCV::VetoMBCV()
   : Veto125MHz("MBCV",20180601)
{
    Init();
}

VetoMBCV::VetoMBCV( const Int_t userflag )
   : Veto125MHz("MBCV",userflag)
{
    Init();
}

VetoMBCV::~VetoMBCV()
{
   ;
}

void VetoMBCV::Init()
{
   const Double_t t0 = 19.4;
   SetT0(t0);

   SetSuppEneThreshold(0.5);

   //SetLooseVetoEneThreshold(9999.);
   SetVetoEneThreshold(1.);
   //SetTightVetoEneThreshold(9999.);

   const Double_t cand_w = 40.;
   const Double_t veto_w = 15.;
   SetCandidateWindow( t0 - cand_w, t0 + cand_w);
   SetVetoWindow( t0 - veto_w, t0 + veto_w);
}

void VetoMBCV::Reset()
{
   CommonReset();
}

bool VetoMBCV::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoMBCV::AddBranches( TTree *tr )
{
   return AddCommonBranches(tr);
}


void VetoMBCV::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t z_surface = MTBP::CSIZPosition + 20.;
   const Double_t radius    = MTBP::CBARInnerInnerR;

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<.5) continue;
      
      const Double_t delz   = TMath::Abs(kl.v().z() - z_surface ); 
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

bool VetoMBCV::IsVeto()
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

bool VetoMBCV::IsLooseVeto()
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

bool VetoMBCV::IsTightVeto()
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
