#include "g5anaKL2pi0g/VetoCC03.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoCC03::VetoCC03()
   : Veto125MHz("CC03",20180601)
{
    Init();
}

VetoCC03::VetoCC03( const Int_t userflag )
   : Veto125MHz("CC03",userflag)
{
    Init();
}

VetoCC03::~VetoCC03()
{
   ;
}

void VetoCC03::Init()
{
   const Double_t t0 = 32.2;
   SetT0(t0);

   SetSuppEneThreshold(1.);
//   SetLooseVetoEneThreshold(5.);
   SetVetoEneThreshold(3.);
//   SetTightVetoEneThreshold(3.);

   const Double_t cand_w = 30.;
   const Double_t veto_w = 15.;
   SetCandidateWindow( t0 - cand_w, t0 + cand_w);
   SetVetoWindow( t0 - veto_w, t0 + veto_w);
}

void VetoCC03::Reset()
{
   CommonReset();
}

bool VetoCC03::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoCC03::AddBranches( TTree *tr )
{
   return AddCommonBranches(tr);
}


void VetoCC03::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t z_surface = MTBP::CC03ZPosition + 20.;   

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<GetSuppEneThreshold() ) continue;

      const Double_t dist   = TMath::Abs(kl.v().z() - z_surface ); 
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

bool VetoCC03::IsVeto()
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

bool VetoCC03::IsLooseVeto()
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

bool VetoCC03::IsTightVeto()
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
