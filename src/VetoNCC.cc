#include "g5anaKL2pi0g/VetoNCC.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoNCC::VetoNCC()
   : Veto125MHz("NCC",20180601)
{
    Init();
}

VetoNCC::VetoNCC( const Int_t userflag )
   : Veto125MHz("NCC",userflag)
{
    Init();
}

VetoNCC::~VetoNCC()
{
   ;
}

void VetoNCC::Init()
{
   const Double_t t0 = 10.9;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
   //SetLooseVetoEneThreshold(5.0);
   SetVetoEneThreshold(1.0);
   //SetTightVetoEneThreshold(1.0);

   SetCandidateWindow( t0-35., t0+70. );
   SetVetoWindow( t0-15., t0+45. );
}

void VetoNCC::Reset()
{
   CommonReset();
}

bool VetoNCC::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoNCC::AddBranches( TTree *tr )
{
   return AddCommonBranches(tr);
}

void VetoNCC::UpdateVars( const KL2pi0g &kl )
{
   Reset();

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( (m_modId[imod]%10)!=0 || m_modId[imod]>=600
                                || m_ene[imod]<GetSuppEneThreshold() ) continue;

      //const Double_t dist   = TMath::Abs(kl.v().z() - MTBP::NCCZPosition ); 
      //const Double_t tof    = dist / (TMath::C()/1.E6);
      //const Double_t vtime  = m_time[imod] - tof - kl.t();
      const Double_t vtime  = m_time[imod] - kl.csi_t();

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

bool VetoNCC::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candEne[imod]>GetVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) ){
         return true;
      } 
   }
   return false;
}

bool VetoNCC::IsLooseVeto()
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

bool VetoNCC::IsTightVeto()
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
