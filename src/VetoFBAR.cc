#include "g5anaKL2pi0g/VetoFBAR.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoFBAR::VetoFBAR()
   : Veto125MHz("FBAR",20180601)
{
    Init();
}

VetoFBAR::VetoFBAR( const Int_t userflag )
   : Veto125MHz("FBAR",userflag)
{
    Init();
}

VetoFBAR::~VetoFBAR()
{
   ;
}

void VetoFBAR::Init()
{
   const Double_t t0 = 40.2;
   SetT0(t0);

   SetSuppEneThreshold(.5);
   //SetLooseVetoEneThreshold(5.);
   SetVetoEneThreshold(1.);
   //SetTightVetoEneThreshold(1.);

   m_isDeadChannel = false;
   if( m_userflag>=20180101 && m_userflag<20180501 ) m_isDeadChannel = true;

   SetCandidateWindow( t0-30., t0+55. );
   SetVetoWindow( t0-15., t0+40. );
}

void VetoFBAR::Reset()
{
   CommonReset();
}

bool VetoFBAR::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoFBAR::AddBranches( TTree *tr )
{
   return AddCommonBranches(tr);
}

void VetoFBAR::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t z_end = MTBP::FBARZPosition + MTBP::FBARLength;

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<GetSuppEneThreshold() ) continue;
      if( m_isDeadChannel && m_modId[imod]==5 ) continue; // run78, mod5 is dead //

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

bool VetoFBAR::IsVeto()
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

bool VetoFBAR::IsLooseVeto()
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

bool VetoFBAR::IsTightVeto()
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
