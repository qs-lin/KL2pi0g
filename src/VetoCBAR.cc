#include "g5anaKL2pi0g/VetoCBAR.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoCBAR::VetoCBAR()
   : Veto125MHz("CBAR",20180601)
{
    Init();
}

VetoCBAR::VetoCBAR( const Int_t userflag )
   : Veto125MHz("CBAR",userflag)
{
    Init();
}

VetoCBAR::~VetoCBAR()
{
   ;
}

void VetoCBAR::Init()
{
   const Double_t t0 = 42.6;
   SetT0(t0);

   SetSuppEneThreshold(0.5);
   //SetLooseVetoEneThreshold(5.);
   SetVetoEneThreshold(1.);
   //SetTightVetoEneThreshold(1.);

   SetCandidateWindow( t0-30., t0+70. );
   SetVetoWindow( t0-15., t0+45. );
}

void VetoCBAR::Reset()
{
   CommonReset();
   for( Int_t i=0; i<s_arrSize; ++i )
      m_candHitZ[i] = 0.;
}

bool VetoCBAR::SetBranchAddresses( TTree *tr )
{
   if( !SetCommonBranches(tr) ) return false;
   tr->SetBranchAddress(Form("%sModuleHitZ",m_detname.c_str()), m_hitz);
   return true;
}

bool VetoCBAR::AddBranches( TTree *tr )
{
   if( !AddCommonBranches(tr) ) return false;
   tr->Branch("CBARCandHitZ",m_candHitZ,"CBARCandHitZ[CBARCandNumber]/D");
   return true;
}

void VetoCBAR::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t zcenter = MTBP::CBARZPosition + MTBP::CBARLength/2.;
   const Double_t csiz    = MTBP::CSIZPosition + 20.;

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<1.) continue;

      const bool isInner = (m_modId[imod] < 32 );
      const Double_t radius     = (isInner) ? MTBP::CBARInnerInnerR : MTBP::CBAROuterInnerR;

      /// overflow treatment
      Double_t hitz       = m_hitz[imod];
      if( hitz < -MTBP::CBARLength/2. ) hitz = -MTBP::CBARLength/2.;
      if( hitz >  MTBP::CBARLength/2. ) hitz =  MTBP::CBARLength/2.;

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


/*
void VetoCBAR::UpdateVars( const KL2pi0g &kl )
{
   Reset();
   const Double_t zcenter = MTBP::CBARZPosition + MTBP::CBARLength/2.;

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<1.) continue;

      const bool isInner = (m_modId[imod] < 32 );
      const Double_t radius     = (isInner) ? MTBP::CBARInnerInnerR : MTBP::CBAROuterInnerR;

      /// overflow treatment
      Double_t hitz       = m_hitz[imod];
      if( hitz < -MTBP::CBARLength/2. ) hitz = -MTBP::CBARLength/2.;
      if( hitz >  MTBP::CBARLength/2. ) hitz =  MTBP::CBARLength/2.;

      hitz += zcenter;
      const Double_t delz   = hitz - kl.v().z(); 

      const Double_t dist   = TMath::Hypot(radius,delz); 
      const Double_t tof    = dist / (TMath::C()/1.E6);
      const Double_t vtime  = m_time[imod] - tof - kl.t();

      if( !IsInsideCandidateWindow(vtime) ) continue;

      /// output vars ///
      m_candId[m_ncand]   = m_modId[imod];
      m_candEne[m_ncand]  = m_ene[imod];
      m_candTime[m_ncand] = vtime;
      m_candHitZ[m_ncand] = hitz;

      m_ncand++;
   }
}
*/
bool VetoCBAR::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candEne[imod]>GetVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) ){
         return true;
      }
   }
   return false;
}

bool VetoCBAR::IsLooseVeto()
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

bool VetoCBAR::IsTightVeto()
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
