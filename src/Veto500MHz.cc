#include "g5anaKL2pi0g/Veto500MHz.h"

Veto500MHz::Veto500MHz()
{
   ;
}

Veto500MHz::Veto500MHz( const std::string detname, const Int_t userflag )
   : Veto(detname,userflag)
{
   ;
}

Veto500MHz::~Veto500MHz()
{
   ;
}

bool Veto500MHz::SetCommonBranches( TTree *tr )
{
   if( tr==NULL ) return false;
   tr->SetBranchAddress(Form("%sNumber",m_detname.c_str()),&m_nch);
   tr->SetBranchAddress(Form("%sModID",m_detname.c_str()), m_chId);
   tr->SetBranchAddress(Form("%snHits",m_detname.c_str()), m_nhit);
   tr->SetBranchAddress(Form("%sEne",m_detname.c_str()), m_ene);
   tr->SetBranchAddress(Form("%sTime",m_detname.c_str()), m_time);

   return true; 
}

