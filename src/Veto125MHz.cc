#include "g5anaKL2pi0g/Veto125MHz.h"

Veto125MHz::Veto125MHz()
{
   ;
}

Veto125MHz::Veto125MHz( const std::string detname, const Int_t userflag )
   : Veto(detname,userflag)
{
   ;
}

Veto125MHz::~Veto125MHz()
{
   ;
}

bool Veto125MHz::SetCommonBranches( TTree *tr )
{
   if( tr==NULL ) return false;
   tr->SetBranchAddress(Form("%sModuleNumber",m_detname.c_str()),&m_nmod);
   tr->SetBranchAddress(Form("%sModuleModID",m_detname.c_str()), m_modId);
   tr->SetBranchAddress(Form("%sModuleEne",m_detname.c_str()), m_ene);
   tr->SetBranchAddress(Form("%sModuleHitTime",m_detname.c_str()), m_time);

   return true; 
}

