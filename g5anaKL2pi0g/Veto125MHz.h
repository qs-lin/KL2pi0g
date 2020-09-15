#ifndef VETO125MHZ_H
#define VETO125MHZ_H

#include <string>

#include "TObject.h"

#include "g5anaKL2pi0g/Veto.h"

class Veto125MHz : public Veto
{
 public:
   Veto125MHz();
   Veto125MHz(const std::string detname, const Int_t m_userflag );

   virtual ~Veto125MHz();

   virtual void Reset() = 0;
   virtual bool SetBranchAddresses( TTree *tr ) = 0;
   virtual bool AddBranches( TTree *tr ) = 0;
   virtual void UpdateVars( const KL2pi0g &kl ) = 0;
   virtual bool IsVeto() = 0;
   virtual bool IsLooseVeto() = 0;
   virtual bool IsTightVeto() = 0;

 protected:
   bool SetCommonBranches( TTree *tr );

 protected:
   /// input vars ///
   Int_t     m_nmod;
   Int_t     m_modId[s_arrSize];
   Float_t   m_ene[s_arrSize];
   Float_t   m_time[s_arrSize];

};

#endif
