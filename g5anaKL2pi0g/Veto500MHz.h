#ifndef VETO500MHZ_H
#define VETO500MHZ_H

#include <string>

#include "TObject.h"

#include "g5anaKL2pi0g/Veto.h"

class Veto500MHz : public Veto
{
 public:
   Veto500MHz();
   Veto500MHz(const std::string detname, const Int_t m_userflag );

   virtual ~Veto500MHz();

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
   static const Int_t s_hitSize = 20;
   Int_t     m_nch;
   Int_t     m_chId[s_arrSize];
   Short_t   m_nhit[s_arrSize];
   Float_t   m_ene[s_arrSize][s_hitSize];
   Float_t   m_time[s_arrSize][s_hitSize];

};

#endif
