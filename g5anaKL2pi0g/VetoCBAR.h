#ifndef VETOCBAR_H
#define VETOCBAR_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCBAR : public Veto125MHz
{
 public:
   VetoCBAR();
   VetoCBAR( const Int_t userflag);

   ~VetoCBAR();

  void Init();
 
  virtual void Reset();
  virtual bool SetBranchAddresses( TTree *tr );
  virtual bool AddBranches( TTree *tr );
  virtual void UpdateVars( const KL2pi0g &kl );
  virtual bool IsVeto();
  virtual bool IsLooseVeto();
  virtual bool IsTightVeto();

 private:
   /// input branch (additional)
   Float_t m_hitz[s_arrSize]; 

   /// output branch (additional)
   Double_t m_candHitZ[s_arrSize]; 
};

#endif
