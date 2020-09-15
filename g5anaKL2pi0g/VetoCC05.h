#ifndef VETOCC05_H
#define VETOCC05_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCC05 : public Veto125MHz
{
 public:
   VetoCC05();
   VetoCC05( const Int_t userflag);

   ~VetoCC05();

  void Init();
 
  virtual void Reset();
  virtual bool SetBranchAddresses( TTree *tr );
  virtual bool AddBranches( TTree *tr );
  virtual void UpdateVars( const KL2pi0g &kl );
  virtual bool IsVeto();
  virtual bool IsLooseVeto();
  virtual bool IsTightVeto();

 private:
};

#endif
