#ifndef VETOMBCV_H
#define VETOMBCV_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoMBCV : public Veto125MHz
{
 public:
   VetoMBCV();
   VetoMBCV( const Int_t userflag);

   ~VetoMBCV();

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
