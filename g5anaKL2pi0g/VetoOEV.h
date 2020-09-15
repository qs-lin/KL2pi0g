#ifndef VETOOEV_H
#define VETOOEV_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoOEV : public Veto125MHz
{
 public:
   VetoOEV();
   VetoOEV( const Int_t userflag);

   ~VetoOEV();

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
