#ifndef VETOCC06_H
#define VETOCC06_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCC06 : public Veto125MHz
{
 public:
   VetoCC06();
   VetoCC06( const Int_t userflag);

   ~VetoCC06();

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
