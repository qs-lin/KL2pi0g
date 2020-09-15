#ifndef VETOCC04_H
#define VETOCC04_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCC04 : public Veto125MHz
{
 public:
   VetoCC04();
   VetoCC04( const Int_t userflag);

   ~VetoCC04();

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
