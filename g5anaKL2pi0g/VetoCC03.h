#ifndef VETOCC03_H
#define VETOCC03_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCC03 : public Veto125MHz
{
 public:
   VetoCC03();
   VetoCC03( const Int_t userflag);

   ~VetoCC03();

   void Init();
 
   virtual void Reset();
   virtual bool SetBranchAddresses( TTree *tr );
   virtual bool AddBranches( TTree *tr );
   virtual void UpdateVars( const KL2pi0g &kl );
   virtual bool IsVeto();
   virtual bool IsLooseVeto();
   virtual bool IsTightVeto();

 private:
   Double_t m_properTime;
   
};

#endif
