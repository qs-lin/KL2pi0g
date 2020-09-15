#ifndef VETOFBAR_H
#define VETOFBAR_H

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoFBAR : public Veto125MHz
{
 public:
   VetoFBAR();
   VetoFBAR( const Int_t userflag);

   ~VetoFBAR();

   void Init();
 
   virtual void Reset();
   virtual bool SetBranchAddresses( TTree *tr );
   virtual bool AddBranches( TTree *tr );
   virtual void UpdateVars( const KL2pi0g &kl );
   virtual bool IsVeto();
   virtual bool IsLooseVeto();
   virtual bool IsTightVeto();

 private:

   bool  m_isDeadChannel; 

};

#endif
