#ifndef VETOBHGC_H
#define VETOBHGC_H

#include <map>

#include "g5anaKL2pi0g/Veto500MHz.h"

class VetoBHGC : public Veto500MHz
{
 public:
   VetoBHGC();
   VetoBHGC( const Int_t userflag);

    ~VetoBHGC();

   void Init();
 
   virtual void Reset();
   virtual bool SetBranchAddresses( TTree *tr );
   virtual bool AddBranches( TTree *tr );
   virtual void UpdateVars( const KL2pi0g &kl );
   virtual bool IsVeto();
   virtual bool IsLooseVeto();
   virtual bool IsTightVeto();

   std::map<Int_t, std::pair<Int_t, Int_t> > GetModIdMap() const;

   void      SetCoinWindowWidth( const Double_t w ){ m_coinWindowWidth = w; }
   Double_t  GetCoinWindowWidth() const { return m_coinWindowWidth; }
   bool      IsCoinHitPair( const Double_t t1, const Double_t t2 ) const;

 private:
   Double_t      m_coinWindowWidth;

};

#endif
