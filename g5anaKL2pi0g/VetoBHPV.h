#ifndef VETOBHPV_H
#define VETOBHPV_H

#include <map>
#include <vector>
#include <list>

#include "g5anaKL2pi0g/Veto500MHz.h"

class VetoBHPV : public Veto500MHz
{
 public:
   VetoBHPV();
   VetoBHPV( const Int_t userflag);

    ~VetoBHPV();

   void Init();
 
   virtual void Reset();
   virtual bool SetBranchAddresses( TTree *tr );
   virtual bool AddBranches( TTree *tr );
   virtual void UpdateVars( const KL2pi0g &kl );
   virtual bool IsVeto();
   virtual bool IsLooseVeto();
   virtual bool IsTightVeto();

   void         SetNcoinThreshold( const Int_t ncoin ){ m_nCoinThreshold = ncoin; }
   void         SetNcoinTightThreshold( const Int_t ncoin ){ m_nCoinTightThreshold = ncoin; }

   Int_t        GetNcoinThreshold() const { return m_nCoinThreshold; }
   Int_t        GetNcoinTightThreshold() const { return m_nCoinTightThreshold; }

   class ModHit{
    public:
      ModHit() : id(-1), e(0.), t(0.) { ; } 
      Int_t    id;
      Double_t e;
      Double_t t;
      
      bool operator<( const ModHit& hit ) const
      {
         return ( t < hit.t );
      } 
   };

   typedef std::pair<Int_t, std::list<Int_t> > HitId_t;
   typedef std::map<Int_t, HitId_t> ChMap_t;
   std::map<Int_t, std::list<ModHit> > GetModHitMap() const;
   std::map<Int_t, std::list<ModHit> > GetModHitMap( ChMap_t &nmap, ChMap_t &smap) const;

 private:
   Int_t        m_nCoinThreshold;
   Int_t        m_nCoinTightThreshold;

   std::vector<Double_t> m_tofvec;
   void SetTofVec();

 private:
   Int_t    m_candNcoin[s_arrSize]; 
   Double_t m_candTimeSpread[s_arrSize]; 

   Double_t m_proper2CoinTime;
   Double_t m_proper3CoinTime; 

};

#endif
