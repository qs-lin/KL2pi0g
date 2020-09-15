#ifndef VETOIB_H
#define VETOIB_H

#include "g5anaKL2pi0g/Veto500MHz.h"

class VetoIB : public Veto500MHz
{
 public:
   VetoIB();
   VetoIB( const Int_t userflag);

    ~VetoIB();

   void Init();
 
   virtual void Reset();
   virtual bool SetBranchAddresses( TTree *tr );
   virtual bool AddBranches( TTree *tr );
   virtual void UpdateVars( const KL2pi0g &kl );
   virtual bool IsVeto();
   virtual bool IsLooseVeto();
   virtual bool IsTightVeto();

   struct ModHit{
      Int_t    id;
      Double_t e;
      Double_t t;
      Double_t z;
      Bool_t   single;
   };

   std::vector<ModHit> GetModHitVec() const;

 private:
   Double_t GetHitTime( const Double_t utime, const Double_t dtime ) const;
   Double_t GetHitZ   ( const Double_t utime, const Double_t dtime ) const;
   Double_t GetHitEne ( const Double_t uene , const Double_t dene  ,
                        const Double_t utime, const Double_t dtime ) const;
   Double_t GetVetoTime( const Double_t hitz, const Double_t hit_time, 
                         const Double_t csi_time ) const;

   bool m_isDeadChannel;
   Double_t m_zcenter;

   void EvalSingleProperTime();

 private:
   /// input branch (additional)
   //Float_t m_hitz[s_arrSize]; 

   /// output branch (additional)
   Double_t m_candHitZ[s_arrSize]; 
   Double_t m_candBackTime[s_arrSize];

   static const Int_t s_hitSize = 300;
   Int_t    m_nsingle;
   Int_t    m_singleId[s_hitSize];
   Double_t m_singleEne[s_hitSize];
   Double_t m_singleTime[s_hitSize];

   Double_t m_singleProperTime;
   Int_t    m_singleProperArrId;
};

#endif
