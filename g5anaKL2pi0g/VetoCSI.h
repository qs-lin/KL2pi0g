#ifndef VETOCSI_H
#define VETOCSI_H

#include <map>

#include "LcDetPar/LcDetParCsiEneSuppThreshold.h"

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCSI : public Veto125MHz
{
 public:
   VetoCSI();
   VetoCSI( const Int_t userflag );

   ~VetoCSI();

  void Init();
 
  virtual void Reset();
  virtual bool SetBranchAddresses( TTree *tr );
  virtual bool AddBranches( TTree *tr );
  virtual void UpdateVars( const KL2pi0g &kl );
  virtual bool IsVeto();
  virtual bool IsLooseVeto();
  virtual bool IsTightVeto();

  Double_t    Get3SigmaSuppEneThreshold( const Int_t modId ) const;

  bool        EvalCsiVeto( const Double_t dist,
                           const Double_t ene,
                           const Double_t ene_thr = 3. ) const;

 private:
   LcLib::LcDetParCsiEneSuppThreshold m_esupp;

   std::map<Int_t, Int_t> GetCsiVetoMap( const KL2pi0g &kl ) const;
   void ExcludeCsiInCluster( std::map<Int_t, Int_t> &idmap, 
                             std::vector<Int_t> const& cls_vec ) const;

   Double_t  CalcDistFromCsi( const Double_t x, const Double_t y, const Int_t id ) const;
 
   void GetClusterInfo( const KL2pi0g &kl, Double_t *x, Double_t *y, Double_t *t ) const;

 private:
   /// input branch (additional)
   Float_t m_dist[s_arrSize]; 

   /// output branch (additional)
   Double_t m_candDist[s_arrSize]; 
};

#endif
