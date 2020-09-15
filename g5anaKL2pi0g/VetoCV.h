#ifndef VETOCV_H
#define VETOCV_H

#include <map>

#include "g5anaKL2pi0g/Veto125MHz.h"

class VetoCV : public Veto125MHz
{
 public:
   VetoCV();
   VetoCV( const Int_t userflag);

   ~VetoCV();

   void Init();
 
   virtual void Reset();
   virtual bool SetBranchAddresses( TTree *tr );
   virtual bool AddBranches( TTree *tr );
   virtual void UpdateVars( const KL2pi0g &kl );
   virtual bool IsVeto();
   virtual bool IsLooseVeto();
   virtual bool IsTightVeto();

   Double_t     GetStripX( const Int_t id ) const;
   Double_t     GetStripY( const Int_t id ) const;

 private:
   typedef std::pair<Double_t, Double_t> metric_t;
   void         LoadCVGeometry();
   std::map<Int_t, metric_t> m_xymap;

};

#endif
