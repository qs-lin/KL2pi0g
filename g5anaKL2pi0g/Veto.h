#ifndef VETO_H
#define VETO_H

#include <string>

#include "TObject.h"
#include "TTree.h"

#include "g5anaKL2pi0g/KL2pi0g.h"

class Veto
{
 public:
   /// constructor
   Veto();
   Veto(const std::string detname, const Int_t userflag );

   /// destructor
   virtual ~Veto();

   /// virtual method
   virtual void Reset() = 0;
   virtual bool SetBranchAddresses( TTree *tr ) = 0; 
   virtual bool AddBranches( TTree *tr ) = 0; 
   virtual void UpdateVars( const KL2pi0g &kl ) = 0;
   virtual bool IsVeto() = 0;
   virtual bool IsLooseVeto() = 0;
   virtual bool IsTightVeto() = 0;

   ///
   void Init();

   ///
   Int_t        GetId() const { return m_id; }
   Double_t     GetT0() const { return m_t0; }
   Double_t     GetSuppEneThreshold() const { return m_eSuppThreshold; }
   Double_t     GetVetoEneThreshold() const { return m_eVetoThreshold; }
   Double_t     GetLooseVetoEneThreshold() const { return m_eLooseVetoThreshold; }
   Double_t     GetTightVetoEneThreshold() const { return m_eTightVetoThreshold; }
   std::string  GetDetectorName() const { return m_detname; }
   bool         GetIsSaveProperTime() const { return m_isSaveProperTime; }

   /// 
   void         SetId( const Int_t id ){ m_id = id; }
   void         SetT0( const Double_t t0 ){ m_t0 = t0; }
   void         SetSuppEneThreshold( const Double_t e ){ m_eSuppThreshold = e; }
   void         SetVetoEneThreshold( const Double_t e ){ m_eVetoThreshold = e; }
   void         SetLooseVetoEneThreshold( const Double_t e ){ m_eLooseVetoThreshold = e; }
   void         SetTightVetoEneThreshold( const Double_t e ){ m_eTightVetoThreshold = e; }
   void         SetSaveProperTime( const bool isSaveProperTime = true )
                { m_isSaveProperTime = isSaveProperTime; }
   void         SetSaveMaxEneHit( const bool isSaveMaxEneHit = true )
                { m_isSaveMaxEneHit = isSaveMaxEneHit; }

   void         SetCandidateWindow( const Double_t t1, const Double_t t2 )
                { m_candT1 = t1; m_candT2 = t2; }
   void         SetVetoWindow( const Double_t t1, const Double_t t2 )
                { m_vetoT1 = t1; m_vetoT2 = t2; }

   ///
   bool         IsInsideCandidateWindow( const Double_t t ) const
                { return ( t > m_candT1 && t < m_candT2 ); }
   bool         IsInsideVetoWindow( const Double_t t ) const
                { return ( t > m_vetoT1 && t < m_vetoT2 ); }
   bool         IsSaveProperTime() const { return m_isSaveProperTime; }
   bool         IsSaveMaxEneHit() const { return m_isSaveMaxEneHit; }

 private:
   Int_t        m_id;
   Double_t     m_t0;
   Double_t     m_eSuppThreshold;
   Double_t     m_eVetoThreshold;
   Double_t     m_eLooseVetoThreshold;
   Double_t     m_eTightVetoThreshold;
   bool         m_isSaveProperTime;
   bool         m_isSaveMaxEneHit;

   Double_t     m_candT1;
   Double_t     m_candT2;
   Double_t     m_vetoT1;
   Double_t     m_vetoT2;

 protected:
   std::string m_detname;
   Int_t       m_userflag; 

   void         CommonReset();
   bool         AddCommonBranches( TTree* tr );  
   void         EvalProperTime();
   void         EvalMaxEneHit();

   /// output vars ///
   static const Int_t s_arrSize = 2716;
   Int_t     m_ncand;
   Int_t     m_candId[s_arrSize];
   Double_t  m_candEne[s_arrSize];
   Double_t  m_candTime[s_arrSize];

   Int_t     m_properArrId;
   Double_t  m_properTime;
   Int_t     m_maxEneArrId;
   Double_t  m_maxEne;

};

#endif
