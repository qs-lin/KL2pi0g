#include "g5anaKL2pi0g/Veto.h"

Veto::Veto()
{
   Init();
}
   
Veto::Veto(const std::string detname, const Int_t userflag )
   : m_detname(detname),
     m_userflag(userflag)
{
   Init();
}

void Veto::Init()
{
   m_isSaveProperTime = false;
   m_isSaveMaxEneHit = false;
}

Veto::~Veto()
{
   ;
}

bool Veto::AddCommonBranches( TTree *tr )
{
   if( tr==NULL ) return false;
   tr->Branch(Form("%sCandNumber",m_detname.c_str()), &m_ncand,
              Form("%sCandNumber/I",m_detname.c_str()));
   tr->Branch(Form("%sCandModId",m_detname.c_str()), m_candId,
              Form("%sCandModId[%sCandNumber]/I", m_detname.c_str(), m_detname.c_str()) );
   tr->Branch(Form("%sCandEne",m_detname.c_str()), m_candEne,
              Form("%sCandEne[%sCandNumber]/D", m_detname.c_str(), m_detname.c_str()) );
   tr->Branch(Form("%sCandTime",m_detname.c_str()), m_candTime,
              Form("%sCandTime[%sCandNumber]/D", m_detname.c_str(), m_detname.c_str()) );

   if( m_isSaveProperTime ){
      tr->Branch(Form("%sProperArrId",m_detname.c_str()), &m_properArrId,
                 Form("%sProperArrId/I",m_detname.c_str()));   
      tr->Branch(Form("%sProperTime",m_detname.c_str()), &m_properTime,
                 Form("%sProperTime/D",m_detname.c_str()));
   }

   if( m_isSaveMaxEneHit ){
      tr->Branch(Form("%sMaxEneHitArrId",m_detname.c_str()), &m_maxEneArrId,
                 Form("%sMaxEneHitArrId/I",m_detname.c_str()));
      tr->Branch(Form("%sMaxEne",m_detname.c_str()), &m_maxEne,
                 Form("%sMaxEne/D",m_detname.c_str()));
   }

   return true;
}

void Veto::CommonReset()
{
   m_ncand = 0;

   for( Int_t i=0; i<s_arrSize; ++i )
   {
      m_candId[i] = -1;
      m_candEne[i] = m_candTime[i] = 0.;
   }

   m_properArrId = m_maxEneArrId = -1;
   m_properTime = -9999.;
   m_maxEne = -1.;

}

void Veto::EvalProperTime()
{
   Double_t min_tdiff = 1000.;

   m_properTime = -9999.;
   m_properArrId = -1;

   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candEne[imod] < GetVetoEneThreshold() ) continue;
      Double_t tdiff = TMath::Abs(m_candTime[imod] - GetT0() );
      if( tdiff < min_tdiff ){
         m_properArrId = imod;
         m_properTime = m_candTime[imod] - GetT0();
         min_tdiff = tdiff;
      }
   }
}

void Veto::EvalMaxEneHit()
{
   m_maxEne = -1.;
   m_maxEneArrId = -1;

   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( !IsInsideVetoWindow(m_candTime[imod]) ) continue;
      if( m_candEne[imod] > m_maxEne ){
         m_maxEneArrId = imod;
         m_maxEne = m_candEne[imod];
      }
   }
}
