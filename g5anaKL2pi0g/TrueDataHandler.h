#ifndef TRUEDATAHANDLER_H
#define TRUEDATAHANDLER_H

#include "TObject.h"
#include "TTree.h"

class TrueDataHandler
{
 public:
   TrueDataHandler();
   ~TrueDataHandler();

   void Init();

   bool SetBranchAddresses( TTree *tr );
   bool AddBranches( TTree *Tr );

   void Reset();
   void UpdateVars();

 private:
   void CountEndPos();

 private:
   static const Int_t s_arrSize = 100;
   Int_t    m_num;
   Int_t    m_track[s_arrSize];
   Int_t    m_mother[s_arrSize];
   Int_t    m_pid[s_arrSize];
   Double_t m_mom[s_arrSize][3];
   Double_t m_ek[s_arrSize];
   Double_t m_mass[s_arrSize];
   Double_t m_time[s_arrSize];
   Double_t m_pos[s_arrSize][3];
   Double_t m_endTime[s_arrSize];
   Double_t m_endPos[s_arrSize][3];

   Double_t m_decayPos[3];
   Int_t    m_ngOnCsi;
   Int_t    m_ngInHole;
   Int_t    m_ng12m;

};

#endif
