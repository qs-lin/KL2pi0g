
#include "g5anaKL2pi0g/VetoFBAR.h"
#include "g5anaKL2pi0g/VetoNCC.h"
#include "g5anaKL2pi0g/VetoCBAR.h"
#include "g5anaKL2pi0g/VetoIB.h"
#include "g5anaKL2pi0g/VetoCV.h"
#include "g5anaKL2pi0g/VetoOEV.h"
#include "g5anaKL2pi0g/VetoIBCV.h"
#include "g5anaKL2pi0g/VetoMBCV.h"
#include "g5anaKL2pi0g/VetoBHPV.h"
#include "g5anaKL2pi0g/VetoCC03.h"
#include "g5anaKL2pi0g/VetoCC04.h"
#include "g5anaKL2pi0g/VetoCC05.h"
#include "g5anaKL2pi0g/VetoCC06.h"
#include "g5anaKL2pi0g/VetoBHGC.h"
#include "g5anaKL2pi0g/VetoCSI.h"

void user_veto(std::vector<Veto*> &vetovec, const Int_t userflag )
{
   enum{ k_FBAR=0, k_NCC,   k_CBAR,   k_IB, 
         k_CV,     k_OEV,   k_IBCV,   k_MBCV,
         k_BHPV,   k_CC03,  k_CC04,   k_CC05, 
         k_CC06,   k_BHGC,  k_CSI };

   /// FBAR ///
   Veto *FBAR = new VetoFBAR(userflag);
   FBAR->SetId(k_FBAR);
   vetovec.push_back(FBAR);

   /// NCC ///
   Veto *NCC = new VetoNCC(userflag);
   NCC->SetId(k_NCC);
   vetovec.push_back(NCC);

   /// CBAR ///
   Veto *CBAR = new VetoCBAR(userflag);
   CBAR->SetId(k_CBAR);
   vetovec.push_back(CBAR);

   /// IB ///
   Veto *IB = new VetoIB(userflag);
   IB->SetId(k_IB);
   vetovec.push_back(IB);

   /// CV ///
   Veto *CV = new VetoCV(userflag);
   CV->SetId(k_CV);
   vetovec.push_back(CV);

   /// OEV ///
   Veto *OEV = new VetoOEV(userflag);
   OEV->SetId(k_OEV);
   vetovec.push_back(OEV);

   /// IBCV ///
   Veto *IBCV = new VetoIBCV(userflag);
   IBCV->SetId(k_IBCV);
   vetovec.push_back(IBCV);

   /// MBCV ///
   Veto *MBCV = new VetoMBCV(userflag);
   MBCV->SetId(k_MBCV);
   vetovec.push_back(MBCV);

   /// BHPV ///
   Veto *BHPV = new VetoBHPV(userflag);
   BHPV->SetId(k_BHPV);
   vetovec.push_back(BHPV);

   /// CC03 ///
   Veto *CC03 = new VetoCC03(userflag);
   CC03->SetId(k_CC03);
   vetovec.push_back(CC03);

   /// CC04 ///
   Veto *CC04 = new VetoCC04(userflag);
   CC04->SetId(k_CC04);
   vetovec.push_back(CC04);

   /// CC05 ///
   Veto *CC05 = new VetoCC05(userflag);
   CC05->SetId(k_CC05);
   vetovec.push_back(CC05);

   /// CC06 ///
   Veto *CC06 = new VetoCC06(userflag);
   CC06->SetId(k_CC06);
   vetovec.push_back(CC06);

   /// BHGC ///
   Veto *BHGC = new VetoBHGC(userflag);
   BHGC->SetId(k_BHGC);
   vetovec.push_back(BHGC);

   /// CSI ///
   Veto *CSI = new VetoCSI(userflag);
   CSI->SetId(k_CSI);
   vetovec.push_back(CSI);

   /// common setting ///
   for(std::vector<Veto*>::iterator it = vetovec.begin(); it!=vetovec.end(); ++it )
   {
      if( (*it)->GetId()==k_BHPV || (*it)->GetId()==k_CSI ) continue;
      (*it)->SetSaveProperTime();
      (*it)->SetSaveMaxEneHit();
   }


}
