#ifndef DATACONTAINER_H
#define DATACONTAINER_H

#include <vector>

#include "TObject.h"
#include "TTree.h"

#include "gnana/E14GNAnaDataContainer.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"

//#include "g5anaKL2pi0g/KLpi0gg.h"
#include "g5anaKL2pi0g/TrueDataHandler.h"

class DataContainer
{
 public:
   DataContainer( const Int_t userflag = 20180601 );
   ~DataContainer();

   void Init();
   void Reset();

   bool SetBranchAddresses( TTree *tr );
   bool AddBranches( TTree *tr );

   bool IsRealData() const { return m_isRealData; }
   bool IsPhysTrig() const;
   bool IsGoodQuality() const; 

   void GetData( std::list<Cluster> &clist ){ m_e14data.getData(clist); }

   //void SetData( const std::vector<KLpi0gg> &klvec );
   //void SetData( const KLpi0gg &kl, const Int_t arr_id );
   void SetData( const Pi0 &pi0, const Int_t arr_id );
   void SetData( const Gamma &g, const Int_t arr_id );
   void SetData( const Cluster &clus, const Int_t arr_id );

   void SetNKlong( const UInt_t nkl ){ m_nkl = nkl; }

   Int_t userflag() const{ return m_userflag;}
   UInt_t nkl() const{ return m_nkl;}

 private:
   Int_t  m_userflag;
   UInt_t m_nkl;

   E14GNAnaDataContainer m_e14data;
   TrueDataHandler       m_trueData;

   bool   m_isRealData;
   
 public: 
  /// tree vars ///

  /// common ///
  Int_t           RunID;
  Int_t           NodeID;
  Int_t           FileID;
  Int_t           DstEntryID;
  Int_t           ClusterEntryID;

  Short_t         SpillID;
  Int_t           TimeStamp;
  Short_t         Error;
  Bool_t          isGoodRun;
  Bool_t          isGoodSpill;
  UInt_t          DetectorBit;
  UInt_t          ScaledTrigBit;
  Int_t           CDTNum;

  /// kinematic ///
  static int const s_arrSize = 120;

  Int_t           GamClusNumber;
  Int_t           GamClusId[s_arrSize];   //[GamClusNumber]
  Double_t        GamClusDepE[s_arrSize];   //[GamClusNumber]
  Double_t        GamClusCoePos[s_arrSize][3];   //[GamClusNumber]
  Double_t        GamClusTime[s_arrSize];   //[GamClusNumber]
  Double_t        GamClusRMS[s_arrSize];  //[GamClusNumber]
  Int_t           GamClusSize[s_arrSize]; //[GamClusNumber]
  Int_t           GamClusCsiId[s_arrSize][s_arrSize];
  Double_t        GamClusCsiE[s_arrSize][s_arrSize];
  Double_t        GamClusCsiTime[s_arrSize][s_arrSize];

  Int_t           GammaNumber;
  Int_t           GammaId[s_arrSize];   //[GammaNumber]
  Double_t        GammaE[s_arrSize];   //[GammaNumber]
  Double_t        GammaPos[s_arrSize][3];   //[GammaNumber]
  Double_t        GammaTime[s_arrSize];   //[GammaNumber]
  Double_t        GammaMom[s_arrSize][3];   //[GammaNumber]
  Double_t        GammaSigmaE[s_arrSize];   //[GammaNumber]
  Double_t        GammaSigmaPos[s_arrSize][3];   //[GammaNumber]
  Double_t        GammaChi2[s_arrSize];  //[GammaNumber]

  Int_t           Pi0Id[s_arrSize]; //[KlongNumber]
  Double_t        Pi0E[s_arrSize]; //[KlongNumber]
  Double_t        Pi0Pos[s_arrSize][3];   //[KlongNumber][3]
  Double_t        Pi0Mom[s_arrSize][3];   //[KlongNumber][3]
  Double_t        Pi0Pt[s_arrSize]; //[KlongNumber]
  Double_t        Pi0Mass[s_arrSize]; //[KlongNumber]
  Double_t        Pi0SigmaMass[s_arrSize];
  Int_t           Pi0GammaId[s_arrSize][2]; //[KlongNumber]

  Int_t           KlongNumber;
  Int_t           KlongId[s_arrSize]; //[KlongNumber]
  Double_t        KlongE[s_arrSize]; //[KlongNumber]
  Double_t        KlongPos[s_arrSize][3];   //[KlongNumber]
  Double_t        KlongMom[s_arrSize][3];  //[KlongNumber]
  Double_t        KlongPt[s_arrSize]; //[KlongNumber]
  Double_t        KlongMass[s_arrSize]; //[KlongNumber]
  Double_t        KlongTime[s_arrSize];
  ////Double_t        K2pi0MassChisq[s_arrSize];

/*
  Double_t        M34[s_arrSize]; //[KlongNumber]
  Double_t        SigmaM34[s_arrSize];
  Int_t           M34GammaId[s_arrSize][2]; //[KlongNumber]
*/
  Double_t        AverageClusterTime;

  Double_t        MinGammaE;
  Double_t        MinFiducialXY;
  Double_t        MaxFiducialR;
  Double_t        MinClusterDistance;
  Double_t        Coe[2];
  Double_t        MaxDeltaVertexTime;

  UInt_t          CutCondition;
  UInt_t          VetoCondition;

};

#endif
