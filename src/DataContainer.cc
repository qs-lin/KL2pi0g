#include "g5anaKL2pi0g/DataContainer.h"

DataContainer::DataContainer( const Int_t userflag )
   : m_nkl(3), // default
     m_isRealData(false)
{
   m_userflag = userflag;
   Init();
}

DataContainer::~DataContainer()
{
   ;
}

void DataContainer::Init()
{
   /// brought from clustering file ///
   RunID = NodeID = FileID = DstEntryID = ClusterEntryID = -1;
   SpillID = -1;
   TimeStamp = -1;;
   Error = 0;
   isGoodRun = isGoodSpill = false;
   DetectorBit = ScaledTrigBit = 0;
   CDTNum = -1;

   /// save in user_cut ///
   CutCondition = VetoCondition = 0;
   MinGammaE = MinFiducialXY =  MaxFiducialR = MinClusterDistance = MaxDeltaVertexTime = 0.;

   /// save in SetData() ///
   Reset();
}

void DataContainer::Reset()
{
  GamClusNumber = 0;
  for( Int_t i=0; i<s_arrSize; ++i )
  {
     GamClusId[i] = -1;
     GamClusSize[i] = 0;
     GamClusDepE[i] = GamClusTime[i] = GamClusRMS[i] = 0.;
     for( Int_t j=0; j<3; ++j ) GamClusCoePos[i][j] = 0.;
     for( Int_t j=0; j<s_arrSize; ++j )
     {
        GamClusCsiId[i][j] = 0;
        GamClusCsiE[i][j] = 0.;
        GamClusCsiTime[i][j] = 0.;
     }
  }

  GammaNumber = 0;
  for( Int_t i=0; i<s_arrSize; ++i )
  {
     GammaId[i] = -1;
     GammaE[i] = GammaTime[i] = GammaChi2[i] = 0.;
     for( Int_t j=0; j<3; ++j ) GammaPos[i][j] = GammaMom[i][j] = 0.; 
  }

  // pi0 and m34 loop
  for( Int_t i=0; i<s_arrSize; ++i )
  {
     Pi0Id[i] = -1;
     Pi0GammaId[i][0] = Pi0GammaId[i][1] = -1;
     Pi0E[i] = Pi0Pt[i] = Pi0Mass[i] = 0.;
     for( Int_t j=0; j<3; ++j ) Pi0Pos[i][j] = Pi0Mom[i][j] = 0.;

     //M34[i] = 0.;
     //M34GammaId[i][0] = M34GammaId[i][1] = -1;
  }

  KlongNumber = 0;
  for( Int_t i=0; i<s_arrSize; ++i )
  {
     KlongId[i] = -1;
     KlongE[i] = KlongPt[i] = KlongMass[i] = 0.;
     for( Int_t j=0; j<3; ++j ) KlongPos[i][j] = KlongMom[i][j] = 0.; 
  }

  Coe[0] = Coe[1] = 0.;
  AverageClusterTime = 0.;
}

bool DataContainer::SetBranchAddresses( TTree *tr )
{
   if( tr==NULL ) return false;

   m_e14data.setBranchAddress(tr);

   tr->SetBranchAddress( "RunID",     &RunID );
   tr->SetBranchAddress( "NodeID",    &NodeID );
   tr->SetBranchAddress( "FileID",    &FileID );
   tr->SetBranchAddress( "DstEntryID",&DstEntryID );

   /// for real data
   if( tr->GetBranch("DetectorBit") ){
      std::cout<<"[Info] This is real data. " << std::endl;
      m_isRealData = true;
      tr->SetBranchAddress("SpillID",        &SpillID );
      tr->SetBranchAddress("TimeStamp",      &TimeStamp);
      tr->SetBranchAddress("Error",          &Error);
      tr->SetBranchAddress("DetectorBit",    &DetectorBit);
      tr->SetBranchAddress("ScaledTrigBit",  &ScaledTrigBit);
      tr->SetBranchAddress("isGoodRun",      &isGoodRun);
      tr->SetBranchAddress("isGoodSpill",    &isGoodSpill);
   }

   if( tr->GetBranch("CDTNum") ){
      tr->SetBranchAddress("CDTNum",&CDTNum);
   }else{
      std::cout<<" Warning: Do not get CDTNum branch! " << std::endl;
   }

   if( tr->GetBranch("GenParticleNumber") ){
      std::cout<<"[Info] This is MC simulation. " << std::endl;
      m_isRealData = false;
      m_trueData.SetBranchAddresses(tr);
   }

   return true;
}
 
bool DataContainer::AddBranches( TTree *tr )
{
   if( tr==NULL ) return false;

   /// common //
   tr->Branch("RunID",&RunID,"RunID/I" );
   tr->Branch("NodeID",&NodeID,"NodeID/I");
   tr->Branch("FileID",&FileID,"FileID/I" );
   tr->Branch("DstEntryID",&DstEntryID,"DstEntryID/I" );
   //tr->Branch("ClusterEntryID",&ClusterEntryID,"ClusterEntryID/I" );
   tr->Branch("CDTNum",&CDTNum,"CDTNum/I");

   if( IsRealData() ){
      tr->Branch("ScaledTrigBit",&ScaledTrigBit,"ScaledTrigBit/i");
      tr->Branch("TimeStamp",&TimeStamp,"TimeStamp/I");
      tr->Branch("SpillID",&SpillID,"SpillID/S");
      tr->Branch("DetectorBit",&DetectorBit,"DetectorBit/i"); 
   }else{
      m_trueData.AddBranches(tr); 
   }

   /// kinematic ///
   tr->Branch("GamClusNumber",&GamClusNumber,"GamClusNumber/I");
   tr->Branch("GamClusId",GamClusId,"GamClusId[GamClusNumber]/I");
   tr->Branch("GamClusDepE",GamClusDepE,"GamClusDepE[GamClusNumber]/D");
   tr->Branch("GamClusCoePos",GamClusCoePos,"GamClusCoePos[GamClusNumber][3]/D");
   tr->Branch("GamClusTime",GamClusTime,"GamClusTime[GamClusNumber]/D");
   tr->Branch("GamClusRMS",GamClusRMS,"GamClusRMS[GamClusNumber]/D");
   tr->Branch("GamClusSize",GamClusSize,"GamClusSize[GamClusNumber]/I");
   tr->Branch("GamClusCsiId",GamClusCsiId,Form("GamClusCsiId[GamClusNumber][%d]/I",s_arrSize));
   tr->Branch("GamClusCsiE",GamClusCsiE,Form("GamClusCsiE[GamClusNumber][%d]/D",s_arrSize));
   tr->Branch("GamClusCsiTime",GamClusCsiTime,Form("GamClusCsiTime[GamClusNumber][%d]/D",s_arrSize));

   tr->Branch("GammaNumber",&GammaNumber,"GammaNumber/I");
   tr->Branch("GammaId",GammaId,"GammaId[GammaNumber]/I");
   tr->Branch("GammaE",GammaE,"GammaE[GammaNumber]/D");
   tr->Branch("GammaPos",GammaPos,"GammaPos[GammaNumber][3]/D");
   tr->Branch("GammaTime",GammaTime,"GammaTime[GammaNumber]/D");
   tr->Branch("GammaMom",GammaMom,"GammaMom[GammaNumber][3]/D");
   tr->Branch("GammaSigmaE",GammaSigmaE,"GammaSigmaE[GammaNumber]/D");
   tr->Branch("GammaSigmaPos",GammaSigmaPos,"GammaSigmaPos[GammaNumber][3]/D");
   tr->Branch("GammaChi2",GammaChi2,"GammaChi2[GammaNumber]/D");

   tr->Branch("KlongNumber",&KlongNumber,"KlongNumber/I");
   tr->Branch("KlongId",KlongId,"KlongId[KlongNumber]/I");
   tr->Branch("KlongE",KlongE,"KlongE[KlongNumber]/D");
   tr->Branch("KlongPos",KlongPos,"KlongPos[KlongNumber][3]/D");
   tr->Branch("KlongMom",KlongMom,"KlongMom[KlongNumber][3]/D");
   tr->Branch("KlongPt",KlongPt,"KlongPt[KlongNumber]/D");
   tr->Branch("KlongMass",KlongMass,"KlongMass[KlongNumber]/D");
   tr->Branch("KlongTime",KlongTime,"KlongTime[KlongNumber]/D");
   //tr->Branch("K2pi0MassChisq",K2pi0MassChisq,"K2pi0MassChisq[KlongNumber]/D");

   tr->Branch("Pi0Id",Pi0Id,"Pi0Id[KlongNumber]/I");
   tr->Branch("Pi0E",Pi0E,"Pi0E[KlongNumber]/D");
   tr->Branch("Pi0Pos",Pi0Pos,"Pi0Pos[KlongNumber][3]/D");
   tr->Branch("Pi0Mom",Pi0Mom,"Pi0Mom[KlongNumber][3]/D");
   tr->Branch("Pi0Pt",Pi0Pt,"Pi0Pt[KlongNumber]/D");
   tr->Branch("Pi0Mass",Pi0Mass,"Pi0Mass[KlongNumber]/D");
   tr->Branch("Pi0SigmaMass",Pi0SigmaMass,"Pi0SigmaMass[KlongNumber]/D");
   tr->Branch("Pi0GammaId",Pi0GammaId,"Pi0GammaId[KlongNumber][2]/I");

/*
   tr->Branch("M34",M34,"M34[KlongNumber]/D");
   tr->Branch("SigmaM34",SigmaM34,"SigmaM34[KlongNumber]/D");
   tr->Branch("M34GammaId",M34GammaId,"M34GammaId[KlongNumber][2]/I");
*/
   tr->Branch("AverageClusterTime",&AverageClusterTime,"AverageClusterTime/D");

   /// cut vars ///
   tr->Branch("MinGammaE",&MinGammaE,"MinGammaE/D");
   tr->Branch("MinFiducialXY",&MinFiducialXY,"MinFiducialXY/D");
   tr->Branch("MaxFiducialR",&MaxFiducialR,"MaxFiducialR/D");
   tr->Branch("MinClusterDistance",&MinClusterDistance,"MinClusterDistance/D");
   tr->Branch("Coe",Coe,"Coe[2]/D");
   tr->Branch("MaxDeltaVertexTime",&MaxDeltaVertexTime,"MaxDeltaVertexTime/D");

   tr->Branch("CutCondition",&CutCondition,"CutCondition/i");
   tr->Branch("VetoCondition",&VetoCondition,"VetoCondition/i");

   return true;
}

bool DataContainer::IsPhysTrig() const
{
   return (m_userflag<20180601) ? (ScaledTrigBit&0x101)==0x101 : (ScaledTrigBit&0x1)==0x1;
}

bool DataContainer::IsGoodQuality() const
{
   return (isGoodRun && isGoodSpill && Error==0);
}

/*
void DataContainer::SetData( const std::vector<KLpi0gg> &klvec )
{
   const UInt_t n_save = (m_nkl > klvec.size()) ? klvec.size() : m_nkl;
   Reset();
   KlongNumber = n_save;
   GamClusNumber = GammaNumber = 4;
   for( UInt_t i=0; i<n_save; ++i )
   {
      SetData(klvec.at(i), i); 
   }
   
   if( !IsRealData() ) m_trueData.UpdateVars();
 
}

void DataContainer::SetData( const KLpi0gg &kl, const Int_t arr_id )
{
   KlongId[arr_id] = kl.id();
   KlongMass[arr_id] = kl.m();
   KlongE[arr_id] = kl.e();
   KlongPt[arr_id] = kl.p3().perp();
   for(int i=0;i<3;i++){
      KlongMom[arr_id][i] = kl.p3()(i);
      KlongPos[arr_id][i] = kl.v()(i);
  }
  KlongTime[arr_id] = kl.t();
  K2pi0MassChisq[arr_id] = kl.k2pi0Chisq();

  SetData(kl.pi0(),arr_id);
  M34[arr_id] = kl.m34();
  M34GammaId[arr_id][0] = kl.g3().id();
  M34GammaId[arr_id][1] = kl.g4().id();

  Coe[0] = kl.coe().x();
  Coe[1] = kl.coe().y();

  Pi0SigmaMass[arr_id] = kl.sigmaM12();
  SigmaM34[arr_id] = kl.sigmaM34();

  /// for gammas (clusters), only the best one is saved.
  if( arr_id == 0 ){ 
     SetData(kl.pi0().g1(),0);
     SetData(kl.pi0().g2(),1);
     SetData(kl.g3(),2);
     SetData(kl.g4(),3);
     AverageClusterTime = kl.csi_t();
  }

}
*/

void DataContainer::SetData( const Pi0 &pi0, const Int_t arr_id )
{
   Pi0Id[arr_id] = pi0.id();
   Pi0Mass[arr_id] = pi0.m();
   Pi0E[arr_id] = pi0.e();
   Pi0Pt[arr_id] = pi0.p3().perp();

   for(int i=0;i<3;i++){
      Pi0Mom[arr_id][i]=pi0.p3()(i);
      Pi0Pos[arr_id][i]=pi0.v()(i);
   }

   Pi0GammaId[arr_id][0] = pi0.g1().id();
   Pi0GammaId[arr_id][1] = pi0.g2().id();

}

void DataContainer::SetData( const Gamma &g, const Int_t arr_id )
{
   GammaId[arr_id] = g.id();
   GammaE[arr_id] = g.e();
   GammaTime[arr_id] = g.t();
   GammaSigmaE[arr_id] = g.sigmaE();

   for( Int_t i=0; i<3; ++i )
   {
      GammaPos[arr_id][i] = g.pos()(i);
      GammaMom[arr_id][i] = g.p3()(i);
      GammaSigmaPos[arr_id][i] = g.sigmaPos()(i);
   }

   GammaChi2[arr_id] = g.chisq();
   SetData(g.cluster(), arr_id );
}

void DataContainer::SetData( const Cluster &clus, const Int_t arr_id )
{
   GamClusId[arr_id] = clus.id();

   GamClusDepE[arr_id] = clus.e();
   GamClusTime[arr_id] = clus.t();

   for(int i=0; i<3; ++i)
   {
      GamClusCoePos[arr_id][i] = clus.pos()(i);
   }

   GamClusRMS[arr_id] = clus.rms();
   GamClusSize[arr_id] = clus.clusterIdVec().size();
    
   if(GamClusSize[arr_id]>s_arrSize){
      GsimMessage::getInstance()->
               report("warning",Form("more than %d crystals in the GammaCluster.",s_arrSize));
      GamClusSize[arr_id] = s_arrSize;
   }

   for(int i=0;i<s_arrSize;i++){
      if( i<GamClusSize[arr_id] ){
         GamClusCsiId[arr_id][i] = clus.clusterIdVec()[i];
         GamClusCsiE[arr_id][i] = clus.clusterEVec()[i];
         GamClusCsiTime[arr_id][i] = clus.clusterTimeVec()[i];
      }else{
        GamClusCsiId[arr_id][i] = 0;
        GamClusCsiE[arr_id][i] = 0;
        GamClusCsiTime[arr_id][i] = 0;
      }
    }
}  
