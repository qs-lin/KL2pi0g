#include "g5anaKL2pi0g/TrueDataHandler.h"

#include "TMath.h"

TrueDataHandler::TrueDataHandler()
{
   Init();
}

TrueDataHandler::~TrueDataHandler()
{
   ;
}

void TrueDataHandler::Init()
{
   Reset();
}

bool TrueDataHandler::SetBranchAddresses( TTree *tr )
{
   if( tr==NULL ) return false;
   tr->SetBranchAddress("GenParticleNumber"  ,&m_num     );
   tr->SetBranchAddress("GenParticleTrack"   , m_track   );
   tr->SetBranchAddress("GenParticleMother"  , m_mother  );
   tr->SetBranchAddress("GenParticlePid"     , m_pid     );
   tr->SetBranchAddress("GenParticleMom"     , m_mom     );
   tr->SetBranchAddress("GenParticleEk"      , m_ek      );
   tr->SetBranchAddress("GenParticleMass"    , m_mass    );
   tr->SetBranchAddress("GenParticleTime"    , m_time    );
   tr->SetBranchAddress("GenParticlePos"     , m_pos     );
   tr->SetBranchAddress("GenParticleEndTime" , m_endTime );
   tr->SetBranchAddress("GenParticleEndPos"  , m_endPos  );
   return true;
}

bool TrueDataHandler::AddBranches( TTree *tr )
{
   if( tr==NULL ) return false;
   tr->Branch("GenParticleNumber"  ,&m_num    ,"GenParticleNumber/I"                       );
   tr->Branch("GenParticleTrack"   , m_track  ,"GenParticleTrack[GenParticleNumber]/I"     );
   tr->Branch("GenParticleMother"  , m_mother ,"GenParticleMother[GenParticleNumber]/I"    );
   tr->Branch("GenParticlePid"     , m_pid    ,"GenParticlePid[GenParticleNumber]/I"       );
   tr->Branch("GenParticleMom"     , m_mom    ,"GenParticleMom[GenParticleNumber][3]/D"    );
   tr->Branch("GenParticleEk"      , m_ek     ,"GenParticleEk[GenParticleNumber]/D"        );
   tr->Branch("GenParticleMass"    , m_mass   ,"GenParticleMass[GenParticleNumber]/D"      );
   tr->Branch("GenParticleTime"    , m_time   ,"GenParticleTime[GenParticleNumber]/D"      );
   tr->Branch("GenParticlePos"     , m_pos    ,"GenParticlePos[GenParticleNumber][3]/D"    );
   tr->Branch("GenParticleEndTime" , m_endTime,"GenParticleEndTime[GenParticleNumber]/D"   );
   tr->Branch("GenParticleEndPos"  , m_endPos ,"GenParticleEndPos[GenParticleNumber][3]/D" );

   tr->Branch("TrueDecayPos"       , m_decayPos  ,"TrueDecayPos[3]/D" );
   tr->Branch("TrueNgammaOnCsi"    ,&m_ngOnCsi   ,"TrueNgammaOnCsi/I" );
   tr->Branch("TrueNgammaInHole"   ,&m_ngInHole  ,"TrueNgammaInHole/I");
  
   //tr->Branch("TrueNgamma12m"      ,&m_ng12m     ,"TrueNgamma12m/I"   );
 
   return true;
}

void TrueDataHandler::Reset()
{
   m_ngOnCsi = m_ngInHole = 0;
   for( Int_t i=0; i<3; ++i ) m_decayPos[i] = 0.;
 
   m_ng12m = 0;
}

void TrueDataHandler::UpdateVars()
{
   Reset();
   Int_t ngamma = 0;
   Double_t  gp[10][3];

   for( Int_t igp=0; igp<m_num; ++igp )
   {
      if( m_pid[igp]==22 ){ // gamma PID 
         if( ngamma>6 ) break;
         for(Int_t i=0; i<3; ++i ) gp[ngamma][i] = m_mom[igp][i];

         ngamma++;
         if( ngamma==1 ){
            for(Int_t i=0; i<3; ++i ) m_decayPos[i] = m_pos[igp][i];
         }     
      }
   }

   if( ngamma>6 ){
      m_ngOnCsi = -1;
      return;
   }

   const Double_t csi_z  = 6168.;
   const Double_t csi_xy = 100. ;
   for( Int_t ig=0; ig<ngamma; ++ig )
   {
      /// fly foward requirement
      if( gp[ig][2] < 0. ) continue;

      Double_t scale = (csi_z - m_decayPos[2]) / gp[ig][2];
      Double_t x     = m_decayPos[0] + gp[ig][0] * scale;
      Double_t y     = m_decayPos[1] + gp[ig][1] * scale;
      Double_t r     = TMath::Hypot(x,y);

      if( r > 950 ) continue;
      if( TMath::Abs(x) < csi_xy && TMath::Abs(y) < csi_xy ){
         m_ngInHole++;
         continue;
      }
      m_ngOnCsi++;
   }  
}

void TrueDataHandler::CountEndPos()
{
   m_ng12m = 0;
  for( Int_t igp=0; igp<m_num; ++igp )
   {
      if( m_pid[igp]==22 ){ // gamma PID
         const Double_t z = m_endPos[igp][2];
         const Double_t r = TMath::Hypot(m_endPos[igp][0],m_endPos[igp][1]);
         if( z > 12050 && r < 500  ){ // BHPV region
            m_ng12m++;
         }
      }
   }
  
}
