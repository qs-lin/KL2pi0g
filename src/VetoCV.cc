#include "g5anaKL2pi0g/VetoCV.h"

#include "MTAnalysisLibrary/MTBasicParameters.h"

VetoCV::VetoCV()
   : Veto125MHz("CV",20180601)
{
    Init();
}

VetoCV::VetoCV( const Int_t userflag )
   : Veto125MHz("CV",userflag)
{
    Init();
}

VetoCV::~VetoCV()
{
   ;
}

void VetoCV::Init()
{
   const Double_t t0 = 54.1;
   SetT0(t0);

   SetSuppEneThreshold(.1);
   //SetLooseVetoEneThreshold(.4);
   SetVetoEneThreshold(.2);
   //SetTightVetoEneThreshold(.3);

   SetCandidateWindow( t0 - 30., t0 + 30.);
   SetVetoWindow( t0 - 8., t0 + 8.);
   LoadCVGeometry();
}

void VetoCV::Reset()
{
   CommonReset();
}

bool VetoCV::SetBranchAddresses( TTree *tr )
{
   return SetCommonBranches(tr);
}

bool VetoCV::AddBranches( TTree *tr )
{
   return AddCommonBranches(tr);
}


void VetoCV::UpdateVars( const KL2pi0g &kl )
{
   Reset();

   const Double_t front_z = MTBP::CVFZPosition + 20.;
   const Double_t rear_z  = MTBP::CVRZPosition + 20.;

   for( Int_t imod=0; imod<m_nmod; ++imod )
   {
      if( m_ene[imod]<GetSuppEneThreshold() ) continue;

      const Int_t id = m_modId[imod];
      const Double_t cv_z = (id<100) ? front_z : rear_z;

      const Double_t radius = TMath::Hypot(kl.v().x()-GetStripX(id), kl.v().y()-GetStripY(id));
      const Double_t dist   = TMath::Hypot(radius, kl.v().z() - cv_z);
      const Double_t tof    = dist / (TMath::C()/1.E6);
      const Double_t vtime  = m_time[imod] - tof - kl.t();

      if( !IsInsideCandidateWindow(vtime) ) continue;

      /// output vars ///
      m_candId[m_ncand]   = m_modId[imod];
      m_candEne[m_ncand]  = m_ene[imod];
      m_candTime[m_ncand] = vtime;

      m_ncand++;
   }

   EvalProperTime();
   EvalMaxEneHit();
}

bool VetoCV::IsVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if( m_candEne[imod]>GetVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) ){
         return true;
      }
   }
   return false;
}

bool VetoCV::IsLooseVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if(   m_candEne[imod]>GetLooseVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
      {
         return true;
      }
   }
   return false;
}

bool VetoCV::IsTightVeto()
{
   for( Int_t imod=0; imod<m_ncand; ++imod )
   {
      if(   m_candEne[imod]>GetTightVetoEneThreshold() && IsInsideVetoWindow(m_candTime[imod]) )
      {
         return true;
      }
   }
   return false;
}

Double_t VetoCV::GetStripX( const Int_t id ) const
{
   std::map<Int_t, metric_t>::const_iterator it = m_xymap.find(id);
   if( it!=m_xymap.end() ) return (it->second).first;
   
   std::cout<<" <Warning> VetoCV::GetStripX: unknown id = " << id << std::endl;
   return 0.;
}
   
Double_t VetoCV::GetStripY( const Int_t id ) const
{
   std::map<Int_t, metric_t>::const_iterator it = m_xymap.find(id);
   if( it!=m_xymap.end() ) return (it->second).second;

   std::cout<<" <Warning> VetoCV::GetStripX: unknown id = " << id << std::endl;
   return 0.;
}

///
void VetoCV::LoadCVGeometry()
{
   const std::string dirName = std::string( std::getenv(MTBP::EnvName.c_str()) ) 
                                + "/AnalysisLibrary/LegacyProjects/MTAnalysisLibrary";
   const std::string mapFile = dirName + "/data/CV_map.csv";

  if( access( mapFile.c_str(), R_OK) != 0 ){
    std::cout << mapFile.c_str() << " is not found." << std::endl;
    exit(1);
  }

  std::string header;
  std::ifstream ifs( mapFile.c_str() );
  std::getline( ifs, header );

  int dummyI, id, quad;
  double dummyD, x_min, x_max, y_min, y_max, z_min, z_max;
  
   while(  ifs >> id >> dummyI >> dummyI >> x_min >> x_max
               >> y_min >> y_max >> z_min >> z_max
               >> dummyD >> dummyD >> dummyD >> dummyD >> quad )
   {
      Double_t x_center = ( x_min + x_max ) / 2.;
      Double_t y_center = ( y_min + y_max ) / 2.; 
      metric_t xy(x_center, y_center);
      m_xymap.insert( std::make_pair(id, xy) );
  }

  ifs.close();
}
