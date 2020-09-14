#include "g5anaKL2pi0g/KL2pi0g.h"

#include "TMath.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "MTAnalysisLibrary/MTBasicParameters.h"

KL2pi0g::KL2pi0g()
   : m_id(0),
     m_status(0),
     m_sortFlag( 2 ),
     m_coe( CLHEP::Hep3Vector( 0.,0.,0. ) ),
     m_v( CLHEP::Hep3Vector( 0.,0.,0. ) ),
     m_p3( CLHEP::Hep3Vector( 0.,0.,0. ) ),
     m_energy( 0.0 ),
     m_mass( 0.0 ),
     m_t( 0.0 ),
     m_csi_t( 0.0 ),
     m_chisqZ( 9999.0 ),
     m_deltaZ( 9999.0 ),
     m_chisqshape( 9999.0 ),
     m_userflag( 20180601 ),
     m_debugLevel(0)
{
   ;
}

KL2pi0g::~KL2pi0g()
{
   ;
}

std::ostream& operator<<( std::ostream& out, const KL2pi0g& kl2pi0g )
{

  out << "Klong::dump :       " << std::endl
      << "               ID = " << kl2pi0g.m_id          << std::endl
      << "              Vtx = " << kl2pi0g.m_v.x()
      << ", "                   << kl2pi0g.m_v.y()
      << ", "                   << kl2pi0g.m_v.z()       << std::endl
      << "              Coe = " << kl2pi0g.m_coe.x()
      << ", "                   << kl2pi0g.m_coe.y()
      << ", "                   << kl2pi0g.m_coe.z()     << std::endl
      << "               P3 = " << kl2pi0g.m_p3.x()
      << ", "                   << kl2pi0g.m_p3.y()
      << ", "                   << kl2pi0g.m_p3.z()      << std::endl
      << "           E,M,PP = " << kl2pi0g.m_energy
      << ", "                   << kl2pi0g.m_mass
      << ", "                   << kl2pi0g.m_p3.mag()    << std::endl
      << "           status = " << kl2pi0g.m_status      << std::endl
      << "           deltaZ = " << kl2pi0g.m_deltaZ      << std::endl
      << "           chisqZ = " << kl2pi0g.m_chisqZ      << std::endl
      << "         sortFlag = " << kl2pi0g.m_sortFlag    << std::endl
      << "         vertex t = " << kl2pi0g.m_t           << std::endl
      << "            csi t = " << kl2pi0g.m_csi_t       << std::endl
      << "             npi0 = " << kl2pi0g.m_pi0.size()  << std::endl;


    //the following are used to test purpose only. To compare with the old code
/*
  out << "Klong::dump :       " << std::endl
      << "               ID = " << kl2pi0g.m_id          << std::endl
      << "              Vtx = " << kl2pi0g.m_v.x()
      << ", "                   << kl2pi0g.m_v.y()
      << ", "                   << kl2pi0g.m_v.z()       << std::endl
      << "               P3 = " << kl2pi0g.m_p3.x()
      << ", "                   << kl2pi0g.m_p3.y()
      << ", "                   << kl2pi0g.m_p3.z()      << std::endl
      << "           E,M,PP = " << kl2pi0g.m_energy
      << ", "                   << kl2pi0g.m_mass
      << ", "                   << kl2pi0g.m_p3.mag()    << std::endl
      << "           status = " << kl2pi0g.m_status      << std::endl
      << "           deltaZ = " << kl2pi0g.m_deltaZ      << std::endl
      << "           chisqZ = " << kl2pi0g.m_chisqZ      << std::endl
      //<< "         sortFlag = " << kl2pi0g.m_sortFlag    << std::endl
      << "             npi0 = " << kl2pi0g.m_pi0.size()  << std::endl;
*/

      for(std::vector<Pi0>::const_iterator i=kl2pi0g.m_pi0.begin();i!=kl2pi0g.m_pi0.end(); i++) {
        out << *i << std::endl;
      }

      for(std::vector<Pi0>::const_iterator i=kl2pi0g.m_pi0.begin();i!=kl2pi0g.m_pi0.end(); i++) {
        out << i->g1() << std::endl;   
        out << i->g2() << std::endl;   
      }
      
      out << kl2pi0g.m_g5  << std::endl
      << std::endl;

   return( out );
}

bool KL2pi0g::operator<( const KL2pi0g& kl2pi0g ) const
{
  if( compare( kl2pi0g ) < 0 ) {
    return( true );
  }
  return( false );
}

bool KL2pi0g::operator==( const KL2pi0g& kl2pi0g) const
{
  if( compare( kl2pi0g) == 0 ) {
    return( true );
  }
  return( false );

}


KL2pi0g& KL2pi0g::operator=(const KL2pi0g& kl2pi0g)
{
  m_id = kl2pi0g.m_id;
  m_status = kl2pi0g.m_status;
  m_sortFlag = kl2pi0g.m_sortFlag;

  m_v  = kl2pi0g.m_v;
  m_p3 = kl2pi0g.m_p3;
  m_coe = kl2pi0g.m_coe;

  m_energy = kl2pi0g.m_energy;
  m_mass   = kl2pi0g.m_mass;
  m_chisqZ = kl2pi0g.m_chisqZ;
  m_t = kl2pi0g.m_t;
  m_csi_t = kl2pi0g.m_csi_t;
  m_deltaZ = kl2pi0g.m_deltaZ;
  m_chisqshape = kl2pi0g.m_chisqshape;

  m_pi0.clear();
  for(unsigned int i=0;i<kl2pi0g.m_pi0.size();i++){
    m_pi0.push_back(kl2pi0g.m_pi0[i]);
  }

  m_g5 = kl2pi0g.m_g5;  
  return *this;
}

Int_t KL2pi0g::compare( KL2pi0g const& kl2pi0g ) const
{
  switch ( m_sortFlag )
  {

    case SORT_BY_CHISQ :
      if( m_chisqZ < kl2pi0g.m_chisqZ ) {
        return( -1 );
      }
      else if( m_chisqZ > kl2pi0g.m_chisqZ ) {
        return( 1 );
      }
      break;

    case SORT_BY_MASS :
    default :

      const Double_t mass = MTBP::KL_MASS;
      Double_t thisDeltaM2 = (m_mass - mass)*(m_mass - mass);
      Double_t thatDeltaM2 = (kl2pi0g.m_mass - mass)*(kl2pi0g.m_mass - mass);

      if( thisDeltaM2 < thatDeltaM2 ) {
        return( -1 );
      }else if( thisDeltaM2 > thatDeltaM2 ) {
        return( 1 );
      }
      break;
  }
  return 0;
}

Int_t KL2pi0g::UpdateVars()
{
  if(m_pi0.size()!=2) return -1;
  if(m_g5.e()==0) return -1;

  // reconstruct Vtx
  Double_t Etot = 0.0;
  Double_t avrX = 0.0;
  Double_t avrY = 0.0;
  Double_t avrZ = 0.0;
  Double_t sig2tot = 0.0;

  for( std::vector<Pi0>::const_iterator p=m_pi0.begin();
       p!=m_pi0.end(); p++ ) {
    // weighted average
    avrZ    += p->recZ()/p->recZsig2();
    sig2tot += 1./p->recZsig2();

    // center of energy
    avrX += p->g1().e()*p->g1().x() + p->g2().e()*p->g2().x();
    avrY += p->g1().e()*p->g1().y() + p->g2().e()*p->g2().y();

    // Etotal
    Etot += p->g1().e() + p->g2().e();
  }
  // 5th gamma (g5) contribution 
  avrX += m_g5.e()*m_g5.x();
  avrY += m_g5.e()*m_g5.y();
  Etot += m_g5.e();


  avrX = avrX/Etot;
  avrY = avrY/Etot;
  avrZ = avrZ/sig2tot;

  Double_t TargetZ      = -21000.; //mm
  Double_t CalorimatorZ =   6148.; //mm
  if( m_userflag>20160101 ) CalorimatorZ =   6148. + 20.; 

  // the coe on CSI
  CLHEP::Hep3Vector coe_vec( avrX, avrY, CalorimatorZ);
  m_coe = coe_vec;

  Double_t scaleFactor  = (avrZ-TargetZ)/(CalorimatorZ-TargetZ);

  avrX = avrX*scaleFactor;
  avrY = avrY*scaleFactor;


  // update Pi0 vars
  for( std::vector<Pi0>::iterator p=m_pi0.begin();
       p!=m_pi0.end(); p++ ) {
    //
    p->setVtx( avrX,avrY,avrZ );  // set Pi0 vertex
    p->updateVars();              // and update Pi0 vars
  }

  CLHEP::Hep3Vector g5_p3  =  CLHEP::Hep3Vector( m_g5.pos().x() - avrX, 
                         m_g5.pos().y() - avrY, 
                         m_g5.pos().z() - avrZ );
  //std::cout << "g5 z check : " << m_g5.pos().z() << std::endl;
  if(m_g5.pos().z() != 6168 && m_userflag>20160601){
    std::cout << "csi Z error" << std::endl;
    std::cout << "z= " << m_g5.pos().z() << std::endl;
    return -1;
  }
  g5_p3.setMag( m_g5.e() );
  m_g5.setP3( g5_p3 );


  // set klong variables
  m_chisqZ = 0.;
  CLHEP::HepLorentzVector p_kl = CLHEP::HepLorentzVector( 0.,0.,0.,0. );
  const Pi0& pi0 = m_pi0.front();
  double    maxZ = pi0.recZ();
  double    minZ = pi0.recZ();
  for( std::vector<Pi0>::iterator p=m_pi0.begin();
       p!=m_pi0.end(); p++ ) {
    // 4 momentum
    CLHEP::HepLorentzVector p_pi0( p->p3(), p->e() );
    p_kl += p_pi0;

    // z chi square
    m_chisqZ += (p->recZ()-avrZ)*(p->recZ()-avrZ)/p->recZsig2();

    //
    if( maxZ < p->recZ() )
      maxZ = p->recZ();
    //
    if( minZ > p->recZ() )
      minZ = p->recZ();
  }
  CLHEP::HepLorentzVector p_g5 ( m_g5.p3() , m_g5.e()  );
  p_kl += p_g5;

  m_p3 = p_kl.vect();
  m_energy = p_kl.e();
  m_mass   = p_kl.m();
  m_deltaZ = maxZ - minZ;
  m_v      = CLHEP::Hep3Vector( avrX, avrY, avrZ );

  if(m_debugLevel==1){
    std::cout << "x: " << avrX << std::endl;
    std::cout << "y: " << avrY << std::endl;
    std::cout << "z: " << avrZ << std::endl;
    std::cout << "chisqz: " << m_chisqZ << std::endl;
  }


  return 1;
}
