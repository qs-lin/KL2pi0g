#ifndef KLPI0GG_H
#define KLPI0GG_H

#include "TObject.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "pi0/Pi0.h"

class KL2pi0g
{

public: 
  // constructor
  KL2pi0g();

  // desctructor
  ~KL2pi0g();

  // operator
  bool operator<( const KL2pi0g& ) const;
  bool operator==( const KL2pi0g& ) const;

  enum{
    SORT_BY_MASS=1, SORT_BY_CHISQ
  };

  // extractor
  Int_t                           id()          const { return m_id; }
  const CLHEP::Hep3Vector&        v()           const { return m_v; }
  CLHEP::Hep3Vector&              v()                 { return m_v; }
  const CLHEP::Hep3Vector&        p3()          const { return m_p3; }
  CLHEP::Hep3Vector&              p3()                { return m_p3; }
  Double_t                        e()           const { return m_energy; }
  Double_t                        m()           const { return m_mass; }
  Double_t                        t()           const { return m_t; }
  Int_t                           status()      const { return m_status; }

  const std::vector<Pi0>&         pi0()         const { return m_pi0; }
  std::vector<Pi0>&               pi0()               { return m_pi0; }

  const Gamma&                    g5()          const { return m_g5; }
  Gamma&                          g5()                { return m_g5; }

  Double_t                        chisqZ()      const { return m_chisqZ; }
  const CLHEP::Hep3Vector&        coe()         const { return m_coe; }
  CLHEP::Hep3Vector&              coe()               { return m_coe; }
  Double_t                        csi_t()       const { return m_csi_t; }

  Double_t                        deltaZ()      const { return m_deltaZ; }
  Int_t                           sortFlag()    const { return m_sortFlag; }

  // method
  void SetId( const int id ){ m_id = id; }
  void SetStatus( const int status){ m_status = status; }
  void SetEnergy( const Double_t energy ) { m_energy = energy; }
  void SetMass( const Double_t mass ) { m_mass = mass; }
  void SetCoe( const CLHEP::Hep3Vector coe ){ m_coe = coe; }
  void SetKlongVertex( const CLHEP::Hep3Vector vtx ){ m_v = vtx; }
  void SetPi0( const Pi0 &pi0, const Pi0 &pi1 ){ m_pi0.clear(); m_pi0.push_back(pi0); m_pi0.push_back(pi1); } 
  void SetGamma( const Gamma &g5 ){ m_g5 = g5; }
  void SetChisqZ( const Double_t chisqZ ){ m_chisqZ = chisqZ; }
  void SetDeltaZ( const Double_t deltaZ){ m_deltaZ = deltaZ;}
  void SetShapeChisq( const Double_t chisq){ m_chisqshape = chisq;}
  void SetSortFlag( const Int_t sortFlag ){ m_sortFlag = sortFlag; }
  void SetTime( const Double_t t ){ m_t = t; }
  void SetCsiTime( const Double_t csi_t ){ m_csi_t = csi_t; }
  void SetUserFlag( const Int_t userflag){ m_userflag = userflag;}
  void SetDebugLevel( const Int_t debug ){m_debugLevel = debug;}

  void  SetVtx( const CLHEP::Hep3Vector& v ){m_v = v;}
  void  SetVtx( double x, double y, double z ) { m_v = CLHEP::Hep3Vector( x,y,z ); }
  void  SetP3( const CLHEP::Hep3Vector& p3 ) { m_p3 = p3; }
  void  SetP3( double px, double py, double pz ) { m_p3 = CLHEP::Hep3Vector( px,py,pz ); }
  void  SetCoe( const CLHEP::Hep3Vector& coe ){m_coe = coe;}
  void  SetCoe( double x, double y, double z ) { m_coe = CLHEP::Hep3Vector( x,y,z ); }

  void PrintId() const{
    std::cout << m_pi0[0].g1().id() << std::endl;
    std::cout << m_pi0[0].g2().id() << std::endl;
    std::cout << m_pi0[1].g1().id() << std::endl;
    std::cout << m_pi0[1].g2().id() << std::endl;
    std::cout << m_g5.id() << std::endl;
  }

  Int_t UpdateVars();

  // friend 
  friend std::ostream& operator<<( std::ostream& out, const KL2pi0g& kl2pi0g );

  KL2pi0g& operator=(const KL2pi0g& );

 private:
  Int_t    m_id;
  Int_t    m_status;
  Int_t    m_sortFlag;

  CLHEP::Hep3Vector m_coe;
  CLHEP::Hep3Vector m_v;
  CLHEP::Hep3Vector m_p3;

  Double_t         m_energy;
  Double_t         m_mass;
  Double_t         m_t;
  Double_t         m_csi_t;
  Double_t         m_chisqZ;
  Double_t         m_deltaZ;
  Double_t         m_chisqshape;
  Int_t            m_userflag;

  Gamma            m_g5;
  std::vector<Pi0> m_pi0;

  Int_t compare( const KL2pi0g& ) const;

  Int_t            m_debugLevel;




};

#endif
