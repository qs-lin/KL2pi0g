#include <iostream>
#include <vector>
#include <list>
#include "klong/RecKlong.h"

#include "TObject.h"

#include "gamma/Gamma.h"
#include "gnana/E14GNAnaFunction.h"
#include "ShapeChi2/ShapeChi2New.h"

#include "g5anaKL2pi0g/KL2pi0g.h"
#include "g5anaKL2pi0g/RecKL2pi0g.h"

bool user_rec(std::list<Gamma> const &glist, std::vector<KL2pi0g> &klvec, Int_t userflag )
{

  RecKL2pi0g recK5g;
  //recK5g.SetDebugLevel(1);
  klvec = recK5g.recK2pi0g(glist,userflag);

  if( klvec.size()==0 ) return false;
  //std::cout << klvec[0] << std::endl;

  /// gamma position correction with angle dependence
  E14GNAnaFunction::getFunction()->correctPosition( klvec[0].pi0()[0] );
  E14GNAnaFunction::getFunction()->correctPosition( klvec[0].pi0()[1] );
  E14GNAnaFunction::getFunction()->correctPosition( klvec[0].g5()  );

  /// gamma energy correction with angle dependence
  E14GNAnaFunction::getFunction()->correctEnergyWithAngle( klvec[0].pi0()[0] );
  E14GNAnaFunction::getFunction()->correctEnergyWithAngle( klvec[0].pi0()[1] );
  E14GNAnaFunction::getFunction()->correctEnergyWithAngle( klvec[0].g5()  );


  /// re-reconstruction with corrected gamma
  std::list<Gamma> glist2;
  glist2.push_back( klvec[0].pi0()[0].g1() );
  glist2.push_back( klvec[0].pi0()[0].g2() );
  glist2.push_back( klvec[0].pi0()[1].g1() );
  glist2.push_back( klvec[0].pi0()[1].g2() );
  glist2.push_back( klvec[0].g5() );

  klvec.clear();
  klvec = recK5g.recK2pi0g(glist2,userflag);
  if(klvec.size()==0) return false;

  static ShapeChi2New chi2Newcalc;
  chi2Newcalc.shapeChi2( klvec[0].pi0()[0].g1() );
  chi2Newcalc.shapeChi2( klvec[0].pi0()[0].g2() );
  chi2Newcalc.shapeChi2( klvec[0].pi0()[1].g1() );
  chi2Newcalc.shapeChi2( klvec[0].pi0()[1].g2() );
  chi2Newcalc.shapeChi2( klvec[0].g5() );


  return true;
}
