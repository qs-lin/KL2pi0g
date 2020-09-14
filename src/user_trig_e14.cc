#include <list>
#include <numeric>
#include <vector>

#include "TObject.h"
#include "TMath.h"

#include "gamma/Gamma.h"

std::list<Gamma> user_trig_e14(std::list<Gamma> const &glist0)
{


  std::vector< double > clusterTimeVec;
  std::list< Gamma > glist;
  
  if(glist0.size() < 5) return glist;
  
  // If there are more than 5 clusters, I select most near time 5 clusters
  for( std::list<Gamma>::const_iterator it=glist0.begin(); it!=glist0.end(); it++ ){
    clusterTimeVec.push_back( it->t() );
    glist.push_back( *it );
  }

  while( glist.size() != 5 ){
    std::sort( clusterTimeVec.begin(), clusterTimeVec.end() );
    while( clusterTimeVec.size() != 5 ){
      double averageTime = std::accumulate( clusterTimeVec.begin(), clusterTimeVec.end(), 0.0);
      averageTime /= clusterTimeVec.size();
      if( averageTime-clusterTimeVec[0] > clusterTimeVec[clusterTimeVec.size()-1]-averageTime ){
        clusterTimeVec.erase( clusterTimeVec.begin() );
      }else{
        clusterTimeVec.pop_back();
      }
    }

    for( std::list<Gamma>::iterator it=glist.begin(); it!=glist.end(); ){
      bool rejectFlag = true;
      for( unsigned int index=0; index<clusterTimeVec.size(); index++ ){
        if( it->t() == clusterTimeVec[index] ){
          rejectFlag = false;
          break;
        }
      }
      if( rejectFlag ){
        it = glist.erase( it );
        continue;
      }
      it++;
    }
  }

  return glist;
}


