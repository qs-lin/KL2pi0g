#include <list>
#include <vector>

#include "TObject.h"
#include "TMath.h"

#include "gamma/Gamma.h"

Double_t GetTrigTime(const Int_t userflag, const bool isMC);

std::list<Gamma> user_trig(std::list<Gamma> const &glist0, const Int_t userflag, const bool isMC )
{
   if( glist0.size() < 5 ) return std::list<Gamma>();

   //std::cout << "test" << std::endl;
   std::list<Gamma> glist;
   const Double_t t0 = GetTrigTime(userflag, isMC );
   const Double_t width = 20.;
   std::vector<Double_t> tvec;

   for(std::list<Gamma>::const_iterator g = glist0.begin(); g!=glist0.end(); ++g )
   {
      if( TMath::Abs(g->t() - t0) < width ) tvec.push_back( g->t() );
   } 

   if( tvec.size()<5 ) return std::list<Gamma>();

   std::sort( tvec.begin(), tvec.end() );
   const Double_t MaxTimeDiff = 5.;
   std::vector<Double_t> tvec2;
   for(std::vector<Double_t>::const_iterator t = tvec.begin(); t!=tvec.end(); ++t )
   {
      if( tvec2.empty() ){
         tvec2.push_back(*t);
         continue;
      }
         
      const Double_t t_mean = TMath::Mean( tvec2.size(), &tvec2[0] );
      if( TMath::Abs( *t - t_mean ) < MaxTimeDiff ){
         tvec2.push_back(*t);
      }else{
         tvec2.clear();
      }
   }

   if( tvec2.size()!=5 ) return std::list<Gamma>();

   for(std::list<Gamma>::const_iterator g = glist0.begin(); g!=glist0.end(); ++g )
   {
      for(std::vector<Double_t>::const_iterator t = tvec2.begin(); t!=tvec2.end(); ++t )
      {
         if( TMath::Abs( g->t() - *t ) < 1.e-3 ) glist.push_back(*g);
      }
   }
   return glist;
}

Double_t GetTrigTime(const Int_t userflag, const bool isMC)
{
   if( isMC ) return 215.;
   if( userflag < 20180501 ) return 204.;
   if( userflag < 20181231 ) return 200.;
   return 200.;
}
