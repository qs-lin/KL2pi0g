#include <iostream>
#include <numeric>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <fstream>

#include "TObject.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

//#include "g5anaKL2pi0g/DataContainer.h"
#include "g5anaKL2pi0g/E14GNAnaDataContainer.h"
#include "g5anaKL2pi0g/KL2pi0g.h"
#include "g5anaKL2pi0g/RecKL2pi0g.h"
#include "gamma/GammaFinder.h"
#include "gamma/Gamma.h"

std::list<Gamma> user_trig(std::list<Gamma> const &glist0, const Int_t userflag, const bool isMC );
std::list<Gamma> user_trig_e14(std::list<Gamma> const &glist0);
bool user_rec(std::list<Gamma> const &glist, std::vector<KL2pi0g> &klvec, Int_t userflag);

int main( int argc, char** argv)
{
  // read argument
  if( argc!=4 ){
    std::cout<<" arg err : " << argv[0] << " (input) (output) (userflag)" << std::endl;
    return -1;
  }

  const std::string ifname = argv[1];
  const std::string ofname = argv[2];
  const Int_t     userflag = std::atoi(argv[3]);

  std::cout<<" input    : " << ifname << std::endl;
  std::cout<<" output   : " << ofname << std::endl;
  std::cout<<" userflag : " << userflag << std::endl;

  /// read input ///
  TChain *itree = new TChain("clusteringTree");
  itree->Add( ifname.c_str() );
  const Long64_t nentry = itree->GetEntries();
  std::cout<<" #entry = " << nentry << std::endl;
  if( nentry==0 ){
    return -1;
  }

  //DataContainer data;
  E14GNAnaDataContainer data;
  data.setBranchAddress(itree);

  /// prepare output ///

  TFile *ofile = new TFile( ofname.c_str(), "RECREATE");
  TTree *otree = new TTree("RecTree","output from g5anaKL2pi0g");
  //data.AddBranches(otree);
  data.branchOfKlong(otree);

  /// loop ///
  GammaFinder gFinder;

  std::ofstream ofiletxt("data.txt",std::ofstream::app);

  //for( Long64_t ientry=0; ientry<100; ++ientry )
  //for( Long64_t ientry=0; ientry<1000; ++ientry )
  for( Long64_t ientry=0; ientry<nentry; ++ientry )
  //for( Long64_t ientry=1303; ientry<1304; ++ientry )
  //for( Long64_t ientry=1303; ientry<1304; ++ientry )
  //for( Long64_t ientry=1555; ientry<1556; ++ientry )
  {

    if( nentry>100 && (ientry%(nentry/10)==0) )
      std::cout<<" " << 10*ientry/(nentry/10) << "%" << std::endl;

    itree->GetEntry(ientry);

    /// pre-selection ///
    //if( data.IsRealData() && !data.IsGoodQuality() ) continue;

    /// load clustering info ///
    std::list<Cluster> clist;
    data.getData(clist);

    std::list<Gamma> glist0;
    gFinder.findGamma(clist, glist0);

    /// Jay's trigger timing requirement ///
    //std::list<Gamma> glist = user_trig(glist0, userflag, !data.IsRealData() );
    /// e14's cluster creteria///
    std::list<Gamma> glist = user_trig_e14(glist0);

    if( glist.size()!=5 ) continue;

    /// reconstruction ///
    std::vector<KL2pi0g> klvec;
    if( !user_rec(glist,klvec,userflag) ) continue;


    //for testing purpose only. Check setData and getData for KL2pi0g

    //std::cout << ientry << std::endl;
    //std::cout << klvec.size() << std::endl;

    data.setData(klvec);
    for(int i=0; i<data.KlongNumber; ++i){
      ofiletxt << klvec[i] << std::endl;
    }

    //std::cout << data.KlongNumber << "  " << klvec.size() << std::endl; 
    otree->Fill();


  }//end of event loop

  // END
  
  ofiletxt.close();
  ofile->cd();
  otree->Write();
  ofile->Close();

  std::cout<<" END of KL2pi0g analysis. " << std::endl;
  return 1;

}
