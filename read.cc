#include <iostream>
#include <numeric>
#include <cstdlib>
#include <string>
#include <algorithm>

#include "TObject.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include "g5anaKL2pi0g/E14GNAnaDataContainer.h"
#include "g5anaKL2pi0g/KL2pi0g.h"
#include "g5anaKL2pi0g/RecKL2pi0g.h"
#include "gamma/GammaFinder.h"
#include "gamma/Gamma.h"

int main( int argc, char** argv){

  TFile *infile = new TFile("test.root");
  //TFile *infile = new TFile("/gpfs/group/had/koto/ps3/klea/work/nobuhiro/data/g5ana_2pi0g_ovl/2pi0g_cls_g5ana_2pi0g_29_ovl.root");
  TTree *tr = (TTree*)infile->Get("RecTree");

  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );

  std::ofstream ofiletxt("data2.txt",std::ofstream::app);

  std::cout << "# of entries = " << tr->GetEntries() << std::endl;

  for( Long64_t ientry=0; ientry<tr->GetEntries(); ++ientry )
  {

    tr->GetEntry(ientry);
    std::vector<KL2pi0g> klist;
    data.getData(klist);


    for(unsigned int i=0; i<klist.size(); ++i){
      ofiletxt << klist[i] << std::endl;
    }




  }//end of event loop


  ofiletxt.close();

return 0;
}





