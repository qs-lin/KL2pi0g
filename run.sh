#I choose 23866 becuase this one will cause segmentation fault during calculation of shapechisq
#make my;
#./bin/g5anaKL2pi0g  /home/had/jaylin/workspace/myana/kpi0gg/pro/run79/mc/51kW/KL3pi0/v0/root/cluster/KL3pi0_0.root test_0.root 20180601
#./bin/g5anaKL2pi0g  /home/had/jaylin/workspace/myana/kpi0gg/pro/run79/mc/51kW/KL3pi0/v0/root/cluster/KL3pi0_23866.root test_23866.root 20180601

rm -f data*.txt
rm -f test.root

./bin/g5anaKL2pi0g  /gpfs/home/old/had/jaylin/workspace/myana/kpi0gg/pro/run79/mc/51kW/KL3pi0/v0/root/cluster/KL3pi0_999.root test.root 20180601 

#./bin/read 
