 #!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately
./AnalyseHiCap_perChr ../../Experiments.HiCap.hicup.txt 1 1000 mE.HiCap.hicup.thr1.$chrname\_ $chrname >mE.HiCap.hicup.thr1.$chrname\.out &
#./AnalyseHiCap_perChr ../../Experiments_HiC.hicup.txt 1 1000 mE.HiC.hicup.thr1.$chrname\_ $chrname >mE.HiC.hicup.thr1.$chrname\.out &


