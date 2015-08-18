 #!/bin/bash -l

# Arguments:  AnalyseHiCap_perChr chrname

chrname=$1

#Detect Promoter Interactions, process each chr separately

./AnalyseHiCap_perChr ../../Experiments_HiC.hicup.txt 0 1000 mE.HiC.hicup.thr0.$chrname\_ $chrname >mE.HiC.hicup.thr0.$chrname\.out &
#echo ./AnalyseHiCap_perChr ../../Experiments_HiC.hicup.txt 3 1000 mE_HiC_merged_hicup_$chrname\_ $chrname 
#./AnalyseHiCap_perChr ../../Experiments_HiC.hicup.txt 3 1000 mE_HiC_merged_hicup_$chrname\_ $chrname >mE_HiC_merged_hicup_$chrname\.out &


