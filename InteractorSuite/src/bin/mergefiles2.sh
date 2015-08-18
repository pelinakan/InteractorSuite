#! /bin/bash

basefilename=$1;
append=$2;

echo $basefilename;
space=" ";
star="*";
tempfile=".noheader";

	for file in `ls $basefilename$star$append`; 
		do filelist[file]=$file$space; 
	done;
	
	echo ${filelist[@]}
	asize=${#filelist[@]}
	echo $asize

	for ((j=0;j<asize;j++))
		do awk <${filelist[j]} 'NR!=1 {print}' >${filelist[j]}$tempfile;
	done;

	for file2 in `ls $basefilename$star$append.noheader`; 
		do filelist2+=($file2$space); 
	done;
	
	echo ${filelist2[@]}
	
	cat ${filelist2[@]} >merged.file

	awk <${filelist[1]} 'NR==1' >header.file
	cat  header.file merged.file >merged.file.wheader

	mv merged.file.wheader $basefilename$append

	rm merged.file
	rm header.file
	rm *.noheader
	unset file
	unset filelist
	unset asize
	unset file2
	unset filelist2

done;
########################
