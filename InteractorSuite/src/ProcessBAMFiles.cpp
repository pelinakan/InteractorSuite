#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <cstring>
#include <string>
#include <time.h>
#include <vector>
#include <sstream>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
using namespace std;

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

//ALGLIB PACKAGE HEADERS
#include "alglibmisc.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "dataanalysis.h"
#include "specialfunctions.h"
#include "solvers.h"
#include "optimization.h"
#include "diffequations.h"
#include "fasttransforms.h"
#include "integration.h"
#include "interpolation.h"
using namespace alglib_impl;


#ifdef _CHAR16T //redefined definition problem
#define CHAR16_T
#endif

#define BOOST_DISABLE_ASSERTS
#define NDEBUG

// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
#pragma warning(disable: 4018) // possible loss of data
#pragma warning(disable: 4244) // possible loss of data
#pragma warning(disable: 4715) // not all control paths return a value

//#define UNIX  //If run on UPPMAX
//#define WINDOWS // If run on WINDOWS
//#define GraphGC
#define CMB  //If run on DNA @ CMB

string dirname;
string ExpFileName;

int MinimumJunctionDistance; // To be entered by the user
int MinNumberofSupportingPairs; // To be entered by the user
int CellType; // 0:mES, 1:XEN, 2:TS // read from the experiments file
const int coreprom_upstream = 1000;
const int coreprom_downstream = 1000;
const int BinSize  = 1000; // Only To Calculate Background Interaction Frequencies
const int NOFEXPERIMENTS = 2; // Number of Experiments
int padding = 350; //For Sequence Capture Probes
double ExpressionThr = 2.0;
double Mappability_Threshold = 0.5;
string whichchr;

#include "linear.h"
#include "Data_Structures.h"
#include "GetOptions.h"
#include "SupplementaryFunctions.h"
#include "RESitesCount_MemAccess.h"
#include "Mappability.h"
#include "Promoters.h"
#include "Probes.h"
#include "Proximities.h"
#include "ProcessBAMFiles.h"
#include "BackgroundInteractionFrequency.h"
#include "Find_Interactions.h"

int main (int argc,char* argv[]){

// Take the parameters from command line
    string BaseFileName;
	if (argc < 6) {
		print_usage();
		return -1;
	}
	ExpFileName = argv[1];
	cout << ExpFileName << endl;
	MinNumberofSupportingPairs = atoi(argv[2]);
	MinNumberofSupportingPairs = double(MinNumberofSupportingPairs);
	MinimumJunctionDistance = atoi(argv[3]);
	BaseFileName = argv[4];
	whichchr = argv[5];

	cout << "Min Number of Supporting Pairs      " << MinNumberofSupportingPairs << endl;
	cout << "Min Junction Distance               " << MinimumJunctionDistance << endl;


//   --        INITIALISE CLASSES   --

	DetermineBackgroundLevels background;
	ProcessBAM bamfile;
    ProximityClass Proximities;

	DetectEnrichedBins EnrichedBins;
    
    ifstream ExpFile(ExpFileName.c_str());
    
    string BAMFILENAME, ExperimentName;
    vector < string > ExperimentNames;
    int ExperimentNo = 0;
    
    //First read the directory where the digested genome and annotation and probe files are
    ExpFile >> dirname;

	MappabilityClass mapp;
	RESitesClass dpnIIsites;
	dpnIIsites.InitialiseVars();
//	mapp.InitialiseVars();
//-------------------//------------------------------
    PromoterClass Promoters;
    Promoters.InitialiseData();
    Promoters.ReadPromoterAnnotation(dpnIIsites, mapp);
	cout << "Promoters Annotated" << endl;
    

	//ReadMetaPeakFile(); // If peaks are already processed.
    ProbeSet Probe_Set;
	Probe_Set.ReadProbeCoordinates();

    //Open the Experiments.txt file to read the bamfile names and other relevant information
    
	
    ExpFile >> BAMFILENAME >> ExperimentName >> CellType;
	do{ // Reads all the pairs in each experiment and fills the interaction maps
		ExperimentNames.push_back(ExperimentName);
		cout << BAMFILENAME << "     will be read" << endl;
        
        bamfile.ProcessTheBAMFile(Promoters, Probe_Set, Proximities, dpnIIsites, BAMFILENAME, ExperimentNo, whichchr);
		
		cout << "Detecting Interactions";
        background.CalculateMeanandStdRegress(Probe_Set, ExperimentName, ExperimentNo);
		
		EnrichedBins.DetectInteractions(Probe_Set,dpnIIsites,mapp,background,ExperimentName,ExperimentNo);
		
		cout << "  finished" << endl;
		
		cout << "Cleaning the background" << endl;
		background.bglevels.mean.clear();
		background.bglevels.stdev.clear();
		++ExperimentNo;
		cout << BAMFILENAME << "     finished" << endl;
		ExpFile >> BAMFILENAME >> ExperimentName >> CellType;
	}while(BAMFILENAME!="END");

	int NofofExperiments = ExperimentNo;

	EnrichedBins.PrintAllInteractions(Probe_Set, BaseFileName, NofofExperiments, ExperimentNames); //Print all same type of interactions

}
