struct InteractionStruct{
	string chr;
    bool dir;
	int resites[2]; // RE fragment start end
	int pos; //Actual read start
	int distance; // Distance from TSS
	int supp_pairs[2]; // the number of reads 
	char type; // U: upstream , D: downstream , X: inter-chromosomal
	bool peakoverlap; // if there is a peak in that RE fragment
	string peakprofile; // Peak profile in binary format
	double mappability;
	double p_val[2];
    double p_val2[2];
};

struct SignalStruct{ // This struct keeps the numbers of restriction fragments where there is at least one pair coming from a feature, one struct per chr
	boost::unordered::unordered_map<int, int* > signals; //key:REsite, [1..n]: number of times observed for experiments 1..n, before: frag starts before
};

struct SignalStruct_CTX{ // This struct keeps the numbers restriction fragments where there is at least one pair coming from a feature, one struct per chr
	boost::unordered::unordered_map<int, int* > signal;
	string maptochrname;
};
struct PeakMap{
	boost::unordered::unordered_set<int> peaks;
}; // each struct is a chromosome, chromosome names are mapped to struct indexes

struct MetaPeakMap{ // Key is the closest REsite to the peak, if there is a particular peak 1 otherwise 0
	boost::unordered::unordered_map<int, string > metapeaks;
};
boost::unordered::unordered_map<string, int> MetaPeakChrMap;


struct FeattoFeatSignalStruct{
	int feat_index; // within proms struct
	double normsignal[NOFEXPERIMENTS]; // The number of reads within the core promoter of the interactor promoter
};

struct FeatureStruct{
	vector <int> ProbeIDs;
	string FeatureType; //Promoter, NegCtrl, GWAS etc..
	string chr;
	int nofRESites; // how many restriction sites the core promoter contains, will be used to normalise the signal
	double featmappability;
	int *closestREsitenums; // offset of the closest RE sites on each side of the negctrl this will be used to index interactions

    //Used for Genes
    string RefSeqName;
    vector < string > TranscriptNames;
    vector < int > isoformpromotercoords;
    int TSS; //Transcription Start Site
    string strand;
    vector < double > expression;
    bool sharedpromoter;
    vector < string > genes_sharingproms;

	SignalStruct Signals; // Intra-chromosomal interactions
    vector < SignalStruct_CTX > Signals_CTX; // Each element of this vector will represent a chromosome
    vector < FeattoFeatSignalStruct > Inter_feature_ints; // Each element represent a promoter

    vector < InteractionStruct > interactions;
};


//For RE Sites
struct REindexes{
	vector <int> binstart;
	vector <int> binend;
	vector <int> offset;
	vector <int> count;
}; // each index struct will keep a chromosome

struct Mappindexes{
	vector <int> binstart;
	vector <int> binend;
	vector <int> offset;
	vector <int> count;
}; // each index struct will keep a chromosome


struct BG_signals{

	boost::unordered::unordered_map< int, double > mean;
	boost::unordered::unordered_map< int, double > stdev;


	double a[2]; // a[0] is upstream a[1] is downstream (power law fit)
	double b[2]; // b[0] is upstream b[1] is downstream (power law fit)
};