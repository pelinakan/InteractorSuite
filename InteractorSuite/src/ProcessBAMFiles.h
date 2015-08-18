class ProcessBAM{
	friend class PromoterClass;
	friend class NegCtrlClass;
	friend class ProbeSet;
    friend class ProximityClass;
public:
    int NofPairs_One_on_Probe;
    int NofPairs_Both_on_Probe;
    int NofPairsNoAnn;
    void ProcessTheBAMFile(PromoterClass&, ProbeSet&, ProximityClass proximities, RESitesClass, string, int, string);
private:
    boost::unordered::unordered_map<int, string> RefIDtoChrNames;
	void ProcessPair(BamTools::BamAlignment&, BamTools::BamAlignment&, string, ProbeSet, PromoterClass, RESitesClass, ProximityClass proximities, int);
	void ProcessPair_CTX(BamTools::BamAlignment& al, BamTools::BamAlignment& almate, ProbeSet, PromoterClass, RESitesClass, ProximityClass proximities, int);
    bool AnnotatePair(string, string, int, int, int, int, ProbeSet&, PromoterClass&, RESitesClass, ProximityClass proximities, int);
};

void ProcessBAM::ProcessPair(BamTools::BamAlignment& al, BamTools::BamAlignment& almate, string whichchr, ProbeSet prs, PromoterClass promoters, RESitesClass dpnII, ProximityClass proximities, int ExperimentNo){
    bool onchr_1 = 0, onchr_2 = 0;
    boost::unordered::unordered_map<int, string>::const_iterator it1 = RefIDtoChrNames.find(al.RefID);
    if(it1 != RefIDtoChrNames.end()){
        if(it1->second == whichchr)
            onchr_1 = 1;
    }
    
    boost::unordered::unordered_map<int, string>::const_iterator it2 = RefIDtoChrNames.find(almate.RefID);
    if(it2 != RefIDtoChrNames.end()){
        if(it2->second == whichchr)
            onchr_2 = 1;
    }
    if(onchr_1 && onchr_2){
       // string chr1 = "chr", chr2 = "chr";
       // chr1.append(it1->second);
       // chr2.append(it2->second);
 //       cout << it1->second << "   " << it2->second << "   " << al.Position << "   " << almate.Position << " on chr" << endl;
        AnnotatePair(it1->second, it2->second, al.Position, almate.Position, al.Length, almate.Length, prs, promoters, dpnII, proximities, ExperimentNo);
    }
  //  cout << it1->second << "   " << it2->second << "   " << al.Position << "   " << almate.Position << endl;
}
void ProcessBAM::ProcessPair_CTX(BamTools::BamAlignment& al, BamTools::BamAlignment& almate, ProbeSet prs, PromoterClass promoters, RESitesClass dpnII, ProximityClass proximities, int ExperimentNo){
	
    string chr_1, chr_2;
	boost::unordered::unordered_map<int, string>::const_iterator it1 = RefIDtoChrNames.find(al.RefID);
	if(it1 != RefIDtoChrNames.end())
		chr_1 = it1->second;

	boost::unordered::unordered_map<int, string>::const_iterator it2 = RefIDtoChrNames.find(almate.RefID);
	if(it2 != RefIDtoChrNames.end())
		chr_2 = it2->second;

	if(chr_1 != chr_2){
        AnnotatePair(it1->second, it2->second, al.Position, almate.Position, al.Length, almate.Length, prs, promoters, dpnII, proximities, ExperimentNo);
        
	}

//		cout << it1->second << "   " << it2->second << "   " << al.Position << "   " << almate.Position << endl;
}

void ProcessBAM::ProcessTheBAMFile(PromoterClass& promoters, ProbeSet& prs, ProximityClass proximities, RESitesClass dpnII, string BAMFILENAME, int ExperimentNo, string whichchr){

	NofPairs_One_on_Probe = 0;
	NofPairs_Both_on_Probe = 0;
	NofPairsNoAnn = 0;
    

	BamReader reader;
	if ( !reader.Open(BAMFILENAME.c_str()) ) {
		cerr << "could not open input BAM files." << endl;
	}
 
	// retrieve 'metadata' from BAM files, these are required by BamWriter
	const SamHeader header = reader.GetHeader();
	const RefVector references = reader.GetReferenceData();

	// Make a map of chr names to RefIDs
	RefVector::const_iterator chrit;
	for (chrit = references.begin(); chrit != references.end(); ++chrit){
		int key = reader.GetReferenceID(chrit->RefName);
		RefIDtoChrNames[key] = chrit->RefName;
		//	cout << chrit->RefName << '\t' << key << endl;
	}
    BamAlignment al, almate;
    unsigned long int pcount = 0;
    if(whichchr == "CTX" || whichchr == "ctx" || whichchr == "Ctx"){ // Only process inter-chr interactions
		cout << "Reading BAM file for  inter-chromosomal pairs..." << endl;
       while ( reader.GetNextAlignmentCore(al) ){
			reader.GetNextAlignmentCore(almate);
			ProcessPair_CTX(al, almate, prs, promoters, dpnII, proximities, ExperimentNo);
        }
    }
	else{ // PROCESS THE CHROMOSOME GIVEN
        cout << "Reading BAM file for " << whichchr << " ..." << endl;
        
        if(reader.LocateIndex()){
            string indexfilename = BAMFILENAME;
            indexfilename.append(".bai");
            cout << indexfilename<< endl;
            reader.OpenIndex(indexfilename.c_str());
            int chrindex;
            boost::unordered::unordered_map<int, string>::const_iterator it1;
            for(it1 = RefIDtoChrNames.begin(); it1 != RefIDtoChrNames.end(); ++it1){
                if(it1->second == whichchr){
                    chrindex = it1->first;
                    break;
                }
            }
            
            reader.Jump(chrindex);
            cout << "jumped to " << whichchr << "  " << atoi(whichchr.c_str()) << endl;
        }
        while ( reader.GetNextAlignmentCore(al) ){
			reader.GetNextAlignmentCore(almate);
            ++pcount;
            if(pcount %50000 == 0)
                cout << pcount << " pairs read" << endl;
			ProcessPair(al, almate, whichchr, prs, promoters, dpnII, proximities, ExperimentNo);
	    }
    }
    cout << "Number of Pairs Both Reads on Probe             " << NofPairs_One_on_Probe << endl;
    cout << "Number of Pairs One Read on Probe               " << NofPairs_Both_on_Probe << endl;
    cout << "Number of Pairs None on Probe                   " << NofPairsNoAnn << endl;
    
}

bool ProcessBAM::AnnotatePair(string chr_1, string chr_2, int startcoord, int endcoord, int firstinpair_len, int secondinpair_len, ProbeSet& Probe_Set, PromoterClass &proms, RESitesClass dpnII, ProximityClass proximities, int ExperimentNo){
    
    bool annprom_firstread = 0;
    bool annprom_secondread = 0;
    
    int *renums1, *renums2;
    int resite_firstinpair, resite_secondinpair;
    renums1 = new int [2];	renums2 = new int [2];
    
    bool re1found = dpnII.GettheREPositions(chr_1, startcoord,renums1);
    bool re2found = dpnII.GettheREPositions(chr_2, (endcoord + 1),renums2);
    
  //  cout << "re coords " << chr_1 << "  " << startcoord << "  " << renums1[0] << "   " << chr_2 << "   " << endcoord << "  " << renums2[0] << endl;
    if(re1found)
        resite_firstinpair = renums1[0];
    else
        resite_firstinpair = startcoord;
    if(re2found)
        resite_secondinpair = renums2[0];
    else
        resite_secondinpair = endcoord;
    
    int probeindex_firstread = Probe_Set.AssociateReadwithProbes(chr_1, startcoord, (startcoord + firstinpair_len), proms); // First read start and end
    if (probeindex_firstread == -1)
        annprom_firstread = proms.AnnotatewithPromoters(chr_1, startcoord, firstinpair_len); // Check if it is close to a promoter without a probe
    
    int probeindex_secondread = Probe_Set.AssociateReadwithProbes(chr_2, endcoord, (endcoord + secondinpair_len), proms); //Second read start and end
    if (probeindex_secondread == -1)
        annprom_secondread = proms.AnnotatewithPromoters(chr_2, endcoord, (endcoord + secondinpair_len)); // Check if it is close to a promoter without a probe
    
    //cout << probeindex_firstread << "  " << probeindex_secondread << endl;
    
    if (probeindex_firstread == -1 && probeindex_secondread == -1){ // If none on probes
        ++NofPairsNoAnn;
        return 0;
    }
    if ((probeindex_firstread != -1) && (probeindex_secondread != -1)) { // If both on probes
        proximities.AnnotateFeatFeatInteraction(Probe_Set, probeindex_firstread, probeindex_secondread, ExperimentNo);
        proximities.AnnotateFeatFeatInteraction(Probe_Set, probeindex_secondread, probeindex_firstread, ExperimentNo);
        ++NofPairs_Both_on_Probe;
        return 1;
    }
    
    if (probeindex_firstread != -1 ) // If first read is annotated with a probe
        proximities.AnnotateDistalInteractor(Probe_Set, probeindex_firstread, chr_1, chr_2, resite_firstinpair, resite_secondinpair, ExperimentNo);
    else
        proximities.AnnotateDistalInteractor(Probe_Set, probeindex_secondread, chr_2, chr_1, resite_secondinpair, resite_firstinpair, ExperimentNo);
    
    ++NofPairs_One_on_Probe;
    return 1;
}


