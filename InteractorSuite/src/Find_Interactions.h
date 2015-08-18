class DetectEnrichedBins{
public:
	void DetectInteractions(ProbeSet&, RESitesClass&, MappabilityClass&, DetermineBackgroundLevels, string, int); // Find bins that have above background level interactions
	void PrintAllInteractions(ProbeSet, string, int,vector < string >);
private:
	bool CheckSupportingPairs(int*);
	double CalculatepVal(boost::unordered::unordered_map< int, double >,boost::unordered::unordered_map< int, double >, int,int);
	void FillInteractionStruct(RESitesClass&, MappabilityClass&, vector<InteractionStruct>&,string, int, int, int, int*,double, double, char,bool);
    int FindClosestTranscriptTSS(int, vector<int>);
};

bool DetectEnrichedBins::CheckSupportingPairs(int* supppairs){
	bool recordit = 0;
	for (int t = 0; t < NOFEXPERIMENTS; ++t){
		if (supppairs[t] >= MinNumberofSupportingPairs){
			recordit = 1;
			break;
		}
	}
	return recordit;
}
double DetectEnrichedBins::CalculatepVal(boost::unordered::unordered_map< int, double > mean,boost::unordered::unordered_map< int, double > std, int bin,int sig){
	double q = 0, p_value = 1;
	if(mean.find(bin) == mean.end()){
		mean[bin] = 0;
		std[bin] = 1;
		p_value = 1;
	}
	else{
        if (std[bin] == 0.0)
            q = 0;
        else
            q = alglib::normaldistribution(( sig - mean[bin]) / std[bin]); //z-score
		p_value = 1- q;
	}
	return p_value;
}

void DetectEnrichedBins::FillInteractionStruct(RESitesClass& dpnIIsites, MappabilityClass& mapp, vector<InteractionStruct>& interactionstr,string chr, int dist, int ExperimentNo, int resite, int* values, double p_value, double p_value2, char inttype, bool direction){
	vector<InteractionStruct>::iterator interactorit;
	bool addedbefore = 0;
	for(interactorit = interactionstr.begin(); interactorit < interactionstr.end(); ++interactorit){
        if(interactorit->resites[0] == resite){ //Added before
            addedbefore = 1;
            interactorit->supp_pairs[ExperimentNo] = values[ExperimentNo];
            interactorit->p_val[ExperimentNo] = p_value;
            interactorit->p_val2[ExperimentNo] = p_value2;
        }
	}
	if(!addedbefore){
		interactionstr.push_back(InteractionStruct());
        interactionstr.back().resites[0] = resite;
		interactionstr.back().peakprofile.clear();
		interactionstr.back().p_val[ExperimentNo] = p_value;
        interactionstr.back().p_val2[ExperimentNo] = p_value2;
		interactionstr.back().supp_pairs[ExperimentNo] = values[ExperimentNo];
		interactionstr.back().chr = chr;
		interactionstr.back().distance = dist;
        int *renums;
        renums = new int [2];
        bool refound = dpnIIsites.GettheREPositions(chr,(resite + 1),renums); // Interactor RE fragment
        if(refound)
            interactionstr.back().resites[1] = renums[1];
        else{
            cout << "could not find resite  " << resite << endl;
            interactionstr.back().resites[1] = resite;
        }
//        interactionstr.back().mappability = mapp.GetMappability(chr,interactionstr.back().resites[0],interactionstr.back().resites[1]);
	}

}
void DetectEnrichedBins::DetectInteractions(ProbeSet& prs, RESitesClass& dpnIIsites, MappabilityClass& mapp, DetermineBackgroundLevels background, string BaseFileName,int ExperimentNo){

	for(int i = 0; i < prs.Probes.size(); ++i){
		boost::unordered::unordered_map<int, int* >::const_iterator it; //key:REpos
        for (it = prs.Probes[i].feature.Signals.signals.begin(); it != prs.Probes[i].feature.Signals.signals.end(); ++it){
			bool enoughpairs = 0;
			int dist = 0;
			enoughpairs = CheckSupportingPairs(it->second);
            
			if (enoughpairs){
                dist = it->first - prs.Probes[i].end;
                int bin = abs(dist) / BinSize;
                double p_value = CalculatepVal(background.bglevels.mean,background.bglevels.stdev, bin,it->second[ExperimentNo]);
                double p_value2 = 1;
                if ((it->second[ExperimentNo] - 1) >= 0 )
                    p_value2 = CalculatepVal(background.bglevels.mean,background.bglevels.stdev,bin,(it->second[ExperimentNo] - 1));
                FillInteractionStruct(dpnIIsites, mapp, prs.Probes[i].feature.interactions, prs.Probes[i].chr, dist, ExperimentNo, it->first, it->second, p_value, p_value2, 'U',0);
            }
        }
	// Inter-chromosomal
        for( int j = 0; j < prs.Probes[i].feature.Signals_CTX.size(); ++j){ // For each chromosome
            for (it = prs.Probes[i].feature.Signals_CTX[j].signal.begin(); it != prs.Probes[i].feature.Signals_CTX[j].signal.end(); ++it){
                bool enoughpairs = 0;
                int dist = 0;
                enoughpairs = CheckSupportingPairs(it->second);
                if (enoughpairs){
                    FillInteractionStruct(dpnIIsites, mapp, prs.Probes[i].feature.interactions,prs.Probes[i].feature.Signals_CTX[j].maptochrname, -1, ExperimentNo, it->first, it->second, -1, -1, 'X',0);
                }
            }
        }
    }
}

void DetectEnrichedBins::PrintAllInteractions(ProbeSet prs, string BaseFileName, int NumberofExperiments, vector < string > ExperimentNames){

//Probe to Distal Interactions
	string FileName;
    
	FileName.append(BaseFileName);
	FileName.append("SignificantInteractions");
	FileName.append("_AllProbes");
	FileName.append(".txt");
	ofstream outf1(FileName.c_str());

	outf1 << "Probe Chr" << '\t' << "Probe Start" << '\t' << "Probe End"    << '\t' << "Probe ID" << '\t'
          << "Target ID" << '\t' << "Probe Name"  << '\t' << "Annotation"   << '\t'
          << "Closest Transcript ID" << '\t' << "TSS of Closest Transcript" << '\t';
    
    outf1 << "Interactor Chr" << '\t' << "Interactor Start" << '\t' << "Interactor End" << '\t' << "abs(distance)" << '\t';
    
	for (int e = 0; e < NumberofExperiments; ++e)
		outf1 << ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_pval" << '\t' << "_pval2" << '\t';
    outf1 << endl;

    vector< CaptureProbes >::iterator probeit;
	for(probeit = prs.Probes.begin() ; probeit != prs.Probes.end(); ++probeit){
        for(int j = 0; j < probeit->feature.interactions.size(); ++j){
            outf1 << probeit->chr << '\t' << probeit->start << '\t' << probeit->end << '\t' << probeit->probeid << '\t'
                  << probeit->target_id << '\t' << probeit->name << '\t' << probeit->annotated << '\t';
            
            if(probeit->annotated == 1) // refseq
            outf1 << probeit->feature.TranscriptNames[probeit->closestTSS_index] << '\t'
                  << probeit->feature.isoformpromotercoords[probeit->closestTSS_index] << '\t';
            else
                outf1 << "N/A" << '\t' << 0 << '\t';

            outf1 << probeit->feature.interactions[j].chr << '\t' << probeit->feature.interactions[j].resites[0] << '\t'
                  << probeit->feature.interactions[j].resites[1] << '\t' << probeit->feature.interactions[j].distance << '\t';
            
			for (int e = 0; e < NumberofExperiments; ++e)
				outf1 << probeit->feature.interactions[j].supp_pairs[e] << '\t' << probeit->feature.interactions[j].p_val[e] << '\t'
                      << probeit->feature.interactions[j].p_val2[e] << '\t';
			outf1 << endl;
		}
	}
	outf1.close();

//Promoter-promoter Interactions
	string FileName3;
	FileName3.append(BaseFileName);
	FileName3.append("Interactions_BetweenProbes");
	FileName3.append(".txt");
	ofstream outf3(FileName3.c_str());

    outf3 << "Probe Chr" << '\t' << "Probe Start" << '\t' << "Probe End"    << '\t' << "Probe ID" << '\t'
          << "Target ID" << '\t' << "Probe Name"  << '\t' << "Annotation"   << '\t';
    outf3 << "Interacting Probe Chr" << '\t' << "Interacting Probe Start" << '\t' << "Interacting Probe End"    << '\t' << "Interacting Probe ID" << '\t'
          << "Interacting Target ID" << '\t' << "Interacting Probe Name"  << '\t' << "Interacting Annotation"   << '\t';

    for (int e = 0; e < NumberofExperiments; ++e)
        outf3 << ExperimentNames[e] << "_SuppPairs" << '\t';
    outf3 << endl;

	for(probeit = prs.Probes.begin() ; probeit != prs.Probes.end(); ++probeit){
		vector< FeattoFeatSignalStruct >::const_iterator it; //first: REpos, second: signal 
		for (it = probeit->feature.Inter_feature_ints.begin(); it != probeit->feature.Inter_feature_ints.end(); ++it){
			int pindex = it->feat_index;
			outf3 << probeit->chr << '\t' << probeit->start << '\t' << probeit->end << '\t' << probeit->probeid << '\t'
                  << probeit->target_id << '\t' << probeit->name << '\t' << probeit->annotated << '\t';

            outf3 << prs.Probes[pindex].chr << '\t' << prs.Probes[pindex].start << '\t' << prs.Probes[pindex].end << '\t' << prs.Probes[pindex].probeid << '\t'
                  << prs.Probes[pindex].target_id << '\t' << prs.Probes[pindex].name << '\t' << prs.Probes[pindex].annotated << '\t';
            
            for (int e = 0; e < NumberofExperiments; ++e)
                outf3 << '\t' << it->normsignal[e];
            outf3 << endl;
		}
	}
	outf3.close();

}

